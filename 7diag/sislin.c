#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "sislin.h"

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline double genRandomA(uint i, uint j, uint k){
    static double invRandMax = 1.0 / (double)RAND_MAX;
    //Eh sempre diagonal dominante
    return ( (i==j) ? (double)(k<<1) : 1.0 )  * (double)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline double genRandomB(uint k){
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return (double)(k<<2) * (double)random() * invRandMax;
}

static inline void initDiag(struct diagMat* A, uint k, uint n){

    A->Diags = aligned_alloc(32, sizeof(double*)*k*n);
    A->n = n;

}

static double sqrVector(double* x, double* y, uint n){
    double sqrVector = 0;
    for (uint i = 0; i < n; i++){
        sqrVector += x[i]*y[i];
    }
    return sqrVector;
}



void genKDiagonal(struct LinearSis *SL, uint k, uint n){
   
    uint half = k/2;
    SL->k = k;
    SL->A = aligned_alloc(32, sizeof(struct diagMat));
    SL->A->Diags = aligned_alloc(32,sizeof(double*)*k*n);
    SL->A->n = n;
    SL->b = aligned_alloc(32,sizeof(double)*n);


    for (uint i = 0; i < k; i++){
        uint diagonal = (i*n); 
        int offset = (i-half);
        uint j;
        for (uint i = 0; i < k; i++){
            uint diagonal = (i*n); 
            int offset = (i-half);
            uint j;
            for (j = 0; j < n - n%UNROLL ; j += UNROLL){
                int j0 = (int)j + offset;
                int j1 = (int)j + offset + 1;
                int j2 = (int)j + offset + 2;
                int j3 = (int)j + offset + 3;
                int j4 = (int)j + offset + 4;
                int j5 = (int)j + offset + 5;
                int j6 = (int)j + offset + 6;
                int j7 = (int)j + offset + 7;

                SL->A->Diags[diagonal+j  ] = (j0 >= 0 && j0 < (int)n) ? genRandomA(j0, j, k) : 0.0;
                SL->A->Diags[diagonal+j+1] = (j1 >= 0 && j1 < (int)n) ? genRandomA(j1, j, k) : 0.0;
                SL->A->Diags[diagonal+j+2] = (j2 >= 0 && j2 < (int)n) ? genRandomA(j2, j, k) : 0.0;
                SL->A->Diags[diagonal+j+3] = (j3 >= 0 && j3 < (int)n) ? genRandomA(j3, j, k) : 0.0;
                SL->A->Diags[diagonal+j+4] = (j4 >= 0 && j4 < (int)n) ? genRandomA(j4, j, k) : 0.0;
                SL->A->Diags[diagonal+j+5] = (j5 >= 0 && j5 < (int)n) ? genRandomA(j5, j, k) : 0.0;
                SL->A->Diags[diagonal+j+6] = (j6 >= 0 && j6 < (int)n) ? genRandomA(j6, j, k) : 0.0;
                SL->A->Diags[diagonal+j+7] = (j7 >= 0 && j7 < (int)n) ? genRandomA(j7, j, k) : 0.0;
            } 
            for (; j < n; j++) {
                int jj = (int)j + offset;
                SL->A->Diags[diagonal + j] = (jj >= 0 && jj < (int)n) ? genRandomA(jj, j, k) : 0.0;
            }
        }
    }
    uint j = 0;
    for(j; j<n%UNROLL; j+=UNROLL){
        SL->b[j] = genRandomB(k);
        SL->b[j+1] = genRandomB(k);
        SL->b[j+2] = genRandomB(k);
        SL->b[j+3] = genRandomB(k);
        SL->b[j+4] = genRandomB(k);
        SL->b[j+5] = genRandomB(k);
        SL->b[j+6] = genRandomB(k);
        SL->b[j+7] = genRandomB(k);
    }

    for (j; j < n ; j++){
        SL->b[j] = genRandomB(k);
    }
   /*int o = 3;
    for (; j < 4; j++){
        double soma = o;
        SL->b[j] = o;
        for (int i = 0; i < o; i++){
            soma += n-1;
            SL->b[j] += soma;
        }
        o++;
    }
    for(; j < n - 3; j++){
        SL->b[j] = SL->b[j-1]+7;
    }

    double sub = 7;
    for( ;j< n; j++){
        SL->b[j] = SL->b[j-1]- sub;
        sub += n+1;
    }*/

}

int genSymmetricPositive(struct LinearSis *restrict SL, struct diagMat *restrict ASP, double *restrict bsp, double *restrict time)
{
    struct diagMat* A = SL->A;
    double* b = SL->b;
    uint n = A->n;

    uint linha = 0;
    double linhaV[7];
    int offset = (linha-3);
    double soma = 0.0;
    initDiag(ASP,13,n);

    for(linha; linha <= 3; linha++){
        bsp[linha] = 0;
        for(uint i = 0; i <= 3+linha; i++){
            linhaV[i] = A->Diags[(3-linha+i)*n+linha];
            bsp[linha] += A->Diags[(3-linha+i)*n+linha]*b[i];
        }
        for(uint col = 0; col < n; col++){
            soma = 0.0;
            for(uint i = 0; i <= 3+linha; i++){
                if(3+i-col >= 0 && 3+i-col < 7 ){
                    soma += linhaV[i] * A->Diags[(3+i-col)*n+col];
                }
            }
            if(6-col+linha >= 0 && 6-col+linha < n ){
                ASP->Diags[(6-col+linha)*n+col] = soma;
            }
        }
    }
    for(linha; linha < n-3; linha++){
        bsp[linha] = 0.0;   
        uint teste = 2*(linha - 3);
        for(uint i = 0; i < 7; i++){
            linhaV[i] = A->Diags[i*n+linha];
            bsp[linha] += A->Diags[i*n+linha]*b[i+linha-3];
        }
        for(uint col = 0; col < n; col++){
            soma = 0.0;
            for(uint i = 0; i <= 3+linha; i++){
                if(linha+i-col >= 0 && linha+i-col < 7 && i < 7 && linha-3 < n){
                    soma += linhaV[i] * A->Diags[(linha+i-col)*n+col];
                }
            }
            if(6-col+linha >= 0 && 6-col+linha < 13 && col < n){
                ASP->Diags[(6-col+linha)*n+col] = soma;
            }
        }
    }

    for(linha; linha < n; linha++){
        bsp[linha] = 0.0;
        for(uint i = 0; i < n+3-linha; i++){
            linhaV[i] = A->Diags[(i)*n+linha];
            bsp[linha] += A->Diags[(i)*n+linha]*b[linha-3+i];
        }
        for(uint col = 0; col <= n; col++){
            soma = 0.0;
            for(uint i = 0; i <= n+2-linha; i++){
                if(6-col+i >= 0 && 6-col+i < 7){
                    soma += linhaV[i] * A->Diags[(6-col+i)*n+(linha-6+col)];
                }
            }
            if(12-col >= 0 && 12-col < 13 && linha+col-6 < n){
                ASP->Diags[(12-col)*n+(linha+col-6)] = soma;
            }
        }
        
    }
}   


/*Um pre condicionamento melhora um SL simetrico, positivo, definido e mal condicionado.*/
int genPreCond(struct diagMat *restrict A, double w, uint n,
        double *restrict M, double *restrict time)
{   
    if (w == -1)
        for(uint i = 0; i < n; i++){
            if(A->Diags[6*n+i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[6)*n+%d] = %f -> A não é definida positiva\n", i, A->Diags[4*n+i]);
                return -2; 
            }
            else
                M[i] = 1.0;
        }
    
    if(time)
        *time = timestamp();
    if (w == 0)
        for(uint i = 0; i < n; i++){
            if(A->Diags[6*n+i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[6)*n+%d] = %f -> A não é definida positiva\n", i, A->Diags[4*n+i]);
                return -2; 
            }
            else
                M[i] = 1.0/A->Diags[6*n+i];
        }
    
    if (time)
        *time = timestamp() - *time;
    return 0;
}

int conjGradientPre(struct diagMat *restrict A, double *restrict B, double *restrict x, double *restrict r, double *restrict M, double *restrict time){
    //d = v1; c = v2
    calcResidue(A, B, 13, x, r, NULL);

    uint n = A->n;

    // Y para calcular o SL com condicionador
    double *Yv = aligned_alloc(32,n * sizeof(double)); 

    
    //Criando a matriz d e c usados para calculos
    double *d = aligned_alloc(32,n * sizeof(double));
    double *c = aligned_alloc(32,n * sizeof(double));

    if(!d || !Yv || !c){
        free(d);
        free(Yv);
        free(c);
        fprintf(stderr, "Falha na alocação de memória\n");
        return -1;
    }
    uint i = 0;
    for (; i < n- n%UNROLL; i+=UNROLL){
        d[i  ] = 0.0;
        d[i+1] = 0.0;
        d[i+2] = 0.0;
        d[i+3] = 0.0;
        d[i+4] = 0.0;
        d[i+5] = 0.0;
        d[i+6] = 0.0;
        d[i+7] = 0.0;
    }
    for(;i<n;i++){
        d[i] = 0.0;
    }

    // Y para calcular o SL com condicionador
    for (i = 0; i < n - n%UNROLL; i+=UNROLL)
        Yv[i  ] = M[i  ] * r[i  ]; // y = M^-1 * r
        Yv[i+1] = M[i+1] * r[i+1];
        Yv[i+2] = M[i+2] * r[i+2];
        Yv[i+3] = M[i+3] * r[i+3];
        Yv[i+4] = M[i+4] * r[i+4];
        Yv[i+5] = M[i+5] * r[i+5];
        Yv[i+6] = M[i+6] * r[i+6];
        Yv[i+7] = M[i+7] * r[i+7];

    for (; i < n; i++)
        Yv[i] = M[i] * r[i]; // y = M^-1 * r


    //Criando a matriz d e c usados para calculos
    for (i = 0; i < n-n%UNROLL; i+=UNROLL){
        d[i  ] = Yv[i  ];    
        d[i+1] = Yv[i+1];   
        d[i+2] = Yv[i+2];   
        d[i+3] = Yv[i+3];   
        d[i+4] = Yv[i+4];   
        d[i+5] = Yv[i+5];   
        d[i+6] = Yv[i+6];   
        d[i+7] = Yv[i+7];    
    }
    for (; i < n; i++)
        d[i] = Yv[i];

    

    double cAd = 0.0; //dkt * Adk
    double alpha; // ak
    double deltaOld = 0.0;
    double deltaNew = 0.0;
    double beta = 0.0;
    double tIter = timestamp();
    uint it = 1;
    uint j = 0;
    deltaOld = sqrVector(r, Yv, n); //Como rk * rkt eh o quadrado nao precisamos multplicar matrizes

    do {
        //Calculando ak = rk * rkt / dkt * A * dk 
        multMatVet(A, d, c, 13); //Precisamos multiplicar matrizes pois A e d nao sao quadrados

       
        //Fazemos a multiplicacao de matrizes manual, ja q n vou criar outro vetor
        cAd = 0.0;
        for (j = 0; j < n - n%UNROLL; j+=UNROLL){
            cAd += c[j]*d[j];
            cAd += c[j+1]*d[j+1];
            cAd += c[j+2]*d[j+2];
            cAd += c[j+3]*d[j+3];
            cAd += c[j+4]*d[j+4];
            cAd += c[j+5]*d[j+5];
            cAd += c[j+6]*d[j+6];
            cAd += c[j+7]*d[j+7];
        }
        for (; j < n; j++){
            cAd += c[j]*d[j];
        }

        if(cAd == 0){
            free(Yv);
            free(d);
            free(c);
            fprintf(stderr,"cAD Divisão por zero\n");
            return -1;
        }
        alpha = deltaOld / cAd; //Calculando ak

        deltaNew = 0.0;
        for(j = 0;j < n%UNROLL; j+=8){
            
            x[j] += alpha * d[j]; //Xk+1 = Xk + akdk
            r[j] -= alpha * c[j]; //rk+1 = rk - akAdk
            Yv[j] = M[j] * r[j];  // y = M⁻¹ * r
            deltaNew += Yv[j] * r[j];

            x[j+1] += alpha * d[j+1]; //Xk+1 = Xk + akdk
            r[j+1] -= alpha * c[j+1]; //rk+1 = rk - akAdk
            Yv[j+1] = M[j+1] * r[j+1];  // y = M⁻¹ * r
            deltaNew += Yv[j+1] * r[j+1];

            x[j+2] += alpha * d[j+2]; //Xk+1 = Xk + akdk
            r[j+2] -= alpha * c[j+2]; //rk+1 = rk - akAdk
            Yv[j+2] = M[j+2] * r[j+2];  // y = M⁻¹ * r
            deltaNew += Yv[j+2] * r[j+2];

            x[j+3] += alpha * d[j+3]; //Xk+1 = Xk + akdk
            r[j+3] -= alpha * c[j+3]; //rk+1 = rk - akAdk
            Yv[j+3] = M[j+3] * r[j+3];  // y = M⁻¹ * r
            deltaNew += Yv[j+3] * r[j+3];

            x[j+4] += alpha * d[j+4]; //Xk+1 = Xk + akdk
            r[j+4] -= alpha * c[j+4]; //rk+1 = rk - akAdk
            Yv[j+4] = M[j+4] * r[j+4];  // y = M⁻¹ * r
            deltaNew += Yv[j+4] * r[j+4];

            x[j+5] += alpha * d[j+5]; //Xk+1 = Xk + akdk
            r[j+5] -= alpha * c[j+5]; //rk+1 = rk - akAdk
            Yv[j+5] = M[j+5] * r[j+5];  // y = M⁻¹ * r
            deltaNew += Yv[j+5] * r[j+5];

            x[j+6] += alpha * d[j+6]; //Xk+1 = Xk + akdk
            r[j+6] -= alpha * c[j+6]; //rk+1 = rk - akAdk
            Yv[j+6] = M[j+6] * r[j+6];  // y = M⁻¹ * r
            deltaNew += Yv[j+6] * r[j+6];

            x[j+7] += alpha * d[j+7]; //Xk+1 = Xk + akdk
            r[j+7] -= alpha * c[j+7]; //rk+1 = rk - akAdk
            Yv[j+7] = M[j+7] * r[j+7];  // y = M⁻¹ * r
            deltaNew += Yv[j+7] * r[j+7];

        }
        for (; j < n; j++) {

            //Xk+1 = Xk + akdk
            x[j] += alpha * d[j];
            //rk+1 = rk - akAdk
            r[j] -= alpha * c[j];
            Yv[j] = M[j] * r[j]; // y = M⁻¹ * r

            deltaNew += Yv[j] * r[j];
        }

        if(deltaOld == 0){
            fprintf(stderr,"DELTA Divisão por zero\n");
            return -1;
        }
        beta = deltaNew / deltaOld;
        for(j = 0; j < n-n%UNROLL; j+=UNROLL){
            d[j  ] = Yv[j  ] + beta * d[j  ];
            d[j+1] = Yv[j+1] + beta * d[j+1];
            d[j+2] = Yv[j+2] + beta * d[j+2];
            d[j+3] = Yv[j+3] + beta * d[j+3];
            d[j+4] = Yv[j+4] + beta * d[j+4];
            d[j+5] = Yv[j+5] + beta * d[j+5];
            d[j+6] = Yv[j+6] + beta * d[j+6];
            d[j+7] = Yv[j+7] + beta * d[j+7];
        }
        for(;j<n;j++){
            d[j] = Yv[j] + beta * d[j];
        }
        deltaOld = deltaNew;

        it++;
        tIter = timestamp() - tIter;
    }while (it < 25);

    *time = tIter/it;

    free(Yv);
    free(d);
    free(c);

    return 0;
}

double calcNormaMax(double *x,double* y, uint n){
    double max = 0.0;
    double aux = 0.0;

    for (uint i = 0; i < n; i++){
        aux = fabs(x[i] - y[i]);
        if (max < aux)
            max = aux;
    }
    return max;
}

double calcNormaEuclidiana(double *x, uint n){
    double aux = 0.0;

    for (uint i = 0; i < n; i++){
        aux += x[i]*x[i];
        //printf("aux(%f) = (%f)*(%f)\n",aux, x[i], x[i]);
    }

    return sqrt(aux);
}

void calcResidue(struct diagMat *restrict A,double *restrict B,uint k, double *restrict x, double *restrict r, double *restrict time)
{ 
    uint n = A->n;

    if (time)
        *time = timestamp();


    double* sum = aligned_alloc(32,sizeof(double)*n);
    multMatVet(A, x, sum, k);

    uint i = 0;
    for (; i < n - n%UNROLL; i+=UNROLL){
        r[i]=B[i]-sum[i];
        r[i+1]=B[i+1]-sum[i+1];
        r[i+2]=B[i+2]-sum[i+2];
        r[i+3]=B[i+3]-sum[i+3];
        r[i+4]=B[i+4]-sum[i+4];
        r[i+5]=B[i+5]-sum[i+5];
        r[i+6]=B[i+6]-sum[i+6];
        r[i+7]=B[i+7]-sum[i+7];
    }
    for(;i<n;i++){
        r[i]=B[i]-sum[i];
    }

    if (time)
        *time = timestamp() - *time;
}

void printVetor(double* vet, uint n){
    for(uint i = 0; i < n; i++){
        printf("%.8g  ", vet[i]);
    }

    printf("\n");
}

void printSis(struct LinearSis *SL){
    uint n = SL->n;
    struct diagMat* A = SL->A;
    uint half = SL->k/2;
    
    for(uint i = 2; i < n-2; i++){
        for(uint j = 0; j < SL->k; j++){
            printf("%f", A->Diags[j*n+(i+(j-half))]);
        }
        printf("\n");
    }
}

void print7Diag(struct diagMat *SL, uint k){

    uint n = SL->n;
    uint half = k/2;
    int j =0;
    uint i = 0;
    
    for(i = 0; i < half; i++){
        for(j = k-1; j >= 0; j--){
            if(i-(j-half) >= 0 && i-(j-half) < n){
                printf("  %2d|%.8g",j, SL->Diags[j*n+(i-(j-half))]);
            } 
        }
        for (j = 0 ;j < n-1-i-half; j++){
            printf("       ");
        }
        printf("  ]\n");
    }

    for(i; i<n-half; i++){
        for (j = 0 ;j < i-half; j++){
            printf("       ");
        }
        for(j = k-1; j >= 0; j--){
            if(i-(j-half) >= 0 && i-(j-half) < n){
                printf("  %2d|%.8g",j, SL->Diags[j*n+(i-(j-half))]);
            }
        }
        for (j = 0 ;j < n-1-i-half; j++){
            printf("       ");
        }
        printf("  ]\n");
    }
    for(i; i<n; i++){
        for (j = 0 ;j < i-half; j++){
            printf("       ");
        }
        for(j = k-1; j >= 0; j--){
            if(i-(j-half) >= 0 && i-(j-half) < n){
                printf("  %2d|%.8g",j, SL->Diags[j*n+(i-(j-half))]);
            } 
        }
    printf("  ]\n");
    }

    /*for(uint i = 0; i < k; i++){
        printf("diag %d:   ",i);
        for(uint j = 0; j < n; j++){
            printf("%f  ",SL->Diags[i*n+j]);
        }
        printf("\n");
    }*/
}


/*ESSA FUNCAO MAIS GERAL FUNCIONA PARA MATRIZES E VETORES, MAS DEPENDE DE COMO O VETOR ESTA ORGANIZADO */
void static multMatVet(struct diagMat *A, double *B, double *C, uint k) {
    uint n = A->n;
    uint linha = 0;
    uint half = k/2;

    for(linha; linha <= half; linha++){
        C[linha] = 0.0;
        for(uint i = 0; i <= half+linha; i++){
            C[linha] += A->Diags[(half+linha-i)*n+i]*B[i];
        }
    }

    for(linha; linha < n-half; linha++){
        C[linha] = 0.0;
        for(uint i = 0; i < k; i++){
            C[linha] += A->Diags[(k-1-i)*n+(i+linha-half)] * B[i+linha-half];
        }
    }
    for(linha; linha < n; linha++){
        C[linha] = 0.0;
        for(uint i = 0; i < n+half-linha; i++){
            C[linha] += A->Diags[(k-1-i)*n+linha-half+i] * B[linha-half+i];  
        }
    }
}
