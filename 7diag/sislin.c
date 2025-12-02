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

void inline static initDiag(struct diagMat* A, uint k, uint n){

    A->Diags = malloc(sizeof(double*)*k*n);
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
    SL->A =malloc(sizeof(struct diagMat));
    initDiag(SL->A, k, n);
    SL->b = malloc(sizeof(double)*n);


    for (uint i = 0; i < k; i++){
        int offset = (i-half);
        for (uint j = 0; j < n ; j++){
            if (j+offset >= 0 && j+offset < n){
                //SL->A->Diags[i*n+j] = i*n+j;
                SL->A->Diags[i*n+j] = genRandomA(j+offset,j,k);
            } else{
                SL->A->Diags[i*n+j] = 0.0;
            }
        } 
    }
    for (uint j = 0; j < n ; j++){
        SL->b[j] = 1; 
        SL->b[j] = genRandomB(k);
    }
    /*SL->b[0] = 90;
    SL->b[1] = 150;
    SL->b[2] = 225;
    SL->b[3] = 315;
    SL->b[4] = 322;
    SL->b[5] = 329;
    SL->b[6] = 336;
    SL->b[7] = 343;
    SL->b[8] = 350;
    SL->b[9] = 357;
    SL->b[10] = 364;
    SL->b[11] = 357;
    SL->b[12] = 335;
    SL->b[13] = 298;*/

}

int genSymmetricPositive(struct LinearSis *SL, struct diagMat *ASP, double *bsp, double *time)
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
            //printf("linha[%d](%f) = A[%d][%d](%f)\n",i,linhaV[i],3-linha+i,linha,A->Diags[(3-linha+i)*n+linha]);
        }
        //printf("\n");
        //printf("linha %d:  ", linha);
        //printVetor(linhaV, 7);
        for(uint col = 0; col <= n; col++){
            soma = 0.0;
            //printf("col: %d\n", col);
            for(uint i = 0; i <= 3+linha; i++){
                //printf("    i: %d\n", i);
                if(3+i-col >= 0 && 3+i-col < 7 ){
                    soma += linhaV[i] * A->Diags[(3+i-col)*n+col];
                    //printf("soma: %f = linhav[%d](%f) * diag[%d][%d](%f)\n",soma,i,linhaV[i], 3+i-col,col, A->Diags[(3+i-col)*n+col]);
                }
            }
            if(6-col+linha >= 0 && 6-col+linha < n ){
                //printf("diag[%d][%d] = %f\n",6-col+linha,col,soma);
                ASP->Diags[(6-col+linha)*n+col] = soma;
            }
            //printf("\n");
        }
    }
    for(linha; linha < n-3; linha++){
        bsp[linha] = 0.0;   
        uint teste = 2*(linha - 3);
        for(uint i = 0; i < 7; i++){
            linhaV[i] = A->Diags[i*n+linha];
            bsp[linha] += A->Diags[i*n+linha]*b[i+linha-3];
            //printf("bsp[%d](%f) = A[%d][%d](%f) * b[%d](%f)\n",linha,bsp[linha],i,linha,A->Diags[i*n+linha], i+linha-3,b[i+linha-3]);
        }
        //printf("linha %d:  ", linha);
        //printVetor(linhaV, 7);
        //printf("bsp[%d](%f)\n",linha,bsp[linha]);
        for(uint col = 0; col <= n; col++){
            soma = 0.0;
            for(uint i = 0; i <= 3+linha; i++){
                if(linha+i-col >= 0 && linha+i-col < 7 && i < 7 && linha-3 < n){
                    soma += linhaV[i] * A->Diags[(linha+i-col)*n+col];
                    //printf("linha %d soma: %f = linhav[%d](%f) * diag[%d][%d](%f)\n",linha,soma,i,linhaV[i], linha+i-col, col, A->Diags[(linha+i-col)*n+col]);
                }
            }
            if(6-col+linha >= 0 && 6-col+linha < 13 && col < n){
                //printf("diagASP[%d][%d] = %f\n",6-col+linha,col,soma);
                ASP->Diags[(6-col+linha)*n+col] = soma;
            }
            //printf("\n");
        }
    }

    for(linha; linha < n; linha++){
        bsp[linha] = 0.0;
        for(uint i = 0; i < n+3-linha; i++){
            linhaV[i] = A->Diags[(i)*n+linha];
            bsp[linha] += A->Diags[(i)*n+linha]*b[linha-3+i];
            //printf("bsp[%d](%f) = A[%d][%d](%f) * b[%d](%f)\n",linha,bsp[linha],i,linha-i,A->Diags[(i)*n+linha-i], i+linha-3,b[i+linha-3]);
        }
        //printf("\n");
        //printf("linha %d:  ", linha);
        //printVetor(linhaV, 7);
        for(uint col = 0; col <= n; col++){
            soma = 0.0;
            for(uint i = 0; i <= n+2-linha; i++){
                //printf("linha %d = linhav[%d] * diag[%d][%d]\n", linha, i, col-i, linha-3+i);
                if(6-col+i >= 0 && 6-col+i < 7){
                    soma += linhaV[i] * A->Diags[(6-col+i)*n+(linha-6+col)];
                    //printf("linha %d soma(%f)= linhav[%d](%f) * diag[%d][%d](%f)\n", linha,soma, i,linhaV[i], 6-col+i, linha-6+col,A->Diags[(6-col+i)*n+linha-6+col]);
                }
            }
            if(12-col >= 0 && 12-col < 13 && linha+col-6 < n){
               //printf("diag[%d][%d] = %f\n",12-col,linha+col-6,soma);
                ASP->Diags[(12-col)*n+(linha+col-6)] = soma;
            }
            //printf("\n");
        }
        
    }
}   


/*Um pre condicionamento melhora um SL simetrico, positivo, definido e mal condicionado.*/
int genPreCond(struct diagMat *A, double w, uint n,
        double* M, double *time)
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

int conjGradientPre(struct diagMat *A, double *B, double *x, double *r,double *norma, double *M, double *time){
    //d = v1; c = v2
    calcResidue(A, B, 13, x, r, NULL);

    uint n = A->n;

    // Y para calcular o SL com condicionador
    double *Yv = malloc(n * sizeof(double)); 

    
    //Criando a matriz d e c usados para calculos
    double *d = calloc(n,sizeof(double));
    double *c = malloc(n * sizeof(double));

    if(!d || !Yv || !c){
        free(d);
        free(Yv);
        free(c);
        fprintf(stderr, "Falha na alocação de memória\n");
        return -1;
    }

    // Y para calcular o SL com condicionador
    for (uint i = 0; i < n; i++)
        Yv[i] = M[i] * r[i]; // y = M^-1 * r

    //Criando a matriz d e c usados para calculos
    for (uint i = 0; i < n; i++)
        d[i] = Yv[i];

    

    double diff = 0.0;
    double cAd = 0.0; //dkt * Adk
    double alpha; // ak
    double deltaOld = 0.0;
    double deltaNew = 0.0;
    double valueNew = 0.0;
    double beta = 0.0;
    double tIter = timestamp();
    uint it = 1;
    deltaOld = sqrVector(r, Yv, n); //Como rk * rkt eh o quadrado nao precisamos multplicar matrizes

    do {
        //Calculando ak = rk * rkt / dkt * A * dk 
        multMatVet(A, d, c, 13); //Precisamos multiplicar matrizes pois A e d nao sao quadrados
       
        //Fazemos a multiplicacao de matrizes manual, ja q n vou criar outro vetor
        cAd = 0.0;
        for (uint i = 0; i < n; i++) 
            cAd += c[i]*d[i];

        if(cAd == 0){
            free(Yv);
            free(d);
            free(c);
            fprintf(stderr,"cAD Divisão por zero\n");
            return -1;
        }
        alpha = deltaOld / cAd; //Calculando ak
       
        deltaNew = 0.0;
        for (uint i = 0; i < n; i++) {

            //Xk+1 = Xk + akdk
            x[i] += alpha * d[i];
            //rk+1 = rk - akAdk
            r[i] -= alpha * c[i];
            Yv[i] = M[i] * r[i]; // y = M⁻¹ * r

            deltaNew += Yv[i] * r[i];
            valueNew += r[i] * r[i];
        }

        if(deltaOld == 0){
            fprintf(stderr,"DELTA Divisão por zero\n");
            return -1;
        }
        beta = deltaNew / deltaOld;
        for(uint i = 0; i < n; i++)
            d[i] = Yv[i] + beta * d[i];

        deltaOld = deltaNew;

        it++;
        tIter = timestamp() - tIter;
        //printf("it: %d", it);
        //printVetor(r,A->n);
    }while (it < 25);

    *norma = diff;
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
    }

    return sqrt(aux);
}

void calcResidue(struct diagMat *A,double *B,uint k, double *x, double *r, double *time)
{ 
    uint n = A->n;

    if (time)
        *time = timestamp();


    double* sum = malloc(sizeof(double)*n);
    //printf("\n\nRESIDUO\n");
    multMatVet(A, x, sum, k);

    for (uint i = 0; i<n; i++){
        r[i]=B[i]-sum[i];
        //printf("r[%d](%f) = B[%d](%f) - r[%d](%f)\n",i,r[i],i,B[i],i,sum[i]);
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
            printf("      ");
        }
        printf("  ]\n");
    }

    for(i; i<n-half; i++){
        for (j = 0 ;j < i-half; j++){
            printf("      ");
        }
        for(j = k-1; j >= 0; j--){
            if(i-(j-half) >= 0 && i-(j-half) < n){
                printf("  %2d|%.8g",j, SL->Diags[j*n+(i-(j-half))]);
            }
        }
        for (j = 0 ;j < n-1-i-half; j++){
            printf("      ");
        }
        printf("  ]\n");
    }
    for(i; i<n; i++){
        for (j = 0 ;j < i-half; j++){
            printf("      ");
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
            //printf("C[%d](%.10g) = A[%d][%d](%.10g)*X[%d](%.10g)\n",linha,C[linha],(half+linha-i),+i, A->Diags[(half+linha-i)*n+i], i ,B[i]);
        }
        //printf("\n");
    }

    for(linha; linha < n-half; linha++){
        C[linha] = 0.0;
        for(uint i = 0; i < k; i++){
            C[linha] += A->Diags[(k-1-i)*n+(i+linha-half)] * B[i+linha-half];
            //printf("C[%d](%.10g) = A[%d][%d](%.10g)*X[%d](%.10g)\n",linha,C[linha],(k-1-i),(i+linha-half), A->Diags[(k-1-i)*n+(i+linha-half)], i+linha-half,B[i+linha-half]);
        }
        //printf("\n");
    }
    for(linha; linha < n; linha++){
        C[linha] = 0.0;
        for(uint i = 0; i < n+half-linha; i++){
            C[linha] += A->Diags[(k-1-i)*n+linha-half+i] * B[linha-half+i];  
            //printf("C[%d](%.10g) = A[%d][%d](%.10g)*X[%d](%.10g)\n",linha,C[linha],(k-1-i),linha-half+i, A->Diags[(k-1-i)*n+linha-half+i], linha-half+i,B[linha-half+i]);
        }
        //printf("\n");
    }
}
