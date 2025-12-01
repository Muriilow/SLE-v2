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

static double sqrVector(double* x, double* y, uint n){
    double sqrVector = 0;
    for (uint i = 0; i < n; i++){
        sqrVector += x[i]*y[i];
    }
    return sqrVector;
}

void genKDiagonal(struct LinearSis *SL){
    int k = SL->k;
    int n = SL->n;
    int half = k/2;
    SL->A = malloc(sizeof (struct diagMat));
    SL->A->n = n;

    SL->A->Diags = malloc(sizeof(double*)*k);
    for(int i = 0; i < k; i++)
        SL->A->Diags[i] = malloc(sizeof(double)*n);
    SL->b = malloc(sizeof(double)*n);


    for (int i = 0; i < k; i++){
        int offset = (i-half);
        for (int j = 0; j < n ; j++){
            if (j+offset >= 0 && j+offset < n){
                SL->A->Diags[i][j] = i+1;
                //SL->A->Diags[i][j] = genRandomA(j+offset,j,k);
            } else{
                SL->A->Diags[i][j] = 0.0;
            }
        } 
    }
    for (int j = 0; j < n ; j++){
        SL->b[j] = j+1; 
        //SL->b[j] = genRandomB(k);
    }
}

int genSymmetricPositive(struct LinearSis *SL, struct diagMat *ASP, double *bsp, double *time)
{
    struct diagMat* A = SL->A;
    double* b = SL->b;
    int n = A->n;

    int half = 3;
    int linha = 0;
    double linhaV[7];
    int offset = (linha-half);
    double soma = 0.0;
    initDiag(ASP,13,n);

    for(linha; linha <= half; linha++){
        bsp[linha] = 0;
        for(int i = 0; i <= half+linha; i++){
            linhaV[i] = A->Diags[half+linha-i][i];
            bsp[linha] += linhaV[i]*b[i];
        }
        //printf("bsp[%d](%f)\n",linha,bsp[linha]);
        //printf("linha %d:  ", linha);
        //printVetor(linhaV, 7);
        for(int col = 0; col <= n; col++){

            soma = 0.0;
            for(int i = 0; i <= half+linha; i++){
                if(half+col-i >= 0 && half+col-i < 7 ){
                    soma += linhaV[i] * A->Diags[half+col-i][i];
                    //printf("soma: %f = linhav[%d](%f) * diag[%d][%d](%f)\n",soma,i,linhaV[i], half+col-i, i, A->Diags[half+col-i][i]);
                }
            }
            if(6-col+linha >= 0 && 6-col+linha < n ){
                //printf("diag[%d][%d] = %f\n",6-col+linha,col,soma);
                ASP->Diags[6-col+linha][col] = soma;
            }
        }
    }
    for(linha; linha < n-half; linha++){
        int teste = 2*(linha - half);
        for(int i = 0; i < 7; i++){
            linhaV[i] = A->Diags[6-i][i+linha-half];
            bsp[linha] += linhaV[i]*b[i+linha-half];
        }
        
        //printf("bsp[%d](%f)\n",linha,bsp[linha]);
        for(int col = 0; col <= n; col++){
            soma = 0.0;
            for(int i = 0; i <= half+linha; i++){
                if(linha-teste+col-i >= 0 && linha-teste+col-i < 7 && i < 7 && linha-half+i < n){
                    soma += linhaV[i] * A->Diags[linha-teste+col-i][linha-3+i];
                    //printf("linha %d soma: %f = linhav[%d](%f) * diag[%d][%d](%f)\n",linha,soma,i,linhaV[i], linha-teste+col-i, linha-half+i, A->Diags[linha-teste+col-i][linha-half+i]);
                }
            }
            if(6-col+linha >= 0 && 6-col+linha < 13 && col < n){
                //printf("diagASP[%d][%d] = %f\n",6-col+linha,col,soma);
                ASP->Diags[6-col+linha][col] = soma;
            }
        }
    }

    for(linha; linha < n; linha++){
        for(int i = 0; i < n+half-linha; i++){
            linhaV[i] = A->Diags[6-i][i];
            bsp[linha] += linhaV[i]*b[linha-half+i];
            //printf("bsp[%d]: %f",linha, bsp[linha]);
        }
        //printf("bsp[%d](%f)\n",linha,bsp[linha]);
        //printf("linha %d:  ", linha);
        //printVetor(linhaV, 7);
        for(int col = 0; col <= n; col++){
            soma = 0.0;
            for(int i = 0; i <= n-1+half-linha; i++){
                //printf("linha %d = linhav[%d] * diag[%d][%d]\n", linha, i, col-i, linha-half+i);
                if(col-i >= 0 && col-i < 7){
                    soma += linhaV[i] * A->Diags[col-i][linha-half+i];
                    //printf("linha %d soma(%f)= linhav[%d](%f) * diag[%d][%d](%f)\n", linha,soma, i,linhaV[i], col-i, linha-half+i,A->Diags[col-i][linha-half+i]);
                }
            }
            //printf("diag[%d][%d] = %f\n",12-col,linha+col-6,soma);
            if(12-col >= 0 && 12-col < 13 && linha+col-6 < n){
                ASP->Diags[12-col][linha+col-6] = soma;
            }
        }
    }
}   


/*Um pre condicionamento melhora um SL simetrico, positivo, definido e mal condicionado.*/
int genPreCond(struct diagMat *A, double w, int n,
        double* M, double *time)
{   
    if (w == -1)
        for(int i = 0; i < n; i++){
            if(A->Diags[6][i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[6][%d] = %f -> A não é definida positiva\n", i, A->Diags[4][i]);
                return -2; 
            }
            else
                M[i] = 1.0;
        }
    
    if(time)
        *time = timestamp();
    if (w == 0)
        for(int i = 0; i < n; i++){
            if(A->Diags[6][i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[6][%d] = %f -> A não é definida positiva\n", i, A->Diags[4][i]);
                return -2; 
            }
            else
                M[i] = 1.0/A->Diags[6][i];
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
    for (int i = 0; i < n; i++)
        Yv[i] = M[i] * r[i]; // y = M^-1 * r

    //Criando a matriz d e c usados para calculos
    for (int i = 0; i < n; i++)
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
        for (int i = 0; i < n; i++) 
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
    }while (it < 25);

    *norma = diff;
    *time = tIter/it;

    free(Yv);
    free(d);
    free(c);

    return 0;
}

double calcNormaMax(double *x,double* y, int n){
    double max = 0.0;
    double aux = 0.0;

    for (int i = 0; i < n; i++){
        aux = fabs(x[i] - y[i]);
        if (max < aux)
            max = aux;
    }
    return max;
}

double calcNormaEuclidiana(double *x, int n){
    double aux = 0.0;

    for (int i = 0; i < n; i++){
        aux += x[i]*x[i];
    }

    return sqrt(aux);
}

void calcResidue(struct diagMat *A,double *B,int k, double *x, double *r, double *time)
{ 
    int n = A->n;

    if (time)
        *time = timestamp();

    multMatVet(A, x, r, k);

    for (int i = 0; i<n; i++)
        r[i]=B[i]-r[i];

    if (time)
        *time = timestamp() - *time;
}

void printVetor(double* vet, int n){
    for(int i = 0; i < n; i++){
        printf("%.8g  ", vet[i]);
    }

    printf("\n");
}

void printSis(struct LinearSis *SL){
    uint n = SL->n;
    int half = SL->k/2;
    
    for(int i = 2; i < n-2; i++){
        for(int j = 0; j < SL->k; j++){
            printf("%f", SL->A->Diags[j][i+(j-half)]);
        }
        printf("\n");
    }
}

void print7Diag(struct diagMat *SL, int k){

    uint n = SL->n;
    int half = k/2;
    int j = k-1;
    int i = 0;
    
    for(i = 0; i < half; i++){
        for(j = k-1; j >= 0; j--){
            if(i-(j-half) >= 0 && i-(j-half) < n){
                printf("  %2d|%.8g",j, SL->Diags[j][i-(j-half)]);
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
                printf("  %2d|%.8g",j, SL->Diags[j][i-(j-half)]);
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
                printf("  %2d|%.8g",j, SL->Diags[j][i-(j-half)]);
            } 
        }
    printf("  ]\n");
    }

    /*for(int i = 0; i < k; i++){
        printf("diag %d:   ",i);
        for(int j = 0; j < n; j++){
            printf("%f  ",SL->Diags[i][j]);
        }
        printf("\n");
    }*/
}


/*ESSA FUNCAO MAIS GERAL FUNCIONA PARA MATRIZES E VETORES, MAS DEPENDE DE COMO O VETOR ESTA ORGANIZADO */
void multMatVet(struct diagMat *A, double *B, double *C, int k) {
    int n = A->n;
    int linha = 0;
    int half = k/2;

    for(linha; linha <= half; linha++){
        C[linha] = 0.0;
        for(int i = 0; i <= half+linha; i++){
            C[linha] += A->Diags[half+linha-i][i]*B[i];
        }
    }

    for(linha; linha < n-half; linha++){
        C[linha] = 0.0;
        for(int i = 0; i < k; i++){
            C[linha] += A->Diags[k-1-i][i+linha-half] * B[i+linha-half];
        }
    }
    for(linha; linha < n; linha++){
        C[linha] = 0.0;
        for(int i = 0; i < n+half-linha; i++){
            C[linha] += A->Diags[k-1-i][i] * B[linha-half+i];
        }
    }
}


void initDiag(struct diagMat* A, int k, int n){

    A->Diags = malloc(sizeof(double*)*k);
    A->n = n;
    for(int i = 0; i < k; i++)
        A->Diags[i] = malloc(sizeof(double)*n);

}