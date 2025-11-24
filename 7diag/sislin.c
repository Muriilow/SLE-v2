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

void gen7Diagonal(struct diag7 *SL){
    int n = SL->n;
    SL->D1 = malloc(sizeof(double)*n);
    SL->D2 = malloc(sizeof(double)*n);
    SL->D3 = malloc(sizeof(double)*n);
    SL->D4 = malloc(sizeof(double)*n);
    SL->D5 = malloc(sizeof(double)*n);
    SL->D6 = malloc(sizeof(double)*n);
    SL->D7 = malloc(sizeof(double)*n);
    SL->B = malloc(sizeof(double)*n);

    SL->D1[0] = 0.0;
    SL->D1[1] = 0.0;
    SL->D1[2] = 0.0;

    SL->D2[0] = 0.0;
    SL->D2[1] = 0.0;
    SL->D2[2] = genRandomA(0,2,7);
    
    SL->D3[0] = 0.0;
    SL->D3[1] = genRandomA(0,1,7);
    SL->D3[2] = genRandomA(1,2,7);

    SL->D4[0] = genRandomA(0,0,7);
    SL->D4[1] = genRandomA(1,1,7);
    SL->D4[2] = genRandomA(2,2,7);

    SL->D5[0] = genRandomA(1,0,7);
    SL->D5[1] = genRandomA(2,1,7);
    SL->D5[2] = genRandomA(3,2,7);

    SL->D6[0] = genRandomA(2,0,7);
    SL->D6[1] = genRandomA(3,1,7);
    SL->D6[2] = genRandomA(4,2,7);

    SL->D7[0] = genRandomA(3,0,7);
    SL->D7[1] = genRandomA(4,1,7);
    SL->D7[2] = genRandomA(5,2,7);

    SL->B[0] = genRandomB(7);
    SL->B[1] = genRandomB(7);
    SL->B[2] = genRandomB(7);

    for(int i = 3; i < n-3; i++){
        SL->B[i] = genRandomB(7);
        SL->D1[i] = genRandomA(i-3,i,7);
        SL->D2[i] = genRandomA(i-2,i,7);
        SL->D3[i] = genRandomA(i-1,i,7);
        SL->D4[i] = genRandomA(i,i,7);
        SL->D5[i] = genRandomA(i+1,i,7);
        SL->D6[i] = genRandomA(i+2,i,7);
        SL->D7[i] = genRandomA(i+3,i,7);
    }

    SL->D1[n-3] = genRandomA(n-6,n-3,7);
    SL->D1[n-2] = genRandomA(n-5,n-2,7);
    SL->D1[n-1] = genRandomA(n-4,n-1,7);

    SL->D2[n-3] = genRandomA(n-5,n-3,7);
    SL->D2[n-2] = genRandomA(n-4,n-2,7);
    SL->D2[n-1] = genRandomA(n-3,n-1,7);
    
    SL->D3[n-3] = genRandomA(n-4,n-3,7);
    SL->D3[n-2] = genRandomA(n-3,n-2,7);
    SL->D3[n-1] = genRandomA(n-2,n-1,7);

    SL->D4[n-3] = genRandomA(n-3,n-3,7);
    SL->D4[n-2] = genRandomA(n-2,n-2,7);
    SL->D4[n-1] = genRandomA(n-1,n-1,7);

    SL->D5[n-3] = genRandomA(n-2,n-3,7);
    SL->D5[n-2] = genRandomA(n-1,n-2,7);
    SL->D5[n-1] = 0.;

    SL->D6[n-3] = genRandomA(n,n-2,7);
    SL->D6[n-2] = 0.0;
    SL->D6[n-1] = 0.0;

    SL->D7[n-3] = 0.0;
    SL->D7[n-2] = 0.0;
    SL->D7[n-1] = 0.0;

    SL->B[n-3] = genRandomB(7);
    SL->B[n-2] = genRandomB(7);
    SL->B[n-1] = genRandomB(7);
    
}

int genSymmetricPositive(struct LinearSis *SL, struct Matrix *ASP, struct Matrix *bsp, double *time)
{
    struct Matrix *A = SL->A;
    struct Matrix *b = SL->b;

    *time = timestamp();

    double *Atv = malloc(A->row*A->column*sizeof(double));
    struct Matrix AT ={Atv, A->row, A->column, A->k};
    if(!Atv){
        free(Atv);
        fprintf(stderr, "Falha na alocação de memória\n");
        return -1;
    }

    //genTranspose(&AT, A);
    
    multMatrix(&AT, A, ASP);
    multMatrix(&AT, b, bsp);

    *time = timestamp() - *time;
    free(Atv);
    return 0;
}

/*Um pre condicionamento melhora um SL simetrico, positivo, definido e mal condicionado.*/
int genPreCond(struct Matrix *A, double w, int n, int k,
        struct Matrix *M, double *time)
{   
    if (w == -1)
        for(int i = 0; i < M->row; i++){
            if(A->v[i*n + i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[%d][%d] = %f -> A não é definida positiva\n", i, i, A->v[i*n + i]);
                return -2; 
            }
            else
                M->v[i] = 1.0;
        }

    *time = timestamp();
    if (w == 0){
        for(int i = 0; i < M->row; i++){
            if(A->v[i*n + i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[%d][%d] = %f -> A não é definida positiva\n", i, i, A->v[i*n + i]);
                return -2; 
            }
            else
                M->v[i] = 1.0/A->v[i*n + i];
        }
    }

    *time = timestamp() - *time;
    return 0;
}

void genTranspose(struct diag7 *A, struct diag7 *T)
{   
    *T->D1 = *A->D7;
    *T->D2 = *A->D6;
    *T->D3 = *A->D5;
}

int conjGradientPre(struct LinearSis *SL, double *x, double *r,double *norma, struct Matrix *M, double *time){

    calcResidue(SL, x, r, NULL);
    uint n = SL->n;

    // Y para calcular o SL com condicionador
    double *Yv = malloc(n * sizeof(double)); 
    struct Matrix y = {Yv, 1, SL->n, 0};
    struct Matrix rMatrix = {r, 1, SL->n, 0};

    
    //Criando a matriz d e c usados para calculos
    double *v1 = calloc(n,sizeof(double));
    struct Matrix d = {v1, n, 1, 0};
    double *v2 = malloc(n * sizeof(double));
    struct Matrix c = {v2, n, 1, 0};

    if(!v1 || !Yv || !v2){
        free(v1);
        free(Yv);
        free(v2);
        fprintf(stderr, "Falha na alocação de memória\n");
        return -1;
    }

    // Y para calcular o SL com condicionador
    for (int i = 0; i < n; i++)
        y.v[i] = M->v[i] * r[i]; // y = M^-1 * r

    //Criando a matriz d e c usados para calculos
    for (int i = 0; i < n; i++)
        d.v[i] = y.v[i];

    

    double diff = 0.0;
    double cAd = 0.0; //dkt * Adk
    double alpha; // ak
    double deltaOld = 0.0;
    double deltaNew = 0.0;
    double valueNew = 0.0;
    double beta = 0.0;
    double tIter = timestamp();
    uint it = 1;
    deltaOld = sqrVector(r, y.v, n); //Como rk * rkt eh o quadrado nao precisamos multplicar matrizes

    do {
        /*Calculando ak = rk * rkt / dkt * A * dk */
        multMatrix(SL->A, &d, &c); //Precisamos multiplicar matrizes pois A e d nao sao quadrados
       
        //Fazemos a multiplicacao de matrizes manual, ja q n vou criar outro vetor
        cAd = 0.0;
        for (int i = 0; i < n; i++) 
            cAd += c.v[i]*d.v[i];

        if(cAd == 0){
            free(Yv);
            free(v1);
            free(v2);
            fprintf(stderr,"cAD Divisão por zero\n");
            return -1;
        }
        alpha = deltaOld / cAd; //Calculando ak
       
        deltaNew = 0.0;
        for (uint i = 0; i < n; i++) {

            /*Xk+1 = Xk + akdk*/
            x[i] += alpha * d.v[i];
            /*rk+1 = rk - akAdk*/
            r[i] -= alpha * c.v[i];
            y.v[i] = M->v[i] * r[i]; // y = M⁻¹ * r

            deltaNew += y.v[i] * r[i];
            valueNew += r[i] * r[i];
        }

        if(deltaOld == 0){
            fprintf(stderr,"DELTA Divisão por zero\n");
            return -1;
        }
        beta = deltaNew / deltaOld;
        for(uint i = 0; i < n; i++)
            d.v[i] = y.v[i] + beta *d.v[i];

        deltaOld = deltaNew;

        it++;
        tIter = timestamp() - tIter;
    }while (it < 25);

    *norma = diff;
    *time = tIter/it;

    free(Yv);
    free(v1);
    free(v2);

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

void calcResidue(struct LinearSis *SL, double *x, double *r, double *time)
{
    if (time)
        *time = timestamp();

    uint n = SL->n;
    double sum = 0.0;

    for (uint i = 0; i < n; i++) {
        sum = 0.0;
        for (uint j = 0; j < n; j++)
            sum += SL->A->v[n*i + j] * x[j];
    
        r[i] = SL->b->v[i] - sum;
    }
    if (time)
        *time = timestamp() - *time;
}

void printVetor(double* vet, int n){
    for(int i = 0; i < n; i++){
        printf("%.16g   ", vet[i]);
    }

    printf("\n");
}
void printSis(struct LinearSis *SL){
    uint n = SL->n;

    for (int i = 0; i < n; i++)
    {
        printf("[  ");
        for(int j = 0; j < n; j++){
            if (SL->A->v[i*n+j] == 0)
                printf("            ");
            else
                printf("%.8f  ", SL->A->v[i*n+j]);
        }
        printf("]   [ %.8f ]\n", SL->b->v[i]);
    }
}

void print7Diag(struct diag7 *SL){

    int n = SL->n;

    printf("  [  %10.6f", SL->D4[0]);
    printf("  %10.6f", SL->D3[1]);
    printf("  %10.6f", SL->D2[2]);
    printf("  %10.6f", SL->D1[3]);
    for(int i = 4; i < n; i++){
        printf("            ");
    }
    printf("  ] [ %10.6f ]\n", SL->B[0]);


    printf("  [  %10.6f", SL->D5[0]);
    printf("  %10.6f", SL->D4[1]);
    printf("  %10.6f", SL->D3[2]);
    printf("  %10.6f", SL->D2[3]);
    printf("  %10.6f", SL->D1[4]);
    for(int i = 5; i < n; i++){
        printf("            ");
    }
    printf("  ] [ %10.6f ]\n", SL->B[1]);


    printf("  [  %10.6f", SL->D6[0]);
    printf("  %10.6f", SL->D5[1]);
    printf("  %10.6f", SL->D4[2]);
    printf("  %10.6f", SL->D3[3]);
    printf("  %10.6f", SL->D2[4]);
    printf("  %10.6f", SL->D1[5]);
    for(int i = 6; i < n; i++){
        printf("            ");
    }
    printf("  ] [ %10.6f ]\n", SL->B[2]);

    for(int i = 3; i < n-3; i++){
        printf("  [");
        for(int j = 3; j < i; j++){
            printf("            ");
        }
        printf("  %10.6f", SL->D7[i-3]);
        printf("  %10.6f", SL->D6[i-2]);
        printf("  %10.6f", SL->D5[i-1]);
        printf("  %10.6f", SL->D4[i]);
        printf("  %10.6f", SL->D3[i+1]);
        printf("  %10.6f", SL->D2[i+2]);
        printf("  %10.6f", SL->D1[i+3]);
        for(int j = i+4; j < n; j++){
            printf("            ");
        }
        printf("  ] [ %10.6f ]\n", SL->B[i]);
    }

    printf("  [");
    for(int i = 0; i < n-6; i++){
        printf("            ");
    }
    printf("  %10.6f", SL->D7[n-6]);
    printf("  %10.6f", SL->D6[n-5]);
    printf("  %10.6f", SL->D5[n-4]);
    printf("  %10.6f", SL->D4[n-3]);
    printf("  %10.6f", SL->D3[n-2]);
    printf("  %10.6f", SL->D2[n-1]);
    printf("  ] [ %10.6f ]\n", SL->B[n-3]);

    printf("  [");
    for(int i = 0; i < n-5; i++){
        printf("            ");
    }
    printf("  %10.6f", SL->D7[n-5]);
    printf("  %10.6f", SL->D6[n-4]);
    printf("  %10.6f", SL->D5[n-3]);
    printf("  %10.6f", SL->D4[n-2]);
    printf("  %10.6f", SL->D3[n-1]);
    printf("  ] [ %10.6f ]\n", SL->B[n-2]);

    printf("  [");
    for(int i = 0; i < n-4; i++){
        printf("            ");
    }
    printf("  %10.6f", SL->D7[n-4]);
    printf("  %10.6f", SL->D6[n-3]);
    printf("  %10.6f", SL->D5[n-2]);
    printf("  %10.6f", SL->D4[n-1]);
    printf("  ] [ %10.6f ]\n", SL->B[n-1]);

}
/*ESSA FUNCAO MAIS GERAL FUNCIONA PARA MATRIZES E VETORES, MAS DEPENDE DE COMO O VETOR ESTA ORGANIZADO */
void multMatrix(struct Matrix *A, struct Matrix *B, struct Matrix *C) {
    if(A->column != B->row)
        return; 

    double sum = 0.0;
    uint aSize = A->column;
    uint bSize = B->column;

    for (uint i = 0; i < A->row; i++) {
        for (uint j = 0; j < B->column ; j++) {
            sum = 0.0;
            for (uint k = 0; k < A->column; k++)
                sum += A->v[i*aSize + k] * B->v[k*bSize + j];
            C->v[i*bSize + j] = sum;
        }
    }
}
