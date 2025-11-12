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

    for(int i = 0; i < n; i++){
        SL->b->v[i] = genRandomB(k);

        for(int j = 0; j < n; j++){
            if ((j > i + k/2 || j < i - k/2)) {
                SL->A->v[n*i+j] = 0.0;
                continue;
            }
            SL->A->v[n*i+j] = genRandomA(i,j,k);
        }
    }
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

    genTranspose(&AT, A);
    
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
                fprintf(stderr, "ERRO: A[%d][%d] = %.6f -> A não é definida positiva\n", i, i, A->v[i*n + i]);
                return -2; 
            }
            else
                M->v[i] = 1.0;
        }

    *time = timestamp();
    if (w == 0){
        for(int i = 0; i < M->row; i++){
            if(A->v[i*n + i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[%d][%d] = %.6f -> A não é definida positiva\n", i, i, A->v[i*n + i]);
                return -2; 
            }
            else
                M->v[i] = 1.0/A->v[i*n + i];
        }
    }

    *time = timestamp() - *time;
    return 0;
}

void genTranspose(struct Matrix *A, struct Matrix *T)
{
    uint n = A->column;
    uint m = A->row;

    for (uint i = 0; i < n; i++)
        for(uint j = 0; j < m; j++)
            A->v[j*n+i] = T->v[i*n+j];
}

int conjGradientPre(struct LinearSis *SL, double *x, double *r,double *norma, struct Matrix *M, uint maxit, double eps, double *time){

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

    //vetor usado para comparação de norma
    double *prevx = calloc(n,sizeof(double));

    if(!v1 || !Yv || !v2 || !prevx){
        free(v1);
        free(Yv);
        free(v2);
        free(prevx);
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
            free(prevx);
            free(Yv);
            free(v1);
            free(v2);
            fprintf(stderr,"cAD Divisão por zero\n");
            return -1;
        }
        alpha = deltaOld / cAd; //Calculando ak
       
        deltaNew = 0.0;
        for (uint i = 0; i < n; i++) {

            prevx[i] = x[i];
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
        diff = calcNormaMax(x, prevx, n); 
        tIter = timestamp() - tIter;
    }while (it < maxit && diff >= eps);

    fprintf(stderr,"diff: %.8g\n",diff);
    fprintf(stderr,"it:%d\n",it);

    *norma = diff;
    *time = tIter/it;

    free(prevx);
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
