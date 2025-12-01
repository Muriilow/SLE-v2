#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include <likwid.h>
int main(){
    LIKWID_MARKER_INIT;
    srandom(20252);
    int status = 0;
    int n;
    int k;
    uint maxit;
    float w;
    double eps;

    scanf("%d", &n);
    scanf("%d", &k);
    scanf("%f", &w);
    if(w != 0 && w !=-1){
        fprintf(stderr,"Opção w inválida\n");
        return -2;
    }
    scanf("%d", &maxit);
    scanf("%lf", &eps);

    /*Criando o sistema Linear, A: Matriz, b: Vetor Indep*/
    double *av = malloc(n*n*sizeof(double));
    struct Matrix A = {av, n, n, k};

    double *bv = malloc(n*sizeof(double));
    struct Matrix b = {bv, n, 1, 0};

    struct LinearSis SL = {&A, &b, n, k};

    genKDiagonal(&SL);


    /*Gerando simetrica positiva, criando Matriz A2, b2 para guardar os novos valores*/
    double *av1 = malloc(n*n*sizeof(double));
    struct Matrix A2 = {av1, n, n, k};

    double *bv1 = malloc(n*sizeof(double));
    struct Matrix b2 = {bv1, n, 1, 0};
    double timePC; //Arrumar esse time e os outros do trabalho

    if(!av || !bv || !av1 || !bv1){
        free(av);
        free(bv);
        free(av1);
        free(bv1);
        fprintf(stderr, "Falha na alocação de memória\n");
        return -1;
    }
    status = genSymmetricPositive(&SL, &A2, &b2, &timePC); 

    struct LinearSis SLNew = {&A2, &b2, n, k};

    free(av);
    free(bv); //Liberando o vetor A e b pois nao preciso mais deles
    
    if(status < 0){
        free(av1);
        free(bv1);
        return -1;
    }

    SL.A = &A2;
    SL.b = &b2; //Agora os valores a matriz e vetor independente sao simetricos

    double* X = (double*) calloc(n, sizeof(double));
    double* r = (double*) malloc(n * sizeof(double));
    double* norma = (double*) malloc(sizeof(double));
    double timeM; 
    double timeGrad;
    double timeRes;
    double *Mv = calloc(n, sizeof(double));
    struct Matrix M = {Mv, n, 1, SL.k};

    if(!Mv || !X || !r || !norma){
        free(av1);
        free(bv1);
        free(Mv);
        free(X);
        free(r);
        free(norma);
        fprintf(stderr, "Falha na alocação de memória\n");
        return -1;
    }
    status = genPreCond(SLNew.A, w, SLNew.n, SLNew.k, &M, &timeM);

    LIKWID_MARKER_START("EXEC_1");
    status += conjGradientPre(&SLNew, X, r, norma,&M, maxit, eps, &timeGrad);
    LIKWID_MARKER_STOP("EXEC_1");

    free(Mv);

    if (status < 0){
        free(norma);
        free(X);
        free(av1);
        free(bv1);
        free(r);
        return -1;
    }

    LIKWID_MARKER_START("RES_1");
    calcResidue(&SL, X, r, &timeRes);
    LIKWID_MARKER_STOP("RES_1");
    
    printf("%d\n",n);
    printVetor(X,n);
    double normaR = calcNormaEuclidiana(r, n);
    printf("%.8g\n", *norma);
    printf("normaR %.8g\n", normaR);
    printf("%.8g\n", timePC + timeM);
    printf("%.8g\n", timeGrad);
    printf("%.8g\n", timeRes);
    
    double timeAv = (timeGrad + timeRes) / 2;

    FILE* time;

    time = fopen("time.txt", "w");
    if (time == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        return 1;
    }
    
    // Escrever no arquivo (similar ao printf)
    fprintf(time, "%.8g\n", timeAv);
    free(norma);
    free(X);
    free(av1);
    free(bv1);
    free(r);
    
    LIKWID_MARKER_CLOSE;
    return 0;
}
