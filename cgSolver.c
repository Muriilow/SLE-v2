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
    uint n;
    uint k;
    uint maxit;
    double eps;
    float w;

    scanf("%d", &n);
    scanf("%d", &k);
    scanf("%f", &w);

    if(w != 0 && w !=-1){
        fprintf(stderr,"Opção w inválida\n");
        return -2;
    }

    scanf("%d", &maxit);
    scanf("%lf", &eps);

    struct LinearSis* newdiag = aligned_alloc(32,sizeof(struct LinearSis));
    struct diagMat* sympos = malloc(sizeof(struct diagMat));
    double* spb = calloc(n,sizeof(double));
    double* M = calloc(n,sizeof(double));

    if(!newdiag || !sympos || !spb || !M){
        fprintf(stderr, "Error alocating memory\n");
        return -1;
    }

    newdiag->n = n;
    newdiag->k = k;

    genKDiagonal(newdiag, k, n);

    genSymmetricPositive(newdiag, sympos, spb, NULL);

    genPreCond(sympos, w, n, M, NULL);
    printVetor(M,n);

    double* x = calloc(n,sizeof(double));
    double* r = calloc(n,sizeof(double));
    double timeGrad;
    double timeRes;

    LIKWID_MARKER_START("EXEC_1");
    status += conjGradientPre(sympos, spb, x, r, M, &timeGrad);
    LIKWID_MARKER_STOP("EXEC_1");

    if (status < 0)
    {
       free(newdiag);
       free(sympos);
       free(spb);
       free(M);
    }

    LIKWID_MARKER_START("RES_1");
    calcResidue(newdiag->A, newdiag->b, k, x, r, &timeRes);
    LIKWID_MARKER_STOP("RES_1");

    double NormaR = calcNormaEuclidiana(r, n);

    printf("%d\n",n);

    printf("X: ");
    printVetor(x,n);

    printf("R: ");
    printVetor(r,n);
    
    printf("normaR: %.8g\n", NormaR);

    FILE* timeExec = fopen("timeExec2.txt", "w");
    FILE* timeResi = fopen("timeRes2.txt", "w");

    if (!timeExec || !timeRes) {
        printf("Erro ao abrir o arquivo!\n");
        return 1;
    }
    
    // Escrever no arquivo (similar ao printf)
    fprintf(timeExec, "%.8g\n", timeGrad);
    fprintf(timeResi, "%.8g\n", timeRes);

    LIKWID_MARKER_CLOSE;
    
    fclose(timeExec);
    fclose(timeResi);

    free(x);
    free(r);
    free(sympos);
    free(spb);
    free(M);
    free(newdiag);

    return 0;
}
