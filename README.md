Modificações feitas - Murilo Bob 

1. Modificar todos os int para unsigned ints quando possivel, seguindo (TODO: LINK)

# Solução de sistemas lineares esparsos usando Pré-condicionadores

Programa que utilizado o método dos Gradientes Conjugados para achar uma solução **x̄** dentro de uma certa tolerância **ε** ou limite de iterações **maxit**.

### Entradas: 
Lê da entrada padrão cinco valores:

    n: (n>10)  a dimensão do Sistema Linear.
    k: (k>1 e k ímpar)  o número de diagonais da matriz A.
    ω:  o pré-condicionador a ser utilizado:

        ω=-1: sem pré-condicionador
        ω=0.0 : pré-condicionador de Jacobi

    maxit:  o número máximo de iterações a serem executadas.
    ε: o erro aproximado absoluto máximo da norma do vetor X.


## Autoria
Guilherme Vitoriano Santana de OLiveira - GRR20245396  
Murilo de Paula Bob - GRR20242184

## Estruturas de Dados:

Struct que guarda uma matriz e informações importantes par trabalhar sobre ela:

```
struct Matrix {
    double *v;
    uint row; 
    uint column;
    uint k;
};
```
   
> v: Matriz diagonal  
> row: Quantidade de linhas   
> column: Quantidade de coluna  
> k: Quatidade de diagonais da matriz, 0 para matriz não diagonal

Struct que define um sistema Linear:

```
struct LinearSis {
    struct Matrix *A;
    struct Matrix *b;
    uint n; 
    uint k; 
};

```
> A: Matriz k-diagonal  
> b: Vetor de valores independentes do sistema  
> n: Ordem da matriz  
> k: Valor que define quantas diagonais não zero


## Módulos

### utils.h

**MACROS**

- #define ABS(num)  ((num) < 0.0 ? -(num) : (num))
- #define numDigits(n)  6  // ( (int) log10(n) + 1 )
- #define isPot2(n) (n && !(n & (n - 1)))

**FUNÇÕES**

- double timestamp(void);
- char* markerName(char* baseName, int n);

### sislin.h

**FUNÇÕES**

- void genKDiagonal(struct LinearSis *SL);
- void genSymmetricPositive(struct LinearSis *SL, struct Matrix *ASP, struct Matrix *bsp, double *time);
- int genPreCond(struct Matrix *A, double w, int n, int k, struct Matrix *M, double *time);
- void genTranspose(struct Matrix *A, struct Matrix *AT);
- int conjGradientPre(struct LinearSis *SL, double *x, double *r, double *norma, struct Matrix *M, uint maxit, double eps, double* time)
- double calcNormaMax(double *x, double *y, int n);
- double calcNormaEuclidiana(double *x, int n);
- void calcResidue(struct LinearSis *SL, double *x, double *r, double* time);
- void printVetor(double* vet, int n);
- void printSis(struct LinearSis *SL);
- void multMatrix(struct Matrix *A, struct Matrix *B, struct Matrix *C);




