#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"


/*ATENCAO: Funcoes inline tem sua implementacao no .h*/

/*
 * Coloca valores x aleatorios no SL.
 * @param SL: Sistema linear a ser modificado
 * */
void genKDiagonal(struct LinearSis *SL, uint k, uint n);

/*
 * Gera um sistema linear simetrico aplicando o metodo CGNE "A^t * Ax = A^t * b"
 * @param SL: Sistema linear a ser corrigido
 * @param ASP: MAtriz A corrigida 
 * @param bsp: Vetor b corrigido
 * @param Variavel para calcular o tempo
 * */
int genSymmetricPositive(struct LinearSis *SL, struct diagMat *ASP, double *bsp, double *time);

/*
 * Gera os pre condicionamentos dependendo do valor de w 
 * @param A: Matriz A que gerara os preCond
 * @param w: Valor que define qual tipo de pre condicionamento sera usado
 * @param n: Ordem da matriz A 
 * @param k: Ordem da diagonal da matriz A 
 * @param M: Matriz preCond
 * @param time: Variavel para calcular o tempo
 * */
int genPreCond(struct diagMat *A, double w, uint n, double* M, double *time);

/*
 * Gera a transposta da Matriz
 * @param A: Matriz A que gerara a transposta
 * @param AT: A matriz transposta
 * */
void genTranspose(struct diagMat *A, struct diagMat *AT, uint k);

/*
 * Algoritmo que resolve Ax = b com o uso de pre condicionamento
 * @param SL: Sistema Linear a ser resolvido 
 * @param x: Vetor solucao 
 * @param r: Vetor do residuo. R = b - Ax
 * @param r: Vetor de diferenças em x entre iterações
 * @param M: Matriz preCond
 * @param maxit: Valor maximo para convergencia
 * @param eps: Valor de parada
 * @param time: Variavel para calcular o tempo
 * */
int conjGradientPre(struct diagMat *A, double *B, double *x, double *r, double *M, double *time);

/*
 * Algoritmo que retorna a maior diferença entre os elementos de dois vetores de tamanho n
 * @param x: Vetor para calcular a norma
 * @param y: Vetor para calcular a norma
 * @param n: Tamanho do Vetor
 * */
double calcNormaMax(double *x, double *y, uint n);

/*
 * Algoritmo que retorna a raiz da soma dos quadrados de um vetor
 * @param x: Vetor para calcular a norma
 * @param n: Tamanho do Vetor
 * */
double calcNormaEuclidiana(double *x, uint n);

/*
 * Algoritmo que calcula r = b - Ax
 * @param SL: Sistema Linear que o residuo ira calcular
 * @param x: Vetor solucao 
 * @param r: Vetor do residuo
 * @param time: Variavel para calcular o tempo
 * */
void calcResidue(struct diagMat *A,double *B, uint k, double *x, double *r, double *time);

/*
 * Printa o vetor, debug
 * @param vet: Vetor que sera printado
 * @param n: Quantidade de itens no vetor 
 * */
void printVetor(double* vet, uint n);

/*
 * Printa o SL
 * @param SL: Sistema linear a ser imprimido
 * */
void printSis(struct LinearSis *SL);

/*
 * A * B = C, funciona para vetores
 * @param A: Matriz A que sera multiplicada
 * @param B: Matriz B que sera multiplicada
 * @param C: Matriz resultado
 * 
void multMatrix(struct Matrix *A, struct Matrix *B, struct Matrix *C);
*/
void print7Diag(struct diagMat *SL, uint k);

void square7diag(struct LinearSis *SL, struct diagMat *C, double* newB);

void inline static initDiag(struct diagMat* A, uint k, uint n);

void static multMatVet(struct diagMat *A, double *B, double *C, uint k);

#endif // __SISLIN_H__;

