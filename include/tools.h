#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

/*!
*
* \file tools.h
*
*/



#include "./struct.h"
#include <stdlib.h>
#include <stdio.h>



float My_drand48(int *initialise);

int Bin2GF(int *U,int GF,int logGF,int **BINGF);

void RandomBinaryGenerator (int N, int M,int GF,int logGF,int **KBIN,int *KSYMB,int **BINGF, int *init_rand);

void GaussianElimination (code_t *code, table_t *table);

int Encoding (code_t *code, table_t *table, int *CodeWord, int **NBIN, int *NSYMB);

int Syndrom (code_t *code, int *decide, table_t *tableGF);

void Decision( int *decision,float **APP,int N,int GF);

#endif
