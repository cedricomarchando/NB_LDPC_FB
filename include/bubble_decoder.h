#ifndef BUBBLE_DECODER_H_INCLUDED
#define BUBBLE_DECODER_H_INCLUDED

/*!
 * \file bubble_decoder.h
 * \brief header for bubble_decoder functions
 */


#include <math.h>
#include "./struct.h"
#include "./init.h"


int maximum(float *Y, int nl);
int minimum(float *Y, int nl);
void CheckPassLogEMS (int node,decoder_t *decoder, code_t *code, table_t *table, int NbOper, float offset);
void CheckPassLogEMS_dc3 (int node,decoder_t *decoder, code_t *code, table_t *table, int NbOper, float offset);



int ElementaryStep(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int GF,int nbMax,int nbOper);
int ElementaryStep_nm(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int GF,int nmU,int nmV,int nmS,int nbOper);



#endif // SYNDROME_DECODER_H_INCLUDED
