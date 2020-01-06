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
void CheckPassLogEMS_bayes (int node,decoder_t *decoder, code_t *code, table_t *table, float offset);
void CheckPassLogEMS_presorting (int node,decoder_t *decoder, code_t *code, table_t *table, int NbOper, float offset,int * stat_bubble,int stat_on);
void CheckPassLogEMS_tree_presorting (int node,decoder_t *decoder, code_t *code, table_t *table,int NbOper, float offset, int* stat_bubble);
void CheckPassLogEMS_tree(int node,decoder_t *decoder, code_t *code, table_t *table, float offset, int* stat_bubble);
void CheckPassLogEMS_presorting_synd (int node,decoder_t *decoder, code_t *code, table_t *table,int NbOper, float offset, int* stat_bubble,int stat_on);



int ElementaryStep(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int GF,int nbMax,int nbOper);
int ElementaryStep_bayes(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int nbMax);

int ElementaryStep_nm(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int GF,int nbMax,int nbOper, int check_number,int * stat_bubble,int stat_on, int nm_out);
int ElementaryStep_nm_rect(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int GF,int nbOper,int nmU,int nmV,int nmS);


int ElementaryStep_tree(float *Input_LLR1, float *Input_LLR2, int *Input_GF1, int *Input_GF2,int *Input_deviation1,int *Input_deviation2,float *Output_LLR, int *Output_GF, int*Output_deviation,
                        int **ADDGF,int GF,int nm_in,int nm_out,int nb_bubble);
int ElementaryStep_synd_rect(float *Input_LLR1, float *Input_LLR2, int *Input_GF1, int *Input_GF2,int *Input_deviation1,int *Input_deviation2,float *Output_LLR,
                         int *Output_GF, int*Output_deviation,int **ADDGF,int GF,int nm_in1, int nm_in2,int nm_out,int check_number,int * stat_bubble,int stat_on);

int ElementaryStep_synd_rect_bayes(float *Input_LLR1, float *Input_LLR2, int *Input_GF1, int *Input_GF2,int *Input_deviation1,int *Input_deviation2,float *Output_LLR,
                         int *Output_GF, int*Output_deviation,int **ADDGF,int GF,int nm_in1, int nm_in2,int nm_out,int check_number,int * stat_bubble,int stat_on);


int ElementaryStep_synd_triangle(float *Input_LLR1, float *Input_LLR2, int *Input_GF1, int *Input_GF2,int *Input_deviation1,int *Input_deviation2,float *Output_LLR,
                         int *Output_GF, int*Output_deviation,int **ADDGF,int GF,int nm_in1, int nm_in2,int nm_out, int check_number, int * stat_bubble, int stat_on);


void update_Mcv_bayes( int dc,decoder_t *decoder , code_t *code, int *GF, float *LLR, int offset);


#endif // SYNDROME_DECODER_H_INCLUDED
