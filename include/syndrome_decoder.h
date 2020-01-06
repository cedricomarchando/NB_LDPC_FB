#ifndef SYNDROME_DECODER_H_INCLUDED
#define SYNDROME_DECODER_H_INCLUDED

/*!
 * \file syndrome_decoder.h
 */

#include "./struct.h"
#include "./init.h"

int syndrome_ems3(int node, decoder_t *decoder, const code_t *code, const table_t *table, int ** config_table, int config_table_size , int dc_max, float offset);
int syndrome_ems2(int node, decoder_t *decoder, const code_t *code, const table_t *table, int ** config_table, int config_table_size , int dc_max, float offset);

int syndrome_ems(int node, decoder_t *decoder, const code_t *code, const table_t *table, int ** config_table, int config_table_size , int dc_max, float offset, int n_cv);





int factorial( int p);
int combin( int n, int r);

int sorting(syndrome_type* syndrome_set, int syndrome_table_size);

float select_sort(float* array_float,int size_array, int order);

float median_median64(float* array_float);
float median_median256(float* array_float);
float median_median_x(float* array_float, int array_float_size, int median_size, int median_select);

int** AllocateArray(int line, int column);



int** build_config_table(int* config_table_size_p,int dc_max,int d_1,int d_2, int d_3);

int sort_config_table(int ** config_table,int config_table_size,int dc_max);

int compute_config_table_size(int dc_max, int d_1, int d_2, int d_3);

int gen_config_table(int ** config_table,int d_c,int d_1,int d_2, int d_3);
int gen_config_table2(int ** config_table,int d_c,int d_1,int d_2, int d_3);
int gen_config_table3(int ** config_table,int d_c,int d_1);


int presorting_mvc(decoder_t *decoder,const code_t *code,int dc_max,int border, int *index_order);

double bayes(double M1 , double M2  );
double f(double x);

int check_deviation(int ** config_table,int dc_max,int index1,int index2,int index3,int *result);

#endif // SYNDROME_DECODER_H_INCLUDED
