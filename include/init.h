#ifndef INIT_H_INCLUDED
#define INIT_H_INCLUDED

/*!
 * \file init.h
 */

#include "./struct.h"


void Table_Add_GF(table_t *table, int GF, int logGF);

void Table_Mul_GF(int **MULGF, int GF);

void Table_Div_GF(int **DIVGF, int GF);

void Table_dec_GF(table_t *table, int GF, int logGF);

void Table_Mul_DEC(table_t *table, int GF);

void Table_Div_DEC(table_t *table, int GF);


void LoadCode (char *FileMatrix, code_t *code);

void FreeCode(code_t *code);

void AllocateDecoder (code_t *code, decoder_t *decoder);

void FreeDecoder (decoder_t *decoder);

void LoadTables (table_t *table, int GF, int logGF);

void FreeTable(table_t *table);


#endif




