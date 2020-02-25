#ifndef CSK_H_INCLUDED
#define CSK_H_INCLUDED


/*
Name: stuct csk_t
Brief : Structure containing parameters of the csk modulation
*/
typedef struct
{
    int PNsize; // size of the PN sequence
    int *PN;    // the PN sequence
    int **CSK_arr; //CSK_arr[GF][PNsize] : table of CSK symbols
} csk_t;


void PNGenerator( csk_t *csk);

void build_natural_csk_mapping(int GF, csk_t *csk);

void build_punctured_csk_mapping(int GF,int logGF, csk_t *csk, int **BINGF);

void build_CSK_map(code_t *code, csk_t *csk);

void GF_to_CSK(int* CodeWord, int N, csk_t *csk, int** NCSK);

void allocate_csk(csk_t *csk, int GF);

void free_csk(csk_t *csk);

void ModelChannel_AWGN_BPSK_CSK (csk_t *csk, code_t *code, decoder_t *decoder, table_t *table, int *CodeWord, float EbN,int *init_rand, float quantif_range_int,float quantif_range_float);




#endif // CSK_H_INCLUDED
