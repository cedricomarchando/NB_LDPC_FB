/*!
 * \file init.c
 * \brief initialization of Non-binary LDPC decoder
 * \author C.Marchand, A. Al-Ghouwahel, Oussama Abassi, L. Conde-Canencia, A. abdmoulah, E. Boutillon
 * \copyright BSD copyright
 * \date 03/03/2015


 */


// initialization
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"


#define STR_MAXSIZE 350    //la taille maximale permise pour les chaines de caractères. Utilisée dans la fonction calloc().

#define KN_matrix

/*!
 * \fn Table_Add_GF
 * \brief Compute the addition table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Add_GF(table_t *table, int GF, int logGF)
{
    int i,j,k;
    int temp[8];

    for(j=0; j<GF; j++)
    {
        for(k=0; k<GF; k++)
        {
            for(i=0; i<logGF; i++)
            {
                temp[i] = (table->BINGF[j][i])^(table->BINGF[k][i]);
            }
            table->ADDGF[j][k] = Bin2GF(temp,GF,logGF,table->BINGF);
        }
    }
}

/*!
 * \fn Table_Mul_GF
 * \brief compute the multiplication table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Mul_GF(int **MULGF, int GF)
{
    int i,j,temp;
    for(i=0; i<GF; i++)
    {
        for(j=0; j<GF; j++)
        {
            if (i==0 || j==0)
                MULGF[i][j] = 0;
            else if (i==1)
                MULGF[i][j] = j;
            else if (j==1)
                MULGF[i][j] = i;
            else
            {
                temp=i+j-2;
                if(temp<GF-1)
                    MULGF[i][j] = temp+1;
                else
                    MULGF[i][j]=(temp%(GF-1))+1;
            }
        }
    }
}

/*!
 * \fn Table_Div_GF
 * \brief compute the division table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Div_GF(int **DIVGF, int GF)
{
    int i,j,nb;
    nb=GF-1;
    for(i=0; i<GF; i++)
    {
        for(j=0; j<GF; j++)
        {
            if(j==0)
            {
                DIVGF[i][j]=-1;
            }
            else if(i==0)
            {
                DIVGF[i][j]=0;
            }
            else if(j==1)
            {
                DIVGF[i][j]=i;
            }
            else
            {
                DIVGF[i][j]=nb--;
            };
            if(nb<1)
            {
                nb=GF-1;
            }
        }
    }
}



/**
 * \fn LoadCode
 * \brief Open a parity-check matrix file. Allocate and initialize the code structures
 * Parameters    :
 * Inputs        :
 * 	- FileMatrix  : Name of the parity-check matrix file
 * 	- code_t code : Structure that describes the code
 * Outputs       :
 */
void LoadCode (char *FileMatrix, code_t *code)
{
    int M,N,GF,logGF;
    int n,m,k;
    char *FileName = malloc(STR_MAXSIZE);
    double temp;

    FILE *f;
    /*
     * Load the files corresponding to code (graph, size, GF)
     */
    strcpy(FileName,FileMatrix);
    f=fopen(FileName,"r");

    fscanf(f,"%d",&N);
    fscanf(f,"%d",&M);
    fscanf(f,"%d",&GF);
    temp = log((double)(GF));
    temp = temp/log((double)2.0);
    temp = rint(temp);
    logGF = (int)temp;
    code->N = N;
    code->M = M;
    code->K = N-M;
    code->rate=(float)(N-M)/N;
    code->GF = GF;
    code->logGF = logGF;


    code->columnDegree = calloc(N,sizeof(int));
    //if ( code->columnDegree == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (n=0; n<N; n++) fscanf(f,"%d",&code->columnDegree[n]);

    code->rowDegree = calloc(M,sizeof(int));
    //if ( code->rowDegree == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (m=0; m<M; m++) fscanf(f,"%d",&code->rowDegree[m]);

    code->mat = calloc(M,sizeof(int *));
    //if ( code->mat == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (m=0; m<M; m++)
    {
        code->mat[m] = calloc(code->rowDegree[m],sizeof(int));
        //if ( code->mat[m] == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    }
    code->matValue = calloc(M,sizeof(int *));
    for (m=0; m<M; m++)
    {
        code->matValue[m] = calloc(code->rowDegree[m],sizeof(int));
        //if ( code->matValue [m] == NULL ) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    }


#ifndef KN_matrix
printf(" \n UBS alist format is used! \n");
    for (m=0; m<M; m++)
        for (k=0; k<code->rowDegree[m]; k++)
            fscanf(f,"%d",&code->mat[m][k]);



    for (m=0; m<M; m++)
        for (k=0; k<code->rowDegree[m]; k++)
            fscanf(f,"%d",&code->matValue[m][k]);

#endif



#ifdef KN_matrix
printf(" \n Normal alist format is used! \n");
    int temp_int;
    for (m=0; m<M; m++)
    {
        for (k=0; k<code->rowDegree[m]; k++)
        {
            fscanf(f,"%d",&temp_int);
            code->mat[m][k]=temp_int-1;
            fscanf(f,"%d",&temp_int);
            code->matValue[m][k]=temp_int+1;
        }

    }


#endif







    fclose(f);

    code->nbBranch=0;
    for (m=0; m<M; m++) code->nbBranch += code->rowDegree[m];

//printf(" \n ");
//for (m=0; m<M; m++)
//{
//        for (k=0; k<code->rowDegree[m]; k++)
//        {
//            printf(" %d ",code->mat[m][k]);
//        }
//        printf(" \n ");
//
//}
//getchar();
//
//printf(" \n ");
//for (m=0; m<M; m++)
//{
//        for (k=0; k<code->rowDegree[m]; k++)
//        {
//            printf(" %d ",code->matValue[m][k]);
//        }
//        printf(" \n ");
//
//}
//getchar();




    printf("LDPC code parameters: \n");
    printf(" \t N \t:%d \n \t K \t:%d \n \t M\t:%d \n \t CR\t:%g \n \t GF \t:%d \n \t logGF \t:%d\n",N,N-M,M,code->rate,GF,logGF);
    fflush(stdout);
    free(FileName);

}









/*!
 * \fn void FreeCode(code_t *code)
 * \brief Free pointers in a code_t structure
 * Inputs        :
 * 	- code_t code : Structure that describes the code
 * Outputs       :
 */
void FreeCode(code_t *code)
{
    int m;
    for (m=0; m<code->M; m++)
    {
        free(code->mat[m]);
        free(code->matValue[m]);
    }
    free(code->matUT[0]);
    free(code->matUT);
    free(code->mat );
    free(code->matValue );
    free(code->rowDegree );
    free(code->columnDegree );
    free(code->Perm );
}

/*!
 * \fn AllocateDecoder
 * \brief Memory allocations for the decoder
 */
void AllocateDecoder (code_t *code, decoder_t *decoder)
{
    const int N = code->N;
    int nbRow, nbCol, k, nbRow_arr;

    decoder->nbBranch 	= code->nbBranch;
    decoder->N 		= code->N;

    nbRow = code->nbBranch;
    nbCol = decoder->nbMax;
    nbRow_arr=code->rowDegree[0];


    /* CtoV [nbBranch][nbMax] */
    decoder->CtoV =calloc((size_t)nbRow,sizeof(softdata_t *));
    //if (decoder->CtoV  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->CtoV [0] = calloc((size_t)nbRow*code->GF,sizeof(softdata_t));
    //if (decoder->CtoV [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) decoder->CtoV[k] = decoder->CtoV[0] + k*code->GF;

    /* M_CtoV_LLR [nbBranch][nbMax] */
    decoder->M_CtoV_LLR =calloc((size_t)nbRow_arr,sizeof(softdata_t *));
    //if (decoder->M_CtoV_LLR  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->M_CtoV_LLR [0] = calloc((size_t)nbRow_arr*code->GF,sizeof(softdata_t));
    //if (decoder->M_CtoV_LLR [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow_arr; k++) decoder->M_CtoV_LLR[k] = decoder->M_CtoV_LLR[0] + k*code->GF;

    /* M_VtoC_LLR [nbBranch][nbMax] */
    decoder->M_VtoC_LLR =calloc((size_t)nbRow_arr,sizeof(softdata_t *));
    //if (decoder->M_VtoC_LLR  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->M_VtoC_LLR [0] = calloc((size_t)nbRow_arr*nbCol,sizeof(softdata_t));
    //if (decoder->M_VtoC_LLR [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow_arr; k++) decoder->M_VtoC_LLR[k] = decoder->M_VtoC_LLR[0] + k*nbCol;

    /* M_CtoV_GF [rowDegree[0]][nbMax] */
    decoder->M_CtoV_GF =calloc((size_t)nbRow_arr,sizeof(int *));
    //if (decoder->M_CtoV_GF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->M_CtoV_GF [0] = calloc((size_t)nbRow_arr*code->GF,sizeof(int));
    //if (decoder->M_CtoV_GF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow_arr; k++) decoder->M_CtoV_GF[k] = decoder->M_CtoV_GF[0] + k*code->GF;
    /* M_VtoC_GF [rowDegree[0]][nbMax] */
    decoder->M_VtoC_GF =calloc((size_t)nbRow_arr,sizeof(int *));
    //if (decoder->M_VtoC_GF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->M_VtoC_GF [0] = calloc((size_t)nbRow_arr*nbCol,sizeof(int));
    //if (decoder->M_VtoC_GF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow_arr; k++) decoder->M_VtoC_GF[k] = decoder->M_VtoC_GF[0] + k*nbCol;



    nbRow = N;
    nbCol = code->GF;
    /* APP [N][GF] */
    decoder->APP =calloc((size_t)nbRow,sizeof(softdata_t *));
    //if (decoder->APP  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->APP [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
    //if (decoder->APP [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) decoder->APP[k] = decoder->APP[0] + k*nbCol;

    nbRow = N;
    nbCol = code->GF;
    /* intrinsic_LLR [N][GF] */   /* VN modified */
    decoder->intrinsic_LLR =calloc((size_t)nbRow,sizeof(softdata_t *));
    //if (decoder->intrinsic_LLR  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->intrinsic_LLR [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
    //if (decoder->intrinsic_LLR [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) decoder->intrinsic_LLR[k] = decoder->intrinsic_LLR[0] + k*nbCol;

    /* intrinsic_GF [N][GF] */  /* VN modified */
    decoder->intrinsic_GF 	= calloc((size_t)nbRow,sizeof(int *));
    //if (decoder->intrinsic_GF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    decoder->intrinsic_GF [0] 	= calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (decoder->intrinsic_GF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) decoder->intrinsic_GF[k] = decoder->intrinsic_GF[0] + k*nbCol;

}



/**
 * \fn void FreeDecoder (decoder_t *decoder)
 * \brief Free pointers in a decoder_t structure
 * Inputs        :
 * 	- decoder_t decoder : Structure that describes the decoder (edges in the decoding graph)
 * Outputs       :
 */
void FreeDecoder (decoder_t *decoder)
{
    free(decoder->CtoV[0] );
    free(decoder->CtoV);
    free(decoder->M_CtoV_LLR[0] );
    free(decoder->M_CtoV_LLR);
    free(decoder->M_VtoC_LLR[0] );
    free(decoder->M_VtoC_LLR);
    free(decoder->M_CtoV_GF[0] );
    free(decoder->M_CtoV_GF);
    free(decoder->M_VtoC_GF[0] );
    free(decoder->M_VtoC_GF);

    free(decoder->APP[0] );
    free(decoder->APP);
    free(decoder->intrinsic_LLR[0] );
    free(decoder->intrinsic_LLR);
    free(decoder->intrinsic_GF[0] );
    free(decoder->intrinsic_GF);

}


/**
 * \fn void LoadTables (table_t *table, int GF, int logGF)
 * \brief Memory allocation for the tables and Initialization of the tables.
 * Inputs        :
 * 	- table_t table : Structure that contains pointers to the tables
 * 	- int GF    : order of the field
 * 	- int logGF : logGF = log2(GF)
 * Outputs       :
 */
void LoadTables (table_t *table, int GF, int logGF)
{
    int nbRow, nbCol, g,k,l;

    if(GF!=16 && GF!=64 && GF!=256)
    {
        printf("The binary image of GF(%d) is not available in this version of the program. Please try GF(64) or GF(256)\n",GF);
        exit(EXIT_FAILURE);
    }

    nbRow = GF;
    nbCol = logGF;
    /* BINGF [GF][logGF] */
    table->BINGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->BINGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->BINGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->BINGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->BINGF[k] = table->BINGF[0] + k*nbCol;

    nbRow = GF;
    nbCol = GF;
    /* ADDGF [GF][logGF] */
    table->ADDGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->ADDGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->ADDGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->ADDGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->ADDGF[k] = table->ADDGF[0] + k*nbCol;

    /* MULGF [GF][logGF] */
    table->MULGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->MULGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->MULGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->MULGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->MULGF[k] = table->MULGF[0] + k*nbCol;

    /* DIVGF [GF][logGF] */
    table->DIVGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->DIVGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->DIVGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->DIVGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->DIVGF[k] = table->DIVGF[0] + k*nbCol;


    if(GF==16)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_16[g][l];
        //printf("Loading of the binary image of GF(64): Success\n");
        //fflush(stdout);
    }


    if(GF==64)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_64[g][l];
        //printf("Loading of the binary image of GF(64): Success\n");
        //fflush(stdout);
    }

    if(GF==256)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_256[g][l];
        //printf("Loading of the binary image of GF(256): Success\n");
        //fflush(stdout);
    }


    /*
     * Build the addition, multiplication and division tables (corresponding to GF[q])
     */
    Table_Add_GF (table,GF,logGF);
    Table_Mul_GF (table->MULGF,GF);
    Table_Div_GF (table->DIVGF,GF);

}


/*!
 * \fn FreeTable(table_t *table)
 * \brief Free the memory in a table_t structure
 * Inputs        :
 * 	- table_t table : Structure that contains pointers to the tables
 * Outputs       :
 */
void FreeTable(table_t *table)
{
    free(table->ADDGF[0]);
    free(table->MULGF[0]);
    free(table->DIVGF[0]);
    free(table->BINGF[0]);
    free(table->ADDGF);
    free(table->MULGF);
    free(table->DIVGF);
    free(table->BINGF);
}




