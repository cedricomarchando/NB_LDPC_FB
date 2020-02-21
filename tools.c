/*!
 * \file tools.c
 * \brief tools for NB decoder
 * \author C.Marchand, A. Al-Ghouwahel, Oussama Abassi, L. Conde-Canencia, A. abdmoulah, E. Boutillon
 * \copyright BSD copyright
 * \date 03/03/2015


 */

#include <string.h>
#include <math.h>

#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"








//// random generator: uniform distribution on [0,1]
//float My_drand48 (void)
//{
//    return((float)(rand())/(float)(RAND_MAX+1.0));
//}


// random generator: uniform distribution on [0,1]

/*!
* \fn My_drand48
* \brief improved randum number generator working for linux and windows
* return an uniformely distributed random number, for windows to replace drand48 of linux
* \author C. Marchand
*******************************************************************************/
float My_drand48(int *initialise)
{

#ifdef _WIN32

    static int s1, s2;
    int k, Z;

    if ( *initialise == -1 )
    {
        s1 = (int)(( rand() % 2147483562 ) + 1);
        s2 = (int)(( rand() % 2147483398 ) + 1);
        *initialise = 0;
    }

    k = s1/53668;
    s1 = 40014*(s1 - k*53668) - k*12211;
    if (s1 < 0)
        s1 += 2147483563;

    k = s2/52774;
    s2 = 40692*(s2 - k*52774) - k*3791;
    if (s2 < 0)
        s2 += 2147483399;

    Z = s1 - s2;
    if (Z < 1)
        Z += 2147483562;

    return(Z/2147483563.0);

#elif __linux__

    return(drand48());

#endif


}



/*!
 * \fn Bin2GF
 * \brief compute a GF(q) symbol corresponding to a frame of log_2(GF) bits
 * Parameters    :
 * Inputs        :
 * 	- int *U    : array representing logGF bits
 * 	- int logGF : size of the array U. logGF = log2 (GF)
 * 	- int GF    : order of the field
 * 	- int ** BINGF: binary mapping table
 * Outputs       :
 *      - index of the non-binary symbol
 */

int Bin2GF(int *U,int GF,int logGF,int **BINGF)
{
    int k;

    for (k=0; k<GF; k++)
    {
        if (memcmp(U,BINGF[k],sizeof(int)*logGF)==0) break;
    }
    return(k);
}




/**
 * \fn RandomBinaryGenerator
 * \brief Uniform random binary generator (generate the information bits of the code (KBIN))
 * 		   The redundancy symbols are initialized to 0.
 * Inputs
 * 	- N 	: Code length
 * 	- M 	: Number of parity non-binary symbols
 * 	- GF    : Order of the field
 * 	- logGF : logGF = log2(GF)
 * Outputs       :
 *      - NBIN  : Binary representation of the codeword
 *      - KIN   : Binary representation of the information symbols
 *      - NSYMB : Non-binary symbols of the codeword
 *      - KSYMB : Information non-binary symbols
 */
void RandomBinaryGenerator (int N, int M,int GF,int logGF, int **KBIN, int *KSYMB,int **BINGF, int *init_rand)
{
    int k,q;

    /* Random and (bitwise) uniformly distributed information symbols */
    for (k=0; k<N-M; k++)
    {
        for (q=0; q<logGF; q++)
            KBIN[k][q]=floor(My_drand48(init_rand)*1.9999); // avoid the case 2

        KSYMB[k]=Bin2GF(KBIN[k],GF,logGF,BINGF);
     }
}


/**
 * \fn GaussianElimination
 * \brief Perform a Gaussian elimination on the parity-check matrix.
 * 		   The procedure stops when the
 * 		   The redundancy symbols are initialized to 0.
 * Inputs
 * 	- code->mat   : parity-check matrix
 * 	- table       : lookup tables for computations in GF(q)
 * Outputs       :
 *      - code->matUT : Upper triangular matrix used for encoding
 *      - code->Perm : Column permutation
 */
void GaussianElimination (code_t *code, table_t *table)
{
    const int N = code->N;
    const int M = code->M;
    int n,m,k,m1,ind, buf;
    int temp[12];
    int i;

    code->matUT = calloc((size_t)M,sizeof(int *));
    //if (code->matUT == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    code->matUT [0] = calloc((size_t)M*N,sizeof(int));
    //if (code->matUT[0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<M; k++) code->matUT[k] = code->matUT[0] + k*N;

    code->Perm 	= calloc(N,sizeof(int));
    //if (code->Perm == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);

    for (n=0; n<N; n++) code->Perm[n] = n;

    for (m=0; m<M; m++)
    {
        for (k=0; k<code->rowDegree[m]; k++)
        {
            code->matUT[m][code->mat[m][k]] = code->matValue[m][k];
        }
    }
    for (m=0; m<M; m++)
    {
        for (ind=m; ind<N; ind++)
        {
            if (code->matUT[m][ind]!=0) break;
        }
        if (ind==N)
        {
            printf("The matrix is not full rank (%d,%d)\n",m,ind);
            exit(EXIT_FAILURE);
        }
        buf = code->Perm[ind];
        code->Perm[ind] = code->Perm[m];
        code->Perm[m] = buf;
        for (m1=0; m1<M; m1++)
        {
            buf = code->matUT[m1][m];
            code->matUT[m1][m] = code->matUT[m1][ind];
            code->matUT[m1][ind] = buf;
        }

        for (m1=m+1; m1<M; m1++)
        {
            if (code->matUT[m1][m]!=0)
            {
                buf=code->matUT[m1][m];
                for (n=m; n<N; n++)
                {
                    if (code->matUT[m1][n]!=0)
                        code->matUT[m1][n] = table->DIVGF[code->matUT[m1][n]][buf];
                }
                for (n=m; n<N; n++)
                {
                    if (code->matUT[m1][n]!=0)
                        code->matUT[m1][n] = table->MULGF[code->matUT[m1][n]][code->matUT[m][m]];
                }
                for (n=m; n<N; n++)
                {
                                for(i=0; i<code->logGF; i++)
                                {
                                    temp[i] = (table->BINGF[code->matUT[m1][n]][i])^(table->BINGF[code->matUT[m][n]][i]);
                                }
                                code->matUT[m1][n] = Bin2GF(temp,code->GF,code->logGF,table->BINGF);

                    //code->matUT[m1][n] = table->ADDGF[code->matUT[m1][n]][code->matUT[m][n]];
                }
            }
        }
    }
}


/**
 * \fn Encoding
 * \brief Encode the information bits into a codeword.
 * 		   matUT beeing upper triangular, the backsubstitution method is used.
 * 		   The M first symbols in NSYMB are redundancy symbols (before deinterleaving)
 * Inputs
 * 	- KSYMB  ( KSYMB are information symbols)
 * Outputs
 *      - Codeword
 *      - NBIN : binary copy of the codeword
 */
int Encoding(code_t *code, table_t *table, int *CodeWord, int **NBIN, int *KSYMB)
{
    const int N = code->N;
    const int M = code->M;
    const int logGF = code->logGF;
    int k,n,m,q,buf;
    int NSYMB[N];
    int temp[12];
    int i;

    for (k=0 ; k<N-M; k++)
        NSYMB[M+k]=KSYMB[k];

    /* Backsubstitution */
    for (m=M-1; m>=0; m--)
    {
        buf=0;
        for (n=m+1; n<N; n++)
        {
            if (code->matUT[m][n]!=0)
            {
                for(i=0; i<code->logGF; i++)
                {
                    temp[i] = (table->BINGF[buf][i])^(table->BINGF[table->MULGF[code->matUT[m][n]][NSYMB[n]]][i]);
                }
                buf = Bin2GF(temp,code->GF,code->logGF,table->BINGF);
              //  buf = table->ADDGF[buf][table->MULGF[code->matUT[m][n]][NSYMB[n]]];
            }
        }
        /* Systematic codeword (interleaved) */
        NSYMB[m] = table->DIVGF[buf][code->matUT[m][m]];
        }

    /* De-interleaving */
    for (n=0; n<N; n++)
        CodeWord[code->Perm[n]] = NSYMB[n];

    /* Binary copy of the codeword: */
    for (n=0; n<N; n++)
    {
        for (q=0; q<logGF; q++)
            NBIN[n][q] = table->BINGF[CodeWord[n]][q];
                   }

    return(0);
}




/**
 * \fn Syndrom
 * \brief Compute the syndom of a message
 * Inputs
 * 	- code structure code_t
 * 	- table_t tableGF : lookup table
 * 	- message
 * Outputs
 * 	- synd is 0 iff. the decided message is a codeword
 * 	(the value of synd is not meaningful if synd != 0 )
 */
int Syndrom (code_t *code, int *decide, table_t *tableGF)
{
    int k,l;
    int synd;
    int temp[12];
    int i;

    synd = 0;
    for (k=0; k<code->M; k++)
    {
        for (l=0; l<code->rowDegree[k]; l++)
        {

            for(i=0; i<code->logGF; i++)
            {
                temp[i] = (tableGF->BINGF[synd][i])^(tableGF->BINGF[tableGF->MULGF[code->matValue[k][l]][decide[code->mat[k][l]]]][i]);
            }
            synd = Bin2GF(temp,code->GF,code->logGF,tableGF->BINGF);
            //synd = tableGF->ADDGF[synd][tableGF->MULGF[code->matValue[k][l]][decide[code->mat[k][l]]]];
        }

        if (synd != 0)
            break;
    }

    return (synd);
}


/**
 * \fn void Decision( int *decision,float **APP,int N,int GF)
 * \brief Make a hard decision given the APP
 * Inputs
 * 	- APP
 * 	- N
 * 	- GF
 * Outputs
 * 	- decision
 */
void Decision( int *decision,float **APP,int N,int GF)
{
    int n,g,ind;
    float min;

    for (n=0; n<N; n++)
    {
        //max = APP[n][0];
        min=+1e5;
        ind=0;
        for (g=0; g<GF; g++)
            if (APP[n][g]<min)
            {
                min=APP[n][g];
                ind=g;
            }
        decision[n]=ind;
    }
}


