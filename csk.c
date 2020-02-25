
/*!
---------------------------------------------------------------

 *\file csk.h
 * \brief csk = code shift keying
 * \author Author : Oussama ABASSI & Cedric MARCHAND
 * \copyright BSD copyright
Description :
    - This file includes definition of all tables and functions associated to csk modulation

Created : Thursday February 21, 2011
Modified : 2020

---------------------------------------------------------------
*/

#include <math.h>
#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"
#include "./include/channel.h"
#include "./include/csk.h"

#define PI 3.1415926536

/*
Name: void PNGenerator( int seed_init, csk_t *csk)
Brief: This function generates a random PN sequence
Arguments:
    - csk : structure describing csk modulation parameters
Return: void
*/
void PNGenerator( csk_t *csk)
{
    int i;
    int a;
    int LowBit;
        // for GF(64)
        if(csk->PNsize == 64)
        {
            a=32;
            for(i=0; i<csk->PNsize; i++)
            {
                //*******************************************
                //** primitive polynomial x**6+x**5+x**4+x+1
                //*******************************************
                csk->PN[i]=BPSK((a >> 5));//the output
                //feedback and shift
                LowBit = a >> 5; //x**6 term
                LowBit = LowBit ^ (a >> 4); //x**5 term
                LowBit = LowBit ^ (a >> 3); //x**4 term
                LowBit = LowBit ^ a; //x term
                a = ((a << 1) + (LowBit & 1)) & 0x03F; //return shifted 6 bit value
                //printf(" %d ",csk->PN[i] );
            }
        //getchar();
        }
        // For GF(256) with polynomial P(x)=X^8+ X^4+X^3+X^2+1
        if(csk->PNsize == 256)
        {
                    a=128;
            for(i=0; i<csk->PNsize; i++)
            {
                csk->PN[i]=BPSK((a >> 7));//the output
                //feedback and shift
                LowBit = a >> 7; //x**8 term

                LowBit = LowBit ^ (a >> 3); //x**4 term
                LowBit = LowBit ^ (a >> 2); //x**3 term
                LowBit = LowBit ^ (a >> 1); //x**2 term
//
////                LowBit = LowBit ^ (a >> 1); //x**2 term
////                LowBit = LowBit ^ a; //x term
//
//                LowBit = LowBit ^ (a >> 6); //x**7 term
//                LowBit = LowBit ^ (a >> 5); //x**6 term
//                LowBit = LowBit ^ (a >> 4); //x**5 term
//                LowBit = LowBit ^ (a >> 3); //x**4 term
//                LowBit = LowBit ^ (a >> 1); //x**2 term


                a = ((a << 1) + (LowBit & 1)) & 0x0FF; //return shifted 8 bit value

            //printf(" %d ",csk->PN[i] );
            }
            //getchar();







//            for(i=1; i<csk->PNsize-1; i++)
//            {
//                csk->PN[i-1]=BPSK(BinGF_256[i][0]);
//
//            //printf(" %d ",csk->PN[i-1] );
//            }
//            csk->PN[csk->PNsize-1]=BPSK(1);
//            //printf(" %d ",csk->PN[256-1] );
//            //getchar();

        }




// For GF(512) with polynomial P(x)=X^9+X^4+1
        if(csk->PNsize == 512)
        {
                    a=256;
            for(i=0; i<csk->PNsize; i++)
            {
                csk->PN[i]=BPSK((a >> 8));//the output
                //feedback and shift
                LowBit = a >> 8; //x**9 term
                LowBit = LowBit ^ (a >> 3); //x**4 term



                a = ((a << 1) + (LowBit & 1)) & 0x1FF; //return shifted 9 bit value
            }
        }







// For GF(1024) with polynomial P(x)=X^10+X^3+1
        if(csk->PNsize == 1024)
        {
                    a=512;
            for(i=0; i<csk->PNsize; i++)
            {
                csk->PN[i]=BPSK((a >> 9));//the output
                //feedback and shift
                LowBit = a >> 9; //x**10 term
                LowBit = LowBit ^ (a >> 2); //x**3 term



                a = ((a << 1) + (LowBit & 1)) & 0x3FF; //return shifted 10 bit value
            }
        }


// For GF(2048) with polynomial P(x)=X^11+X^2+1
        if(csk->PNsize == 2048)
        {
                    a=1024;
            for(i=0; i<csk->PNsize; i++)
            {
                csk->PN[i]=BPSK((a >> 10));//the output
                //feedback and shift
                LowBit = a >> 10; //x**11 term
                LowBit = LowBit ^ (a >> 1); //x**2 term



                a = ((a << 1) + (LowBit & 1)) & 0x07FF; //return shifted 11 bit value
            }
        }




// For GF(4096) with polynomial P(x)=X^12+X^9+X^7+X^6+1
        if(csk->PNsize == 4096)
        {
                    a=2048;
            for(i=0; i<csk->PNsize; i++)
            {
                csk->PN[i]=BPSK((a >> 11));//the output
                //feedback and shift
                LowBit = a >>11; //x**12 term

                LowBit = LowBit ^ (a >> 8); //x**9 term
                LowBit = LowBit ^ (a >> 6); //x**7 term
                LowBit = LowBit ^ (a >> 5); //x**6 term
                LowBit = LowBit ^ (a); //x term





                a = ((a << 1) + (LowBit & 1)) & 0x0FFF; //return shifted 12 bit value
            }
        }

    printf(" \n PN generation: Success\n");

}

/*
Name: shiftPN(int GF, csk_t *csk)
Brief: This function generates CSK_arr by shifting PN vector
Arguments:
    - GF : galois field order
    - csk : structure describing csk modulation parameters
Return: void
*/
void build_natural_csk_mapping(int GF, csk_t *csk)
{
    int i, j;

    for (i=0; i<csk->PNsize; i++)
    {
        csk->CSK_arr[0][i]=csk->PN[i];
        //printf("%d ",csk->PN[i] );

    }

//getchar();


    for (i=1; i<GF; i++)
    { //printf(" \n %d:",i);

        for (j=0; j<csk->PNsize; j++)
        {
            csk->CSK_arr[i][j]=csk->CSK_arr[i-1][(j+1)%csk->PNsize];
            //printf("%d ",(csk->CSK_arr[i][j]+1)/2 );

        }
        //printf(" \n ");
    }


    //getchar();
    printf("Filling CSK_arr: Success\n");
    fflush(stdout);
}


// search in the csk array the GF indices giving first logGF bits that are different
void build_punctured_csk_mapping(int GF,int logGF, csk_t *csk, int **BINGF)
{
    int i, j;
    int bin_temp[12];
    int GF_temp;
    int GF_selected[GF];
    int PN_index =0;

    for (j=0; j<GF; j++) {GF_selected[j]= 0;}


    for (i=0; i<csk->PNsize; i++)
    { //printf(" \n %d:",i);

        for (j=0; j<logGF; j++)
        {
            bin_temp[j]=(csk->PN[(j+i)%csk->PNsize] +1)/2;
            //printf("%d ",bin_temp[j]);
        }
        GF_temp = Bin2GF(bin_temp,GF,logGF,BINGF);
        //printf(" GF_temp: %d", GF_temp ); getchar();
        if (GF_selected[GF_temp]==0)
        {
            for (j=0; j<csk->PNsize; j++)
            {
                csk->CSK_arr[PN_index][j] = csk->PN[(j+i)%csk->PNsize];
            }
            GF_selected[GF_temp]=1;
            //printf("CSK_index: %d \t GF_temp:%d \t PN_index:%d ",i,GF_temp,PN_index); getchar();
            PN_index = PN_index + 1;
        }
        //printf(" \n ");
    }


    //getchar();
    printf("build_punctured_mapping: Success\n");
    fflush(stdout);
}





// construction of a specific CSK mapping, not based on a PN
void build_CSK_map(code_t *code, csk_t *csk)
{
    //construction of the table
    int tmp1,i,j,k;

// int tmp1,tmp2,tmp3,tmp4
//    for (i=0; i<code->GF; i++)
//    {
//        for (j=0; j<csk->PNsize; j++)
//        {
//            //matrix multiplication
//            tmp4=0;
//            for (k=0;k<code->logGF;k++ )
//            {
//
//                //printf(" GF:%d k:%d shiftk:%d  mask:%d  \n", i, k, i>>k, (i>>k) & 0x01 );
//
//                tmp1=(i>>k) & 0x01;
//                tmp2=(j>>k) & 0x01;
//                tmp3=tmp1 & tmp2;
//                tmp4=tmp4 ^ tmp3;
//                //printf(" tmp1:%d  tmp2:%d  tmp3:%d tmp4:%d \n",tmp1, tmp2, tmp3, tmp4 );
//            }
//            csk->CSK_arr[i][j]= BPSK(tmp4);
//        //printf("tmp4:%d, map[%d][%d]=%d", tmp4,i,j, csk->CSK_arr[i][j] );
//        //getchar();
//        }
//    }


       // repeat construction
    for (i=0; i<code->GF; i++)
    {
        for (j=0; j<csk->PNsize / code->logGF; j++)
        {
            for (k=0;k<code->logGF;k++ )
            {
                tmp1=(i>>k) & 0x01;
                //printf("CSK[%d][%d] = %d \n ",i,j*code->logGF + k, tmp1 );
                csk->CSK_arr[i][j*code->logGF + k]= BPSK(tmp1);
            }
        }
        //getchar();
    }




            printf("Filling CSK_arr with matrix: Success\n");
    fflush(stdout);


}

/*
Name : GF_to_CSK
Brief: This function converts a GF symbols of a given codeword into its corresponding CCSK sequences
Argument:
    - CodeWord: the codeword to be converted
    - N: the length of the codeword
    - csk: structure describing csk modulation parameters
    - NCSK : the resulting CCSK sequences
Return: void
*/
void GF_to_CSK(int* CodeWord, int N, csk_t *csk, int** NCSK)
{
    int i,j; //counters
    for (i=0; i<N; i++)
    {

        //printf("GF= %d \t:",CodeWord[i]);
        for (j=0; j<csk->PNsize; j++)
        {
            NCSK[i][j] = csk->CSK_arr[CodeWord[i]][j];
            //printf("%d ",NCSK[i][j]);
        }
        //getchar();
    //printf(" \n");

    }
    //getchar();

}

/*
Name: allocate_csk(csk_t csk, int GF, int N)
Brief: This function allocates memory to a csk_t structure
Arguments:
    - csk    :  csk structure
    - GF     :  order of the galois field
    - N      :  the code length
Return: void
*/
void allocate_csk(csk_t *csk, int GF)
{
    int i;

    //allocate space of memory for CSK_arr[GF][csk_bitsize]
    csk->CSK_arr = calloc(GF,sizeof(int *));
    csk->CSK_arr[0] = calloc (GF*csk->PNsize,sizeof(int));
    for ( i=1; i<GF; i++ )
    {
        csk->CSK_arr[i] = csk->CSK_arr[0] + i*csk->PNsize;
    }

    //allocate space of memory for PN[PNsize]
    csk->PN = calloc( csk->PNsize, sizeof(int));
    printf("Allocate memory space for CSK structure: Success\n");
    fflush(stdout);

}

/*
Name: free_csk(csk_t *csk, int GF, int N)
Brief: free allocated memory
Arguments:
    - csk    :  csk structure
    - GF     :  order of the galois field
    - N      :  the code length
Return: void
*/
void free_csk(csk_t *csk)
{
    free(csk->PN);
    free(csk->CSK_arr[0]); free(csk->CSK_arr);
}


/*!
 * \fn void ModelChannel_AWGN_BPSK (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN,int *init_rand)
 * \brief BPSK modulation + AWGN noise on the codeword.
 *                 The function computes the intrinsic_LLRs corresponding to the noisy observations.
 * Inputs
 * 	- NBIN : Binary copy of the codeword
 * 	- EbN  : Signal to noise ratio in terms of Eb/No (dB)
 * Outputs
 *      - decoder->intrinsic_LLR
 */
void ModelChannel_AWGN_BPSK_CSK (csk_t *csk, code_t *code, decoder_t *decoder, table_t *table, int *CodeWord, float EbN,int *init_rand, float quantif_range_int,float quantif_range_float)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    float temp;

    float **NoisyBin = calloc(N,sizeof(float *));
    for (n=0; n<N; n++) {NoisyBin[n] = calloc(csk->PNsize,sizeof(float));}

    int **NCSK= calloc(N,sizeof(int *));
    for (n=0; n<N; n++) {NCSK[n] = calloc(csk->PNsize,sizeof(int));}


    //Map symbols to CSK
    for (n=0; n<N; n++)
    {

        //printf("GF= %d \t:",CodeWord[n]);
        for (q=0; q<csk->PNsize; q++)
        {
            NCSK[n][q] = csk->CSK_arr[CodeWord[n]][q];
            //printf("%d ",NCSK[n][q]);
        }
       // getchar();
    //printf(" \n");
    }





    /* Binary-input AWGN channel : */
    //sigma = sqrt(1.0/(2.0*code->rate*pow(10,EbN/10.0))); // for EbNo
    sigma = sqrt(1.0/(pow(10,EbN/10.0))); // for EsNo

    for (n=0; n<N; n++)
    {
        for (q=0; q<csk->PNsize; q++)
        {
            u=My_drand48(init_rand);
            while(u==0.0)
			{
				u=My_drand48(init_rand);
			}
            v=My_drand48(init_rand);
            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
            NoisyBin[n][q] = NCSK[n][q] + sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v);
        }

    }


    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n=0; n<N; n++)
    {
        for (g=0; g<code->GF; g++)
        {
            TMP[g]=0.0;
            for (q=0; q<csk->PNsize; q++)
            {
                //TMP[g] = TMP[g] + SQR(NoisyBin[n][q]-BPSK(table->BINGF[g][q]))/(2.0*SQR(sigma));

                //TMP[g] = TMP[g] + SQR(NoisyBin[n][q]-   csk->CSK_arr[g][q]    )/(2.0*SQR(sigma));


                temp = NoisyBin[n][q] *   csk->CSK_arr[g][q];
                //temp = BPSK_quantif(temp, quantif_range_int, quantif_range_float);

                TMP[g] = TMP[g] - temp;


            }

        }

        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }



    }



    for (n=0; n<N; n++) free(NCSK[n]);
    free(NCSK);

    for (n=0; n<N; n++) free(NoisyBin[n]);
    free(NoisyBin);



}

