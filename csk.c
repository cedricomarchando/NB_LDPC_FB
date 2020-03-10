
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

            for(i=0; i<csk->PNsize; i++)
            {
                csk->PN[i] =PN1024[i];
                printf(" %d ",csk->PN[i]);
            }
            getchar();
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

            for(i=0; i<csk->PNsize; i++)
            {
                csk->PN[i] =PN4096[i];
                printf(" %d ",csk->PN[i]);
            }
            getchar();


        }

    printf(" \n PN generation: Success\n");

}




void CHU_Generator( float *chu_real,float *chu_imag,int N)
{
    int i;
    int R=1;
if(N == 64) R=1;
//if(N == 128) R=3;
    if(N == 256) R=1;
//    if(N == 512) R=7;
    if(N == 1024) R=1;

        for (i=0; i<N; i++)
        {
            chu_real[i] = cos( i*R*(i+2)*PI /N);
            chu_imag[i] = sin( i*R*(i+2)*PI /N);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }
    printf(" \n CHU generation: Success\n");
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
            //printf("CSK_index: %d \t GF_temp:%d \t PN_index:%d \n",i,GF_temp,PN_index); //getchar();
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
                temp = NoisyBin[n][q] *   csk->CSK_arr[g][q];
                //temp = SQR(NoisyBin[n][q]+ csk->CSK_arr[g][q])/(2.0*SQR(sigma));
                TMP[g] = TMP[g] - temp;
            }
        }
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }
    }

//    // !!!!!
//    for (n=0; n<N/2; n++)
//    {
//        for (g=0; g<code->GF; g++)
//        {
//            TMP[g]=0.0;
//            for (q=0; q<csk->PNsize-1; q++)
//            {
//                temp = NoisyBin[n][q] *   csk->CSK_arr[g][q];
//                TMP[g] = TMP[g] - temp;
//            }
//        }
//        for(k=0; k<code->GF; k++)
//        {
//            decoder->intrinsic_LLR[n][k] = TMP[k];
//        }
//    }




    for (n=0; n<N; n++) free(NCSK[n]);
    free(NCSK);

    for (n=0; n<N; n++) free(NoisyBin[n]);
    free(NoisyBin);



}








/*!
 * \fn ModelChannel_AWGN_64 (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN, int *init_rand)
 * \brief 64 modulation + AWGN noise on the codeword.
 *                 The function computes the intrinsic_LLRs corresponding to the noisy observations.
 * Inputs
 * 	- NBIN : Binary copy of the codeword
 * 	- EbN  : Signal to noise ratio in terms of Eb/No (dB)
 * Outputs
 *      - decoder->intrinsic_LLR
 */
void ModelChannel_AWGN_64_CSK(csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;


    float **NoisyBin = calloc(csk->PNsize,sizeof(float *));
    for (q=0; q<csk->PNsize; q++) NoisyBin[q] = calloc(2,sizeof(float));
    /* Binary-input AWGN channel : */

    int i;

    float modulation[code->GF][2];

    const int CCSK64[128]={

//        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
//    41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,

            2,57,33,29,26,4,58,10,27,43,37,7,19,20,0,63,6,51,8,36,9,62,53,15,18,23,5,34,13,55,45,16,46,35,48,21,61,22,40,12,56,47,24,1,39,
    41,42,44,38,59,32,49,30,17,50,52,11,54,14,3,60,31,25,28,

    0,21,61,48,2,8,57,49,19,14,62,53,16,6,20,45,22,24,26,54,28,30,1,5,59,10,60,43,55,32,11,46,18,34,52,36,9,38,12,58,40,
    23,25,50,27,42,13,51,29,31,56,4,33,17,44,3,63,35,37,39,41,15,47,7,


//    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
//    41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63


    };


    ////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < code->GF ; i++)
    {
        // norm_factor = table_64QAM[i][0]*table_64QAM[i][0] + table_64QAM[i][1]*table_64QAM[i][1]+norm_factor;//compute sum
        norm_factor = table_64APSK[i][0]*table_64APSK[i][0] + table_64APSK[i][1]*table_64APSK[i][1]+norm_factor;//compute sum
    }
    norm_factor = sqrt( code->GF / norm_factor);
// printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< code->GF ; i++)
    {
        modulation[i][0]=norm_factor*table_64APSK[i][0];
        modulation[i][1]=norm_factor*table_64APSK[i][1];
    }


    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {
        som=0;

        for (q=0; q<6; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        //printf("\n %d \n",som );

        for (q=0; q<csk->PNsize; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                NoisyBin[q][i] = modulation[CCSK64[som+q]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }

         for(k=0; k<code->GF; k++)
        {
            TMP[k] =0;
        }

        for (g=0; g<csk->PNsize; g++)
        {
            for(k=0; k<code->GF; k++)
            {
                som=0;

                for (q=0; q<6; q++)
                {
                    som = som + BinGF_64[k][q]*pow(2,q);
                }
                TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[som+g]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[som+g]][1])/(2.0*SQR(sigma));
                //printf("%d %f \n",k, TMP[k]);
            }
        }


        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }
    }

    for (q=0; q<csk->PNsize; q++) free(NoisyBin[q]);
    free(NoisyBin);
}





void ModelChannel_CHU_CSK(float *chu_real,float *chu_imag, csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;
    int transmited = csk->PNsize;
    //int transmited = 12;

    float **NoisyBin = calloc(transmited,sizeof(float *));
    for (q=0; q<transmited; q++) NoisyBin[q] = calloc(2,sizeof(float));
    int i;

    //sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));
    sigma = sqrt(1.0/pow(10,EbN/10.0));

    for (n=0; n<N; n++)
    {
        som=0;
        for (q=0; q<6; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        for (q=0; q<transmited; q++)
        {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                NoisyBin[q][0] = chu_real[(som+q)%csk->PNsize]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                NoisyBin[q][1] = chu_imag[(som+q)%csk->PNsize]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
        }

         for(k=0; k<code->GF; k++)
        {
            TMP[k] =0;
        }

        for (g=0; g<transmited; g++)
        {
            for(k=0; k<code->GF; k++)
            {
                som=0;

                for (q=0; q<6; q++)
                {
                    som = som + BinGF_64[k][q]*pow(2,q);
                }
                TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-chu_real[(som+g)%csk->PNsize])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-chu_imag[(som+g)%csk->PNsize])/(2.0*SQR(sigma));
                //printf("%d %f \n",k, TMP[k]);
            }
        }

        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }
    }

    for (q=0; q<transmited; q++) free(NoisyBin[q]);
    free(NoisyBin);
}












void ModelChannel_AWGN_256_CSK(csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;


    float **NoisyBin = calloc(csk->PNsize,sizeof(float *));
    for (q=0; q<csk->PNsize; q++) NoisyBin[q] = calloc(2,sizeof(float));
    /* Binary-input AWGN channel : */

    int i;

    float modulation[code->GF][2];

    // test0
    const int CCSK256[256]={56, 107, 10, 45, 108, 177, 227, 194, 22, 98, 236, 36, 66, 99, 234, 82, 62, 237,
     161, 198, 83, 252, 119, 239, 104, 151, 6, 202, 64, 157, 144, 78, 43, 253, 152, 84, 243, 168, 229, 129, 77, 224, 213, 187, 254, 29,
     175, 75, 154, 28, 170, 135, 199, 55, 225, 127, 32, 128, 5, 102, 65, 126, 114, 74, 35, 54, 130, 200, 46, 146, 211, 49, 110, 172, 94,
     191, 121, 72, 141, 38, 14, 118, 91, 87, 183, 52, 134, 116, 180, 142, 18, 31, 188, 81, 51, 122, 53, 186, 164, 23, 255, 185, 206, 79, 217,
      181, 133, 25, 250, 16, 233, 136, 228, 249, 21, 92, 232, 12, 189, 167, 89, 50, 100, 140, 238, 148, 70, 245, 37, 223, 8, 179, 85, 44, 219,
       113, 220, 201, 58, 30, 93, 176, 57, 158, 155, 88, 17, 203, 212, 208, 147, 24, 150, 131, 138, 90, 159, 197, 0, 153, 247, 244, 207, 156,
       69, 231, 178, 117, 105, 80, 41, 34, 137, 195, 26, 242, 86, 60, 184, 222, 173, 95, 165, 143, 3, 7, 125, 210, 19, 209, 221, 193, 160, 2,
       48, 47, 139, 241, 218, 248, 163, 226, 76, 115, 27, 20, 4, 251, 42, 111, 182, 192, 101, 162, 103, 39, 132, 61, 169, 205, 73, 97, 196, 190,
       71, 171, 230, 214, 68, 13, 112, 145, 124, 1, 174, 215, 123, 67, 106, 96, 246, 240, 40, 33, 109, 15, 149, 235, 59, 63, 216, 11, 204, 9,
       120, 166};

       //test1
//const int CCSK256[256]={130, 74, 86, 16, 234, 15, 87, 80, 118, 223, 97, 27, 198, 9, 153, 103, 62, 60, 2, 210, 196, 197, 184, 66, 139, 26, 53, 145, 151, 208, 119, 96,
//75, 220, 163, 147, 10, 91, 61, 207, 149, 48, 195, 152, 22, 177, 105, 238, 182, 143, 20, 95, 63, 248, 225, 18, 132, 52, 222, 90, 172, 133, 70,
//69, 237, 137, 40, 157, 47, 64, 44, 111, 174, 201, 39, 49, 219, 169, 235, 55, 199, 121, 251, 71, 215, 156, 13, 54, 144, 154, 200, 230, 94, 247,
//1, 242, 30, 244, 32, 56, 236, 193, 167, 240, 214, 46, 155, 136, 146, 178, 6, 253, 212, 218, 99, 173, 25, 164, 43, 129, 188, 7, 59, 190, 24, 38,
//11, 114, 203, 14, 221, 35, 171, 131, 122, 73, 148, 37, 36, 192, 51, 93, 45, 128, 42, 216, 209, 176, 175, 34, 79, 88, 135, 78, 202, 12, 227, 141,
// 231, 134, 106, 217, 0, 58, 232, 82, 72, 161, 92, 252, 102, 109, 98, 150, 185, 81, 23, 31, 160, 110, 124, 3, 65, 83, 213, 50, 246, 243, 186, 159,
//  100, 179, 205, 194, 140, 123, 255, 249, 245, 68, 187, 206, 162, 107, 181, 76, 17, 165, 104, 112, 166, 33, 254, 57, 241, 8, 120, 158, 170, 168,
//   115, 126, 183, 101, 127, 85, 28, 125, 117, 224, 41, 204, 67, 189, 138, 116, 226, 5, 250, 142, 233, 191, 108, 84, 180, 77, 113, 21, 4, 239, 89,
//    228, 229, 211, 29, 19};
//
//
//    };



float table_GF[code->GF][2];

    int ring,angle;
    float constellation_radius[8]= {1, 1.794, 2.409, 2.986, 3.579, 4.045, 4.5, 5.2};

    for (ring=0; ring<8; ring++)
    {
        for (angle=0; angle<32; angle++)
        {
            table_GF[ring*32+angle][0]= constellation_radius[ring]*cos(angle*PI/16);
            table_GF[ring*32+angle][1]= constellation_radius[ring]*sin(angle*PI/16);;
        }

    }







    ////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < code->GF ; i++)
    {
        //norm_factor = table_256_circleQAM[i][0]*table_256_circleQAM[i][0] + table_256_circleQAM[i][1]*table_256_circleQAM[i][1]+norm_factor;//compute sum
        norm_factor = table_GF[i][0]*table_GF[i][0] + table_GF[i][1]*table_GF[i][1]+norm_factor;//compute sum

    }
    norm_factor = sqrt( code->GF / norm_factor);
 //printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< code->GF ; i++)
    {
        //modulation[i][0]=norm_factor*table_256_circleQAM[i][0];
        //modulation[i][1]=norm_factor*table_256_circleQAM[i][1];
        modulation[i][0]=norm_factor*table_GF[i][0];
        modulation[i][1]=norm_factor*table_GF[i][1];
    }


    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));


    for (n=0; n<N; n++)
    {
        som=0;

        for (q=0; q<8; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        //printf("\nGF: %d \n",som );

        for (q=0; q<csk->PNsize; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                NoisyBin[q][i] = modulation[CCSK256[(som+q)%256]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }


         for(k=0; k<code->GF; k++)
        {
            TMP[k] =0;
        }

        for (g=0; g<csk->PNsize; g++)
        {
            for(k=0; k<code->GF; k++)
            {
                som=0;

                for (q=0; q<8; q++)
                {
                    som = som + BinGF_256[k][q]*pow(2,q);
                }
                TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK256[(som+g)%256]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK256[(som+g)%256]][1])/(2.0*SQR(sigma));
                //printf("%d %f \n",k, TMP[k]);
            }
        }


        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }
    }

    for (q=0; q<csk->PNsize; q++) free(NoisyBin[q]);
    free(NoisyBin);
}




