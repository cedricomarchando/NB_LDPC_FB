
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


        if(csk->PNsize == 128)
        {
            a=64;
            for(i=0; i<csk->PNsize; i++)
            {
                //*******************************************
                //** primitive polynomial x**7+x**3+1
                //*******************************************
                csk->PN[i]=BPSK((a >> 6));//the output
                //feedback and shift
                LowBit = a >> 6; //x**7 term
                LowBit = LowBit ^ (a >> 2); //x**3 term
                a = ((a << 1) + (LowBit & 1)) & 0x07F; //return shifted 6 bit value
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
                //printf(" %d ",csk->PN[i]);
            }
            //getchar();
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
                //printf(" %d ",csk->PN[i]);
            }
            //getchar();


        }

    printf(" \n PN generation: Success\n");

}




void CHU_Generator( float *chu_real,float *chu_imag,int N)
{
    int i;
    int R=1;
//if(N == 64) R=1;
////if(N == 128) R=3;
//    if(N == 256) R=1;
////    if(N == 512) R=7;
//    if(N == 1024) R=1;

        for (i=0; i<N; i++)
        {
            chu_real[i] = cos( i*R*(i)*PI /N);
            chu_imag[i] = sin( i*R*(i)*PI /N);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }
    printf(" \n CHU generation: Success\n");
}




void CHU_AM_Generator( float *chu_real,float *chu_imag,int N)
{
    int i;
    int R=1;
    float chu_energie;
    float norm_factor;


        for (i=0; i<N; i++)
        {
            chu_real[i] = cos( i*R*(i+1)*PI /N);
            chu_imag[i] = sin( i*R*(i+1)*PI /N);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }


//        for (i=1; i<N/2 ; i++)
//        {
//            chu_real[i] = chu_real[i]*i;
//            chu_imag[i] = chu_imag[i]*i;
//            chu_real[N-i] = chu_real[N-i]*i;
//            chu_imag[N-i] = chu_imag[N-i]*i;
//        }
//        chu_real[N/2] = chu_real[N/2]*N/2;
//        chu_imag[N/2] = chu_imag[N/2]*N/2;

        for (i=0; i<N ; i++)
        {
            chu_real[i] = chu_real[i]*sin((i+1)*PI/N);
            chu_imag[i] = chu_imag[i]*sin((i+1)*PI/N);
        }


//        for (i=0; i<N; i++)
//        {
//           printf("%d \t %f \t %f \t %f\n",i,chu_real[i],chu_imag[i],sin((i+1)*PI/N));
//        }
//        getchar();



        chu_energie = 0.0;
        for (i=0; i<N; i++)
        {
            chu_energie = chu_energie + chu_real[i]*chu_real[i] + chu_imag[i]* chu_imag[i];
        }
        norm_factor=sqrt(N/chu_energie);

        for (i=0; i<N; i++)
        {
            chu_real[i] = norm_factor * chu_real[i];
            chu_imag[i] = norm_factor *chu_imag[i];
        }

//        for (i=0; i<256; i++)
//        {
//           printf("%d  %f  %f \n",i, chu_real[i],chu_imag[i]);
//        }
//        getchar();


    printf(" \n CHU generation: Success\n");
}




void CHU_Generator_64apsk( float *chu_real,float *chu_imag,int N)
{
    int i;
    int R=1;
    float chu_energie;
    float norm_factor;
    int nb_c1,nb_c2,nb_c3,nb_c4;
    float r_c2,r_c3,r_c4;

    nb_c1=32;
    nb_c2=16;
    nb_c3=8;
    nb_c4=8;

    r_c2=0.8;
    r_c3=0.55;
    r_c4=0.3;

        for (i=0; i<nb_c1; i++)
        {
            chu_real[i] = cos( i*R*(i+1)*PI /nb_c1);
            chu_imag[i] = sin( i*R*(i+1)*PI /nb_c1);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }

        for (i=0; i<nb_c2; i++)
        {
            chu_real[i+nb_c1] = r_c2*cos( i*R*(i+1)*PI /nb_c2 + PI);
            chu_imag[i+nb_c1] = -r_c2*sin( i*R*(i+1)*PI /nb_c2 + PI);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }


        for (i=0; i<nb_c3; i++)
        {
            chu_real[i+nb_c1 +nb_c2] =r_c3* cos( i*R*(i+1)*PI /nb_c3);
            chu_imag[i+nb_c1 + nb_c2] = r_c3*sin( i*R*(i+1)*PI /nb_c3);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }

        for (i=0; i<nb_c4; i++)
        {
            chu_real[i+nb_c1 +nb_c2 + nb_c3] = r_c4* cos( i*R*(i+1)*PI /nb_c4 );
            chu_imag[i+nb_c1 +nb_c2 + nb_c3] = -r_c4*sin( i*R*(i+1)*PI /nb_c4 );
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }



        chu_energie = 0.0;
        for (i=0; i<64; i++)
        {
            chu_energie = chu_energie + chu_real[i]*chu_real[i] + chu_imag[i]* chu_imag[i];
        }
        norm_factor=sqrt(64/chu_energie);

        for (i=0; i<64; i++)
        {
            chu_real[i] = norm_factor * chu_real[i];
            chu_imag[i] = norm_factor *chu_imag[i];
        }

//        for (i=0; i<64; i++)
//        {
//           printf("%f  %f \n",chu_real[i],chu_imag[i]);
//        }
//        getchar();


    printf(" \n CHU generation: Success\n");
}






void CHU_Generator_256apsk( float *chu_real,float *chu_imag,int N)
{
    int i;
    int R=1;
    float chu_energie;
    float norm_factor;
    int nb_c1,nb_c2,nb_c3,nb_c4,nb_c5,nb_c6;
    float r_c1,r_c2,r_c3,r_c4,r_c5,r_c6;

    nb_c1=16;
    nb_c2=16;
    nb_c3=32;
    nb_c4=64;
    nb_c5=64;
    nb_c6=64;

    r_c1=0.15;
    r_c2=0.25;
    r_c3=0.4;
    r_c4=0.6;
    r_c5=0.8;
    r_c6=1.0;

        for (i=0; i<nb_c1; i++)
        {
            chu_real[i] = r_c1*cos( i*R*(i+1)*PI /nb_c1);
            chu_imag[i] = -r_c1*sin( i*R*(i+1)*PI /nb_c1);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }

        for (i=0; i<nb_c2; i++)
        {
            chu_real[i+nb_c1] = r_c2*cos( i*R*(i+1)*PI /nb_c2 + PI);
            chu_imag[i+nb_c1] = r_c2*sin( i*R*(i+1)*PI /nb_c2 + PI);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }


        for (i=0; i<nb_c3; i++)
        {
            chu_real[i+nb_c1 +nb_c2] =r_c3* cos( i*R*(i+1)*PI /nb_c3);
            chu_imag[i+nb_c1 + nb_c2] = r_c3*sin( i*R*(i+1)*PI /nb_c3);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }
R=11;
        for (i=0; i<nb_c4; i++)
        {
            chu_real[i+nb_c1 +nb_c2 + nb_c3] = r_c4* cos( i*R*(i+1)*PI /nb_c4 );
            chu_imag[i+nb_c1 +nb_c2 + nb_c3] = r_c4*sin( i*R*(i+1)*PI /nb_c4 );
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }
R=1;
        for (i=0; i<nb_c5; i++)
        {
            chu_real[i+nb_c1 +nb_c2 + nb_c3+nb_c4] = r_c5* cos( i*R*(i+1)*PI /nb_c5 );
            chu_imag[i+nb_c1 +nb_c2 + nb_c3+nb_c4] = -r_c5*sin( i*R*(i+1)*PI /nb_c5 );
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }

        for (i=0; i<nb_c6; i++)
        {
            chu_real[i+nb_c1 +nb_c2 + nb_c3+nb_c4 +nb_c5] = r_c6* cos( i*R*(i+1)*PI /nb_c6 );
            chu_imag[i+nb_c1 +nb_c2 + nb_c3+nb_c4 +nb_c5] = r_c6*sin( i*R*(i+1)*PI /nb_c6 );
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }

        chu_energie = 0.0;
        for (i=0; i<256; i++)
        {
            chu_energie = chu_energie + chu_real[i]*chu_real[i] + chu_imag[i]* chu_imag[i];
        }
        norm_factor=sqrt(256/chu_energie);

        for (i=0; i<256; i++)
        {
            chu_real[i] = norm_factor * chu_real[i];
            chu_imag[i] = norm_factor *chu_imag[i];
        }

//        for (i=0; i<256; i++)
//        {
//           printf("%d %f  %f \n",i,chu_real[i],chu_imag[i]);
//        }
//        getchar();


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
void ModelChannel_AWGN_BPSK_CSK (csk_t *csk, code_t *code, decoder_t *decoder, table_t *table, int *CodeWord, float EbN,int *init_rand)
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
                //temp = -SQR(NoisyBin[n][q]- csk->CSK_arr[g][q])/(2.0*SQR(sigma)); //distance
                TMP[g] = TMP[g] - temp;
            }
        }
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }
    }

//    // !!!!! for hybrid CCSK
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

  //      0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
  //  41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,


//            2,57,33,29,26,4,58,10,27,43,37,7,19,20,0,63,6,51,8,36,9,62,53,15,18,23,5,34,13,55,45,16,46,35,48,21,61,22,40,12,56,47,24,1,39,
//    41,42,44,38,59,32,49,30,17,50,52,11,54,14,3,60,31,25,28,
//
//    0,21,61,48,2,8,57,49,19,14,62,53,16,6,20,45,22,24,26,54,28,30,1,5,59,10,60,43,55,32,11,46,18,34,52,36,9,38,12,58,40,
//    23,25,50,27,42,13,51,29,31,56,4,33,17,44,3,63,35,37,39,41,15,47,7



//    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
//    41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63


////new max 179
//17, 39, 6, 61, 54, 27, 38, 31, 18, 9, 43, 37, 56, 14, 35, 46, 21, 3, 57, 11, 47, 29, 1, 62, 55, 7, 49, 28, 24, 4, 44, 36, 53, 8,
//40, 41, 26, 42, 15, 30, 52, 59, 5, 25, 10, 0, 45, 58, 2, 32, 33, 34, 19, 13, 48, 50, 51, 22, 63, 12, 23, 60, 20, 16, 49, 7, 40, 56,
//35, 48, 45, 51, 5, 27, 31, 62, 0, 3, 9, 58, 63, 46, 60, 29, 23, 52, 43, 12, 11, 38, 44, 6, 18, 14, 15, 36, 42, 26, 50, 16, 20, 21, 32,
//41, 55, 47, 8, 30, 37, 59, 4, 34, 61, 39, 33, 54, 53, 19, 57, 25, 2, 17, 10, 13, 28, 22, 24, 1

//new2
//42,  10,  51,  24,  0,  35,  48,  57,  56,  45,  18,  26,  19,  5,  30,  17,  40,  31,  38,  36,  29,  47,  23,  43,  32,  15,  44,  33,  8,  52,  54,  49,  9,  61,  20,  22,  4,  55,  1,  34,  14,  50,  16,  2,  59,  53,  39,  21,  60,  12,  11,  41,  13,  3,  37,  7,  62,  58,  6,  46,  25,  28,  27,  63,  38,  18,  12,  33,  56,  44,  6,  46,  61,  32,  17,  51,  40,  34,  10,  49,  4,  43,  55,  57,  20,  39,  58,  1,  28,  16,  31,  45,  35,  41,  27,  7,  26,  22,  11,  25,  13,  59,  24,  9,  29,  21,  50,  15,  0,  2,  63,  36,  8,  42,  48,  23,  3,  60,  37,  30,  52,  54,  19,  62,  5,  47,  53,  14

//new3 585
//50,  0,  30,  8,  44,  15,  12,  5,  2,  47,  6,  29,  49,  11,  40,  28,  16,  32,  19,  62,  33,  63,  41,  3,  53,  24,  55,  13,  48,  31,  38,  56,  1,  61,  46,  20,  22,  10,  35,  18,  9,  14,  60,  27,  21,  54,  36,  26,  45,  59,  37,  4,  57,  39,  51,  58,  23,  52,  42,  17,  43,  34,  25,  7,  27,  40,  26,  24,  53,  1,  50,  49,  35,  2,  31,  41,  42,  39,  33,  8,  56,  29,  16,  45,  0,  5,  48,  37,  11,  60,  28,  18,  6,  62,  46,  9,  23,  12,  36,  51,  25,  30,  44,  52,  55,  22,  21,  47,  19,  32,  17,  7,  58,  20,  59,  61,  14,  10,  54,  3,  57,  38,  4,  43,  15,  63,  34,  13
//new 4 max198
//22,  8,  46,  16,  56,  1,  47,  32,  31,  28,  0,  23,  54,  2,  40,  11,  33,  38,  57,  9,  49,  3,  58,  37,  42,  12,  36,  5,  24,  20,  35,  52,  7,  45,  4,  26,  43,  13,  60,  29,  48,  50,  41,  59,  18,  44,  19,  61,  53,  15,  62,  21,  14,  51,  10,  34,  39,  63,  55,  25,  17,  27,  30,  6,  23,  2,  54,  17,  36,  60,  46,  38,  3,  41,  53,  29,  50,  55,  25,  56,  42,  49,  44,  18,  12,  45,  47,  8,  30,  39,  61,  6,  40,  19,  43,  0,  31,  7,  51,  10,  22,  13,  20,  58,  57,  26,  5,  59,  21,  14,  16,  15,  28,  4,  9,  63,  11,  34,  32,  35,  1,  37,  27,  48,  33,  62,  52,  24
//new 5 sum12149
//18,  55,  44,  4,  36,  27,  26,  35,  52,  24,  40,  51,  2,  48,  22,  7,  42,  58,  9,  34,  6,  14,  41,  53,  1,  61,  8,  0,  62,  3,  43,  16,  60,  17,  33,  49,  10,  5,  39,  59,  11,  29,  21,  47,  56,  54,  25,  46,  20,  23,  31,  32,  63,  28,  45,  38,  13,  30,  57,  15,  50,  19,  37,  12,  2,  41,  49,  33,  48,  29,  10,  55,  0,  12,  21,  38,  5,  24,  52,  7,  62,  50,  53,  40,  9,  57,  11,  1,  47,  23,  32,  8,  28,  56,  58,  39,  36,  14,  60,  37,  16,  44,  59,  46,  34,  51,  13,  19,  3,  22,  63,  45,  27,  30,  31,  43,  6,  26,  25,  20,  17,  42,  61,  15,  35,  54,  18,  4
// new 6 sum(ccorr.^2) 1M48 sum12194 max230
//1,  32,  18,  16,  28,  41,  35,  34,  27,  46,  29,  49,  15,  53,  26,  40,  3,  56,  57,  24,  48,  22,  50,  55,  23,  60,  51,  4,  5,  21,  33,  25,  17,  45,  47,  2,  20,  54,  11,  0,  37,  13,  10,  44,  30,  8,  31,  61,  6,  43,  14,  9,  63,  7,  12,  36,  58,  52,  62,  38,  59,  19,  39,  42,  8,  46,  36,  49,  15,  56,  37,  23,  41,  10,  48,  22,  40,  45,  43,  38,  2,  16,  52,  31,  25,  54,  47,  34,  12,  27,  28,  19,  29,  58,  61,  30,  35,  33,  3,  39,  59,  21,  50,  11,  57,  53,  1,  44,  7,  14,  42,  5,  0,  20,  4,  17,  55,  32,  13,  62,  26,  51,  6,  60,  18,  63,  9,  24
//new7 distance.^2 47M53
//32,  46,  62,  60,  39,  27,  35,  23,  22,  8,  29,  10,  53,  19,  18,  13,  42,  54,  45,  36,  44,  0,  6,  59,  16,  55,  28,  63,  9,  3,  51,  61,  56,  11,  25,  49,  41,  40,  14,  52,  43,  38,  31,  50,  1,  4,  17,  26,  21,  15,  58,  7,  5,  2,  37,  57,  24,  34,  20,  30,  12,  33,  48,  47,  0,  30,  55,  24,  31,  59,  49,  40,  5,  17,  1,  50,  62,  52,  43,  6,  16,  7,  28,  14,  25,  63,  58,  27,  4,  18,  57,  35,  42,  45,  44,  20,  19,  23,  37,  34,  2,  13,  39,  8,  9,  41,  47,  51,  15,  38,  10,  60,  54,  36,  61,  26,  12,  53,  32,  21,  48,  11,  29,  46,  56,  3,  33,  22

//sum(real(ccorr)=6820
20,  34,  46,  39,  25,  40,  60,  3,  10,  38,  32,  37,  28,  29,  43,  2,  5,  61,  0,  62,  49,  7,  12,  1,  41,  18,  14,  56,  51,  9,  53,  8,  55,  26,  15,  17,  23,  30,  11,  22,  48,  45,  63,  52,  47,  13,  57,  27,  4,  6,  35,  19,  24,  50,  21,  33,  54,  44,  58,  16,  36,  31,  42,  59,  15,  41,  13,  44,  28,  59,  31,  30,  23,  14,  25,  55,  63,  62,  48,  39,  40,  60,  46,  6,  38,  2,  10,  0,  22,  56,  26,  27,  34,  16,  51,  45,  8,  5,  3,  7,  47,  24,  21,  32,  43,  1,  9,  20,  18,  17,  36,  33,  12,  58,  53,  52,  35,  29,  50,  42,  61,  57,  19,  49,  11,  4,  37,  54
    };



    ////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < code->GF ; i++)
    {
        //norm_factor = table_64QAM[i][0]*table_64QAM[i][0] + table_64QAM[i][1]*table_64QAM[i][1]+norm_factor;//compute sum
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
                NoisyBin[q][i] = modulation[CCSK64[(som+q)%128]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }

         for(k=0; k<code->GF; k++)
        {
            TMP[k] =0.0;
        }


        if (n<N) // per symbol puncturing
        {

            for(k=0; k<code->GF; k++)
            {
                for (g=0; g<csk->PNsize; g++)
                {

                    som=0;

                    for (q=0; q<6; q++)
                    {
                        som = som + BinGF_64[k][q]*pow(2,q);
                    }
                    //TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%128]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%128]][1])/(2.0*SQR(sigma));
                    //printf("%d %f \n",k, TMP[k]);

                   //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    TMP[k] = TMP[k]-( NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][0] + NoisyBin[g][1]*modulation[CCSK64[(som+g)%128]][1] ); //real part


                }
            }
        }
        else
        {

            for(k=0; k<code->GF; k++)
            {
                for (g=0; g<csk->PNsize-1; g++)
                {

                    som=0;

                    for (q=0; q<6; q++)
                    {
                        som = som + BinGF_64[k][q]*pow(2,q);
                    }
                    //TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%128]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%128]][1])/(2.0*SQR(sigma));
                    //printf("%d %f \n",k, TMP[k]);
                    //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    TMP[k] = TMP[k]-( NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][0] + NoisyBin[g][1]*modulation[CCSK64[(som+g)%128]][1] );// + SQR(  NoisyBin[g][1] * modulation[CCSK64[(som+g)%128]][0] -  NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][1]        );

                }
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




//LLR generation using Zadoff-Chu sequence, with puncturing and hybrid puncturing
void ModelChannel_CHU_CSK(float *chu_real,float *chu_imag, csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;
    //int transmited = csk->PNsize;
    int transmited = 128 ;
    int mapping_step =2;
    int mappinp_start =0;

    float **NoisyBin = calloc(transmited,sizeof(float *));
    for (q=0; q<transmited; q++) NoisyBin[q] = calloc(2,sizeof(float));
    int i;

    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));
    //sigma = sqrt(1.0/pow(10,EbN/10.0));

    for (n=0; n<N; n++)
    {
        som=0;
        for (q=0; q<code->logGF; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        for (q=0; q<transmited; q++)
        {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                NoisyBin[q][0] = chu_real[(mappinp_start+som*mapping_step+q)%csk->PNsize]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                NoisyBin[q][1] = chu_imag[(mappinp_start+som*mapping_step+q)%csk->PNsize]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
        }

        for(k=0; k<code->GF; k++)
        {
            TMP[k] =0;
        }

        if(n< N) // for hybrid puncturing, set for example n<N/2
        {
            for (g=0; g<transmited; g++)
                {
                    for(k=0; k<code->GF; k++)
                    {
                        som=0;

                        for (q=0; q<code->logGF; q++)
                        {
                            //som = som + BinGF_256[k][q]*pow(2,q);
                            som = som + BinGF_64[k][q]*pow(2,q);
                        }
                        //TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-chu_real[(mappinp_start+som*mapping_step+g)%csk->PNsize])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-chu_imag[(mappinp_start+som*mapping_step+g)%csk->PNsize])/(2.0*SQR(sigma));
                        //printf("%d %f \n",k, TMP[k]);

                    //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    TMP[k] = TMP[k]-( NoisyBin[g][0]*chu_real[(mappinp_start+som*mapping_step+g)%csk->PNsize]     + NoisyBin[g][1]*chu_imag[(mappinp_start+som*mapping_step+g)%csk->PNsize]    );// + SQR(  NoisyBin[g][1] * modulation[CCSK64[(som+g)%128]][0] -  NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][1]        );


                    }
                }
        }
        else
        {
                for (g=0; g<transmited -1; g++)
                {
                    for(k=0; k<code->GF; k++)
                    {
                        som=0;

                        for (q=0; q<code->logGF ; q++)
                        {
                            //som = som + BinGF_256[k][q]*pow(2,q);
                            som = som + BinGF_64[k][q]*pow(2,q);
                        }
                        //TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-chu_real[(mappinp_start+som*mapping_step+g)%csk->PNsize])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-chu_imag[(mappinp_start+som*mapping_step+g)%csk->PNsize])/(2.0*SQR(sigma));
                        //printf("%d %f \n",k, TMP[k]);
                        TMP[k] = TMP[k]-( NoisyBin[g][0]*chu_real[(mappinp_start+som*mapping_step+g)%csk->PNsize]     + NoisyBin[g][1]*chu_imag[(mappinp_start+som*mapping_step+g)%csk->PNsize]    );// + SQR(  NoisyBin[g][1] * modulation[CCSK64[(som+g)%128]][0] -  NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][1]        );

                    }
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

//    // test0
//    const int CCSK256[512]={56, 107, 10, 45, 108, 177, 227, 194, 22, 98, 236, 36, 66, 99, 234, 82, 62, 237,
//     161, 198, 83, 252, 119, 239, 104, 151, 6, 202, 64, 157, 144, 78, 43, 253, 152, 84, 243, 168, 229, 129, 77, 224, 213, 187, 254, 29,
//     175, 75, 154, 28, 170, 135, 199, 55, 225, 127, 32, 128, 5, 102, 65, 126, 114, 74, 35, 54, 130, 200, 46, 146, 211, 49, 110, 172, 94,
//     191, 121, 72, 141, 38, 14, 118, 91, 87, 183, 52, 134, 116, 180, 142, 18, 31, 188, 81, 51, 122, 53, 186, 164, 23, 255, 185, 206, 79, 217,
//      181, 133, 25, 250, 16, 233, 136, 228, 249, 21, 92, 232, 12, 189, 167, 89, 50, 100, 140, 238, 148, 70, 245, 37, 223, 8, 179, 85, 44, 219,
//       113, 220, 201, 58, 30, 93, 176, 57, 158, 155, 88, 17, 203, 212, 208, 147, 24, 150, 131, 138, 90, 159, 197, 0, 153, 247, 244, 207, 156,
//       69, 231, 178, 117, 105, 80, 41, 34, 137, 195, 26, 242, 86, 60, 184, 222, 173, 95, 165, 143, 3, 7, 125, 210, 19, 209, 221, 193, 160, 2,
//       48, 47, 139, 241, 218, 248, 163, 226, 76, 115, 27, 20, 4, 251, 42, 111, 182, 192, 101, 162, 103, 39, 132, 61, 169, 205, 73, 97, 196, 190,
//       71, 171, 230, 214, 68, 13, 112, 145, 124, 1, 174, 215, 123, 67, 106, 96, 246, 240, 40, 33, 109, 15, 149, 235, 59, 63, 216, 11, 204, 9,
//       120, 166,
//       //};
//
//       //test1
////const int CCSK256[256]={
//130, 74, 86, 16, 234, 15, 87, 80, 118, 223, 97, 27, 198, 9, 153, 103, 62, 60, 2, 210, 196, 197, 184, 66, 139, 26, 53, 145, 151, 208, 119, 96,
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

    const int CCSK256[512]={56, 107, 10, 45, 108, 177, 227, 194, 22, 98, 236, 36, 66, 99, 234, 82, 62, 237,
     161, 198, 83, 252, 119, 239, 104, 151, 6, 202, 64, 157, 144, 78, 43, 253, 152, 84, 243, 168, 229, 129, 77, 224, 213, 187, 254, 29,
     175, 75, 154, 28, 170, 135, 199, 55, 225, 127, 32, 128, 5, 102, 65, 126, 114, 74, 35, 54, 130, 200, 46, 146, 211, 49, 110, 172, 94,
     191, 121, 72, 141, 38, 14, 118, 91, 87, 183, 52, 134, 116, 180, 142, 18, 31, 188, 81, 51, 122, 53, 186, 164, 23, 255, 185, 206, 79, 217,
      181, 133, 25, 250, 16, 233, 136, 228, 249, 21, 92, 232, 12, 189, 167, 89, 50, 100, 140, 238, 148, 70, 245, 37, 223, 8, 179, 85, 44, 219,
       113, 220, 201, 58, 30, 93, 176, 57, 158, 155, 88, 17, 203, 212, 208, 147, 24, 150, 131, 138, 90, 159, 197, 0, 153, 247, 244, 207, 156,
       69, 231, 178, 117, 105, 80, 41, 34, 137, 195, 26, 242, 86, 60, 184, 222, 173, 95, 165, 143, 3, 7, 125, 210, 19, 209, 221, 193, 160, 2,
       48, 47, 139, 241, 218, 248, 163, 226, 76, 115, 27, 20, 4, 251, 42, 111, 182, 192, 101, 162, 103, 39, 132, 61, 169, 205, 73, 97, 196, 190,
       71, 171, 230, 214, 68, 13, 112, 145, 124, 1, 174, 215, 123, 67, 106, 96, 246, 240, 40, 33, 109, 15, 149, 235, 59, 63, 216, 11, 204, 9,
       120, 166,
       130, 74, 86, 16, 234, 15, 87, 80, 118, 223, 97, 27, 198, 9, 153, 103, 62, 60, 2, 210, 196, 197, 184, 66, 139, 26, 53, 145, 151, 208, 119, 96,
75, 220, 163, 147, 10, 91, 61, 207, 149, 48, 195, 152, 22, 177, 105, 238, 182, 143, 20, 95, 63, 248, 225, 18, 132, 52, 222, 90, 172, 133, 70,
69, 237, 137, 40, 157, 47, 64, 44, 111, 174, 201, 39, 49, 219, 169, 235, 55, 199, 121, 251, 71, 215, 156, 13, 54, 144, 154, 200, 230, 94, 247,
1, 242, 30, 244, 32, 56, 236, 193, 167, 240, 214, 46, 155, 136, 146, 178, 6, 253, 212, 218, 99, 173, 25, 164, 43, 129, 188, 7, 59, 190, 24, 38,
11, 114, 203, 14, 221, 35, 171, 131, 122, 73, 148, 37, 36, 192, 51, 93, 45, 128, 42, 216, 209, 176, 175, 34, 79, 88, 135, 78, 202, 12, 227, 141,
 231, 134, 106, 217, 0, 58, 232, 82, 72, 161, 92, 252, 102, 109, 98, 150, 185, 81, 23, 31, 160, 110, 124, 3, 65, 83, 213, 50, 246, 243, 186, 159,
  100, 179, 205, 194, 140, 123, 255, 249, 245, 68, 187, 206, 162, 107, 181, 76, 17, 165, 104, 112, 166, 33, 254, 57, 241, 8, 120, 158, 170, 168,
   115, 126, 183, 101, 127, 85, 28, 125, 117, 224, 41, 204, 67, 189, 138, 116, 226, 5, 250, 142, 233, 191, 108, 84, 180, 77, 113, 21, 4, 239, 89,
    228, 229, 211, 29, 19

       };





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

                     for(k=0; k<code->GF; k++)
                    {
                        TMP[k] =0;
                    }

        if (n<N)
        {


                    for (q=0; q<csk->PNsize; q++)
                    {
                        for (i=0; i<2; i++)
                        {
                            u=My_drand48(init_rand);
                            v=My_drand48(init_rand);
                            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                            NoisyBin[q][i] = modulation[CCSK256[(som+q)%512]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
                        }
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
                            TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK256[(som+g)%512]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK256[(som+g)%512]][1])/(2.0*SQR(sigma));
                            //printf("%d %f \n",k, TMP[k]);
                        }
                    }
        }
        else
        {

                    for (q=0; q<csk->PNsize; q++)
                    {
                        for (i=0; i<2; i++)
                        {
                            u=My_drand48(init_rand);
                            v=My_drand48(init_rand);
                            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                            NoisyBin[q][i] = modulation[CCSK256[(som+q)%512]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
                        }
                    }

                    for (g=0; g<csk->PNsize-1; g++)
                    {
                        for(k=0; k<code->GF; k++)
                        {
                            som=0;

                            for (q=0; q<8; q++)
                            {
                                som = som + BinGF_256[k][q]*pow(2,q);
                            }
                            TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK256[(som+g)%512]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK256[(som+g)%512]][1])/(2.0*SQR(sigma));
                            //printf("%d %f \n",k, TMP[k]);
                        }
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

// NG(256) is mapped on 64NB-CCSK
void ModelChannel_AWGN_64APSK_CSK256(csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
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

    float modulation[64][2];

    const int CCSK64[256]={

            2,57,33,29,26,4,58,10,27,43,37,7,19,20,0,63,6,51,8,36,9,62,53,15,18,23,5,34,13,55,45,16,46,35,48,21,61,22,40,12,56,47,24,1,39,
    41,42,44,38,59,32,49,30,17,50,52,11,54,14,3,60,31,25,28,

    0,21,61,48,2,8,57,49,19,14,62,53,16,6,20,45,22,24,26,54,28,30,1,5,59,10,60,43,55,32,11,46,18,34,52,36,9,38,12,58,40,
    23,25,50,27,42,13,51,29,31,56,4,33,17,44,3,63,35,37,39,41,15,47,7,

32, 40, 22, 34, 35, 6, 55, 3, 16, 11, 54, 30, 45, 60, 62, 51, 33, 7, 38, 58, 42, 28, 17, 41, 47, 14, 46, 56, 63, 8, 59, 5,
 48, 53, 29, 21, 25, 52, 37, 64, 31, 49, 27, 61, 50, 26, 43, 19, 44, 15, 1, 36, 23, 2, 4, 18, 24, 39, 13, 9, 20, 57, 10, 12,

59, 35, 28, 27, 55, 48, 63, 57, 32, 4, 5, 51, 37, 41, 13, 49, 10, 14, 8, 6, 43, 22, 39, 36, 12, 17, 25, 64, 56, 47, 34, 16,
40, 29, 53, 3, 20, 26, 33, 19, 42, 15, 44, 45, 46, 24, 23, 60, 30, 38, 9, 61, 52, 18, 7, 62, 1, 50, 21, 11, 31, 2, 58, 54

    };



    ////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < 64 ; i++)
    {
        //norm_factor = table_64QAM[i][0]*table_64QAM[i][0] + table_64QAM[i][1]*table_64QAM[i][1]+norm_factor;//compute sum
        norm_factor = table_64APSK[i][0]*table_64APSK[i][0] + table_64APSK[i][1]*table_64APSK[i][1]+norm_factor;//compute sum
    }
    norm_factor = sqrt( 64 / norm_factor);
// printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< 64 ; i++)
    {
        modulation[i][0]=norm_factor*table_64APSK[i][0];
        modulation[i][1]=norm_factor*table_64APSK[i][1];
    }


    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {
        som=0;

        for (q=0; q<8; q++)
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
                NoisyBin[q][i] = modulation[CCSK64[(som+q)%256]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }

         for(k=0; k<code->GF; k++)
        {
            TMP[k] =0;
        }


        if (n<N) // per symbol puncturing
        {

            for(k=0; k<code->GF; k++)
            {
                for (g=0; g<csk->PNsize; g++)
                {

                    som=0;

                    for (q=0; q<8; q++)
                    {
                        som = som + BinGF_256[k][q]*pow(2,q);
                    }
                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%256]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%256]][1])/(2.0*SQR(sigma));
                    //printf("%d %f \n",k, TMP[k]);
                }
            }
        }
//        else
//        {
//
//            for(k=0; k<code->GF; k++)
//            {
//                for (g=0; g<csk->PNsize-1; g++)
//                {
//
//                    som=0;
//
//                    for (q=0; q<6; q++)
//                    {
//                        som = som + BinGF_64[k][q]*pow(2,q);
//                    }
//                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[som+g]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[som+g]][1])/(2.0*SQR(sigma));
//                    //printf("%d %f \n",k, TMP[k]);
//                }
//            }
//        }

        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }
    }

    for (q=0; q<csk->PNsize; q++) free(NoisyBin[q]);
    free(NoisyBin);
}



