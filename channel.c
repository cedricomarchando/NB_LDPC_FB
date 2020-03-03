/*!
 * \file channel.c
 * \brief AWGN and Rayleigh channel
 * \author C.Marchand, A. Al-Ghouwahel, Oussama Abassi, L. Conde-Canencia, A. abdmoulah, E. Boutillon
 * \copyright BSD copyright
 * \date 03/03/2015


 */

#include <math.h>
#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"
#include "./include/channel.h"



#define APSK
//#define QAM
//#define QAM_R //rotated QAM

//#define rayleigh_fading
//#define rayleigh_fading_SSD
//#define erasure

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
void ModelChannel_AWGN_BPSK (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN,int *init_rand)
{
    const int N = code->N;
    const int logGF = code->logGF;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[4096];

    float **NoisyBin = calloc(N,sizeof(float *));
    for (n=0; n<N; n++) NoisyBin[n] = calloc(logGF,sizeof(float));

    /* Binary-input AWGN channel : */

    //sigma = sqrt(1.0/(2.0*code->rate*pow(10,EbN/10.0)));  // considering EbNo
    sigma = sqrt(1.0/(pow(10,EbN/10.0))); // considering SNR
    for (n=0; n<N; n++)
    {
        for (q=0; q<logGF; q++)
        {

            u=My_drand48(init_rand);
            while(u==0.0)
			{
                u=My_drand48(init_rand);
			}
            v=My_drand48(init_rand);
            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
            NoisyBin[n][q] = BPSK(NBIN[n][q]) + sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v);
        }

    }


    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n=0; n<N; n++)
    {
        for (g=0; g<code->GF; g++)
        {
            TMP[g]=0.0;
            for (q=0; q<logGF; q++)
            {
                //TMP[g] = TMP[g] + SQR(NoisyBin[n][q]-BPSK(table->BINGF[g][q]))/(2.0*SQR(sigma));
                TMP[g] = TMP[g] - NoisyBin[n][q]*BPSK(table->BINGF[g][q]);
            }

        }

        for(k=0; k<code->GF; k++)
        {
//            decoder->intrinsic_LLR[n][k] = +1e5;
//            decoder->intrinsic_GF[n][k] = -1;
//            for (g=0; g<code->GF; g++)
//            {
//                if (TMP[g] < decoder->intrinsic_LLR[n][k])
//                {
//                    decoder->intrinsic_LLR[n][k] = TMP[g];
//                    decoder->intrinsic_GF[n][k] = g;
//                }
//            }
//            TMP[decoder->intrinsic_GF[n][k]] = +1e5;


            decoder->intrinsic_LLR[n][k] = TMP[k];
        }


    }

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
void ModelChannel_AWGN_64(code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;
    float **NoisyBin = calloc(N,sizeof(float *));
    for (n=0; n<N; n++) NoisyBin[n] = calloc(2,sizeof(float));
    /* Binary-input AWGN channel : */

    int i;

    float modulation[code->GF][2];




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
//    modulation[i][0]=norm_factor*table_64QAM[i][0];
//    modulation[i][1]=norm_factor*table_64QAM[i][1];

        modulation[i][0]=norm_factor*table_64APSK[i][0];
        modulation[i][1]=norm_factor*table_64APSK[i][1];
    }



#ifdef rayleigh_fading
    float attenuation[N];
    float attenuation_temp;
    float rand_temp;
#endif // rayleigh_fading

    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {

#ifdef  rayleigh_fading

        rand_temp=My_drand48(init_rand);
        attenuation_temp=sqrt(-log(rand_temp));
        attenuation[n]=attenuation_temp;
#endif // rayleigh_fading


        som=0;

        for (q=0; q<6; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        //printf("\n %d \n",som );

        for (q=0; q<2; q++)
        {
            u=My_drand48(init_rand);
            v=My_drand48(init_rand);
            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
#ifdef rayleigh_fading
            NoisyBin[n][q] = modulation[som][q]*attenuation_temp+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
#else
            NoisyBin[n][q] = modulation[som][q]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
#endif
        }
        //printf("%f %f %f %f \n",NoisyBin[n][0], modulation[som][0] , NoisyBin[n][1] , modulation[som][1]);getchar();
    }

    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n=0; n<N; n++)
    {

        //printf("%d \n",code->GF);getchar();
        for(k=0; k<code->GF; k++)
        {
            som=0;
            for (q=0; q<6; q++)
            {
                som = som + BinGF_64[k][q]*pow(2,q);
            }

#ifdef rayleigh_fading
            if(attenuation[n]>0 )
                TMP[k] = SQR(NoisyBin[n][0]-modulation[som][0]*attenuation[n])/(2.0*SQR(sigma))+SQR(NoisyBin[n][1]-attenuation[n]*modulation[som][1])/(2.0*SQR(sigma));
            else
                TMP[k]=0;
#else
            TMP[k] = +SQR(NoisyBin[n][0]-modulation[som][0])/(2.0*SQR(sigma))+SQR(NoisyBin[n][1]-modulation[som][1])/(2.0*SQR(sigma));
#endif // rayleigh_fading



            //printf("%f \n", TMP[k]);
        }
        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }

    }

    for (n=0; n<N; n++) free(NoisyBin[n]);
    free(NoisyBin);
}





/*!
 * \fn ModelChannel_AWGN_256QAM(code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN, int *init_rand)
 * \brief 256QAM modulation + AWGN noise on the codeword.
 *                 The function computes the likelihoods corresponding to the noisy observations.
 * Inputs
 * 	- NBIN : Binary copy of the codeword
 * 	- EbN  : Signal to noise ratio in terms of Eb/No (dB)
 * Outputs
 *      - decoder->intrinsic_LLR
 */
void ModelChannel(code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    const int GF = code->GF;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[GF];
    int som;
    float **NoisyBin = calloc(N,sizeof(float *));
    for (n=0; n<N; n++) NoisyBin[n] = calloc(2,sizeof(float));
    /* Binary-input AWGN channel : */
    int i;

    float table_GF[GF][2];
    int pos_gf_to_bin[GF];

if (GF==256)
{
#ifdef QAM_R
// rotated QAM constellation construction
    float alpha=31.7*PI/180;


    for(i=0; i<GF; i++)
    {
        table_GF[i][0]=table_256QAM[i][0]*cos(alpha)-table_256QAM[i][1]*sin(alpha);
        table_GF[i][1]=table_256QAM[i][1]*cos(alpha)+table_256QAM[i][0]*sin(alpha);
   //       table_GF[i][0]=table_256_circleQAM[i][0]*cos(alpha)-table_256_circleQAM[i][1]*sin(alpha);
   //     table_GF[i][1]=table_256_circleQAM[i][1]*cos(alpha)+table_256_circleQAM[i][0]*sin(alpha);
    }


#endif
#ifdef APSK
// APSK constellation construction
//from DVB-S2x code identifier 135/180
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
#endif
#ifdef QAM
    for (i=0; i<256; i++)
    {
        table_GF[i][0]=table_256QAM[i][0];
        table_GF[i][1]=table_256QAM[i][1];
//        table_GF[i][0]=table_256_circleQAM[i][0];
//        table_GF[i][1]=table_256_circleQAM[i][1];
        pos_gf_to_bin[i] = pos_gf256_to_bin[i];


//printf("%d \t %f \t %f \n",i,table_256APSK[i][0],table_256APSK[i][1] );
    }
#endif


    for (i=0; i<256; i++)
    {
        pos_gf_to_bin[i] = pos_gf256_to_bin[i];
    }


}
//getchar();


if (GF==64)
{
#ifdef QAM_R
// rotated QAM constellation construction
    float alpha=31.7*PI/180;


    for(i=0; i<GF; i++)
    {
        table_GF[i][0]=table_64QAM[i][0]*cos(alpha)-table_64QAM[i][1]*sin(alpha);
        table_GF[i][1]=table_64QAM[i][1]*cos(alpha)+table_64QAM[i][0]*sin(alpha);
    }


#endif
#ifdef APSK
// APSK constellation construction



   // from DVB S2X 8+16+20+20APSK constellation
float table_APSK[64][2]=
{
  {  2.2*cos(PI*25/16), 2.2*sin(PI*25/16)},
{2.2*cos(PI*23/16), 2.2*sin(PI*23/16)},
{2.2*cos(PI*7/16), 2.2*sin(PI*7/16)},
{2.2*cos(PI*9/16), 2.2*sin(PI*9/16)},
{5.2*cos(PI*7/4 ), 5.2*sin(PI*7/4 )},
{5.2*cos(PI*5/4 ), 5.2*sin(PI*5/4 )},
{5.2*cos(PI*1/4 ), 5.2*sin(PI*1/4 )},
{5.2*cos(PI*3/4 ), 5.2*sin(PI*3/4 )},
{2.2*cos(PI*27/16 ),2.2*sin(PI*27/16 )},
{2.2*cos(PI*21/16 ), 2.2*sin(PI*21/16 )},
{2.2*cos(PI*5/16 ), 2.2*sin(PI*5/16 )},
{2.2*cos(PI*11/16 ), 2.2*sin(PI*11/16 )},
{3.6*cos(PI*7/4 ) , 3.6*sin(PI*7/4 )},
{3.6*cos(PI*5/4 ), 3.6*sin(PI*5/4 )},
{3.6*cos(PI*1/4 ), 3.6*sin(PI*1/4 )},
{3.6*cos(PI*3/4 ), 3.6*sin(PI*3/4 )},
{5.2*cos(PI*31/20 ), 5.2*sin(PI*31/20 )},
{5.2*cos(PI*29/20 ),5.2*sin(PI*29/20 )},
{5.2*cos(PI*9/20 ), 5.2*sin(PI*9/20 )},
{5.2*cos(PI*11/20 ), 5.2*sin(PI*11/20 )},
{5.2*cos(PI*33/20 ),5.2*sin(PI*33/20 )},
{5.2*cos(PI*27/20 ), 5.2*sin(PI*27/20 )},
{5.2*cos(PI*7/20 ),5.2*sin(PI*7/20 )},
{5.2*cos(PI*13/20 ), 5.2*sin(PI*13/20 )},
{3.6*cos(PI*31/20 ),3.6*sin(PI*31/20 )},
{3.6*cos(PI*29/20 ),3.6*sin(PI*29/20 )},
{3.6*cos(PI*9/20 ), 3.6*sin(PI*9/20 )},
{3.6*cos(PI*11/20 ), 3.6*sin(PI*11/20 )},
{3.6*cos(PI*33/20 ),  3.6*sin(PI*33/20 )},
{3.6*cos(PI*27/20 ),3.6*sin(PI*27/20 )},
{3.6*cos(PI*7/20 ),3.6*sin(PI*7/20 )},
{3.6*cos(PI*13/20 ), 3.6*sin(PI*13/20 )},
{cos(PI*13/8 ),sin(PI*13/8 )},
{cos(PI*11/8 ), sin(PI*11/8 )},
{cos(PI*3/8 ),sin(PI*3/8 )},
{cos(PI*5/8 ),sin(PI*5/8 )},
{5.2*cos(PI*37/20 ),5.2*sin(PI*37/20 )},
{5.2*cos(PI*23/20 ),  5.2*sin(PI*23/20 )},
{5.2*cos(PI*3/20 ) ,5.2*sin(PI*3/20 )},
{5.2*cos(PI*17/20 ),5.2*sin(PI*17/20 )},
{2.2*cos(PI*29/16 ), 2.2*sin(PI*29/16 )},
{2.2*cos(PI*19/16 ), 2.2*sin(PI*19/16 )},
{2.2*cos(PI*3/16 ), 2.2*sin(PI*3/16 )},
{2.2*cos(PI*13/16 ), 2.2*sin(PI*13/16 )},
{3.6*cos(PI*37/20 ), 3.6*sin(PI*37/20 )},
{3.6*cos(PI*23/20 ),  3.6*sin(PI*23/20 )},
{3.6*cos(PI*3/20 ), 3.6*sin(PI*3/20 )},
{3.6*cos(PI*17/20 ),3.6*sin(PI*17/20 )},
{cos(PI*15/8 ),sin(PI*15/8 )},
{cos(PI*9/8 ), sin(PI*9/8 )},
{cos(PI*1/8 ),sin(PI*1/8 )},
{cos(PI*7/8 ),sin(PI*7/8 )},
{5.2*cos(PI*39/20 ),5.2*sin(PI*39/20 )},
{5.2*cos(PI*21/20 ),5.2*sin(PI*21/20 )},
{5.2*cos(PI*1/20 ),5.2*sin(PI*1/20 )},
{5.2*cos(PI*19/20 ),5.2*sin(PI*19/20 )},
{2.2*cos(PI*31/16 ),2.2*sin(PI*31/16 )},
{2.2*cos(PI*17/16 ),2.2*sin(PI*17/16 )},
{2.2*cos(PI*1/16 ),2.2*sin(PI*1/16 )},
{2.2*cos(PI*15/16 ),2.2*sin(PI*15/16 )},
{3.6*cos(PI*39/20 ),3.6*sin(PI*39/20 )},
{3.6*cos(PI*21/20 ), 3.6*sin(PI*21/20 )},
{3.6*cos(PI*1/20 ),3.6*sin(PI*1/20 )},
{3.6*cos(PI*19/20 ),3.6*sin(PI*19/20 )}
};



        for (i=0; i<GF; i++)
        {
            table_GF[i][0]= table_APSK[i][0];
            table_GF[i][1]= table_APSK[i][1];
        }



#endif
#ifdef QAM
    for (i=0; i<GF; i++)
    {
        table_GF[i][0]=table_64QAM[i][0];
        table_GF[i][1]=table_64QAM[i][1];

//printf("%d \t %f \t %f \n",i,table_GF[i][0],table_GF[i][1] );
    }
#endif
    for (i=0; i<GF; i++)
    {
        pos_gf_to_bin[i] = pos_gf64_to_bin[i];
    }

}
//getchar();




if (GF==16)
{
#ifdef QAM_R
// rotated QAM constellation construction
    float alpha=31.7*PI/180;


    for(i=0; i<GF; i++)
    {
        table_GF[i][0]=table_16QAM[i][0]*cos(alpha)-table_16QAM[i][1]*sin(alpha);
        table_GF[i][1]=table_16QAM[i][1]*cos(alpha)+table_16QAM[i][0]*sin(alpha);
   //       table_GF[i][0]=table_256_circleQAM[i][0]*cos(alpha)-table_256_circleQAM[i][1]*sin(alpha);
   //     table_GF[i][1]=table_256_circleQAM[i][1]*cos(alpha)+table_256_circleQAM[i][0]*sin(alpha);
    }


#endif
#ifdef QAM
    for (i=0; i<GF; i++)
    {
        table_GF[i][0]=table_16QAM[i][0];
        table_GF[i][1]=table_16QAM[i][1];
//        table_GF[i][0]=table_16_circleQAM[i][0];
//        table_GF[i][1]=table_16_circleQAM[i][1];


//printf("%d \t %f \t %f \n",i,table_16APSK[i][0],table_16APSK[i][1] );
    }


#endif

    for (i=0; i<GF; i++)
    {
            pos_gf_to_bin[i] = pos_gf16_to_bin[i];
    }
}






////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < GF ; i++)
    {
        norm_factor = table_GF[i][0]*table_GF[i][0] + table_GF[i][1]*table_GF[i][1]+norm_factor;//compute sum
    }
    norm_factor = sqrt( code->GF / norm_factor);
//printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< GF ; i++)
    {
        table_GF[i][0]=norm_factor*table_GF[i][0];
        table_GF[i][1]=norm_factor*table_GF[i][1];
    }




float attenuation_temp;
#ifdef rayleigh_fading // rayleigh_fading
    float attenuation[N];
    float rand_temp;
#endif
#ifdef rayleigh_fading_SSD // Signal Space Diversity
    float attenuation[N][2];
    float rand_temp;
#endif // rayleigh_fading
#ifdef erasure
float erasure_proba=0.1;
#endif // erasure



    sigma = sqrt(1.0/(2*pow(10,EbN/10.0)));  //SNR

    //sigma = sqrt(1.0/(2*code->rate * code->logGF * pow(10,EbN/10.0)));  //Eb/No

//opfile=fopen("./data/bubble_stat.txt","a");
//    for (k=0; k <code.rowDegree[0]*3; k++ )
//    {
//
//    printf( "bubble check number: %d \n",k );
//
//    for (i=0; i<decoder.nbMax; i++)
//    {
//        for (l=0; l<decoder.nbMax; l++)
//        {
//printf(" %d ", stat_bubble[i + decoder.nbMax * l + k * decoder.nbMax * decoder.nbMax ]);
//fprintf(opfile," %d ", stat_bubble[i + decoder.nbMax * l + k * decoder.nbMax * decoder.nbMax ]);
//        }
//        printf(" \n ");
//        fprintf(opfile," \n ");
//    }


     //printf("sigma = %f \n",sigma);getchar();
attenuation_temp=1.0;

    for (n=0; n<N; n++)
    {
        som=0;

#ifdef  rayleigh_fading
        rand_temp=My_drand48(init_rand);
        attenuation_temp=sqrt(-log(rand_temp));
        attenuation[n]=attenuation_temp;
#endif // rayleigh_fading


        for (q=0; q< code->logGF; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
            //printf(" %d",NBIN[n][q] );
        }

        //printf("\n somme= %d \n",som );
        //getchar();

        for (i=0; i<2; i++)
        {
            u=My_drand48(init_rand);
            v=My_drand48(init_rand);

#ifdef rayleigh_fading_SSD
            rand_temp=My_drand48(init_rand);
            attenuation_temp=sqrt(-log(rand_temp));
attenuation[n][i]=attenuation_temp;
#ifdef erasure

rand_temp = My_drand48(init_rand);

//printf("rand: %f , attenuation1: %f",rand_temp,attenuation_temp);
if (rand_temp< erasure_proba)
{
    attenuation_temp=0.0;
}
else
{
attenuation_temp=attenuation_temp/sqrt(1-erasure_proba);

}

attenuation[n][i]=attenuation_temp;

#endif // erasure


#endif

NoisyBin[n][i] = table_GF[som][i]*attenuation_temp+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;



        }
        //printf("%f %f %f %f \n",NoisyBin[n][0], table_GF[som][0] , NoisyBin[n][1] , table_GF[som][1]);getchar();
    }

    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n=0; n<N; n++)
    {


        //printf("%d \n",code->GF);getchar();
        for(k=0; k<GF; k++)
        {

som=pos_gf_to_bin[k];




#ifdef rayleigh_fading
            if(attenuation[n]>0 )
                TMP[k] = SQR(NoisyBin[n][0]-table_GF[som][0]*attenuation[n])/(2.0*SQR(sigma))+SQR(NoisyBin[n][1]-attenuation[n]*table_GF[som][1])/(2.0*SQR(sigma));
            else
                TMP[k]=0;
#else
#ifdef rayleigh_fading_SSD
            TMP[k]=0;
            if(attenuation[n][0]>0 )
            {
                TMP[k] =         SQR(NoisyBin[n][0]-attenuation[n][0]*table_GF[som][0])/(2.0*SQR(sigma));
            }

            if(attenuation[n][1]>0 )
            {
                TMP[k] = TMP[k] +SQR(NoisyBin[n][1]-attenuation[n][1]*table_GF[som][1])/(2.0*SQR(sigma));
            }


#else
            TMP[k] = +SQR(NoisyBin[n][0]-table_GF[som][0])/(2.0*SQR(sigma))+SQR(NoisyBin[n][1]-table_GF[som][1])/(2.0*SQR(sigma));
#endif
#endif // rayleigh_fading

            //printf("%f \n", TMP[k]);
        }
        //getchar();

        for(k=0; k<code->GF; k++)
        {
        decoder->intrinsic_LLR[n][k] = TMP[k];
        }

    }

    for (n=0; n<N; n++) free(NoisyBin[n]);
    free(NoisyBin);
}


void ModelChannel_AWGN_256QAM_4D(code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[256];
    int som;
    float **NoisyBin = calloc(N,sizeof(float *));
    for (n=0; n<N; n++) NoisyBin[n] = calloc(4,sizeof(float));
    /* Binary-input AWGN channel : */
    int i;

    float table_256[256][4];




// Signal Space Diversity
    float attenuation[N][4];
   float attenuation_temp;
    float rand_temp;
#ifdef erasure
float erasure_proba=0.1;
#endif // erasure


////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < code->GF ; i++)
    {
//        table_256[i][0]=table_256QAM_4D_16QAM[i][0];
//        table_256[i][1]=table_256QAM_4D_16QAM[i][1];
//        table_256[i][2]=table_256QAM_4D_16QAM[i][2];
//        table_256[i][3]=table_256QAM_4D_16QAM[i][3];

        table_256[i][0]=table_256QAM_4D_16QAM_R[i][0];
        table_256[i][1]=table_256QAM_4D_16QAM_R[i][1];
        table_256[i][2]=table_256QAM_4D_16QAM_R[i][2];
        table_256[i][3]=table_256QAM_4D_16QAM_R[i][3];


//        table_256[i][0]=table_256QAM_4D[i][0];
//        table_256[i][1]=table_256QAM_4D[i][1];
//        table_256[i][2]=table_256QAM_4D[i][2];
//        table_256[i][3]=table_256QAM_4D[i][3];

    }

    for(i=0; i < code->GF ; i++)
    {
        norm_factor = table_256[i][0]*table_256[i][0] + table_256[i][1]*table_256[i][1]+norm_factor;//compute sum
        norm_factor = table_256[i][2]*table_256[i][2] + table_256[i][3]*table_256[i][3]+norm_factor;//compute sum
    }
    norm_factor = sqrt(2* code->GF / norm_factor);
//printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< code->GF ; i++)
    {
        table_256[i][0]=norm_factor*table_256[i][0];
        table_256[i][1]=norm_factor*table_256[i][1];
        table_256[i][2]=norm_factor*table_256[i][2];
        table_256[i][3]=norm_factor*table_256[i][3];
    }

    sigma = sqrt(1.0/(2*pow(10,EbN/10.0)));  //SNR

    //sigma = sqrt(1.0/(2*code->rate * code->logGF * pow(10,EbN/10.0)));  //Eb/No

    //printf("sigma = %f \n",sigma);getchar();


    for (n=0; n<N; n++)
    {
        som=0;


        for (q=0; q<8; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        //printf("\n %d \n",som );

        for (q=0; q<4; q++)
        {
            u=My_drand48(init_rand);
            v=My_drand48(init_rand);


            rand_temp=My_drand48(init_rand);
            attenuation_temp=sqrt(-log(rand_temp));
attenuation[n][q]=attenuation_temp;
NoisyBin[n][q] = table_256[som][q]*attenuation_temp+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;



#ifdef erasure

rand_temp = My_drand48(init_rand);

//printf("rand: %f , attenuation1: %f",rand_temp,attenuation_temp);
if (rand_temp< erasure_proba)
{
    attenuation_temp=0.0;
}
else
{
attenuation_temp=attenuation_temp/sqrt(1-erasure_proba);
}

attenuation[n][q]=attenuation_temp;

#endif // erasure





        }
        //printf("%f %f %f %f \n",NoisyBin[n][0], table_256[som][0] , NoisyBin[n][1] , table_256[som][1]);getchar();
    }

    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n=0; n<N; n++)
    {


        //printf("%d \n",code->GF);getchar();
        for(k=0; k<code->GF; k++)
        {
            som=0;
            for (q=0; q<8; q++)
            {
                som = som + BinGF_256[k][q]*pow(2,q);
            }


            TMP[k]=0;
            if(attenuation[n][0]>0 )
            {
                TMP[k] =         SQR(NoisyBin[n][0]-attenuation[n][0]*table_256[som][0])/(2.0*SQR(sigma));
            }

            if(attenuation[n][1]>0 )
            {
                TMP[k] = TMP[k] +SQR(NoisyBin[n][1]-attenuation[n][1]*table_256[som][1])/(2.0*SQR(sigma));
            }

            if(attenuation[n][2]>0 )
            {
                TMP[k] = TMP[k] +SQR(NoisyBin[n][2]-attenuation[n][2]*table_256[som][2])/(2.0*SQR(sigma));
            }

            if(attenuation[n][3]>0 )
            {
                TMP[k] = TMP[k] +SQR(NoisyBin[n][3]-attenuation[n][3]*table_256[som][3])/(2.0*SQR(sigma));
            }


            //printf("%f \n", TMP[k]);
        }
        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = +1e5;
            decoder->intrinsic_GF[n][k] = -1;
            for (g=0; g<code->GF; g++)
            {
                if (TMP[g] < decoder->intrinsic_LLR[n][k])
                {
                    decoder->intrinsic_LLR[n][k] = TMP[g];
                    decoder->intrinsic_GF[n][k] = g;
                }
            }
            TMP[decoder->intrinsic_GF[n][k]] = +1e5;
        }

    }

    for (n=0; n<N; n++) free(NoisyBin[n]);
    free(NoisyBin);
}

void build_QAM_table(int logGF)
{
    int i,j;
    for( i = - 1<<((logGF>>1)-1) ; i< 1<<((logGF>>1)-1); i++)
    {
        for( j = - 1<<((logGF>>1)-1) ; j< 1<<((logGF>>1)-1); j++)
        {

            printf( "%f \t %f \n", i+0.5,j+0.5);
        }
    }
}
