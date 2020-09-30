/*!
 * \file NB_LDPC.c
 * \brief Non-binary LDPC reduced complexity decoder with horizontal scheduling
 * \author C.Marchand, A. Al-Ghouwahel, Oussama Abassi, L. Conde-Canencia, A. abdmoulah, E. Boutillon
 * \copyright BSD copyright
 * \date 03/03/2015
 * \details

   Extended Min-Sum Decoder
  Horizontal Scheduling
    Layered decoder
    syndrome based architecture

 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "./include/NB_LDPC.h"

/// preprocessing directives
#define CCSK // use of Code-shift keying modulation



/*!
 * \fn int main(int argc, char * argv[])
 * \brief main program

 * Inputs
 *
 *		NbMonteCarlo     : # simulated frames
 *		NbIterMax        : # of maximum decoding iteration
 *		FileMatrix       : File name of the parity-check matrix
 *		EbN              : Eb/No (dB)
 *		n_vc             : size of truncated messages from Variable to Check
 *      n_cv             : size of truncated messages from Check to Variable
 *		Offset           : offset correction factor (0.4 -- 1)
 *		NbOper           : Maximum number of operations for sorting
 * Output
 *              Frame Error Rate for the given simulation parameters
 * Input File : 'FileMatrix' is an ASCII file with the parity-check matrix description in aList format.
 * Output File : in ./data file
 *
 *
 under linux
 compile using make
 then Run the executable with the following parameters

./essai 2000 10 ./matrices/KN/N576_K480_GF64.txt 3.5 10 20 0.3 25

giving

<0> FER = 40/751 = 0.053262 BER = 520 / x = 0.001443 avr_it = 2.58

*/
int main(int argc, char * argv[])
{
    int 		k,l,n,iter,i,g;
    int 		**KBIN,*KSYMB,**NBIN,*NSYMB;
    int 		*decide,*CodeWord;
    int 		nb,NbIterMax;
    float 		EbN;
    int 		NbOper,NbMonteCarlo;
    float 	offset;
    char *FileName,*FileMatrix,*name;
    int synd=0, nbErrors, nbErroneousFrames = 0, nbUndetectedErrors = 0;
    int total_errors =0;

    code_t code;
    table_t table;
    decoder_t decoder;

    int node;

    int Idum=-1; // initialization of random generator
    srand(5);

    /*
     * Command line arguments
     */
    if (argc < 8)
    {
        printf("File:\n %s\n ",argv[0]);
        printf(usage_txt);
        return (EXIT_FAILURE);
    }
    FileName 	= malloc(STR_MAXSIZE);
    FileMatrix 	= malloc(STR_MAXSIZE);
    name 		= malloc(STR_MAXSIZE);


    NbMonteCarlo 		= atoi(argv[1]);
    NbIterMax 		= atoi(argv[2]);
    strcpy(FileMatrix,argv[3]);
    EbN 			= atof(argv[4]);
    decoder.n_vc 	= atoi(argv[5]);
    decoder.n_cv 	= atoi(argv[6]);
    offset  		= atof(argv[7]);
    NbOper  		= atoi(argv[8]);
    printf(" Monte-Carlo simulation of Non-Binary LDPC decoder \n\n");
    printf("Simulation parameters:\n");
    printf("\n\t NbMonteCarlo     : %d", NbMonteCarlo);
    printf("\n\t NbIterMax        : %d", NbIterMax);
    printf("\n\t FileMatrix       : %s", FileMatrix);
    printf("\n\t Eb/No (dB)       : %g", EbN);
    printf("\n\t n_vc            : %d", decoder.n_vc);
    printf("\n\t n_cv            : %d", decoder.n_cv);
    printf("\n\t Offset           : %g", offset);
    printf("\n\t NbOper           : %d\n",NbOper);

    printf("Load code  ... ");
    LoadCode (FileMatrix, &code);
    printf(" OK \n Load table ...");
    LoadTables (&table, code.GF, code.logGF);
    printf("OK \n Allocate decoder ... ");
    AllocateDecoder (&code, &decoder);
    printf("OK \n Gaussian Elimination ... ");
    GaussianElimination (&code, &table);
    printf(" OK \n");


// output results in a file
    FILE *opfile;
    char note[40]="K120GF64_CCSK";
    printf("\n\t Note             : %s\n",note);
    char file_name [70];
    time_t start_time;
    time_t end_time;
    double exec_time;
    char* c_time_string;


    #ifdef CCSK
    // CCSK: build CCSK table
    csk_t csk;
    //csk.PNsize =code.GF;
    csk.PNsize = 64;
    printf("\n\t PN is generated using an LFSR \n");
    allocate_csk(&csk, csk.PNsize);
    PNGenerator( &csk ); //generate a PN sequence for csk modulation
    build_natural_csk_mapping(code.GF, &csk, table.BINGF ); //fills the csk_arr with shifted versions of PN sequence
    //build_equidistant_csk_mapping(code.GF, &csk, table.BINGF);//optimized for 1024
    //build_punctured_csk_mapping(code.GF,code.logGF, &csk, table.BINGF);
    //build_CSK_map(&code, &csk); //construction of a mapping without PN sequence









    float chu_real[csk.PNsize];
    float chu_imag[csk.PNsize];
    CHU_Generator(chu_real,chu_imag,csk.PNsize);
    //CHU_AM_Generator(chu_real,chu_imag,csk.PNsize);
    //CHU_Generator_64apsk(chu_real,chu_imag,csk.PNsize);
    //CHU_Generator_256apsk(chu_real,chu_imag,csk.PNsize);





    csk.PNsize = 64;  // for "short" CCSK mapping

    #endif


    sprintf (file_name,"./data/results_N%d_CR%0.2f_GF%d_IT%d_Offset%0.1f_nm%d_%s.txt",code.N,code.rate,code.GF,NbIterMax, offset,decoder.n_cv,note);


    start_time = time(NULL);
    c_time_string = ctime(&start_time);
    printf("Simulation started at time: %s \n", c_time_string);


    /*
     * Memory  allocation
     */
    NBIN=(int **)calloc(code.N,sizeof(int *));
    for (n=0; n<code.N; n++)  	NBIN[n]=(int *)calloc(code.logGF,sizeof(int));
    KBIN=(int **)calloc(code.K,sizeof(int *));
    for (k=0; k<code.K; k++) 	KBIN[k]=(int *)calloc(code.logGF,sizeof(int));

    NSYMB=(int *)calloc(code.N,sizeof(int));
    KSYMB=(int *)calloc(code.K,sizeof(int));
    CodeWord=(int *)calloc(code.N,sizeof(int));

    decide=(int *)calloc(code.N,sizeof(int));



    // check that dc is constant
    int dc_min=100;
    int dc_max=0;
    for (node=0; node<code.M; node++) /* Loop for the M Check nodes */
    {
        if (dc_max < code.rowDegree[node])
        {
            dc_max = code.rowDegree[node];
        }
        if (dc_min > code.rowDegree[node])
        {
            dc_min = code.rowDegree[node];
        }
    }
    if (dc_min != dc_max)
    {
        printf("d_c is not constant: dc_min= %d ; dc_max=%d !!!!!! \n", dc_min,dc_max);
    }


    int sum_it;


//    getchar();

    softdata_t Mvc_temp[dc_max][code.GF];
    softdata_t Mvc_temp2[dc_max][code.GF];
    sum_it=0;




    for (nb=1; nb<=NbMonteCarlo; nb++)
    {
        /* Decoder re-initialization */

        /* Generate uniformly distributed information bits (KBIN)) */
        RandomBinaryGenerator (code.N, code.M, code.GF, code.logGF, KBIN, KSYMB, table.BINGF,&Idum);

        /* Encode the information bits KBIN to a (non binary) codeword NSYMB */
        Encoding (&code, &table, CodeWord, NBIN, KSYMB);

        /* Noisy channel (AWGN)*/
        #ifdef CCSK
        ModelChannel_AWGN_BPSK_CSK (&csk,&code, &decoder, &table, CodeWord, EbN,&Idum);
        //ModelChannel_CHU_CSK(chu_real,chu_imag, &csk,&code, &decoder, NBIN, EbN, &Idum);

        //ModelChannel_AWGN_64_CSK (&csk,&code, &decoder, NBIN, EbN,&Idum);
        //ModelChannel_AWGN_256_CSK (&csk,&code, &decoder, NBIN, EbN,&Idum);
        //ModelChannel_AWGN_64APSK_CSK256(&csk,&code, &decoder, NBIN, EbN, &Idum);
        #endif
        #ifndef CCSK
            ModelChannel_AWGN_BPSK (&code, &decoder, &table,  NBIN, EbN,&Idum);
           // ModelChannel_AWGN_BPSK_repeat (&code, &decoder, &table,  NBIN, EbN,&Idum);
        //ModelChannel_AWGN_64 (&code, &decoder, NBIN, EbN,&Idum);
        //ModelChannel(&code, &decoder,  NBIN, EbN,&Idum);
        #endif


        /***********************************************/
        /* Implementation of the horizontal scheduling */

        // init Mvc with intrinsic
        for (n=0; n<code.N; n++)
        {
            for (k=0; k<code.GF; k++)
            {
                decoder.VtoC[n][k] = decoder.intrinsic_LLR[n][k];
            }
        }



        /* Decoding iterations*/
        for (iter=0; iter < NbIterMax - 1; iter++)
        {

            for (node=0; node<code.M; node++) /* Loop for the M Check nodes */
            {
                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for (k=0; k < code.GF; k++)
                    {
                        Mvc_temp[i][k] = decoder.VtoC[code.mat[node][i]][k];
                        Mvc_temp2[i][k] = Mvc_temp[i][k];
                    }
                }


                // sorting Mvc values
                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for(k=0; k<decoder.n_vc; k++)
                    {
                        decoder.M_VtoC_LLR[i][k]=+1e5;
                        decoder.M_VtoC_GF[i][k]=0;
                        for (g=0; g<code.GF; g++)
                        {
                            if (Mvc_temp2[i][g] < decoder.M_VtoC_LLR[i][k])
                            {
                                decoder.M_VtoC_LLR[i][k]=Mvc_temp2[i][g];
                                decoder.M_VtoC_GF[i][k]=g;
                            }
                        }
                        Mvc_temp2[i][decoder.M_VtoC_GF[i][k]]=+1e5;
                    }

                    // Normalisation  */
                    for(g=1; g<decoder.n_vc; g++)
                        {
                            decoder.M_VtoC_LLR[i][g]= decoder.M_VtoC_LLR[i][g]-decoder.M_VtoC_LLR[i][0];
                            //printf(" GF:%d  LLR:%0.2f \n",decoder.M_VtoC_GF[i][g],decoder.M_VtoC_LLR[i][g] );
                        }
                        //getchar();
                    decoder.M_VtoC_LLR[i][0]=0.0;
                }



                CheckPassLogEMS (node, &decoder, &code, &table,NbOper,offset);
                //CheckPassLogEMS_dc3(node, &decoder, &code, &table,NbOper,offset);


                //compute SO
                for (i=0; i < code.rowDegree[node]; i++)
                {
                    for (k=0; k < code.GF; k++)
                    {
                        decoder.APP[code.mat[node][i]][k] = decoder.M_CtoV_LLR[i][k] + Mvc_temp[i][k];
                        decoder.VtoC[code.mat[node][i]][k]= decoder.M_CtoV_LLR[i][k] + decoder.intrinsic_LLR[code.mat[node][i]][k];// compute Mvc and save RAM
                    }
                }


                /*******************************************************************************************************************/
                /*******************************************************************************************************************/

            } /* End of the node update*/

            Decision (decide, decoder.APP, code.N, code.GF);
            synd = Syndrom (&code, decide, &table);
            if (synd == 0)
                break;
        }

        sum_it= sum_it+ iter +1;



        /* Compute the Bit Error Rate (BER)*/
        nbErrors = 0;
        for (k=0; k<code.K; k++)
        {
            for (l=0; l<code.logGF; l++)
                if (table.BINGF[decide[k]][l] != NBIN[k][l])
                    nbErrors ++ ;
        }


        total_errors = total_errors + nbErrors;
        if (nbErrors != 0)
        {
            nbErroneousFrames ++;
            if (synd == 0)
                nbUndetectedErrors ++;
        }
        if (nb%10 ==0)
        {
                    printf("\r<%d> FER= %d / %d = %f BER= %d / x = %f  avr_it=%.2f",
               nbUndetectedErrors, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb,total_errors,(double)total_errors/(nb*code.K*code.logGF),(double)(sum_it)/nb);
        fflush(stdout);
        }



        if (nbErroneousFrames == 40)
            break;



    }

    printf("\r<%d> FER= %d / %d = %f BER= %d / x = %f  avr_it=%.2f",
    nbUndetectedErrors, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb,total_errors,(double)total_errors/(nb*code.K*code.logGF),(double)(sum_it)/nb);

    printf(" \n results are printed in file %s \n",file_name);



    end_time = time(NULL);
    c_time_string = ctime(&end_time);
    exec_time = difftime(end_time,start_time);
    opfile=fopen(file_name,"a");
//opfile=fopen("./data/results.txt","a");
    if ((opfile)==NULL)
    {
        printf(" \n !! file not found \n ");
    }
    else
    {
        fprintf(opfile," SNR:%.2f: \t FER= %d / %d = %f ", EbN, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb );
        fprintf(opfile," \t BER= %d / x = \t %f  avr_it= \t %.2f \t time: %s",total_errors, (double)total_errors/(double)(nb*code.K*code.logGF),(double)(sum_it)/nb , c_time_string );
    }
    fclose(opfile);
    printf(" \n results printed \n ");


    printf("\n");
    printf("Simulation complete at time: %s", c_time_string);

    printf("execution time:%0.2f",exec_time);

    getchar();

    free(FileName);
    free(FileMatrix);
    free(name);
    free(decide);
    free(CodeWord);
    free(KSYMB);
    free(NSYMB);

    for (n=0; n<code.N; n++) free(NBIN[n]);
    free(NBIN);
    for (k=0; k<code.K; k++) free(KBIN[k]);
    free(KBIN);

    FreeCode(&code);
    FreeDecoder(&decoder);
    FreeTable(&table);
    return (EXIT_SUCCESS);
}
