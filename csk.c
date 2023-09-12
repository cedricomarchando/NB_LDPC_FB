
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
            a=1;
            for(i=0; i<csk->PNsize ; i++)
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

            for(i=0; i<csk->PNsize; i++)
            {
                //csk->PN[i] = PN64[i];
                //printf(" %d ",(csk->PN[i]-1)/-2);
            }
        //getchar();
        }


        if( (csk->PNsize == 128)||(csk->PNsize == 127))
        {
            a=1;
            for(i=0; i<csk->PNsize; i++)
            {
                //*******************************************
                //** primitive polynomial x**7+x**3+1
                //*******************************************
                csk->PN[i]=BPSK((a >> 6));//the output
                //feedback and shift
                LowBit = a >> 6; //x**7 term
                //LowBit = LowBit ^ (a >> 2); //x**3 term
                LowBit = LowBit ^ a; //x term
                a = ((a << 1) + (LowBit & 1)) & 0x07F; //return shifted 6 bit value
                //printf(" %d ",csk->PN[i] );
            }
        //getchar();
        }





        // For GF(256) with polynomial P(x)=X^8+ X^4+X^3+X^2+1
        if(csk->PNsize == 256)
        {
                    a=1;
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
                    a=1;
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
                    a=1;
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
                    a=1;
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
                    a=1;
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
            chu_real[i] = cos( i*R*(i+1)*PI /N);
            chu_imag[i] = sin( i*R*(i+1)*PI /N);
            //printf(" real:%f \t imag:%f ",chu_real[i],chu_imag[i]); getchar();
        }

        // print constellation in a file
    FILE *opfile;
    opfile=fopen("chu.txt","a");
    for(i=0; i< N ; i++)
    {
        fprintf(opfile," %f %f ;", chu_real[i],chu_imag[i]);
    }
    fclose(opfile);
    printf(" \n CHU generation: Success\n");
        getchar();
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
void build_natural_csk_mapping(int GF, csk_t *csk, int **BINGF)
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


        //check that all shifted PN are different
    int bin_temp[12];
    int GF_temp;
    int GF_check[GF];
    int logGF = log2(GF);
    int pi_permut[GF];
    for (i=0; i<GF; i++){GF_check[i]=0;}

    for (i=0; i<GF; i++)
    {

        for (j=0; j<logGF; j++)
        {
            bin_temp[j]=(csk->CSK_arr[i][j] -1)/-2;
          //  printf("%d ",bin_temp[j]);
        }
        GF_temp = Bin2GF(bin_temp,GF,logGF,BINGF);
        GF_check[GF_temp]=GF_check[GF_temp] + 1;
        //printf(" = %d \t", GF_temp);
        pi_permut[i]=GF_temp;

    }
    //for (i=0; i<GF; i++){printf("%d ",GF_check[i]);}
    //printf("\n");
    //for (i=0; i<GF; i++){printf("%d \t",pi_permut[i]);}

    //getchar();


//permut CSK_arr such that punctured p et equivalent to BPSK
     for (i=0; i<GF; i++)
    { //printf(" \n %d:",i);

        for (j=0; j<csk->PNsize; j++)
        {
            csk->CSK_arr[pi_permut[i]][j]=csk->PN[(j+i)%csk->PNsize];;
            //printf("%d ",(csk->CSK_arr[i][j]+1)/2 );

        }
        //printf(" \n ");
    }


//    // check the result
//    for (i=0; i<GF; i++)
//    { printf(" \n %d:",i);
//
//        for (j=0; j<logGF; j++)
//        {
//            printf("%d ",(csk->CSK_arr[i][j]-1)/(-2) );
//
//        }
//        printf(" \n ");
//    }
//
//
//    getchar();


    printf("Filling CSK_arr: Success\n");
    fflush(stdout);
}


void build_equidistant_csk_mapping(int GF, csk_t *csk, int **BINGF)
{
    int i, j,step;

    for (i=0; i<csk->PNsize; i++)
    {
        csk->CSK_arr[0][i]=csk->PN[i];
        //printf("%d ",csk->PN[i] );

    }

//getchar();
step=csk->PNsize/GF;

    for (i=1; i<GF; i++)
    { //printf(" \n %d:",i);

        for (j=0; j<csk->PNsize; j++)
        {
            csk->CSK_arr[i][j]=csk->PN[(j+i*step)%csk->PNsize];
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
            //printf("%d ",i); //getchar();

            PN_index = PN_index + 1;

        }


        //printf(" \n ");
    }

//
//    //check that all shifted PN are different
//    for (i=0; i<GF; i++)
//    {
//
//        for (j=0; j<logGF; j++)
//        {
//            bin_temp[j]=(csk->CSK_arr[i][j] +1)/2;
//            //printf("%d ",bin_temp[j]);
//        }
//        GF_temp = Bin2GF(bin_temp,GF,logGF,BINGF);
//
//        printf(" %d \n", GF_temp);
//
//    }
//    getchar();






        if (PN_index < GF )
        {
            printf("error in build_punctured_csk_mapping: nb_index_found < GF_size !" ); getchar();
        }


    //getchar();
    printf("\n build_punctured_mapping: Success\n");
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

int sequence_tmp[300];
       // repeat construction
    for (i=0; i<code->GF; i++)
    {
        for (j=0; j<csk->PNsize / code->logGF +1; j++)
        {
            for (k=0;k<code->logGF;k++ )
            {
                tmp1=(i>>k) & 0x01;
                //printf("CSK[%d][%d] = %d \n ",i,j*code->logGF + k, tmp1 );
                //csk->CSK_arr[i][j*code->logGF + k]= BPSK(tmp1);
                sequence_tmp[j*code->logGF + k]= BPSK(tmp1);
            }
        }
        for (j=0; j<csk->PNsize; j++)
        {
            csk->CSK_arr[i][j]= sequence_tmp[j];
            //printf("CSK[%d][%d] = %d \n ",i,j, csk->CSK_arr[i][j]);
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
    sigma = sqrt(1.0/(pow(10,EbN/10.0))); // for EsNo considering Bi-AWGN
//sigma = sqrt(1.0/(2*pow(10,EbN/10.0))); // for EsNo considering nose on I and Q

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
                temp = NoisyBin[n][q] *   csk->CSK_arr[g][q]; //convolution
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
//    for (n=N/2; n<N; n++)
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

//    // !!!!! for hybrid CCSK with over-puncturing on redundant symbols
//    for (n=code->K; n<N; n++)
//    {
//        for (g=0; g<code->GF; g++)
//        {
//            TMP[g]=0.0;
//            for (q=0; q<2; q++)
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
 * \fn ModelChannel_AWGN_16 (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN, int *init_rand)
 * \brief 16 modulation + AWGN noise on the codeword.
 *                 The function computes the intrinsic_LLRs corresponding to the noisy observations.
 * Inputs
 * 	- NBIN : Binary copy of the codeword
 * 	- EbN  : Signal to noise ratio in terms of Eb/No (dB)
 * Outputs
 *      - decoder->intrinsic_LLR
 */
void ModelChannel_AWGN_16_CSK(csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
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



int pn_size = 16;


int CCSK16[16]={0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

//    for(i=0; i< 16 ; i++)
//    {
//CCSK16[i]=CCSK16[i]-1;
//   }
    //getchar();

    ////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < code->GF ; i++)
    {
        norm_factor = table_16QAM[i][0]*table_64QAM[i][0] + table_16QAM[i][1]*table_16QAM[i][1]+norm_factor;//compute sum
    }
    norm_factor = sqrt( code->GF / norm_factor);
// printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< code->GF ; i++)
    {
        modulation[i][0]=norm_factor*table_16QAM[i][0];
        modulation[i][1]=norm_factor*table_16QAM[i][1];

    }

// print constellation in a file
    FILE *opfile;
    opfile=fopen("16QAM.txt","a");
    for(i=0; i< code->GF ; i++)
    {
        fprintf(opfile," %f %f ;", modulation[i][0],modulation[i][1]);
    }
    fclose(opfile);
    getchar();




    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {
        som=0;
        for (q=0; q<4; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        //printf("\n %d \n",som );

		for(k=0; k<code->GF; k++)
        {
            TMP[k] =0.0;
        }

        for (q=0; q<csk->PNsize; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                NoisyBin[q][i] = modulation[CCSK16[(som+q)%pn_size]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }




        if (n<N) // per symbol puncturing
        {

            for(k=0; k<code->GF; k++)
            {

				som=0;
				for (q=0; q<4; q++)
				{
					som = som + BinGF_16[k][q]*pow(2,q);
				}

				for (g=0; g<csk->PNsize; g++)
                {

                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK16[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK16[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
                }
            }
        }
        else
        {

            for(k=0; k<code->GF; k++)
            {

				som=0;
				for (q=0; q<4; q++)
				{
					som = som + BinGF_64[k][q]*pow(2,q);
				}

                for (g=0; g<csk->PNsize-1; g++)
                {
                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK16[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK16[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
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



int pn_size = 64;


// test of build_nb_pn2 (increase minimal distance) (it works!)
// int CCSK64[64]={11, 8, 18, 56, 37, 51, 61, 5, 42, 7, 4, 50, 17, 3, 59, 20, 64, 53, 12, 24, 38, 1, 13, 2, 23, 39, 26, 30, 62, 31, 44, 6, 47, 45, 19, 33, 15, 40, 49, 10, 60, 54, 14, 21, 29, 28, 46, 55, 43, 25, 48, 63, 36, 9, 16, 35, 22, 32, 57, 41, 52, 58, 34, 27
//};

// min=3.16
// int CCSK64[64]={4, 25, 14, 26, 15, 16, 43, 1, 28, 39, 38, 58, 19, 59, 48, 6, 54, 27, 17, 3, 32, 46, 21, 13, 20, 50, 41, 61, 12, 7, 55, 9, 22, 56, 49, 23, 33, 52, 18, 5, 24, 40, 10, 11, 47, 51, 53, 57, 34, 62, 30, 42, 37, 36, 35, 8, 63, 31, 29, 44, 2, 45, 64, 60
//};

//test
//  int CCSK64[64]={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64
//  };

// min=3.222
//const int CCSK64[64]={  36, 58, 54, 28, 60, 35, 22, 45, 39, 57, 4, 6, 27, 49, 52, 1, 14, 59, 42, 19, 47, 31, 61, 30, 5, 12, 33, 23, 7, 56, 44, 55, 3, 16, 0, 11, 62, 21, 46, 25, 41, 38, 51, 37, 13, 15, 18, 10, 2, 53, 32, 63, 17, 26, 24, 48, 34, 40, 20, 43, 9, 50, 29, 8
//  };
// min=5.2 considering square distance
//int CCSK64[64]={21, 23, 9, 42, 60, 30, 34, 5, 63, 3, 18, 36, 51, 35, 20, 53, 22, 11, 62, 28, 52, 4, 58, 37, 55, 50, 10, 15, 7, 0, 25, 27, 2, 41, 29, 38, 44, 47, 19, 6, 13, 8, 39, 17, 24, 43, 26, 46, 49, 54, 45, 61, 14, 40, 16, 31, 48, 1, 32, 12, 33, 59, 57, 56   };

// min = 5.26
//const int CCSK64[64]={ 58, 8, 56, 36, 0, 53, 45, 46, 14, 57, 54, 62, 16, 3, 22, 44, 41, 9, 42, 26, 13, 37, 31, 19, 12, 34, 55, 63, 32, 25, 21, 39, 40, 6, 4, 23, 30, 15, 47, 61, 11, 17, 59, 5, 1, 51, 48, 33, 27, 35, 29, 18, 52, 2, 28, 38, 49, 60, 50, 7, 24, 20, 43, 10};

// min = 5.48
//const int CCSK64[64]={44, 54, 11, 46, 17, 0, 43, 21, 48, 56, 3, 20, 15, 63, 39, 27, 52, 18, 26, 4, 61, 13, 29, 59, 32, 24, 33, 22, 25, 37, 57, 62, 14, 30, 53, 40, 6, 31, 51, 55, 58, 7, 5, 42, 8, 23, 36, 49, 41, 19, 45, 50, 38, 9, 10, 16, 34, 47, 35, 12, 28, 2, 60, 1};

// min=7.8 from EB  best!!!!
//int CCSK64[64]={48, 25, 51, 37, 13, 55, 36, 57, 41, 56, 33, 39, 45, 18, 14, 35, 22, 58, 24, 10, 63, 6, 1, 27, 53, 8, 43, 54, 59, 20, 60, 19, 40, 12, 21, 42, 9, 46, 31, 50, 64, 38, 62, 44, 34, 26, 15, 3, 47, 2, 5, 16, 61, 11, 17, 29, 7, 32, 28, 49, 52, 4, 23, 30};

// min=8.6 from EB
//int CCSK64[64]={63, 42, 39, 41, 53, 51, 29, 19, 26, 55, 32, 14, 9, 2, 40, 45, 48, 1, 21, 50, 24, 25, 38, 60, 37, 4, 64, 28, 3, 33, 15, 56, 16, 52, 12, 7, 46, 5, 13, 57, 20, 34, 49, 62, 10, 58, 11, 43, 8, 31, 17, 22, 61, 47, 36, 30, 18, 44, 54, 35, 27, 6, 23, 59};

// min=6.47 from CM
//int CCSK64[64]={28, 17, 24, 40, 5, 51, 3, 16, 39, 63, 1, 18, 0, 62, 31, 52, 26, 27, 46, 59, 7, 2, 29, 44, 36, 61, 55, 48, 30, 6, 19, 34, 47, 45, 54, 41, 32, 20, 43, 21, 15, 38, 50, 23, 33, 13, 58, 9, 25, 49, 42, 60, 22, 35, 4, 56, 12, 53, 11, 37, 14, 10, 8, 57};

//min=9.25
//int CCSK64[64]={34, 36, 21, 16, 45, 46, 50, 9, 20, 12, 55, 2, 10, 48, 6, 61, 3, 19, 8, 27, 15, 24, 51, 32, 1, 18, 14, 29, 52, 54, 22, 49, 39, 43, 38, 11, 30, 28, 62, 40, 58, 47, 44, 60, 23, 26,  5, 17, 25, 63, 13,  7, 42, 53, 64,  4, 35, 37, 33, 31, 41, 57, 56, 59
//min=10.55 alexandru
//int CCSK64[64]={21, 57, 63, 25, 38, 51, 7, 26, 16, 8, 43, 5, 36, 54, 2, 3, 55, 61, 14, 45, 37, 1, 41, 56, 15, 49, 22, 52, 59, 12, 4, 19, 50, 48, 64, 13, 23, 27, 29, 18, 9, 58, 39, 60, 6, 28, 46, 62, 32, 33, 20, 11, 40, 30, 47, 31, 24, 44, 17, 42, 34, 53, 35, 10 };
//10.5525,14.7053,24.3705
//int CCSK64[64]={24,44,17,42,34,53,35,10,21,57,63,29,26,16,8,43,37,1,41,56,15,49,18,9,58,39,60,6,28,46,62,32,33,20,11,40,30,23,47,31,25,38,51,7,22,52,59,12,4,19,50,48,64,13,5,36,54,2,3,55,61,14,45,27};

//max(ccorr(2:64))=27.7
//int CCSK64[64]={35,15,14,62,24,53,34,2,49,32,42,13,57,25,59,51,63,8,61,46,29,5,16,28,47,18,27,0,20,6,36,45,38,23,31,11,39,54,26,55,4,48,17,10,60,30,44,43,12,3,22,50,1,58,9,33,21,40,7,41,52,19,56,37};


//QAM
// 0.47
//int CCSK64[64]={29,42,55,44,35,45,18,13,46,59,41,2,27,56,62,14,6,9,25,39,32,15,40,0,26,54,53,19,47,23,7,61,48,22,43,58,8,60,16,51,4,63,50,1,49,12,31,20,30,34,57,24,33,52,5,11,21,36,17,38,10,3,28,37};

//0.47  1.05  2.1  2.9  4.47 5.3
//int CCSK64[64]={2,22,46,40,35,60,41,8,10,29,38,31,63,54,34,15,43,13,4,55,24,37,5,14,7,18,52,48,1,25,53,26,11,58,6,20,12,62,51,21,9,59,32,30,28,42,56,39,3,50,49,45,19,57,0,44,61,33,23,47,17,16,27,36};






//64PSK
//int CCSK64[64]={17,20,39,59,32,9,8,53,50,54,45,62,10,51,12,63,33,36,55,11,48,25,24,5,2,6,61,14,26,3,28,15,49,52,7,27,64,41,40,21,18,22,13,30,42,19,44,31,1,4,23,43,16,57,56,37,34,38,29,46,58,35,60,47};
//int CCSK64[64]={3,53,48,46,23,42,52,17,44,29,56,54,63,18,59,25,19,5,64,62,39,58,4,33,60,45,8,6,15,34,11,41,35,21,16,14,55,10,20,49,12,61,24,22,31,50,27,57,51,37,32,30,7,26,36,1,28,13,40,38,47,2,43,9};
//epicycloid 4 cusps
//int CCSK64[64]={0,63,58,33,36,51,62,53,8,39,2,9,44,27,6,29,16,15,10,49,52,3,14,5,24,55,18,25,60,43,22,45,32,31,26,1,4,19,30,21,40,7,34,41,12,59,38,61,48,47,42,17,20,35,46,37,56,23,50,57,28,11,54,13};


int CCSK64[64]={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,64};

    for(i=0; i< 64 ; i++)
    {
    CCSK64[i]=CCSK64[i]-1;
   }




   // getchar();
//
//    ////compute normalization factor so that average power of one point of constellation is equal to one
    float norm_factor=0.0;
    for(i=0; i < code->GF ; i++)
    {
        //norm_factor = table_64QAM[i][0]*table_64QAM[i][0] + table_64QAM[i][1]*table_64QAM[i][1]+norm_factor;//compute sum
        //norm_factor = table_64APSK[i][0]*table_64APSK[i][0] + table_64APSK[i][1]*table_64APSK[i][1]+norm_factor;//compute sum
        norm_factor = table_64APSK_EB[i][0]*table_64APSK_EB[i][0] + table_64APSK_EB[i][1]*table_64APSK_EB[i][1]+norm_factor;//compute sum

    }
    norm_factor = sqrt( code->GF / norm_factor);
 //printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< code->GF ; i++)
    {

        //modulation[i][0]=norm_factor*table_64QAM[i][0];
        //modulation[i][1]=norm_factor*table_64QAM[i][1];
        //modulation[i][0]=norm_factor*table_64APSK[i][0];
        //modulation[i][1]=norm_factor*table_64APSK[i][1];
        modulation[i][0]=norm_factor*table_64APSK_EB[i][0];
        modulation[i][1]=norm_factor*table_64APSK_EB[i][1];

    }



//      //  for 64PSK
//    for (i=0;i<code->GF; i++)
//    {
//        modulation[i][0]=sin(i*2*PI/64);
//        modulation[i][1]=cos(i*2*PI/64);
//       // modulation[i][0]=table_64PSK[i-1][0];
//       // modulation[i][1]=table_64PSK[i-1][1];
//       //modulation[i][0]=table_64CC[i][0];
//       //modulation[i][1]=table_64CC[i][1];
//    }


//// print constellation in a file
//    FILE *opfile;
//    opfile=fopen("64QAM.txt","a");
//    for(i=0; i< code->GF ; i++)
//    {
//        fprintf(opfile," %f %f ;", modulation[i][0],modulation[i][1]);
//    }
//    fclose(opfile);
//    getchar();




    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {
        som=0;
        for (q=0; q<6; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        //printf("\n %d \n",som );

		for(k=0; k<code->GF; k++)
        {
            TMP[k] =0.0;
        }

        for (q=0; q<csk->PNsize; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                NoisyBin[q][i] = modulation[CCSK64[(som+q)%pn_size]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }




        if (n < N) // per symbol puncturing
        {

            for(k=0; k<code->GF; k++)
            {

				som=0;
				for (q=0; q<6; q++)
				{
					som = som + BinGF_64[k][q]*pow(2,q);
				}

				for (g=0; g<csk->PNsize; g++)
                {

                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
                    //printf("%d %f \n",k, TMP[k]);

                   //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    //TMP[k] = TMP[k]-( NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][0] + NoisyBin[g][1]*modulation[CCSK64[(som+g)%128]][1] ); //real part


                }
            }
        }
        else
        {

            for(k=0; k<code->GF; k++)
            {

				som=0;
				for (q=0; q<6; q++)
				{
					som = som + BinGF_64[k][q]*pow(2,q);
				}

                for (g=0; g<csk->PNsize-1; g++)
                {
                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
                    //printf("%d %f \n",k, TMP[k]);
                    //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    //TMP[k] = TMP[k]-( NoisyBin[g][0]*modulation[CCSK64[(som+g)%pn_size]][0] + NoisyBin[g][1]*modulation[CCSK64[(som+g)%pn_size]][1] );// + SQR(  NoisyBin[g][1] * modulation[CCSK64[(som+g)%128]][0] -  NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][1]        );

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






 // map GF64 symbol on 2 8ary CCSK
void ModelChannel_AWGN_64_8CSK(csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;
    float **NoisyBin = calloc(csk->PNsize,sizeof(float *));
    for (q=0; q<csk->PNsize; q++) NoisyBin[q] = calloc(2,sizeof(float));
    int i;
    float modulation[code->GF][2];


    int pn_size = 8;
    int CCSK8[8] = {0,1,6,7,4,5,2,3};


    for (i=0;i<code->GF; i++)//for 8PSK
    {
        modulation[i][0]=table_8PSK[i][0];
        modulation[i][1]=table_8PSK[i][1];
    }


    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {


		for(k=0; k<code->GF; k++){TMP[k] =0.0;}

        // first 8CCSK sequence
        som=0;
        for (q=0; q<3; q++){som = som + NBIN[n][q]*pow(2,q);}
        for (q=0; q<csk->PNsize; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                // CCSK modulation + AWGN noise (Box Muller method for Gaussian sampling)
                NoisyBin[q][i] = modulation[CCSK8[(som+q)%pn_size]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }
        for(k=0; k<code->GF; k++)
        {
            som=0;
            for (q=0; q<3; q++){som = som + BinGF_64[k][q]*pow(2,q);}

            for (g=0; g<csk->PNsize; g++)
            {
                TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK8[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK8[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
            }
        }

        //second 8CCSK sequence
        som=0;
        for (q=3; q<6; q++){som = som + NBIN[n][q]*pow(2,q-3);}
        for (q=0; q<csk->PNsize; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                // BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling)
                NoisyBin[q][i] = modulation[CCSK8[(som+q)%pn_size]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
            }
        }
        for(k=0; k<code->GF; k++)
        {
            som=0;
            for (q=3; q<6; q++){som = som + BinGF_64[k][q]*pow(2,q-3);}

            for (g=0; g<csk->PNsize; g++)
            {
                TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK8[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK8[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
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
    int transmited = 1 ;
    int mapping_step =1;
    int mappinp_start =0;

    float **NoisyBin = calloc(transmited,sizeof(float *));
    for (q=0; q<transmited; q++) NoisyBin[q] = calloc(2,sizeof(float));
    int i;

    sigma = sqrt(1.0/(2.0*pow(10,EbN/10.0)));
    //sigma = sqrt(1.0/pow(10,EbN/10.0));

    for (n=0; n<4*N/5; n++)
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
                        TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-chu_real[(mappinp_start+som*mapping_step+g)%csk->PNsize])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-chu_imag[(mappinp_start+som*mapping_step+g)%csk->PNsize])/(2.0*SQR(sigma));
                        //printf("%d %f \n",k, TMP[k]);

                    //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    //TMP[k] = TMP[k]-( NoisyBin[g][0]*chu_real[(mappinp_start+som*mapping_step+g)%csk->PNsize]     + NoisyBin[g][1]*chu_imag[(mappinp_start+som*mapping_step+g)%csk->PNsize]    );// + SQR(  NoisyBin[g][1] * modulation[CCSK64[(som+g)%128]][0] -  NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][1]        );


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

 /*   const int CCSK256[512]={56, 107, 10, 45, 108, 177, 227, 194, 22, 98, 236, 36, 66, 99, 234, 82, 62, 237,
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
       };*/
// dist min = 0.1
//const int CCSK256[256]={84, 32, 206, 96, 192, 112, 161, 63, 70, 163, 7, 159, 79, 22, 188, 72, 135, 23, 251, 122, 146, 36, 195, 221, 131, 125, 69, 202, 28, 183, 101, 148, 59, 107, 61, 223, 211, 33, 230, 141, 162, 29, 117, 120, 25, 216, 166, 129, 97, 255, 154, 144, 74, 152, 226, 170, 45, 217, 124, 19, 200, 67, 83, 186, 254, 121, 73, 239, 198, 143, 204, 173, 114, 49, 172, 149, 5, 57, 94, 196, 66, 106, 160, 156, 244, 60, 126, 212, 246, 238, 142, 199, 133, 128, 222, 134, 164, 17, 243, 240, 116, 190, 87, 191, 34, 91, 176, 37, 147, 65, 132, 241, 31, 123, 182, 119, 18, 77, 165, 215, 64, 247, 55, 6, 95, 13, 44, 250, 89, 46, 136, 48, 110, 145, 139, 137, 127, 210, 235, 50, 231, 40, 232, 24, 85, 43, 93, 194, 75, 181, 174, 171, 12, 62, 26, 245, 10, 169, 208, 71, 115, 108, 158, 1, 233, 35, 104, 105, 86, 82, 11, 98, 41, 209, 16, 130, 138, 90, 0, 100, 197, 109, 38, 51, 118, 229, 207, 76, 236, 249, 111, 234, 179, 21, 47, 58, 187, 227, 78, 140, 225, 180, 157, 52, 178, 184, 242, 224, 9, 155, 102, 2, 185, 42, 81, 203, 88, 253, 201, 15, 220, 168, 151, 80, 150, 99, 27, 3, 167, 218, 54, 56, 175, 103, 189, 113, 237, 39, 30, 53, 92, 14, 20, 193, 4, 248, 153, 228, 219, 214, 68, 8, 205, 177, 213, 252
//};
// dist min = 0.12
//const int CCSK256[256]={81, 67, 54, 139, 8, 158, 70, 97, 152, 163, 165, 23, 156, 120, 172, 29, 225, 86, 157, 210, 178, 227, 219, 32, 58, 104, 213, 216, 254, 68, 140, 34, 77, 207, 209, 30, 241, 146, 49, 246, 12, 105, 16, 116, 48, 190, 186, 168, 233, 217, 235, 199, 127, 96, 122, 63, 124, 143, 229, 155, 52, 53, 175, 94, 193, 83, 60, 223, 76, 123, 188, 14, 33, 113, 103, 167, 196, 82, 247, 133, 228, 153, 248, 191, 10, 47, 253, 142, 42, 41, 208, 73, 84, 169, 7, 119, 100, 27, 95, 204, 51, 244, 64, 137, 230, 243, 121, 71, 179, 24, 203, 154, 5, 136, 106, 249, 101, 222, 232, 3, 176, 206, 87, 174, 138, 197, 144, 149, 11, 62, 226, 211, 9, 69, 99, 184, 39, 218, 110, 55, 92, 74, 0, 252, 131, 224, 93, 31, 151, 135, 35, 161, 194, 88, 20, 17, 160, 237, 61, 65, 221, 147, 177, 126, 15, 21, 118, 38, 46, 195, 242, 4, 245, 198, 56, 6, 240, 130, 236, 90, 166, 201, 1, 181, 107, 145, 43, 234, 112, 231, 109, 28, 18, 89, 182, 22, 44, 80, 50, 141, 251, 255, 180, 128, 214, 239, 205, 85, 238, 36, 159, 162, 2, 25, 102, 183, 13, 220, 111, 78, 108, 75, 45, 26, 192, 150, 170, 129, 19, 79, 215, 115, 66, 202, 171, 57, 148, 72, 173, 134, 59, 250, 40, 114, 37, 189, 164, 125, 91, 185, 98, 187, 132, 117, 200, 212
//};
// dist min=0.13
//const int CCSK256[256]={194, 143, 15, 180, 102, 233, 155, 217, 69, 245, 238, 164, 34, 198, 231, 42, 58, 204, 140, 83, 91, 122, 144, 168, 183, 196, 19, 135, 178, 205, 38, 211, 24, 25, 137, 208, 171, 49, 207, 182, 225, 10, 18, 105, 48, 75, 216, 131, 7, 127, 20, 77, 16, 81, 218, 40, 145, 1, 70, 95, 65, 229, 224, 67, 138, 118, 175, 185, 22, 112, 110, 201, 86, 115, 136, 248, 74, 63, 200, 104, 13, 97, 103, 244, 87, 134, 152, 251, 44, 167, 130, 32, 111, 52, 12, 255, 132, 162, 147, 241, 202, 184, 227, 242, 150, 154, 141, 191, 197, 85, 203, 76, 237, 235, 92, 179, 5, 236, 128, 172, 148, 153, 249, 78, 247, 47, 94, 226, 163, 174, 109, 66, 4, 9, 55, 250, 213, 114, 68, 113, 177, 151, 6, 234, 108, 252, 176, 62, 206, 230, 14, 239, 121, 27, 124, 139, 8, 36, 31, 209, 123, 35, 43, 160, 246, 89, 100, 23, 37, 173, 88, 149, 106, 228, 156, 133, 165, 214, 166, 220, 3, 253, 199, 84, 39, 60, 98, 222, 93, 71, 80, 56, 73, 129, 21, 33, 120, 107, 146, 187, 82, 51, 159, 210, 30, 190, 50, 117, 215, 54, 142, 96, 41, 221, 223, 11, 61, 126, 188, 17, 186, 158, 192, 72, 240, 29, 119, 53, 181, 57, 0, 193, 2, 125, 59, 79, 45, 116, 26, 195, 219, 189, 243, 90, 157, 99, 232, 161, 101, 170, 169, 212, 254, 64, 46, 28
//};

// dist min=0.182
//const int CCSK256[256]={1, 240, 61, 138, 43, 198, 111, 246, 83, 40, 87, 131, 254, 125, 84, 244, 168, 129, 176, 88, 172, 64, 220, 68, 56, 252, 127, 229, 158, 122, 130, 139, 212, 0, 78, 190, 86, 103, 199, 67, 110, 164, 22, 184, 196, 81, 236, 21, 201, 221, 13, 132, 30, 243, 124, 219, 228, 192, 197, 74, 77, 227, 34, 202, 134, 123, 37, 31, 6, 107, 52, 121, 238, 231, 180, 249, 137, 223, 32, 16, 187, 23, 226, 167, 63, 225, 162, 18, 203, 204, 29, 59, 133, 189, 93, 113, 94, 232, 46, 250, 4, 237, 71, 109, 135, 9, 70, 153, 76, 151, 222, 217, 11, 26, 95, 178, 152, 19, 159, 163, 150, 24, 155, 72, 115, 80, 7, 96, 118, 215, 5, 251, 58, 55, 146, 248, 91, 143, 186, 213, 41, 48, 85, 205, 14, 179, 128, 73, 183, 200, 208, 2, 119, 62, 230, 207, 47, 209, 75, 104, 20, 112, 147, 10, 255, 42, 156, 120, 188, 49, 105, 28, 173, 108, 193, 116, 195, 69, 39, 100, 181, 90, 148, 53, 45, 92, 44, 160, 98, 191, 8, 177, 218, 242, 140, 89, 65, 60, 247, 241, 38, 214, 79, 224, 210, 99, 157, 211, 175, 126, 25, 182, 235, 149, 142, 233, 17, 50, 165, 171, 144, 161, 36, 12, 33, 102, 35, 114, 170, 57, 3, 194, 234, 174, 166, 245, 117, 169, 216, 154, 145, 141, 185, 101, 239, 51, 253, 82, 15, 206, 54, 27, 136, 106, 66, 97
//};
// dist_min=0.2
//const int CCSK256[256]={52, 156, 234, 53, 46, 201, 240, 193, 91, 249, 165, 242, 218, 229, 4, 116, 163, 14, 245, 92, 8, 113, 236, 33, 109, 180, 143, 104, 149, 28, 47, 167, 254, 51, 63, 103, 107, 202, 96, 205, 128, 100, 226, 17, 173, 122, 85, 129, 97, 184, 125, 189, 78, 139, 43, 212, 18, 185, 194, 5, 251, 42, 93, 49, 148, 137, 38, 207, 58, 70, 233, 82, 142, 75, 25, 183, 204, 198, 182, 168, 195, 181, 159, 217, 192, 211, 64, 179, 81, 239, 41, 40, 21, 224, 127, 246, 22, 67, 72, 250, 190, 94, 200, 117, 11, 76, 162, 208, 79, 36, 169, 133, 186, 44, 230, 48, 56, 10, 223, 9, 136, 124, 119, 220, 123, 214, 248, 57, 187, 247, 45, 178, 243, 215, 219, 55, 140, 62, 172, 228, 50, 244, 7, 210, 24, 199, 99, 158, 170, 26, 238, 102, 90, 253, 95, 175, 132, 141, 23, 59, 164, 161, 30, 160, 157, 84, 237, 174, 150, 154, 191, 197, 120, 152, 101, 6, 225, 130, 88, 111, 15, 13, 37, 196, 255, 69, 147, 60, 121, 106, 206, 35, 20, 105, 77, 151, 176, 39, 110, 188, 135, 114, 222, 146, 252, 80, 29, 213, 71, 0, 231, 27, 177, 108, 209, 235, 12, 98, 112, 134, 16, 144, 216, 66, 68, 61, 138, 115, 232, 89, 241, 1, 31, 87, 131, 34, 203, 118, 32, 126, 65, 166, 54, 227, 73, 153, 3, 155, 145, 86, 83, 19, 171, 2, 74, 221
//};
// dist_min=0.207
//const int CCSK256[256]={180, 7, 241, 95, 77, 214, 249, 189, 165, 42, 135, 254, 155, 137, 128, 106, 200, 148, 55, 236, 108, 0, 130, 64, 131, 8, 145, 201, 217, 78, 122, 237, 96, 59, 142, 115, 70, 169, 238, 56, 178, 91, 240, 183, 203, 171, 157, 86, 150, 101, 127, 226, 51, 100, 52, 247, 220, 16, 228, 181, 117, 110, 87, 202, 58, 17, 97, 136, 5, 244, 152, 195, 193, 129, 239, 125, 139, 197, 4, 74, 208, 79, 179, 170, 33, 11, 123, 211, 159, 175, 105, 99, 251, 184, 174, 229, 49, 213, 94, 221, 27, 190, 41, 172, 246, 34, 176, 22, 80, 223, 248, 167, 93, 73, 83, 218, 81, 21, 243, 48, 143, 10, 31, 196, 103, 173, 102, 88, 92, 28, 62, 232, 69, 1, 168, 166, 188, 242, 132, 133, 163, 206, 205, 111, 68, 67, 154, 255, 29, 156, 6, 36, 23, 85, 212, 161, 186, 43, 149, 18, 234, 162, 76, 146, 141, 45, 222, 140, 153, 19, 209, 84, 9, 216, 121, 199, 182, 187, 245, 25, 192, 134, 20, 3, 204, 119, 47, 13, 147, 164, 57, 227, 40, 250, 104, 112, 215, 210, 71, 107, 39, 194, 30, 24, 72, 61, 230, 75, 109, 233, 65, 126, 35, 160, 191, 32, 151, 144, 253, 82, 114, 63, 113, 44, 54, 207, 2, 185, 118, 224, 38, 12, 50, 231, 90, 37, 198, 177, 120, 158, 124, 89, 219, 235, 15, 225, 14, 252, 53, 26, 66, 138, 46, 60, 116, 98
//};
// dist_min=0.217 to check
//const int CCSK256[256]={171, 91, 123, 107, 243, 82, 103, 210, 74, 230, 39, 240, 142, 254, 37, 59, 118, 48, 206, 214, 183, 185, 179, 63, 136, 197, 255, 109, 229, 172, 157, 67, 205, 30, 237, 156, 175, 88, 163, 40, 68, 146, 220, 249, 73, 212, 247, 207, 49, 15, 208, 44, 66, 256, 139, 75, 184, 151, 147, 100, 19, 61, 111, 22, 83, 64, 248, 35, 69, 143, 6, 250, 209, 218, 1, 170, 180, 104, 78, 221, 178, 130, 114, 238, 33, 77, 232, 176, 195, 150, 198, 24, 120, 245, 94, 164, 99, 119, 253, 93, 7, 144, 161, 154, 216, 51, 149, 213, 87, 89, 174, 167, 122, 34, 85, 101, 241, 81, 177, 27, 222, 239, 124, 158, 134, 14, 121, 54, 236, 72, 162, 173, 41, 116, 190, 233, 53, 155, 219, 12, 201, 140, 191, 32, 11, 228, 96, 182, 141, 193, 112, 246, 115, 227, 79, 128, 8, 47, 203, 204, 166, 215, 25, 217, 98, 138, 9, 252, 133, 165, 80, 36, 126, 160, 97, 31, 29, 187, 196, 23, 145, 110, 152, 62, 194, 202, 189, 76, 21, 65, 18, 55, 105, 223, 13, 71, 186, 58, 95, 60, 50, 211, 20, 234, 2, 131, 125, 169, 17, 132, 90, 28, 244, 3, 153, 38, 137, 200, 70, 42, 57, 251, 192, 168, 188, 52, 226, 5, 224, 92, 45, 106, 127, 108, 135, 4, 242, 231, 46, 148, 117, 10, 235, 16, 199, 113, 43, 225, 159, 84, 129, 56, 181, 26, 102, 86
//};
//const int CCSK256[256]={149, 139, 254, 108, 59, 214, 117, 41, 180, 132, 201, 220, 146, 245, 196, 95, 134, 159, 221, 162, 243, 105, 30, 248, 33, 119, 116, 178, 144, 157, 137, 231, 190, 152, 203, 212, 232, 92, 229, 16, 224, 115, 84, 191, 193, 250, 56, 198, 68, 124, 227, 24, 4, 128, 236, 62, 98, 89, 223, 186, 74, 163, 69, 195, 177, 76, 1, 77, 147, 0, 237, 192, 43, 255, 216, 106, 225, 78, 110, 244, 176, 136, 46, 197, 54, 202, 143, 226, 189, 125, 251, 48, 206, 219, 28, 104, 170, 72, 166, 209, 130, 107, 233, 65, 111, 246, 57, 112, 55, 18, 75, 239, 79, 96, 32, 181, 64, 187, 184, 9, 199, 122, 129, 151, 174, 135, 175, 185, 100, 215, 82, 169, 23, 142, 3, 234, 165, 204, 45, 241, 61, 63, 158, 51, 167, 183, 240, 121, 208, 26, 161, 131, 31, 91, 27, 242, 120, 87, 210, 19, 230, 7, 249, 114, 172, 67, 15, 113, 66, 168, 228, 8, 86, 150, 160, 5, 12, 34, 39, 99, 83, 13, 235, 6, 247, 21, 252, 22, 222, 71, 58, 205, 38, 141, 138, 118, 35, 207, 97, 217, 47, 11, 94, 42, 182, 123, 49, 2, 29, 52, 88, 200, 179, 188, 37, 153, 20, 101, 211, 127, 164, 238, 10, 25, 156, 148, 93, 145, 90, 36, 133, 44, 70, 213, 173, 85, 253, 80, 218, 171, 126, 109, 17, 194, 50, 60, 14, 154, 155, 103, 81, 102, 40, 73, 53, 140
//};
// dist_min = 0.218
//const int CCSK256[256]={61, 231, 122, 232, 218, 185, 129, 182, 56, 201, 192, 141, 137, 230, 114, 133, 151, 204, 128, 58, 226, 4, 181, 211, 139, 90, 46, 81, 72, 147, 180, 104, 26, 245, 78, 199, 98, 31, 6, 17, 132, 41, 176, 142, 35, 237, 59, 256, 10, 250, 251, 144, 48, 148, 166, 243, 102, 19, 194, 255, 51, 253, 111, 156, 118, 136, 14, 173, 117, 60, 190, 93, 109, 248, 193, 167, 223, 219, 131, 249, 112, 172, 203, 11, 183, 207, 220, 155, 105, 150, 184, 64, 202, 197, 108, 178, 73, 53, 80, 99, 67, 161, 154, 244, 113, 239, 70, 174, 57, 160, 175, 134, 44, 33, 143, 42, 168, 39, 121, 187, 94, 212, 126, 191, 170, 152, 252, 87, 2, 242, 205, 158, 246, 145, 36, 157, 153, 86, 52, 95, 221, 43, 188, 22, 236, 206, 77, 140, 18, 12, 54, 228, 227, 214, 83, 124, 32, 106, 45, 69, 101, 3, 146, 186, 38, 9, 254, 149, 75, 27, 130, 49, 110, 29, 196, 84, 171, 233, 76, 213, 222, 16, 210, 209, 116, 88, 216, 37, 247, 215, 92, 115, 8, 195, 189, 229, 15, 125, 164, 91, 21, 107, 238, 5, 97, 225, 217, 163, 165, 28, 79, 89, 66, 20, 30, 240, 82, 24, 71, 241, 1, 85, 100, 179, 25, 55, 177, 162, 74, 159, 40, 138, 103, 96, 198, 7, 68, 47, 208, 224, 127, 63, 50, 200, 119, 13, 234, 135, 169, 62, 65, 120, 23, 123, 235, 34
//};
//dist_min=0.239
//const int CCSK256[256]={135,8,159,14,197,151,90,143,158,145,78,182,118,217,26,113,72,112,89,0,42,172,149,83,69,25,240,81,125,22,56,99,132,204,184,66,219,29,227,47,106,35,109,101,124,216,53,43,64,222,58,230,108,79,236,88,139,195,28,249,40,201,24,157,167,162,39,187,247,7,104,59,225,21,242,98,95,156,189,76,114,238,9,1,220,213,130,136,128,2,23,203,102,121,103,248,210,153,150,144,30,180,169,214,62,207,218,134,33,50,19,215,75,92,186,241,3,91,196,137,115,55,194,146,34,243,37,193,185,221,63,133,160,250,252,116,166,246,191,164,126,65,181,16,179,171,96,234,36,117,27,70,163,54,212,208,87,200,93,77,12,255,155,138,239,61,209,188,49,253,154,245,140,31,32,161,205,18,67,38,11,46,190,105,13,235,80,177,129,52,229,148,4,100,175,206,223,244,45,86,110,232,147,131,44,176,68,168,73,254,85,119,94,5,178,211,231,20,6,120,224,183,226,122,127,142,97,107,199,165,82,141,228,41,198,123,10,174,71,51,233,15,237,60,202,170,173,251,74,152,48,111,84,192,57,17};
// 0.283
//const int CCSK256[256]={161,51,242,151,249,34,35,190,175,116,36,209,132,158,44,160,255,82,31,197,121,251,124,227,149,178,159,166,2,95,217,42,196,108,115,11,241,33,130,68,73,105,27,102,19,23,147,117,188,213,101,67,134,112,81,131,215,79,208,155,60,59,150,75,154,183,243,89,125,72,252,76,142,156,49,233,186,22,230,126,100,231,54,195,15,91,77,99,86,30,107,199,204,62,12,207,232,40,4,248,18,244,206,133,179,239,10,83,144,205,201,17,216,226,176,138,120,113,198,70,109,106,225,185,189,139,177,20,218,80,153,224,211,78,247,137,114,228,65,48,85,127,246,173,57,245,45,181,52,214,152,63,212,61,187,146,58,202,6,203,229,143,1,24,98,220,165,235,200,110,210,145,180,162,104,69,29,50,171,21,119,148,169,191,219,38,256,240,26,3,182,5,122,103,123,174,46,167,129,64,168,253,7,71,184,87,53,47,13,194,222,93,172,254,16,37,238,28,221,55,90,164,25,223,66,141,237,8,14,74,94,96,88,236,84,193,118,135,43,56,41,250,136,234,170,140,39,97,9,157,32,111,128,92,192,163};
// from alex
//const int CCSK256[256]={140,203,40,78,59,189,244,160,242,227,12,168,84,23,157,229,46,101,85,214,186,122,120,43,130,21,89,27,243,161,231,178,210,139,10,37,212,253,58,170,175,176,114,159,224,61,105,177,49,184,64,128,136,151,63,28,91,145,125,158,26,225,204,191,208,60,44,110,102,194,79,205,252,51,33,90,195,251,221,7,135,41,30,235,206,104,14,127,228,5,211,167,123,82,200,108,62,113,116,169,240,222,47,155,20,24,42,88,226,215,50,190,19,11,188,68,92,97,182,77,163,134,216,111,73,166,245,76,179,223,56,133,198,31,209,241,35,75,143,147,45,86,100,141,81,119,237,38,6,132,197,154,80,238,109,187,39,250,164,126,118,22,207,172,52,15,4,18,115,239,230,150,246,1,213,16,55,201,87,236,106,247,129,171,98,193,3,71,36,83,248,234,70,74,192,107,137,53,199,149,96,66,69,0,93,165,72,183,233,196,146,174,142,152,185,29,131,99,112,217,103,65,25,2,254,255,138,232,95,162,156,57,117,121,220,13,124,153,17,173,94,249,9,181,148,32,180,202,34,8,48,67,219,218,144,54};

// 0.116 0.51  1.2 2.12 3.02 4.4 5.6 6.9 8.1
const int CCSK256[256]={101,166,108,58,159,195,147,142,240,36,78,245,65,190,217,21,207,122,103,110,95,227,69,153,67,157,43,220,161,251,28,8,199,139,232,42,84,221,168,189,83,80,254,15,19,243,109,186,219,17,167,211,61,16,59,99,200,132,79,4,118,77,45,185,236,57,60,25,176,163,112,51,136,11,141,137,228,89,13,253,31,226,164,205,150,18,181,170,2,224,138,20,41,204,144,151,52,12,117,180,165,188,225,7,160,194,76,133,91,212,206,29,182,62,123,174,10,235,90,184,208,230,6,72,249,0,81,104,55,242,56,246,24,111,171,140,73,238,33,22,222,241,53,66,74,233,128,100,155,92,177,214,178,145,172,86,143,5,210,169,96,201,9,154,27,203,255,97,215,192,71,44,54,75,98,14,248,115,252,70,87,105,209,48,179,127,216,40,162,39,148,129,113,234,175,23,193,124,3,173,47,223,114,191,125,213,131,85,152,63,250,120,198,64,88,202,49,229,183,34,68,130,26,218,107,149,37,231,146,196,38,32,239,50,46,156,126,121,247,116,1,94,197,158,35,134,187,119,30,237,102,93,135,106,244,82};




//for 256QAM 0.047, 0.235, 0.8 , 1.176  ccorr= 278, 288
//const int CCSK256[256]={221,54,238,10,223,103,200,34,144,135,145,250,248,189,1,251,48,252,31,90,37,96,181,77,30,115,81,112,108,233,255,196,230,89,247,105,42,94,246,86,169,217,133,206,92,158,184,98,241,209,134,58,124,165,7,243,182,164,122,202,210,55,88,151,159,215,128,227,70,201,222,132,166,198,231,172,24,224,193,3,152,61,49,214,4,197,18,244,188,140,253,160,139,76,60,72,107,22,190,12,41,95,63,229,75,28,203,113,137,245,191,154,68,220,194,157,84,50,216,156,80,208,125,104,129,170,239,213,199,141,218,117,177,240,57,254,228,120,46,99,56,195,114,59,73,26,40,29,27,36,9,211,256,121,6,234,8,16,123,186,32,20,78,142,162,185,136,212,106,5,53,85,180,66,111,183,71,101,74,235,82,147,100,148,97,118,102,179,146,226,126,2,173,13,116,207,79,65,91,204,93,131,62,35,205,110,225,174,149,83,52,192,242,138,232,17,51,167,163,39,25,23,45,44,119,236,219,155,187,47,109,33,176,171,38,67,150,161,11,168,143,43,153,130,69,14,178,127,237,249,21,87,175,15,19,64};
// for 256QAM 0.07  0.54  1.13 2.14  3.1  4.2  5.7   6.7
//const int CCSK256[256]={249,64,209,86,223,38,71,109,23,126,111,98,181,136,26,27,84,194,211,79,197,229,200,89,202,16,199,134,153,132,80,55,73,0,92,76,20,118,254,218,54,245,226,19,99,162,57,219,125,234,244,35,228,188,96,44,8,238,193,180,128,85,124,88,122,236,156,43,28,70,169,97,252,60,61,214,170,155,22,158,179,117,110,6,90,184,166,210,174,203,206,237,100,101,11,123,33,78,190,95,235,2,46,13,82,107,68,93,148,138,74,113,231,172,133,41,102,230,5,173,3,34,205,192,163,185,63,114,144,143,56,151,47,17,119,196,225,227,239,15,14,104,65,49,48,246,147,232,168,201,115,224,183,220,32,67,36,121,241,81,58,137,25,154,105,24,167,116,75,140,243,120,108,69,142,213,221,18,53,21,175,10,171,31,217,253,159,131,182,91,7,208,45,198,141,42,216,77,233,135,187,87,222,66,29,30,160,176,106,251,103,59,139,50,152,177,195,145,12,207,1,186,164,62,189,255,161,9,165,157,39,248,150,240,212,204,4,40,37,247,215,52,127,83,51,149,242,250,130,178,191,112,94,146,129,72};



//const int CCSK256[256]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};


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
        //norm_factor = table_256QAM[i][0]*table_256QAM[i][0] + table_256QAM[i][1]*table_256QAM[i][1]+norm_factor;//compute sum
        norm_factor = table_GF[i][0]*table_GF[i][0] + table_GF[i][1]*table_GF[i][1]+norm_factor;//compute sum

    }
    norm_factor = sqrt( code->GF / norm_factor);
 //printf(" norm_factor = %f ", norm_factor); getchar();

    for(i=0; i< code->GF ; i++)
    {
        //modulation[i][0]=norm_factor*table_256QAM[i][0];
        //modulation[i][1]=norm_factor*table_256QAM[i][1];
        modulation[i][0]=norm_factor*table_GF[i][0];
        modulation[i][1]=norm_factor*table_GF[i][1];
    }



//// print constellation in a file
//    FILE *opfile;
//    opfile=fopen("256QAM.txt","a");
//    for(i=0; i< code->GF ; i++)
//    {
//        fprintf(opfile," %f %f ;", modulation[i][0],modulation[i][1]);
//    }
//    fclose(opfile);
//    getchar();



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

        if (n < N) // no per symbol puncturing: -1 , half:  N/2
        {
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
                        som=0;
                        for (q=0; q<8; q++)
                        {
                            som = som + BinGF_256[k][q]*pow(2,q);
                        }
                        for (g=0; g<csk->PNsize; g++)
                        {

                                TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK256[(som+g)%256]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK256[(som+g)%256]][1])/(2.0*SQR(sigma));
                                //printf("%d %f \n",k, TMP[k]);
                        }
                    }
        }
        else
        {

                    for (q=0; q<csk->PNsize-1; q++)
                    {
                        for (i=0; i<2; i++)
                        {
                            u=My_drand48(init_rand);
                            v=My_drand48(init_rand);
                            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                            NoisyBin[q][i] = modulation[CCSK256[(som+q)%256]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
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
                            TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK256[(som+g)%256]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK256[(som+g)%256]][1])/(2.0*SQR(sigma));
                            //printf("%d %f \n",k, TMP[k]);
                        }
                    }
        }




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



