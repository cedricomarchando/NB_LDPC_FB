/*!
 * \file bubble_decoder.c
 * \brief Check node based on bubbles
 * \author C.Marchand
 * \copyright BSD copyright
 * \date 03/03/2015
 */

#include <stdlib.h>
#include <stdio.h>
#include "./include/syndrome_decoder.h"
#include "./include/bubble_decoder.h"





int maximum(float *Y, int nl)
{

    float aux = 0.0;
    int pos = 0;
    int i;

    for (i=0; i<nl ; i++)
    {
        if (Y[i] > aux)
        {
            aux = Y[i];
            pos = i;
        }
    }

    return pos;

}

int minimum(float *Y, int nl)
{

    float aux = 1e5;
    int pos = 0;
    int i;

    for (i=0; i<nl ; i++)
    {
        if (Y[i] < aux)
        {
            aux = Y[i];
            pos = i;
        }
    }

    return pos;

}







/**
 * \fn CheckPassLogEMS
 * \brief Process the Check to Variable messages in the decoding graph.
 * Inputs
 * 	- decoder->VtoC
 * Outputs
 * 	- decoder->CtoV
 */
void CheckPassLogEMS (int node,decoder_t *decoder, code_t *code, table_t *table,int NbOper, float offset)
{
    const int nbMax = decoder->nbMax;

    int        t,k,g,kk,S,c,k1,Stp;
    int        **MatriceInterIndice,*OutForwardIndice,*OutBackwardIndice,*OutForwardIndice1,*OutBackwardIndice1;
    float      **MatriceInter,*OutForward,*OutBackward,*OutForward1,*OutBackward1;
    float LLR_tmp[code->GF];



    //Temporary buffers (for the F/B recursivity)

    OutForward=(float *)calloc(nbMax,sizeof(float));
    OutForwardIndice=(int *)calloc(nbMax,sizeof(int));
    OutBackward=(float *)calloc(nbMax,sizeof(float));
    OutBackwardIndice=(int *)calloc(nbMax,sizeof(int));

    OutForward1=(float *)calloc(nbMax,sizeof(float));
    OutForwardIndice1=(int *)calloc(nbMax,sizeof(int));
    OutBackward1=(float *)calloc(nbMax,sizeof(float));
    OutBackwardIndice1=(int *)calloc(nbMax,sizeof(int));


    // Number of steps F/B = S
    S=2*(code->rowDegree[node]-2);

// MatriceInter - store the results (reeal values) of the F/B recursivity
    MatriceInter=(float **)calloc(S,sizeof(float*));
    for (k=0; k<S; k++)
    {
        MatriceInter[k]=(float *)calloc(nbMax,sizeof(float));
        for (k1=0; k1<nbMax; k1++)
        {
            MatriceInter[k][k1]=1e5;
        }
    }


//MatriceInterIndice - that store the results (GF(q) values) of the F/B recursivity
    MatriceInterIndice=(int **)calloc(S,sizeof(int*));
    for (k=0; k<S; k++)
    {
        MatriceInterIndice[k]=(int *)calloc(nbMax,sizeof(int));
        for (k1=0; k1<nbMax; k1++)
        {
            MatriceInterIndice[k][k1]=-1;
        }
    }

//Initialization of the temporary buffers
    for (k=0; k<nbMax; k++)
    {
        OutForward[k]=1e5;
        OutBackward[k]=1e5;
        OutForwardIndice[k]=-1;
        OutBackwardIndice[k]=-1;
    }


//Rotation par la valeur non nulle dans la matrice pour VtoC
    for(t=0; t<code->rowDegree[node]; t++)
    {
        //printf("\n node=%d, coef=%d \n",node,code->matValue[node][t] );
        for(k=0; k<nbMax; k++)
        {

            //printf(" %d ->",decoder->M_VtoC_GF[t][k]);


            if (decoder->M_VtoC_GF[t][k]!=-1)
            {
                c=decoder->M_VtoC_GF[t][k];
                decoder->M_VtoC_GF[t][k]=table->MULGF[c][code->matValue[node][t]];
            }
            else printf("M_VtoC_GF[%d][%d]=%d\n",t,k,decoder->M_VtoC_GF[t][k]);

            //printf(" %d ->",decoder->M_VtoC_GF[t][k]);

        }
    }



//Initialisation algorithm
    for(k=0; k<nbMax; k++)
    {
        OutForward[k]=decoder->M_VtoC_LLR[0][k];
        OutForwardIndice[k]=decoder->M_VtoC_GF[0][k];
        OutBackward[k]=decoder->M_VtoC_LLR[code->rowDegree[node]-1][k];
        OutBackwardIndice[k]=decoder->M_VtoC_GF[code->rowDegree[node]-1][k];
    }

// Start of the recusivity
    for(kk=1; kk<(code->rowDegree[node]-1); kk++)
    {


        for(k=0; k<nbMax; k++)
        {
            //forward step
            OutForward1[k]=decoder->M_VtoC_LLR[kk][k];
            OutForwardIndice1[k]=decoder->M_VtoC_GF[kk][k];
            //backward step
            OutBackward1[k]=decoder->M_VtoC_LLR[code->rowDegree[node]-kk-1][k];
            OutBackwardIndice1[k]=decoder-> M_VtoC_GF[code->rowDegree[node]-kk-1][k];
        }

        for(k=0; k<nbMax; k++)
        {
            // Storage of the intermediate vectors
            MatriceInter[kk-1][k]=OutForward[k];
            MatriceInterIndice[kk-1][k]=OutForwardIndice[k];
            MatriceInter[2*(code->rowDegree[node]-2)-kk][k]=OutBackward[k];
            MatriceInterIndice[2*(code->rowDegree[node]-2)-kk][k]=OutBackwardIndice[k];
        }

        //forward step
        //printf(" \n forward step %d \n",kk);
        if(kk<(code->rowDegree[node]-1))
            ElementaryStep(OutForward,OutForward1,OutForwardIndice,OutForwardIndice1,OutForward,OutForwardIndice,table->ADDGF,code->GF,nbMax,NbOper);


        //backward step
        //printf(" \n backward step %d \n",kk);
        if(kk<(code->rowDegree[node]-1))
            ElementaryStep(OutBackward,OutBackward1,OutBackwardIndice,OutBackwardIndice1,OutBackward,OutBackwardIndice,table->ADDGF,code->GF,nbMax,NbOper);

    }


    //Update of vectors CtoV (first and last)
    for(k=0; k<nbMax; k++)
    {
        // last vector M_CtoV_LLR
        decoder->M_CtoV_LLR[code->rowDegree[node]-1][k]=OutForward[k];
        decoder->M_CtoV_GF[code->rowDegree[node]-1][k]=OutForwardIndice[k];
        // first vector M_CtoV_LLR
        decoder->M_CtoV_LLR[0][k]=OutBackward[k];
        decoder->M_CtoV_GF[0][k]=OutBackwardIndice[k];

    }

    //printf("\n merge \n");
    //Update of the others vectors CtoV
    for(k=0; k<(code->rowDegree[node]-2); k++)
    {


        ElementaryStep(MatriceInter[k],MatriceInter[(code->rowDegree[node]-2)+k],MatriceInterIndice[k],MatriceInterIndice[(code->rowDegree[node]-2)  +k],OutForward,OutForwardIndice,table->ADDGF,code->GF,nbMax,NbOper);
        for(g=0; g<nbMax; g++)
        {
            decoder->M_CtoV_LLR[k+1][g]=OutForward[g];
            decoder->M_CtoV_GF[k+1][g]=OutForwardIndice[g];
        }
    }


//Rotation par la valeur non nulle dans la matrice pour CtoV
    for(t=0; t<code->rowDegree[node]; t++)
    {
        for(k=0; k<nbMax; k++)
        {

            if (decoder->M_CtoV_GF[t][k]== -1)
            {
                Stp=k;
                break;
            }
            else
                Stp=nbMax;
        }





        for(k=0; k<Stp; k++)
        {
            decoder->M_CtoV_GF[t][k]=table->DIVGF[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];


        }


/// reorder in GF order and add offset
//printf("sat=%0.3f \n", decoder->M_CtoV_LLR[t][Stp-1] + offset );
//getchar();


        for(k=0; k<code->GF; k++)
        {
            LLR_tmp[k] = decoder->M_CtoV_LLR[t][Stp-1] + offset;
        }

        for(k=0; k<Stp; k++)
        {
            LLR_tmp[decoder->M_CtoV_GF[t][k]]= decoder->M_CtoV_LLR[t][k];
        }



        for(k=0; k<code->GF; k++)
        {
            decoder->M_CtoV_GF[t][k] = k;
            decoder->M_CtoV_LLR[t][k] = LLR_tmp[k];
        }


    }

//    for(g=0; g<code->GF; g++)
//        {
//            printf(" GF=%d, LLR=%0.2f \n",decoder->M_CtoV_GF[0][g],decoder->M_CtoV_LLR[0][g]  );
//        }
//        getchar();



    for (k=0; k<S; k++) free(MatriceInterIndice[k]);
    free(MatriceInterIndice);
    for (k=0; k<S; k++) free(MatriceInter[k]);
    free(MatriceInter);

    free(OutForward);
    free(OutForwardIndice);
    free(OutBackward);
    free(OutBackwardIndice);
    free(OutForward1);
    free(OutForwardIndice1);
    free(OutBackward1);
    free(OutBackwardIndice1);

}






/**
 * \fn ElementaryStep with x bubbles
 * \brief Elementary bubble check
 */
int ElementaryStep(float *Input1,float *Input2,int *IndiceInput1,int *IndiceInput2,float *Output,int *IndiceOut,int **ADDGF,int GF,int nbMax,int nbOper)
{

    int nmU = nbMax;
    int nmV = nbMax;
    int nmS = nbMax;

    float loc_Output[nbMax];
    int loc_IndiceOut[nbMax];

    int  i,j, s, ss, Indice_aux=0;
    int nb_bubble=8;

    int pos;

    float tab_aux[nbMax][nbMax];

    float tab_comp[3][nb_bubble];
    int GFvalues_already_in_Out[GF];

    //int debug = 0;


//            printf("\n GF in 1 \n");
//            for(i=0;i<nmU;i++)
//                {
//                 printf(" %d",IndiceInput1[i] );
//
//                }
//            printf("\n GF in 2 \n");
//            for(i=0;i<nmV;i++)
//                {
//                 printf(" %d",IndiceInput2[i] );
//
//                }
//                   printf(" \n input1 \n");
//    for(i=0; i<nbMax; i++)
//    {
//        printf(" %f  ",Input1[i]);
//    }
//       printf(" \n input2 \n");
//    for(i=0; i<nbMax; i++)
//    {
//        printf(" %f  ",Input2[i]);
//    }
//    printf("\n");


    // init
    for (i=0; i<GF; i++)
    {
        GFvalues_already_in_Out[i]=-1;
    }

    for (i=0; i<nbMax; i++)  // if nbOper < number of operations needed to completely fill de nbMax elements of Output, then -1000 and -1 are provided
    {
        loc_Output[i] = 1e5;
        loc_IndiceOut[i] = -1;
    }
    for (i=0; i<nmU; i++)
    {
        for (j=0; j<nmV; j++)
        {
            tab_aux[i][j] = 0;
        }
    }


// pre-compute the addition of bubble check
    // horizontal
    for (j=0; j< nb_bubble>>1; j++)
    {
        for (i=0; i< nmU; i++)
        {
            tab_aux[i][j] = Input1[i]+Input2[j]; //min-sum algorithm

//            // min-max algorithm
//            if (Input1[i]>Input2[j])
//            {
//             tab_aux[i][j] = Input1[i];
//            }
//            else
//            {
//            tab_aux[i][j] = Input2[j];
//            }




        }
    }
    // vertical
     for (i=0; i<nb_bubble>>1; i++)
    {
        for (j=nb_bubble>>1; j<nmV; j++)
        {
            tab_aux[i][j] = Input1[i]+Input2[j]; // min-sim algorithm
//                        // min-max algorithm
//            if (Input1[i]>Input2[j])
//            {
//             tab_aux[i][j] = Input1[i];
//            }
//            else
//            {
//            tab_aux[i][j] = Input2[j];
//            }
        }
    }


//printf(" nbmax=%d \n",nbMax);
//    for (i=0; i<nbMax; i++)
//        {
//        for (j=0; j<nbMax; j++)
//            {
//            printf(" %0.2f \t",tab_aux[i][j]);
//            }
//            printf(" \n");
//        }
//        getchar();





    // tab_comp, comparator matrix, contains the "competitor elements" and their coordinates in "tab_aux"
    // Initialisation of tab_comp

    //vertical competitors
    for (j=0; j<nb_bubble/2; j++)
    {
        tab_comp[0][j]= tab_aux[j][0];
        //tab_aux[j][0]=-1;   /* for each element getting from tab_aux to tab_comp, its position in tab_aux becomes -1
        // This way, we control the problem described in *modified 12/02/2008*/
        tab_comp[1][j]=j;
        tab_comp[2][j]=0;
    }

    //horizontal competitors
    for (j=0; j<nb_bubble/2; j++)
    {
        tab_comp[0][j + nb_bubble/2]= tab_aux[nb_bubble/2][j];
        tab_comp[1][j + nb_bubble/2]=nb_bubble/2;
        tab_comp[2][j + nb_bubble/2]=j;
    }







    // filling Out and IndiceOut

    s = 0;

    for(ss=0; ss < nbOper; ss++)
    {

//        printf("s=%d; ss= %d \n", s, ss);
        pos = minimum(tab_comp[0], nb_bubble);

        if ((IndiceInput1[(int)(tab_comp[1][pos])]==-1) || (IndiceInput2[(int)(tab_comp[2][pos])]==-1))
        {

            printf("\n \n out of bounds : IndiceInput1 = %d , IndiceInput1 = %d \n", IndiceInput1[(int)(tab_comp[1][pos])], IndiceInput2[(int)(tab_comp[2][pos])]);

            break;
        }

        Indice_aux = ADDGF[IndiceInput1[(int)(tab_comp[1][pos])]][IndiceInput2[(int)(tab_comp[2][pos])]];

        // control redundancy in the output list
        // transfer from tab_comp to output Out, IndiceOut
        if (GFvalues_already_in_Out[Indice_aux] == -1)
        {
            loc_Output[s] = tab_comp[0][pos];
            loc_IndiceOut[s] = Indice_aux;
            GFvalues_already_in_Out[Indice_aux] = 1;
            s++;
        }


//printf("candiate %d, x=%d, y=%d, llr=%f, GF=%d + %d =%d \n",s,(int)(tab_comp[1][pos]),(int)(tab_comp[2][pos]),
//       tab_comp[0][pos],IndiceInput1[(int)(tab_comp[1][pos])],IndiceInput2[(int)(tab_comp[2][pos])],Indice_aux  );

        if (s==nmS) break;


        //control limits of tab_aux
        if ((tab_comp[1][pos]>= nmU-1)||(tab_comp[2][pos]>=nmV-1))
        {


        //    printf("\n\n  out of bounds pos=%d, tab_aux: x=%d, y=%d , candidate=%d, output number=%d \n",pos,(int)(tab_comp[1][pos]),(int)(tab_comp[2][pos]),ss,s );

//            for(i=0;i<GF;i++){
//                    printf("%d",GFvalues_already_in_Out[i]+1);
//            }
//            printf("\n GF in 1 \n");
//            for(i=0;i<nmU;i++)
//                {
//                 printf(" %d",IndiceInput1[i] );
//
//                }
//            printf("\n GF in 2 \n");
//            for(i=0;i<nmV;i++)
//                {
//                 printf(" %d",IndiceInput2[i] );
//
//                }
//                   printf(" \n input1 \n");
//    for(i=0; i<nbMax; i++)
//    {
//        printf(" %f  ",Input1[i]);
//    }
//       printf(" \n input2 \n");
//    for(i=0; i<nbMax; i++)
//    {
//        printf(" %f  ",Input2[i]);
//    }
//
//
//getchar();

break;


        }

        // update tab_comp with next value from tab_aux

        if (pos > nb_bubble/2-1)
        {
            tab_comp[1][pos]=tab_comp[1][pos] + 1;
        }
        else
        {
            tab_comp[2][pos]=tab_comp[2][pos] + 1;
        }

        tab_comp[0][pos]= tab_aux[(int)tab_comp[1][pos]][(int)tab_comp[2][pos]];




    }

    for(i=0; i<nbMax; i++)
    {
        Output[i]=loc_Output[i];
        IndiceOut[i]=loc_IndiceOut[i];
    }



//    printf("\n input2 \n");
//   for(i=0; i<nbMax; i++)
//    {
//        printf(" %f  ",Input2[i]);
//    }
//    printf("\n output \n");
//    for(i=0; i<nbMax; i++)
//    {
//        printf(" %f  ",Output[i]);
//    }
//    for(i=0; i<nbMax; i++)
//    {
//        printf(" %d  ",IndiceOut[i]);
//    }
//    printf(" \n");
//    getchar();




return(0);
}







