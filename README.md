# Welcome to the c code for NB-LDPC simulation

You can use this code to simulation NB-LDPC matrices using the Extented-Min Sum (EMS) algorithm. The Check Node (CN) is processed using Forward Backward(FB) algorithm. The FB algorithm splits CN   in elementary CNs (ECN)

#usage

## input argument
 there are 7 arguments
 1.		NbMonteCarlo     : # simulated frames
 1.		NbIterMax        : # of maximum decoding iteration
 1.		FileMatrix       : File name of the parity-check matrix
 1.		EbN              : Eb/No (dB)
 1.		NbMax            : size of truncated messages
 1.		Offset           : offset correction factor (0.4 -- 1)
 1.		NbOper           : Maximum number of operations for sorting
 
 ## output

Frame Error Rate for the given simulation parameters

## input and output files
 * Input File : 'FileMatrix' is an ASCII file with the parity-check matrix description in aList format.
 * Output File : in ./data file

## Simulation on windows

you may use CodeBlocks IDE
in the Projet-> Set Programs' arguments
2000 10 ./matrices/KN/N576_K480_GF64.txt 3.5 20 0.3 25

## Simulation on Linux

make

./essai 2000 10 ./matrices/KN/N576_K480_GF64.txt 3.5 20 0.3 25

you may use the start.sh shell script to run simulation for multiple snr
sh ./start.sh

## Simulation results

 Monte-Carlo simulation of Non-Binary LDPC decoder

Simulation parameters:

         NbMonteCarlo     : 2000
         NbIterMax        : 10
         FileMatrix       : ./matrices/KN/N576_K480_GF64.txt
         Eb/No (dB)       : 3.5
         NbMax            : 20
         Offset           : 0.3
         NbOper           : 25

 Normal alist format is used!
LDPC code parameters:
         N      :96
         K      :80
         M      :16
         CR     :0.833333
         GF     :64
         logGF  :6

         Note             : FB30
Simulation started at time: Wed Jan 08 17:34:02 2020

<0> FER= 40 / 751 = 0.053262 BER= 520 / x = 0.001443  avr_it=2.58
 results are printed in file ./data/results_N96_CR0.83_GF64_IT10_Offset0.3_nm20_FB30.txt

 results printed

Simulation complete at time: Wed Jan 08 17:34:08 2020

