# la commande d'exécution est:
# ./nom_du_fichier_binaire nombre_de_trame nombre_iteration nom_de_la_matrice EbN nm offset nombre_operation_dans_ecn
# script to lunch multiple SNR
# if edited in window, dos2linux file.sh
#! /bin/bash
echo "Bash version ${BASH_VERSION}..."
nb_frame=2000000000
it=30
matrix="./matrices/KN/N360_K120_GF64.txt " # Mat212_N480_M80  N96_K48_G64 or Mat24_N48_M24
n_vc=20
n_cv=20
offset=1.0
nb_oper=25

LANG="en_US.utf8" #define us language to have floating point with point with a point 

echo "start multiple decoding"

for SNR in $(seq 0.0 0.2 0.4)
do
  echo "simulation with SNR= $SNR "
  xterm -xrm '*hold:true' -e ./essai $nb_frame $it ${matrix} $SNR $n_vc $n_vc $offset $nb_oper &
done

