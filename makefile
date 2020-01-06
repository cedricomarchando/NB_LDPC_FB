#makefile
all: essai

############ Generation des fichiers objets ###############

#init.o
obj/init.o: init.c include/init.h include/tools.h include/struct.h
	gcc -O3 -Wall -Iinclude/ -c init.c -o obj/init.o

#tools.o
obj/tools.o: tools.c include/tools.h include/init.h include/struct.h
	gcc -O3 -Wall -Iinclude/ -c tools.c -o obj/tools.o

#channel.o
obj/channel.o: channel.c include/channel.h include/init.h include/tools.h include/struct.h
	gcc -O3 -Wall -Iinclude/ -c channel.c -o obj/channel.o

#bubble_decoder.o
obj/bubble_decoder.o: bubble_decoder.c include/bubble_decoder.h include/init.h include/tools.h include/struct.h
	gcc -O3 -Wall -Iinclude/ -c bubble_decoder.c -o obj/bubble_decoder.o

#syndrome_decoder.o
obj/syndrome_decoder.o: syndrome_decoder.c include/syndrome_decoder.h include/init.h include/tools.h include/struct.h
	gcc -O3 -Wall -Iinclude/ -c syndrome_decoder.c -o obj/syndrome_decoder.o

#NB_LDPC.o
obj/NB_LDPC.o: NB_LDPC.c include/syndrome_decoder.h include/bubble_decoder.h include/init.h include/tools.h include/channel.h include/struct.h include/NB_LDPC.h
	gcc -O3 -Wall -Iinclude/ -c NB_LDPC.c -o obj/NB_LDPC.o


########## Generation de l'executable #######################
essai: obj/init.o obj/tools.o obj/channel.o obj/bubble_decoder.o obj/syndrome_decoder.o obj/NB_LDPC.o
	gcc  -O3 -Wall -lm -o essai obj/init.o obj/tools.o obj/channel.o obj/bubble_decoder.o obj/syndrome_decoder.o obj/NB_LDPC.o

#fin makefile
