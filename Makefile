CC=g++
OPTIONFLAG=-Wall -std=c++11
EXEC=main

all:compile run

compile:
	@${CC} *.cpp -fopenmp -o ${EXEC} ${OPTIONFLAG}
run:
	@./${EXEC} 	
