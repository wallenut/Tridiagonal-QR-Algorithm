CC=gcc
CFLAGS= -std=c99 -pedantic -Wall -g3

#####
# Instructions to make hw3
#####

hw3: hw3.o
	${CC} ${CFLAGS} -o hw3 hw3.o
