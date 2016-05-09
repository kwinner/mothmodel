INCS=-I$(CURDIR) -I/usr/local/include -I/usr/local/include/flint
LIBS=-L$(CURDIR) -L/usr/local/lib -L/usr/local/lib -lflint -lmpfr -lgmp -lm -lpthread 

# CFLAGS=-fno-common -ansi -pedantic -Wall -O2 -funroll-loops -g -mpopcnt 
CFLAGS=

all:: main

main:
	gcc $(CFLAGS) $(INCS) pgfba.c $(LIBS) -o pgfba.out