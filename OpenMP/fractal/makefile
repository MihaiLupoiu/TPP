CFLAGS  = -O3
CC       = cc
LIBS     = -lm 
OMP		 = -fopenmp

#mandel: mandel.o misc.o
#	$(CC) -o $@ mandel.o misc.o $(LDFLAGS) $(LIBS)


all:
	$(CC) -o mandel mandel.c misc.c $(LDFLAGS) $(LIBS) $(OMP)
	$(CC) -o mandelOriginal mandelOriginal.c misc.c $(LDFLAGS) $(LIBS)
	
clean:
	rm mandel mandelOriginal