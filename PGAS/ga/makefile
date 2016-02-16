CFLAGS  = -I/scratch/jroman/ga-5-4/include -O3 -DMPI
CC      = mpicc
LDFLAGS = -L/scratch/jroman/ga-5-4/lib
LIBS    = -lga -larmci -llapack -lblas

hello_world: hello_world.o
	$(CC) -o $@ $< $(LDFLAGS) $(LIBS)

matrix: matrix.o
	$(CC) -o $@ $< $(LDFLAGS) $(LIBS)

transp1D: transp1D.o
	$(CC) -o $@ $< $(LDFLAGS) $(LIBS)

transp1Dnb: transp1Dnb.o
	$(CC) -o $@ $< $(LDFLAGS) $(LIBS)

