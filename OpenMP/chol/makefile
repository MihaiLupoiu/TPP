CFLAGS  = -O
CC       = cc
LIBS     = -lm -llapack -lblas /labos/alumnos/milu/TPP/OpenMP/chol/quark-0.9.0/libquark.a -lhwloc -lpthread

chol: chol.o
	$(CC) -o $@ chol.o $(LDFLAGS) $(LIBS)

original: chol_origianl.o
	$(CC) -o $@ chol_origianl.o $(LDFLAGS) $(LIBS)

