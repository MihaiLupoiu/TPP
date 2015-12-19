CC=nvcc
CLINKER=${CC}
#CFLAGS=-pg
#LDFLAGS=-pg
LDFLAGS=-L/usr/lib64/nvidia

default: matmult

OBJS = matmult.o errors.o

matmult: ${OBJS}
	${CLINKER} -o $@ ${LDFLAGS} ${OBJS} -lOpenCL
	${RM} ${OBJS}
