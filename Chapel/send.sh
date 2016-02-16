#!/bin/sh
#PBS -l nodes=1,walltime=01:00:00
#PBS -q mcpd
#PBS -d .
#PBS -N jacobi_milu
for t in 1 2 4 8 16 32; do
	export CHPL_RT_NUM_THREADS_PER_LOCALE=$t
	echo "--" $t "threads"
	./a.out --n=600 --epsilon=1e-4
done