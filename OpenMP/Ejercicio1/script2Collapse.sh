#!/bin/sh
#PBS -l nodes=1,walltime=01:00:00
#PBS -q mcpd
#PBS -d .

for (( j=800; j <= 1600; j=j+200 ))
do
	for (( i=2; i <= 4; i=i+1 ))
	do
		OMP_NUM_THREADS=$i ./pmatmat_v4 $j >> tiempoV4.txt
	done
done