#!/bin/sh
#PBS -l nodes=1,walltime=05:00:00
#PBS -q mcpd
#PBS -d .

for (( k=256; k<=8192; k=k*2 ));
do
	for (( i=16; i<=256; i=i*2 ))
	do
	     ./original $k $i >> tiempoOriginal_$k.txt
	    for (( j=2; j<=8; j=j+1 ))
		do
			./chol $k $i $j >> tiempo_$k_$j.txt
		done
	done
done
