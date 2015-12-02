#!/bin/sh
#PBS -l nodes=1,walltime=03:00:00
#PBS -q mcpd
#PBS -d .

## declare an array variable
declare -a array=("static" "dynamic")
# get length of an array
arraylength=${#array[@]}

echo "Secuencial"
./mandelOriginal -width 4096 -height 2169 -ulx -1.37 -uly 0.04 -lly -0.01 >> tiempoSecuencial.txt
./mandelOriginal -width 4096 -height 2169 -ulx -1.37 -uly 0.014 -lly 0.0135 >> tiempoSecuencial.txt
./mandelOriginal -width 4096 -height 2169 -ulx -1.369538 -uly 0.013797 -lly 0.01368 >> tiempoSecuencial.txt

for (( i=0; i<${arraylength}; i++ ));
do
	echo "OpenMP " ${array[$i]}

	for (( j=2; j <= 8; j=j+1 ))
	do
		echo "Threads " $j	
			
		OMP_NUM_THREADS=$j  OMP_SCHEDULE=${array[$i]} ./mandel -width 4096 -height 2169 -ulx -1.37 -uly 0.04 -lly -0.01 >> tiemposOpenMP_v1_${array[$i]}.txt
		OMP_NUM_THREADS=$j  OMP_SCHEDULE=${array[$i]} ./mandel -width 4096 -height 2169 -ulx -1.37 -uly 0.014 -lly 0.0135 >> tiemposOpenMP_v2_${array[$i]}.txt
		OMP_NUM_THREADS=$j  OMP_SCHEDULE=${array[$i]} ./mandel -width 4096 -height 2169 -ulx -1.369538 -uly 0.013797 -lly 0.01368 >> tiemposOpenMP_v3_${array[$i]}.txt
	done

done