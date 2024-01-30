#!/bin/bash
omp="omp"
if [ $1 = $omp ]
then
gcc -Wall -g main.c -o mainOMP -lm -lgomp -fopenmp -std=c11
else
gcc -Wall -g main.c -o mainProc -lm -lgomp -std=c11
fi
