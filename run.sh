#!/bin/bash
module load xl
mpicc -g -qasm -qtm -lpthread main.c -o tm_test
if [ $? -eq 0 ]
  then srun --partition=small --time=15 --nodes=3 ./tm_test
fi
