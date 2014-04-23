#!/bin/bash
module load xl
mpixlc_r -g -qasm -qtm -lpthread main.c -o tm_test
if [ $? -eq 0 ]
  then srun --partition=small --time=15 --nodes=3 ./tm_test
fi
