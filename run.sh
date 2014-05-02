#!/bin/bash
module load xl
mpixlc_r -g -qasm -qtm -lpthread main.c -o tm_test
if [ $? -eq 0 ]
  then
  for NODES in 63 510 1023
	if [ $NODES -lt 64 ]; then
		PARTITION="small"
	elif [ $NODES -lt 512 ]; then
		PARTITION="medium"
	else
		PARTITION="large"
	fi
        srun --partition=$PARTITION --time=15 --overcommit --runjob-opts="--mapping TEDCBA" --nodes=$NODES --ntasks-per-node=16 ./tm_test >> results.txt
	if [ $? -ne 0 ] 
	then
		echo "Error!"
		exit 0
	fi
done
fi
