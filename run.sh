#!/bin/bash
module load xl
mpixlc_r -g -qasm -qtm -lpthread main.c -o tm_test
export TM_REPORT_ENABLE=YES
export TM_REPORT_LOG=SUMMARY
if [ $? -eq 0 ]
  then
  for NODES in 64 128 256; do
    PARTITION="small"
	if [ $NODES -gt 64 ]; then
		PARTITION="medium"
	fi
    rm -r $NODES
    mkdir $NODES
    cd $NODES
    srun --partition=$PARTITION --time=5 --runjob-opts="--mapping TEDCBA" --nodes=$NODES ../tm_test >> results.txt
    cd ..
	if [ $? -ne 0 ]; then
		echo "Error!"
		exit 0
	fi
  done
fi
