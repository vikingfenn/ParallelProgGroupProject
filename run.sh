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
    srun --partition=$PARTITION --time=1 --runjob-opts="--mapping TEDCBA" --nodes=$NODES ../tm_test >> results.txt
    let "LOGS = $NODES / 3 - 1"
    for ((LOG=0; LOG<=$LOGS; LOG++)); do
        let "ROLLBACKS += $(sed -n 451p tm_report.log.$LOG | sed 's/[^0-9]//g')"
        let "TRANSACTIONS += $(sed -n 452p tm_report.log.$LOG | sed 's/[^0-9]//g')"
    done
    let "ROLLBACKS /= $LOGS"
    let "TRANSACTIONS /= $LOGS"
    echo -e "Average rollbacks:\t\t$ROLLBACKS" >> averages.txt
    echo -e "Average transactions:\t$TRANSACTIONS" >> averages.txt
    cd ..
    if [ $? -ne 0 ]; then
        echo "Error!"
        exit 0
    fi
  done
fi
