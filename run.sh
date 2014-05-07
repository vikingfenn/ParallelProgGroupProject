#!/bin/bash

echo -n "Compiling..."

module load xl
mpixlc_r -g -qasm -qtm -lpthread main.c -o tm_test

if [ $? -ne 0 ]
then
    echo "Compilation error!"
    exit 1
fi

echo " Done"

export TM_REPORT_ENABLE=YES
export TM_REPORT_LOG=SUMMARY

if [ $? -eq 0 ]
then

    for THREADS in 2 4 8 16 32; do

        rm -rf $THREADS
        mkdir $THREADS
        cd $THREADS

        echo -n "Running with $THREADS threads..."

        srun --partition=small --time=10 --runjob-opts="--mapping TEDCBA" --nodes=63 ../tm_test $THREADS >> results.txt

        if [ $? -ne 0 ]; then
            echo "Error!"
            exit 1
        fi

        echo " Done"

        for ((LOG=0; LOG<=20; LOG++)); do
            let "TRANSACTIONS += $(sed -n 452p tm_report.log.$LOG | sed 's/[^0-9]//g')"
            let "ROLLBACKS += $(sed -n 451p tm_report.log.$LOG | sed 's/[^0-9]//g')"
        done

        echo -e "Total transactions:\t$TRANSACTIONS" >> totals.txt
        echo -e "Total rollbacks:\t$ROLLBACKS" >> totals.txt

        cd ..

    done
fi
