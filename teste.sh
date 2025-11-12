#!/bin/bash

PROG=cgSolver
CPU=3

DATA_DIR="Dados/"

mkdir -p ${DATA_DIR}

make purge
make

INPUTS=("SEM_1.in" "COM_1.in" "SEM_2.in" "COM_2.in" "SEM_3.in" "COM_3.in" "SEM_4.in" "COM_4.in" "SEM_5.in" "COM_5.in" "SEM_6.in" "COM_6.in")

    
LIKWID_CSV="${DATA_DIR}/${m}.csv"
    rm -f ${TEMPOS}

for input_file in "${INPUTS[@]}"; do
    echo -e "teste com $input_file"

    output_file="${input_file%.*}.out"

    ./"$PROG" <  "$DATA_DIR""$input_file" > "$DATA_DIR""$output_file"
    expected_file="${input_file%.*}"
    #diff <(head -n 2 "$DATA_DIR""$output_file") "$DATA_DIR""$expected_file"

done

