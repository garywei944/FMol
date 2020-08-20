#!/bin/bash

TARGET=$1
TARGET_FILE=$2
MODEL_DIR=$3
OUTPUT_DIR="${TARGET}_out"

echo -e "Saving output to ${OUTPUT_DIR}/\n"

for replica in 0 1 2 3; do
	if [ $replica -eq 0 ]; then mode="D B T"; else mode="D B"; fi
	for m in $mode; do
		echo "Launching model: ${m} ${replica}"
		python alphafold.py -i $TARGET_FILE -o $OUTPUT_DIR -m $MODEL_DIR -r $replica -t $m | cat &
	done
#	wait
done

echo -e "All models running, waiting for them to complete\n"
wait

echo "Ensembling all replica outputs & Pasting contact maps"
python alphafold.py -i $TARGET_FILE -o $OUTPUT_DIR -e
