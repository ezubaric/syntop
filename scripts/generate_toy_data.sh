#!/bin/bash

#dir dataset_name suffix num_topic num_doc num_sen 
DIR=$1
DATASET=$2
SUFFIX=$3
NUM_TOPIC=$4
NUM_DOC=$5
NUM_SEN=$6

#predefined number of words in toy data
NUM_VOCAB=45

#directory for parameters files
PARAMS_DIR=$DIR/$DATASET\_params_$SUFFIX


if [ $# -ne 6 ]
then
	#./generate_toy_data.sh toy 0 5 20 45 10
    echo "USAGE: generate_toy_data.sh dir dataset_name suffix num_topic num_doc num_sen"
    exit
fi

echo "creating dataset directory: $DIR"
echo "creating parameters directory: $PARAMS_DIR"
mkdir -p $PARAMS_DIR

#python toy/toy 20 10
python scripts/create_synthetic.py $DIR/$DATASET $NUM_DOC $NUM_SEN

#python 5 20 45 synthetic_data 
python scripts/seed.py $NUM_TOPIC $NUM_DOC $NUM_VOCAB $PARAMS_DIR/ $DIR/$DATASET.0.prs


cp $PARAMS_DIR/beta.origin $PARAMS_DIR/0.beta
cp $PARAMS_DIR/nu.origin $PARAMS_DIR/0.nu
cp $PARAMS_DIR/tau.origin $PARAMS_DIR/0.tau
cp $PARAMS_DIR/gamma.origin $PARAMS_DIR/0.gamma

cp $PARAMS_DIR/0.beta $PARAMS_DIR/current.beta
cp $PARAMS_DIR/0.nu $PARAMS_DIR/current.nu
cp $PARAMS_DIR/0.tau $PARAMS_DIR/current.tau
cp $PARAMS_DIR/0.gamma $PARAMS_DIR/current.gamma

