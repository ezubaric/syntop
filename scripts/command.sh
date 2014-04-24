#!/bin/bash

if [ $# -ne 7 ]
then
    echo "USAGE: command.sh dataset doc_limit topic_num vocab_size suffix dir online[true/false]"
    exit
fi

# wackypedia doc 983 voc 5050
# toy doc 20 voc 45

DATASET=$1
DOC_LIMIT=$2
TOPIC_NUM=$3
VOCAB_SIZE=$4
SUFFIX=$5
DIR=$6
ONLINE=$7

LOCAL="true"
RESUME="false"

ITERATIONS=50

# ignore trans = latent dirichlet allocation

FINITE="false"
IGNORE_TRANS="false"
IGNORE_DOCS="false"
SHORTCUT_GSL="false"

#make clean

if [ $LOCAL = "true" ]; then
    make simulation
else
	echo "Hadoop version is not implemented"
fi

DATASET_NAME=$DATASET
GLOBAL_DIR=$DIR

for chunk in $DIR/$DATASET_NAME*; 
do	
	echo "========================================================"
	echo $chunk
	echo "========================================================"
	
	DATASET=$(basename $chunk)
	DIR=$GLOBAL_DIR/$DATASET
	PARAMS_DIR=$DIR/$DATASET\_params_$SUFFIX
	
	if [ $ONLINE = 'true' ];
	then
		GLOBAL_PARAMS_DIR=$GLOBAL_DIR/global_${DATASET_NAME}_params_$SUFFIX	
		# copying global parameters to data chunk directory
		echo "cp $GLOBAL_PARAMS_DIR/* $PARAMS_DIR/"
		cp $GLOBAL_PARAMS_DIR/0.nu $PARAMS_DIR/current.nu
		cp $GLOBAL_PARAMS_DIR/0.tau $PARAMS_DIR/current.tau		

		# creating current local paramenters
		cp $PARAMS_DIR/0.beta $PARAMS_DIR/current.beta
		cp $PARAMS_DIR/0.gamma $PARAMS_DIR/current.gamma
	else 
		cp $PARAMS_DIR/0.nu $PARAMS_DIR/current.nu
		cp $PARAMS_DIR/0.tau $PARAMS_DIR/current.tau
		cp $PARAMS_DIR/0.beta $PARAMS_DIR/current.beta
		cp $PARAMS_DIR/0.gamma $PARAMS_DIR/current.gamma						
	fi
		
	TMP_DIR=$DIR/tmp-output_$SUFFIX	
	rm -rf $TMP_DIR
    mkdir -p $TMP_DIR
	
	echo "running dataset $DATASET"
	
	for ((i=1;i<=ITERATIONS;i++));
	do
	    echo "Iteration $i"
        echo "--truncation=$TOPIC_NUM --doc_limit=$DOC_LIMIT --vocab_size=$VOCAB_SIZE --directory=$DIR --dataset=$DATASET --suffix=$SUFFIX --ignore_doc=$IGNORE_DOCS --ignore_trans=$IGNORE_TRANS --shortcut_gsl=$SHORTCUT_GSL --finite=$FINITE > $TMP_DIR/output$i.txt"
		./simulation --truncation=$TOPIC_NUM --doc_limit=$DOC_LIMIT --vocab_size=$VOCAB_SIZE --directory=$DIR --dataset=$DATASET --suffix=$SUFFIX --ignore_doc=$IGNORE_DOCS --ignore_trans=$IGNORE_TRANS --shortcut_gsl=$SHORTCUT_GSL --finite=$FINITE > $TMP_DIR/output$i.txt

        echo "scripts/merge.py $TOPIC_NUM $DOC_LIMIT $VOCAB_SIZE $i $DATASET"
		python scripts/merge.py $TOPIC_NUM $DOC_LIMIT $VOCAB_SIZE $i $DIR/$DATASET $SUFFIX
        cp $PARAMS_DIR/tmp-lhood.txt $PARAMS_DIR/$i.lhood

		echo "========================================================"
		echo "LIKELIHOOD IS " $i $(more $PARAMS_DIR/tmp-lhood.txt)

	    cp $PARAMS_DIR/"$i.nu" $PARAMS_DIR/current.nu
	    cp $PARAMS_DIR/"$i.tau" $PARAMS_DIR/current.tau
	    
	    cp $PARAMS_DIR/"$i.gamma" $PARAMS_DIR/current.gamma
	    cp $PARAMS_DIR/current.beta $PARAMS_DIR/"$i.beta"

		mkdir -p $PARAMS_DIR/temp-$i
		cp $PARAMS_DIR/tmp-* $PARAMS_DIR/temp-$i
		sort $PARAMS_DIR/tmp-tokens.txt > $PARAMS_DIR/temp-$i/tmp-tokens.txt
				
	done
	
	if [ $ONLINE = "true" ]; 
	then
		#interpolation
		echo "python scripts/interpolate.py $GLOBAL_PARAMS_DIR/current.tau $PARAMS_DIR/current.tau 0.1 $GLOBAL_PARAMS_DIR/current.tau true"		
		python scripts/interpolate.py $GLOBAL_PARAMS_DIR/current.tau $PARAMS_DIR/current.tau 0.1 $GLOBAL_PARAMS_DIR/current.tau true
		
		echo "scripts/interpolate.py $GLOBAL_PARAMS_DIR/current.nu $PARAMS_DIR/current.nu 0.1 $GLOBAL_PARAMS_DIR/current.nu false"		
		python scripts/interpolate.py $GLOBAL_PARAMS_DIR/current.nu $PARAMS_DIR/current.nu 0.1 $GLOBAL_PARAMS_DIR/current.nu false

	else
		break
	fi
done

echo "python scripts/map_words.py $PARAMS_DIR/current.tau $DIR\/$DATASET\.voc $PARAMS_DIR/tau.map"
python scripts/map_words.py $PARAMS_DIR/current.tau $DIR\/$DATASET\.voc $PARAMS_DIR/tau.map

echo "python scripts/visualize_topics.py $DATASET\/$DATASET\.voc $DATASET\_params/current $TOPIC_NUM 10 .1"
python scripts/visualize_topics.py $DIR/$DATASET\.voc $PARAMS_DIR/current $TOPIC_NUM 10 .2

#echo "python scripts/visualize_topics.py $DATASET\/$DATASET\.voc $DATASET\_params/current $TOPIC_NUM 10 .2"
#python scripts/visualize_topics.py $DATASET\/$DATASET\.voc $DATASET\_params/current $TOPIC_NUM 10 .2
#python map_words.py $DATASET\_params/current.tau $DATASET\/$DATASET\.voc $DATASET\_params/tau.map

python scripts/merge_lhood.py $PARAMS_DIR
echo "set terminal png; set output \"$PARAMS_DIR/lhood.png\"; plot \"$PARAMS_DIR/current.lhood\" with line" | gnuplot
# global likelihood plot
echo "set terminal png; set output \"$PARAMS_DIR/global_lhood.png\"; plot \"$PARAMS_DIR/current.global_lhood\" with line" | gnuplot
# local likelihood plot
echo "set terminal png; set output \"$PARAMS_DIR/local_lhood.png\"; plot \"$PARAMS_DIR/current.local_lhood\" with line" | gnuplot


#./syntop2 --truncation=5 --doc_limit=20 --vocab_size=45 --directory=toy1 --ignore_doc=false --ignore_trans=false --shortcut_gsl=false --finite=true --dataset=toy1
