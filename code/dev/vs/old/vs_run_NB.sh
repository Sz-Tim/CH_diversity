#!/bin/bash

nWarm=1500
nSamp=20

echo "Starting W_NB"
for i in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $i
	echo " "
	outi="out/vs/W_NB_"
	outi+=$i
	code/vs/vs_W_NB sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/vs/vs_30.Rdump \
      init=0 \
      output file=$outi.csv \
      refresh=1 &
done


echo "Starting Y_NB"
for j in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/vs/Y_NB_"
	outj+=$j
	code/vs/vs_Y_NB sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/vs/vs_30.Rdump \
      init=0 \
      output file=$outj.csv \
      refresh=1 &
done


echo "Starting WY_NB"
for k in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $k
	echo " "
	outk="out/vs/WY_NB_"
	outk+=$k
	code/vs/vs_WY_NB sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/vs/vs_30.Rdump \
      init=0 \
      output file=$outk.csv \
      refresh=1 &
done


wait
