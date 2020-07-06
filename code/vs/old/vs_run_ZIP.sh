#!/bin/bash

nWarm=1500
nSamp=20

echo "Starting W_ZIP"
for i in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $i
	echo " "
	outi="out/vs/W_ZIP_"
	outi+=$i
	code/vs/vs_W_ZIP sample \
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


echo "Starting Y_ZIP"
for j in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/vs/Y_ZIP_"
	outj+=$j
	code/vs/vs_Y_ZIP sample \
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


echo "Starting WY_ZIP"
for k in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $k
	echo " "
	outk="out/vs/WY_ZIP_"
	outk+=$k
	code/vs/vs_WY_ZIP sample \
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
