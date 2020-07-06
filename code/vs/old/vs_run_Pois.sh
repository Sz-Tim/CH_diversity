#!/bin/bash

nWarm=1500
nSamp=20

echo "Starting W_Pois"
for i in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $i
	echo " "
	outi="out/vs/W_Pois_"
	outi+=$i
        code/vs/vs_W_Pois sample \
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


echo "Starting Y_Pois"
for j in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/vs/Y_Pois_"
	outj+=$j
	code/vs/vs_Y_Pois sample \
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


echo "Starting WY_Pois"
for k in {1..8}
do
	echo " "
	echo "||||---------- Starting chain" $k
	echo " "
	outk="out/vs/WY_Pois_"
	outk+=$k
	code/vs/vs_WY_Pois sample \
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
