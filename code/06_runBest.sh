#!/bin/bash

# 05_runBest.sh
# Tim Szewczyk
#
# This script runs models with the best covariates through cmdstan: 
#    code/mods/best_Y 
#    code/mods/best_WY
# Dataset is data/stan_data/pred.Rdump
# Output is stored in out/

# warmup and sampling per chain
nWarm=100
nSamp=20
nChain=12

echo "Starting Y"
for j in {1..$nChain}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/best_Y_"
	outj+=$j
	code/mods/best_Y sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/stan_data/pred.Rdump \
      init=0 \
      output file=$outj.csv \
      refresh=1 &
done


echo "Starting WY"
for j in {1..$nChain}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/best_WY_"
	outj+=$j
	code/mods/best_WY sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/stan_data/pred.Rdump \
      init=0 \
      output file=$outj.csv \
      refresh=1 &
done


wait
