#!/bin/bash

# This script runs models with all covariates through cmdstan: 
#    code/mods/full_Y 
#    code/mods/full_WY
# Dataset is data/stan_data/vs_no_pred.Rdump
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
	outj="out/vs_Y_"
	outj+=$j
	code/mods/full_Y sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/stan_data/vs_no_pred.Rdump \
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
	outj="out/vs_WY_"
	outj+=$j
	code/mods/full_WY sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/stan_data/vs_no_pred.Rdump \
      init=0 \
      output file=$outj.csv \
      refresh=1 &
done


wait
