#!/bin/bash

# 02_runVS_altRE.sh
# Tim Szewczyk
#
# This script runs models with all covariates through cmdstan: 
#    code/mods/full_Y_altRE 
#    code/mods/full_WY_altRE
# Dataset is data/stan_data/vs_no_pred.Rdump
# Output is stored in out/

# warmup and sampling per chain
nWarm=100
nSamp=20

echo "Starting Y"
for j in {1..12}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/vs_Y_altRE_G_"
	outj+=$j
	code/mods/full_Y_altRE_G sample \
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
for j in {1..12}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/vs_WY_altRE_G_"
	outj+=$j
	code/mods/full_WY_altRE_G sample \
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
