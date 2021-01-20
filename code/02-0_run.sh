#!/bin/bash

# 02-0_run.sh
# Tim Szewczyk
#
# This script runs models with no covariates through cmdstan: 
#    code/mods/Y_fwdSearch
#    code/mods/WY_fwdSearch
# Dataset is data/fwdSearch/0__*.Rdump
# Output is stored in out/fwdSearch/

# warmup and sampling per chain
nWarm=2000
nSamp=1000

FILES=data/fwdSearch/0__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))

for f in {0..3} 
do
  df="${dfull[$f]}"
  dname="${df%.Rdump}"
  
  # Y
  out="out/fwdSearch/Y_"
  out+=$dname
  out+="_"
  
  for j in {1..3}
  do
    outj=$out
    outj+=$j
    code/mods/Y_fwdSearch sample \
     	  algorithm=hmc \
     	    engine=nuts \
          metric=diag_e \
        num_samples=$nSamp \
        num_warmup=$nWarm \
        data file="${d[$f]}" \
        init=0 \
        output file=$outj.csv \
        refresh=10 &
  done
  
  
  # WY
  out="out/fwdSearch/WY_"
  out+=$dname
  out+="_"
  
  for j in {1..3}
  do
    outj=$out
    outj+=$j
    code/mods/WY_fwdSearch sample \
     	  algorithm=hmc \
    	    engine=nuts \
          metric=diag_e \
        num_samples=$nSamp \
        num_warmup=$nWarm \
        data file="${d[$f]}" \
        init=0 \
        output file=$outj.csv \
        refresh=10 &
  done
done

wait
