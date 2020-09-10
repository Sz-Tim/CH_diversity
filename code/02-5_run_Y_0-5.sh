#!/bin/bash

# 02-5_run_Y_0-5.sh
# Tim Szewczyk
#
# This script runs models with 5 covariates through cmdstan: 
#    code/mods/Y_fwdSearch
# Dataset is data/fwdSearch/Y_5__*.Rdump
# Output is stored in out/fwdSearch/

# warmup and sampling per chain
nWarm=2000
nSamp=70

FILES=data/fwdSearch/Y_5__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))

for f in {0..5}
do 
  df="${dfull[$f]}"
  dname="${df%.Rdump}"
  out="out/fwdSearch/"
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
done

wait

