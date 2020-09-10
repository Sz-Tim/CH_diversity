#!/bin/bash

# 02-2_run_WY_8-14.sh
# Tim Szewczyk
#
# This script runs models with 2 covariates through cmdstan: 
#    code/mods/WY_fwdSearch
# Dataset is data/fwdSearch/WY_2__*.Rdump
# Output is stored in out/fwdSearch/

# warmup and sampling per chain
nWarm=2000
nSamp=70

FILES=data/fwdSearch/WY_2__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))

for f in {8..14}
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
