#!/bin/bash

# 02-5-1_Y.sh
# Tim Szewczyk
#
# This script runs models with one covariate through cmdstan: 
#    code/mods/Y_fwdSearch
# Dataset is data/fwdSearch/Y_5__*.Rdump
# Output is stored in out/fwdSearch/

# warmup and sampling per chain
nWarm=2000
nSamp=500

FILES=data/fwdSearch/Y_5__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))

for f in {0..7}
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
