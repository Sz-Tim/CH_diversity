#!/bin/bash

# 04_run_opt.sh
# Tim Szewczyk
#
# This script runs models the optimal covariate sets through cmdstan: 
#    code/mods/Y_opt
#    code/mods/WY_opt
# Dataset is data/opt/*_null__opt_var_set.Rdump
# Output is stored in out/opt/

# warmup and sampling per chain
nWarm=2000
nSamp=100

FILES=data/opt/cov_Y_null__opt_var_set.Rdump
d=($FILES)
dfull=($(basename -a $FILES))

f=0
df="${dfull[$f]}"
dname="${df%.Rdump}"
out="out/opt/"
out+=$dname
out+="_"
  
for j in {1..12}
do
  outj=$out
  outj+=$j
  code/mods/Y_opt sample \
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


FILES=data/opt/cov_WY_null__opt_var_set.Rdump
d=($FILES)
dfull=($(basename -a $FILES))

f=0
df="${dfull[$f]}"
dname="${df%.Rdump}"
out="out/opt/"
out+=$dname
out+="_"
  
for j in {1..12}
do
  outj=$out
  outj+=$j
  code/mods/WY_opt sample \
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

wait
