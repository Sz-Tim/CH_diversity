#!/bin/bash

# Tim Szewczyk
#
# This script runs models through cmdstan: 
#    code/mods/*_fwd_LV
# Dataset is data/fwdSearch/4__*.Rdump
# Output is stored in out/fwdSearch/

# warmup and sampling per chain
nWarm=2000
nSamp=500

FILES=data/fwdSearch/LV_Y_4__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))


# Run Y
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
	  code/mods/Y_fwd_LV sample  
   	    algorithm=hmc  
   	      engine=nuts  
          metric=diag_e  
        num_samples=$nSamp  
        num_warmup=$nWarm  
        data file="${d[$f]}"  
        init=0  
        output file=$outj.csv  
        refresh=10 &
  done
done

wait


# Run WY after Y is done
FILES=data/fwdSearch/LV_WY_4__*.Rdump
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
	  code/mods/WY_fwd_LV sample  
   	    algorithm=hmc  
   	      engine=nuts  
          metric=diag_e  
        num_samples=$nSamp  
        num_warmup=$nWarm  
        data file="${d[$f]}"  
        init=0  
        output file=$outj.csv  
        refresh=10 &
  done
done

wait

