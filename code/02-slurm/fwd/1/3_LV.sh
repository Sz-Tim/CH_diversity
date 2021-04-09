#!/bin/bash

# Tim Szewczyk
#
# This script runs models through cmdstan: 
#    code/mods/*_fwd_LV
# Dataset is data/fwdSearch/1__*.Rdump
# Output is stored in out/fwdSearch/

# warmup and sampling per chain
nWarm=2000
nSamp=500

FILES=data/fwdSearch/1__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))


# Run Y
for f in {16..23}
do 
  df="${dfull[$f]}"
  dname="${df%.Rdump}"
  out="out/fwdSearch/LV_Y_"
  out+=$dname
  out+="_"
  
  for j in {1..3}
  do
	  outj=$out
	  outj+=$j
	  code/mods/Y_fwd_LV sample num_samples=$nSamp num_warmup=$nWarm data file="${d[$f]}" init=0 output file=$outj.csv refresh=10 &
  done
done

wait


# Run WY after Y is done
FILES=data/fwdSearch/1__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))


for f in {16..23}
do 
  df="${dfull[$f]}"
  dname="${df%.Rdump}"
  out="out/fwdSearch/LV_WY_"
  out+=$dname
  out+="_"
  
  for j in {1..3}
  do
	  outj=$out
	  outj+=$j
	  code/mods/WY_fwd_LV sample num_samples=$nSamp num_warmup=$nWarm data file="${d[$f]}" init=0 output file=$outj.csv refresh=10 &
  done
done

wait

