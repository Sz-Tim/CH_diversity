#!/bin/bash

# Tim Szewczyk
#
# This script runs models through cmdstan: 
#    code/mods/*_fwd_cov
# Dataset is data/fwdSearch/7__*.Rdump
# Output is stored in out/fwdSearch/

# warmup and sampling per chain
nWarm=2000
nSamp=500

FILES=data/fwdSearch/cov_Y_7__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))


# Run Y
for f in {8..15}
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
	  code/mods/Y_fwd_cov sample num_samples=$nSamp num_warmup=$nWarm data file="${d[$f]}" init=0 output file=$outj.csv refresh=10 &
  done
done

wait


# Run WY after Y is done
FILES=data/fwdSearch/cov_WY_7__*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))


for f in {8..15}
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
	  code/mods/WY_fwd_cov sample num_samples=$nSamp num_warmup=$nWarm data file="${d[$f]}" init=0 output file=$outj.csv refresh=10 &
  done
done

wait

