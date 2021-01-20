#!/bin/bash

nWarm=15
nSamp=10

# echo "Starting W_GP"
# for i in {1..8}
# do
# 	echo " "
# 	echo "||||---------- Starting chain" $i
# 	echo " "
# 	outi="out/vs/W_GP_"
# 	outi+=$i
# 	code/vs/vs_W_GP sample \
#    	  algorithm=hmc \
#    	    engine=nuts \
#         metric=diag_e \
#       num_samples=$nSamp \
#       num_warmup=$nWarm \
#       data file=data/vs/vs_30.Rdump \
#       init=0 \
#       output file=$outi.csv \
#       refresh=1 &
# done


echo "Starting Y_GP"
for j in {1..2}
do
	echo " "
	echo "||||---------- Starting chain" $j
	echo " "
	outj="out/full/Y_GP_"
	outj+=$j
	code/vs/full_Y_GP sample \
   	  algorithm=hmc \
   	    engine=nuts \
        metric=diag_e \
      num_samples=$nSamp \
      num_warmup=$nWarm \
      data file=data/stan_data/full.Rdump \
      init=0 \
      output file=$outj.csv \
      refresh=1 &
done


# echo "Starting WY_GP"
# for k in {1..8}
# do
# 	echo " "
# 	echo "||||---------- Starting chain" $k
# 	echo " "
# 	outk="out/vs/WY_GP_"
# 	outk+=$k
# 	code/vs/vs_WY_GP sample \
#    	  algorithm=hmc \
#    	    engine=nuts \
#         metric=diag_e \
#       num_samples=$nSamp \
#       num_warmup=$nWarm \
#       data file=data/vs/vs_30.Rdump \
#       init=0 \
#       output file=$outk.csv \
#       refresh=1 &
# done


wait
