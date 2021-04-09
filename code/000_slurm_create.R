

sh_i <- tibble(folder=1:14,
               nScripts=rep(7:1, each=2))

for(i in sh_i$folder) {
  
  dir.create('code/02-slurm/fwd')
  # master scripts for submitting jobs
  sink(paste0("code/02-slurm/fwd_", i, "_cov.sh"))
  cat('#!/bin/bash', '\n',
      '', '\n',
      '# Slurm foward search master script', '\n',
      '# 02-', i, '\n',
      '# Covariate ', '\n',
      '', '\n',
      paste0(paste0('sbatch code/02-slurm/fwd/', i, '/'), 
             1:sh_i$nScripts[i], "_cov_slurm.sh\n"),
      sep=""
  )
  sink()
  
  sink(paste0("code/02-slurm/fwd_", i, "_LV.sh"))
  cat('#!/bin/bash', '\n',
      '', '\n',
      '# Slurm foward search master script', '\n',
      '# 02-', i, '\n',
      '# Latent variable ', '\n',
      '', '\n',
      paste0(paste0('sbatch code/02-slurm/fwd/', i, '/'), 
             1:sh_i$nScripts[i], "_LV_slurm.sh\n"),
      sep=""
  )
  sink()
  
  dir.create(paste0("code/02-slurm/fwd/", i), showWarnings=F)
  
  for(j in 1:sh_i$nScripts[i]) {
    
    # set data file indexes
    f_start <- (j-1)*8
    f_end <- j*8 - 1
    # if i is even and this is the end, reduce by 4
    if(!(i %% 2) && j == sh_i$nScripts[i]) {
      f_end <- f_end - 4
    }
    
    # 02-slurm/i/j_cov_slurm.sh
    sink(paste0("code/02-slurm/fwd/", i, "/", j, "_cov_slurm.sh"))
    cat('#!/bin/bash \n',
        '#SBATCH --ntasks=1 \n',
        '#SBATCH --cpus-per-task=24 \n',
        '#SBATCH --job-name="cov_', i, '-', j, '" \n',
        '#SBATCH --output=logs/cov_', i, '-', j, '.output \n',
        '\n',
        'module purge \n',
        'module load linuxbrew/colsa \n',
        '\n',
        'srun bash code/02-slurm/fwd/', i, '/', j, '_cov.sh',
        sep="") 
    sink()
    
    # 02-slurm/i/j_LV_slurm.sh
    sink(paste0("code/02-slurm/fwd/", i, "/", j, "_LV_slurm.sh"))
    cat('#!/bin/bash \n',
        '#SBATCH --ntasks=1 \n',
        '#SBATCH --cpus-per-task=24 \n',
        '#SBATCH --job-name="LV_', i, '-', j, '" \n',
        '#SBATCH --output=logs/LV_', i, '-', j, '.output \n',
        '\n',
        'module purge \n',
        'module load linuxbrew/colsa \n',
        '\n',
        'srun bash code/02-slurm/fwd/', i, '/', j, '_LV.sh',
        sep="") 
    sink()
    
    # 02-slurm/i/j_cov.sh
    sink(paste0("code/02-slurm/fwd/", i, "/", j, "_cov.sh"))
    cat('#!/bin/bash', '\n',
        '', '\n',
        '# Tim Szewczyk', '\n',
        '#', '\n',
        '# This script runs models through cmdstan: ', '\n',
        '#    code/mods/*_fwd_cov', '\n',
        '# Dataset is data/fwdSearch/', i, '__*.Rdump', '\n',
        '# Output is stored in out/fwdSearch/', '\n',
        '', '\n',
        '# warmup and sampling per chain', '\n',
        'nWarm=2000', '\n',
        'nSamp=500', '\n',
        '', '\n',
        'FILES=data/fwdSearch/', ifelse(i==1, '', 'cov_Y_'), i, '__*.Rdump', '\n',
        'd=($FILES)', '\n',
        'dfull=($(basename -a $FILES))', '\n',
        '', '\n',
        '', '\n',
        '# Run Y', '\n',
        'for f in {', f_start, '..', f_end, '}', '\n',
        'do ', '\n',
        '  df="${dfull[$f]}"', '\n',
        '  dname="${df%.Rdump}"', '\n',
        '  out="out/fwdSearch/', ifelse(i==1, 'cov_Y_"', '"'), '\n',
        '  out+=$dname', '\n',
        '  out+="_"', '\n',
        '  ', '\n',
        '  for j in {1..3}', '\n',
        '  do', '\n',
        '	  outj=$out', '\n',
        '	  outj+=$j', '\n',
        '	  code/mods/Y_fwd_cov sample ',
            'num_samples=$nSamp ', 
            'num_warmup=$nWarm ',
            'data file="${d[$f]}" ',
            'init=0 ',
            'output file=$outj.csv ',
            'refresh=10 &', '\n',
        '  done', '\n',
        'done', '\n',
        '', '\n',
        'wait', '\n',
        '', '\n',
        '', '\n',
        '# Run WY after Y is done', '\n',
        'FILES=data/fwdSearch/', ifelse(i==1, '', 'cov_WY_'), i, '__*.Rdump', '\n',
        'd=($FILES)', '\n',
        'dfull=($(basename -a $FILES))', '\n',
        '', '\n',
        '', '\n',
        'for f in {', f_start, '..', f_end, '}', '\n',
        'do ', '\n',
        '  df="${dfull[$f]}"', '\n',
        '  dname="${df%.Rdump}"', '\n',
        '  out="out/fwdSearch/', ifelse(i==1, 'cov_WY_"', '"'), '\n',
        '  out+=$dname', '\n',
        '  out+="_"', '\n',
        '  ', '\n',
        '  for j in {1..3}', '\n',
        '  do', '\n',
        '	  outj=$out', '\n',
        '	  outj+=$j', '\n',
        '	  code/mods/WY_fwd_cov sample ', 
              'num_samples=$nSamp ', 
              'num_warmup=$nWarm ', 
              'data file="${d[$f]}" ',
              'init=0 ', 
              'output file=$outj.csv ', 
              'refresh=10 &', '\n',
        '  done', '\n',
        'done', '\n',
        '', '\n',
        'wait', '\n',
        '', '\n',
        sep="")
    sink()
    
    # 02-slurm/i/j_LV.sh
    sink(paste0("code/02-slurm/fwd/", i, "/", j, "_LV.sh"))
    cat('#!/bin/bash', '\n',
        '', '\n',
        '# Tim Szewczyk', '\n',
        '#', '\n',
        '# This script runs models through cmdstan: ', '\n',
        '#    code/mods/*_fwd_LV', '\n',
        '# Dataset is data/fwdSearch/', i, '__*.Rdump', '\n',
        '# Output is stored in out/fwdSearch/', '\n',
        '', '\n',
        '# warmup and sampling per chain', '\n',
        'nWarm=2000', '\n',
        'nSamp=500', '\n',
        '', '\n',
        'FILES=data/fwdSearch/', ifelse(i==1, '', 'LV_Y_'), i, '__*.Rdump', '\n',
        'd=($FILES)', '\n',
        'dfull=($(basename -a $FILES))', '\n',
        '', '\n',
        '', '\n',
        '# Run Y', '\n',
        'for f in {', f_start, '..', f_end, '}', '\n',
        'do ', '\n',
        '  df="${dfull[$f]}"', '\n',
        '  dname="${df%.Rdump}"', '\n',
        '  out="out/fwdSearch/', ifelse(i==1, 'LV_Y_"', '"'), '\n',
        '  out+=$dname', '\n',
        '  out+="_"', '\n',
        '  ', '\n',
        '  for j in {1..3}', '\n',
        '  do', '\n',
        '	  outj=$out', '\n',
        '	  outj+=$j', '\n',
        '	  code/mods/Y_fwd_LV sample ', 
              'num_samples=$nSamp ',
              'num_warmup=$nWarm ', 
              'data file="${d[$f]}" ',
              'init=0 ', 
              'output file=$outj.csv ', 
              'refresh=10 &', '\n',
        '  done', '\n',
        'done', '\n',
        '', '\n',
        'wait', '\n',
        '', '\n',
        '', '\n',
        '# Run WY after Y is done', '\n',
        'FILES=data/fwdSearch/', ifelse(i==1, '', 'LV_WY_'), i, '__*.Rdump', '\n',
        'd=($FILES)', '\n',
        'dfull=($(basename -a $FILES))', '\n',
        '', '\n',
        '', '\n',
        'for f in {', f_start, '..', f_end, '}', '\n',
        'do ', '\n',
        '  df="${dfull[$f]}"', '\n',
        '  dname="${df%.Rdump}"', '\n',
        '  out="out/fwdSearch/', ifelse(i==1, 'LV_WY_"', '"'), '\n',
        '  out+=$dname', '\n',
        '  out+="_"', '\n',
        '  ', '\n',
        '  for j in {1..3}', '\n',
        '  do', '\n',
        '	  outj=$out', '\n',
        '	  outj+=$j', '\n',
        '	  code/mods/WY_fwd_LV sample ', 
              'num_samples=$nSamp ', 
              'num_warmup=$nWarm ', 
              'data file="${d[$f]}" ',
              'init=0 ',
              'output file=$outj.csv ', 
              'refresh=10 &', '\n',
        '  done', '\n',
        'done', '\n',
        '', '\n',
        'wait', '\n',
        '', '\n',
        sep="")
    sink()
  }
}
