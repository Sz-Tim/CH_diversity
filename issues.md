# Issues to resolve  

- Any method to control for the same colony being sampled twice at different flags from the same plot? 
  - Using data as-is would over-estimate colony density when that happens, but lumping to presence/absence at the plot-level would reduce the alignment of a poisson-distributed variable, as well as underestimate the abundance of species with smaller colonies (e.g., multiple Lepto/Temno in separate branches).  
  - what about linking a presence/absence function to the PPM lambda?
    - essentially turning the PPM into a method for estimating pr(Pres)
    - Y.ijs>0 ~ Binom(1-exp(lambda.ijs))
  - better to assume data represent independent colonies I think
  - SOLVED: Assume colonies are actually separate colonies

- NPP is NA (=maxValue) in a lot of cells where many W observations occur
  - SOLVED: Occurs in cities where NPP is negligible; recoded to 0
  
  



Next:
- how to assess / compare for variable selection
  - SOLVED: projpred
- how to compare W Y WY 
  - I think it makes more sense to only compare Y & WY
  - what does citizen science data bring?
  - COVARIATES
    - included in optimal model using Y, WY as refmodel
    - magnitude of slope estimates
- output quantities
  - Shannon H
  - beta
  - prP calculated each iteration directly as pr(y>0 | lambda, LAMBDA)
  - map: predicted richness at 1km2 scale
  - D
  



# Variables  
Regional:  
- grwnDD0, grwnDD0_sq, AP  
- npp  
- lcH, Forest, Edge  
- bldgPer, rdLen  
- aspctN, aspctE  
  
Local:  
- SoilTSt  
- VegTot  
- CnpyOpn, CnpyMxd  
- Pasture, Crop, Built  
- aspctN, aspctE  


# Script organization:  
- 000_*.R: Exploratory / draft scripts  
- 00_fn.R: Functions for munging, etc  
- 01_makeDatasets.R: Makes 3 datasets: 'pred', 'vs_30', 'vs_no_pred'  
- 02_runVS.sh: Uses cmdstan to run full Y, WY using 'vs_no_pred'  
- 02a_slurmVS.sh: Cluster job submission - 02_runVS.sh
- 03_runShell.R: Fits an rstanarm mixed model  
- 03a_slurmShell.sh: Cluster job submission - 03_runShell.R
- 04_projpredY.R: Updates shell model with full Y, then runs variable selection  
- 05_projpredWY.R: Updates shell model with full WY, then runs variable selection  
- 06_runBest.sh: Uses cmdstan to run optimal Y, WY using 'pred'  
- 06a_slurmBest.sh: Cluster job submission - 06_runBest.R
- 07_aggregateOutput.R: Reads and munges best model outputs  
- 08_analyseOutput.R: Analyses aggregated outputs  
- 09_figures.R: Makes figures  


# Models: use GP  
- full_Y.stan: all covariates, no predictions  
- full_WY.stan: all covariates, no predictions  
- best_Y.stan: select covariates, predict VD  
- best_WY.stan: select covariates, predict VD  

