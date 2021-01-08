
library(tidyverse); theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

agg.ls <- map(paste0("out/agg_opt_", c("WY", "Y"), ".rds"), readRDS) %>%
  map(., ~discard(.x, is.null))
agg <- vector("list", n_distinct(c(names(agg.ls[[1]]), names(agg.ls[[2]])))) %>%
  setNames(unique(names(agg.ls[[1]]), names(agg.ls[[2]])))
for(i in names(agg)) {
  agg[[i]] <- rbind(agg.ls[[1]][[i]], agg.ls[[2]][[i]]) %>%
    mutate(., model=case_when(model=="Y" ~ "Structured", 
                              model=="WY" ~ "Joint"))
}
rm(agg.ls)

d.ls <- readRDS("data/opt/Y__opt_var_set_ls.rds")
d.i <- readRDS("data/opt/Y__opt_var_set_i.rds")

ants <- load_ant_data(str_type="soil", clean_spp=T)
tax_i <- read_csv("data/tax_i.csv") %>%
  mutate(subf=as.factor(case_when(genus=="Apha" ~ "Myrmicinae",
                                  genus=="Camp" ~ "Formicinae",
                                  genus=="Colo" ~ "Formicinae",
                                  genus=="Crem" ~ "Myrmicinae",
                                  genus=="Doli" ~ "Dolichoderinae",
                                  genus=="Form" ~ "Formicinae",
                                  genus=="Formx" ~ "Myrmicinae",
                                  genus=="Harp" ~ "Myrmicinae",
                                  genus=="Lasi" ~ "Formicinae",
                                  genus=="Lept" ~ "Myrmicinae",
                                  genus=="Mani" ~ "Myrmicinae",
                                  genus=="Myrm" ~ "Myrmicinae",
                                  genus=="Myrme" ~ "Myrmicinae",
                                  genus=="Plag" ~ "Formicinae",
                                  genus=="Poly" ~ "Formicinae",
                                  genus=="Pone" ~ "Ponerinae",
                                  genus=="Sole" ~ "Myrmicinae",
                                  genus=="Sten" ~ "Myrmicinae",
                                  genus=="Tapi" ~ "Dolichoderinae",
                                  genus=="Temn" ~ "Myrmicinae",
                                  genus=="Tetr" ~ "Myrmicinae")),
         fam=as.factor("Formicidae"))

det_Y <- which(colSums(d.ls$Y)>0)

site.mns <- ants$str %>% filter(TypeOfSample=="soil") %>%
  st_set_geometry(NULL) %>% group_by(BDM, Plot_id) %>%
  summarise(nTubes=n(), S=n_distinct(SPECIESID), mnt25=mean(mnt25)) %>%
  group_by(BDM) %>% 
  summarise(mnTubes=sum(nTubes)/25, seTubes=sd(nTubes)/5,
            mnS=sum(S)/25, seS=sd(S)/5,
            el=mean(mnt25))

det.base <- ants$all %>% sf::st_set_geometry(NULL) %>%
  rename(sppName=SPECIESID) %>%
  group_by(sppName, source) %>% summarise(nDet=n()) %>%
  group_by(sppName) %>% mutate(pDet=nDet/sum(nDet))




########------------------------------------------------------------------------
## DATA SUMMARIES
########------------------------------------------------------------------------







########------------------------------------------------------------------------
## SLOPES
########------------------------------------------------------------------------
agg$beta







########------------------------------------------------------------------------
## LAMBDAS
########------------------------------------------------------------------------







########------------------------------------------------------------------------
## Richness
########------------------------------------------------------------------------









########------------------------------------------------------------------------
## DIVERSITY
########------------------------------------------------------------------------



