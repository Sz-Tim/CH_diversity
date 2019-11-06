#  II: Leveraging citizen science to assess richness, diversity, and abundance  

--------  
## Background  
What structures ant species richness, diversity, and community structure at different spatial scales? We know that at a regional scale, climate is generally important, with some support for phylogenetically conserved temperature preferences. At a local scale, richness typically decreases with canopy cover. In general, abundance seems to be more idiosyncratic and variable, both temporally and spatially. This is particularly true for small species. In ants, abundance can be measured as either the number of colonies in a particular area (i.e., colony density) or as the number of workers (i.e., worker density). 
Increasingly, ecologists have access to occurrence data collected in various and haphazard ways, typically in the form of online databases or citizen science projects. These data are commonly used for species distribution models, but their use in predicting richness or diversity directly has been somewhat more limited. There are several reasons for this. First, the data do not generally come from communities or assemblages, but rather an aggregation of detections from many different collectors across a variety of time spans. We like to think of diversity and richness as properties of communities, and these are decidedly not samples of communities, obscuring the ability to detect or account for interactions among species. Second, the collections for each species may have differing spatial biases, rendering any simple aggregation methods erroneous. Third, there are biases in the species that are more likely to be detected, such that any estimate of richness or diversity will necessarily be of a subset of the community biased toward larger, more active, or more interesting species.

However, that doesn't mean these data can't be useful. Instead, the geographic breadth and rather indiscriminate collection methods can capture occurrences in unexpected locations or species that may be missed in alternative, more structured sampling methods. Here, we combine species occurrences of ants collected in a citizen science project in western Switzerland with a concurrent structured sampling effort. In a sequential Bayesian framework, we use the citizen science data to estimate the regional variables driving each species' distribution, then use the species-specific, genus-specific, and overall responses as prior distributions in a cross-scale model of the structured samples. We then estimate species diversity in each community, accounting for species that may not have been detected, and predicting the posterior alpha, beta, and gamma diversities, capturing the uncertainty in the community composition. Based on the species composition in each sample from the posterior distribution, we calculate taxonomic and phylogenetic diversity across spatial scales. For comparison, we run the model with uninformed prior distributions as well to assess the impact of the citizen science data.
 


--------  
## Questions  

1. What environmental factors drive ant diversity at different spatial scales?  
  - At a regional scale (i.e., 1 km^2 sites), richness will most likely be determined by temperature and productivity  
  - At a local scale (i.e., soil plots), richness willl most likely be determined by habitat type and relative soil temperature
  - Differences between richness and diversity? A lot of the theory is really about richness...?  
2. Is there density compensation?  
3. How much do results differ based on the inclusion or exclusion of each data source? Does W add anything to predictions about local communities?  


--------  
## Approaches  

For the model structure, there are 5 broad approaches. The goal of each is to predict local communities, including relative abundances, by drawing on information from both datasets. We are making the assumptions that 1) the structured samples are more representative of the local community on average, and 2) the citizen science data can be used to inform broader-scale distributions, given the lack of local environmental information. According to simulations in Pacifici et al (2017), using correlated spatial random effects to link the two data sources is a good way to join information when the two datasets are not expected to be of the same quality. It's unclear to me exactly why that would give the better dataset more influence though... A challenge with the correlated spatial random effects is that this is on a geographic and taxonomic scale that is quite likely prohibitive, even ignoring the complications of potential correlations among species.  
1. **Sequential + covariate layer.** The citizen science data are summarized in some manner to produce a raster(s), which can then be used as covariate(s) in the structured sampling model. Conceptually simple, but in practice it gets complicated for an entire community. This would likely require building separate SDMs for each species. Donc, it seems that this approach is perhaps not the most practical, not the most interesting, and not really aligned with questions about diversity.  
2. **Sequential + informed priors + $\Lambda_W$.** Separate models are fit in sequence, where the posterior distributions generated by the citizen science model are used as prior distributions for the structured sampling model. Both models are point process models focused on the intensity $\lambda$. The citizen science data would be aggregated to 1-km^2 cells, with counts for each species and a total estimate of effort. The counts would be estimates of $\Lambda$ scaled by the effort, where $\Lambda$ represents the latent abundance per unit of sampling effort within each cell, and is predicted by a set of covariates $X_r$ and slopes $\beta_r$. The structured sampling data are used as-is, with regional processes at the 1-km^2 scale, similarly with $\Lambda$ predicted by covariates $X_r$ and slopes $\beta_r$, and local processes at the soil plot scale, with $\lambda$ predicted by covariates $V_l$ and slopes $\alpha_l$. The challenges for this approach involve aligning the definitions of $\Lambda$ such that the two quantities are comparable across the models, and appropriately accounting for sampling effort in the citizen science model (both in terms of the effort across the environmental gradient, and in terms of uneven effort among species).  
3. **Sequential + informed priors + $\Psi_W$.** Separate models are fit in sequence, where the posterior distributions generated by the citizen science model are used as prior distributions for the structured sampling model. The structured sampling model is a point process model with local intensity $\lambda$, but the citizen science model is a presence/absence model with probability of presence $\Psi$. Most likely, the structured sampling model would involve a site-level presence-absence term $\Z$ instead of a site-level intensity $\Lambda$ to align the regressions and allow for the informed prior distributions. The slopes $\beta_r$ then relate to the probability of presence within a 1-km^2 area rather than the intensity. Thus, this approach loses the ability to predict relative abundance. The challenges for this approach involve aligning these probabilities of presence at the 1-km^2 level, and accounting for sampling effort in the citizen science data. This could take a bit of thought, as the examples in the literature (e.g., Miller et al 2019) involve species-specific individual detection probabilities, which we do not have.  
4. **Joint + shared $\beta$ + $\Lambda_W$.** Similar structure as Approach 2, but a single joint model estimates the slopes $\beta_r$ simultaneously. This might be better able to account for undetected species in the local communities, since the slopes will be jointly estimated at each taxonomic level. The taxonomic/phylogenetic hierarchical structure of the slopes makes that difficult if the models are fit sequentially, since there is no direct prior on each $\b_s$ or $\B_g$. As above, the challenge is in aligning the definitions of $\Lambda$, accounting for uneven effort in the citizen science model, but also perhaps needing to weight the two datasets differently. The number of points may be prohibitive, but the correlated spatial random effects in the examples in Miller et al 2019 could provide a solution.  
5. **Joint + shared $\beta$ + $\Psi_W$.** Similar structure as in Approach 3, but a single joint model estimates the slopes $\beta_r$ simultaneously. The benifits are identical to the benefits described for Approach 4, and the challenges are identical to the challenges for Approach 3. As with Approach 4, it may be necessary to weight the datasets, since the sample sizes would naturally give more weight to the citizen science dataset which we are assuming is somewhat lower in quality.  


--------  
## Outcomes  

- Bayesian framework makes the most sense for combining disparate data sources  
- Focus on local communities rather than elevational bins  
- The local component can also include covariance among species to look at associations  
- Focus on aggregate measures (richness, a/b/g diversity) instead of individual species  
- Point Process Model  
- Could compare patterns across different datasets:  
	1. citizen science only  
	2. structured sampling only  
	3. simple regional totals  
	4. hierarchical model  
- Citizen science data will capture some species that were not detected in the structured sampling, so predicted communities will reflect the likely presence of these species.  


--------  
## Covariates  

1. log($\Lambda$), logit($\Psi$) ~ MAP + Tmin/mean/max + DTR + TAR + NPP + H'Hab  
2. log($\lambda$) ~ Veg + SoilTempAnomaly + Habitat/CanopyType  


--------  
## Challenges  

The two datasets were collected in quite different manners. The citizen science data, $W$, are georeferenced, haphazard detections with clear spatial bias clustering around cities and hiking trails. The observations are likely biased toward large, active, conspicuous species. The structured sampling data, $Y$, are stratified random samples of assemblages of soil-dwelling ants, with sites arranged on a grid. The observations are biased toward species nesting in the soil. I have ignored the transect samples so far, but they are not as spatially independent as the soil samples, and intentionally include only mound-building species.


- Different collection methods, which results in:  
	1. Different estimated quantities (density of soil colonies/mounds, frequency of easily detected species[?])  
	2. Different biases in species compositions  
		- using presence/absence rather than intensity would reduce this effect, but lose information  
		- bias should affect the *intercepts* b0 rather than the slopes b1–R  
		- a species-specific offset could be used: $W_ks ~ Poisson(\Lambda_{ks} E_k D_s)$), where $E_k$ is the proportional effort and $D_s$ is the species-specific proportional bias, with $D_s \lt 1$ meaning the species is less likely to be detected, and $D_s \gt 1$ meaning the species is more likely to be detected.  
	3. Different spatial distribution of collections  
	  - I think this doesn't matter if the datasets are only linked by the covariate estimation, so long as the covariates used are the same and are scaled together.  
	4. Different amounts of local environmental data  
	  - 
	5. Different locational error  
	6. Unique species in each dataset (probably)  
- Relative abundance is better represented in structured sampling data  
- Do species with different nesting/activity habits need to be treated differently?  
- Does the model need to account for space? Could be limited by number of samples...
  - Spatial random effects for soil plots within squares?


--------  
## Modelling questions  
- In the joint likelihood implementation, how much do prediction differ based on the inclusion or exclusion of each data source? Does W add anything?  





