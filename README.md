# EHE
The R code used in the article of 2012 on the role of asexual multiplication in poplar rust epidemic

---
  
The article related to this GitHub repository has been published in Molecular Ecology in 2012: [Barrès B, Dutech C, Andrieux A, Halkett F, Frey P (2012) Exploring the role of asexual multiplication in poplar rust epidemics: impact on diversity and genetic structure. Molecular Ecology 21(20): 4996–5008.] (http://dx.doi.org/10.1111/mec.12008)  

Originally, the R code as well as the dataset have been published on Dryad:
[Barrès B, Dutech C, Andrieux A, Halkett F, Frey P (2012) Data from: Exploring the role of asexual multiplication in poplar rust epidemics: impact on diversity and genetic structure. Dryad Digital Repository.] (http://dx.doi.org/10.5061/dryad.415sg) The code has been adapted in order to be compatible with the new version of `adegenet` R package. The functions have been separated from the rest of the R code, in order to make things easier to follow.  




![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/EHEsamplingdesign.png "Cartoon of the sampling design of the survey")



##List of the different scripts

* **EHE_functions.R:** the function used in the script EHE_compute.R. You have to run this code first in order to load the functions used in the other script. 
* **EHE_compute.R:** the code for producing the analyses found in the paper.  


##The dataset

* **EHEdata.txt**: the data used for the paper. This file and the one you can find on the Dryad repository are identical. Details for column names: 
 + there is no column header for the Individual identifier
 + pop_ID: unique population identifier (combination between location and collection date information). "AC" stand for Charrey-A, "AD" for Amance, "AN" for Prelles-A, "AX" for Charrey-B and "CB" for Prelles-B
 + eff_loc: number of individuals in the corresponding population
 + nb_tr: number of trees sampled in the corresponding population
 + tree_ID: unique tree identifier
 + eff_tr: number of individuals sampled in the corresponding tree
 + nb_tw: number of twigs sampled in the corresponding tree
 + twig_ID: unique twig identifier
 + eff_tw: number of individuals sampled in the corresponding twig
 + nb_lv: number of leaves sampled in the corresponding twig
 + leave_ID: unique leave identifier
 + eff_lv: number of individuals sampled in the corresponding leave
 + ind_ID: a simplified individual identifier
 + geno_ID: Genotype identifier based on the multilocus profile (see publication for details)
 + [MLP09 to MLP38]: allele scores at 9 microsatellite loci (see publication and [Barrès et al 2008](http://www.sciencedirect.com/science/article/pii/S1567134808000725) for details). Individuals are diploid and are coded in 6 digits (ie 2 * 3 digits representing the 2 allele sizes in bp). 



