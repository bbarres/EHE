# Exploring the role of asexual multiplication in poplar rust epidemics: impact on diversity and genetic structure
*The R code used in the article of 2012 on the role of asexual multiplication in poplar rust epidemic*

---
  
The article related to this GitHub repository has been published in Molecular Ecology in 2012: [Barrès B, Dutech C, Andrieux A, Halkett F, Frey P (2012) Exploring the role of asexual multiplication in poplar rust epidemics: impact on diversity and genetic structure. Molecular Ecology 21(20): 4996–5008.](http://dx.doi.org/10.1111/mec.12008)  

Originally, the R code, as well as the dataset, have been published on Dryad:
[Barrès B, Dutech C, Andrieux A, Halkett F, Frey P (2012) Data from: Exploring the role of asexual multiplication in poplar rust epidemics: impact on diversity and genetic structure. Dryad Digital Repository.](http://dx.doi.org/10.5061/dryad.415sg) The code has been adapted in order to be compatible with the new version of `adegenet` R package. The functions have been separated from the rest of the R code, in order to make things easier to follow.  

---

## The question
As numerous pathogens, *Melampsora larici-populina* has the ability to reproduce both sexually and asexually. During the epidemic, both local multiplication (via asexual reproduction) and migration of spores from neighbors infected trees (which derive from independent sexual reproduction events) can influence the disease evolution. But too which extent each type of infection is involved in the evolution of the epidemic remains largely unknown.   

![alt text](https://q77bda.db.files.1drv.com/y4mjzmBC5Xg6S19GrhAC6f4c6r4fQa_N0LPqX3vca6P4sF3u8K-oTatD_pAzEuxYPaTnBUBjge8GE_4Vzm_e0R8LTAloGhA9OE9dgnHUluB64gbteVCDt1AwIIcV_umnSNSAbYa3af1Zrv5vauba8x1DzUVodGLpHJfnd8MJoaOU-PbUHcT-uCHJ-JCEvwWLRJOW7VxVHppLejpljHHY00Exg?width=904&height=354&cropmode=none "Cartoon of the different strains that may influence the course of the epidemic")  

---

## The assay
This study tries to quantify the relative importance of the local *vs.* migrant strains during the course of the epidemic. We sampled strains at several hierarchical levels (from leaves to trees and sites) at the beginning and at the end of the epidemic. Those strains were characterized using molecular markers (microsatellites) and we explored the genetic and genotypic diversity variation in space and time.  

![alt text](https://q77zda.db.files.1drv.com/y4mIITEGA_QJgSfCuVZ8lX0Gd5bitpxZWN1-Z0Vq-RBaEdCaA9WmOvDo2k2ucrkMa273oUDXMBQu_UT0jGh1FA2NXl3PD7jJ0PyYkKAF3CbRWVTTCEpRHoAjE6LsIhSADyEvEphgGV4-4ZPtsjfklHcwX6Ag6c04jf9VcsO3EV-RVtzGifakGBTE9hzsTKaBySqSyZFn6QD0eIvoz_z-9S_Yg?width=925&height=255&cropmode=none "Cartoon of the sampling design of the survey")  

---

## The dataset

* **EHEdata.txt**: the data used for the paper. This file and the one you can find on the Dryad repository are identical. Details for column names: 
 + there is no column header for the Individual identifier
 + *pop_ID*: unique population identifier (combination between location and collection date information). "AC" stand for Charrey-A, "AD" for Amance, "AN" for Prelles-A, "AX" for Charrey-B and "CB" for Prelles-B
 + *eff_loc*: number of individuals in the corresponding population
 + *nb_tr*: number of trees sampled in the corresponding population
 + *tree_ID*: unique tree identifier
 + *eff_tr*: number of individuals sampled in the corresponding tree
 + *nb_tw*: number of twigs sampled in the corresponding tree
 + *twig_ID*: unique twig identifier
 + *eff_tw*: number of individuals sampled in the corresponding twig
 + *nb_lv*: number of leaves sampled in the corresponding twig
 + *leave_ID*: unique leave identifier
 + *eff_lv*: number of individuals sampled in the corresponding leave
 + *ind_ID*: a simplified individual identifier
 + *geno_ID*: Genotype identifier based on the multilocus profile (see publication for details)
 + *[MLP09 to MLP38]*: allele scores at 9 microsatellite loci (see publication and [Barrès et al 2008](http://www.sciencedirect.com/science/article/pii/S1567134808000725) for details). Individuals are diploid and are coded in 6 digits (ie 2 * 3 digits representing the 2 allele sizes in bp).  

---

## List of the different scripts

* **EHE_functions.R:** the function used in the script EHE_compute.R. You have to run this code first in order to load the functions used in the other script. 
* **EHE_compute.R:** the code for producing the analyses found in the paper.  

