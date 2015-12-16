# EHE
The R code used in the article of 2012 on the role of asexual multiplication in poplar rust epidemic

---

[Barrès B, Dutech C, Andrieux A, Halkett F, Frey P (2012) Exploring the role of asexual multiplication in poplar rust epidemics: impact on diversity and genetic structure. Molecular Ecology 21(20): 4996–5008.] (http://dx.doi.org/10.1111/mec.12008)

Originally, the R code as well as the dataset have been published on Dryad:
[Barrès B, Dutech C, Andrieux A, Halkett F, Frey P (2012) Data from: Exploring the role of asexual multiplication in poplar rust epidemics: impact on diversity and genetic structure. Dryad Digital Repository.] (http://dx.doi.org/10.5061/dryad.415sg) The code has been adapted in order to be compatible with the new version of `adegenet` R package. The functions have been separated from the rest of the R code, in order to make things easier to follow. 




![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/EHEsamplingdesign.png "Cartoon of the sampling design of the survey")



##List of the different scripts

* **EHE_functions.R:** the function used in the script EHE_compute.R. You have to run this code first in order to load the functions used in the other script. 
* **EHE_compute.R:** the code for producing the analyses found in the paper. 


##The dataset

* **EHEdata.txt**: the data used for the paper. This file is identical to the data file you can find on the Dryad repository


