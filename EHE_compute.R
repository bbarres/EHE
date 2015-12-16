###############################################################################
###############################################################################
#computation of gene and genotypic diversity indices
###############################################################################
###############################################################################

#before using this code, you have to run 'EHE_functions.R' first
#then you have to set the right working directory (where the data file is)
setwd("~/work/Rfichiers/Githuber/EHE")

#loading the needed packages
require(adegenet)
require(combinat)
require(pegas)

#loading the dataset
EHE<-read.table(file="EHEdata.txt", header=T)
#see the structure of the file
head(EHE)

#because the minimum number of individual at the lower level (leaf) is low 
#(N=3), we have restricted the sample at the leaf and the twig levels only 
#to population with 3 individuals (at the leaf level), and 3 leaves (at the 
#twig level). This restriction is not neccessary at the tree and the site 
#levels because there are always three differents twigs per tree or three 
#different tree per site that comprise at least one sample. 
#select individuals on leaves with 3 individuals
EHElv<-EHE[EHE$eff_lv==3,]
#select individuals on twigs with 3 leaves
EHEtw<-EHE[EHE$nb_lv==3,]

###############################################################################
#Resampling individuals at the infralevel
###############################################################################

#creation of a resampling matrix, each column correspond to the list of ind_ID 
#resampled at the level below (here, leaf level)
grTW<-matrix(nrow=0,ncol=1000, byrow=F)
for (i in (1:(length(levels(as.factor(EHEtw$leave_ID)))/3))) {
  gr<-matrix(nrow=9,ncol=1000, byrow=F)
  for (j in (1:3)){
    ifelse(length(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+1],]$ind_ID)>1,
           gr[3*(j-1)+1,]<-sample(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+1],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+1,]<-rep(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+1],]$ind_ID,
                               1000))
    ifelse(length(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+2],]$ind_ID)>1,
           gr[3*(j-1)+2,]<-sample(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+2],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+2,]<-rep(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+2],]$ind_ID,
                               1000))
    ifelse(length(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+3],]$ind_ID)>1,
           gr[3*(j-1)+3,]<-sample(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+3],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+3,]<-rep(EHEtw[EHEtw$leave_ID==levels(as.factor(EHEtw$leave_ID))[3*(i-1)+3],]$ind_ID,
                               1000))
  }
  grTW<-rbind(grTW,gr)
}

#same thing for the tree level (level below = twig)
grTR<-matrix(nrow=0,
             ncol=1000, byrow=F)
for (i in (1:(length(levels(as.factor(EHE$twig_ID)))/3))) {
  gr<-matrix(nrow=27,ncol=1000, byrow=F)
  for (j in (1:9)){
    ifelse(length(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+1],]$ind_ID)>1,
           gr[3*(j-1)+1,]<-sample(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+1],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+1,]<-rep(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+1],]$ind_ID,
                               1000))
    ifelse(length(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+2],]$ind_ID)>1,
           gr[3*(j-1)+2,]<-sample(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+2],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+2,]<-rep(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+2],]$ind_ID,
                               1000))
    ifelse(length(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+3],]$ind_ID)>1,
           gr[3*(j-1)+3,]<-sample(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+3],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+3,]<-rep(EHE[EHE$twig_ID==levels(as.factor(EHE$twig_ID))[3*(i-1)+3],]$ind_ID,
                               1000))
  }
  grTR<-rbind(grTR,gr)
}

#same thing for the site level (level below = tree)
grSI<-matrix(nrow=0,
             ncol=1000, byrow=F)
for (i in (1:(length(levels(as.factor(EHE$tree_ID)))/5))) {
  gr<-matrix(nrow=135,ncol=1000, byrow=F)
  for (j in (1:45)){
    ssArb<-(sample(c(levels(as.factor(EHE$tree_ID))[5*(i-1)+1],
                     levels(as.factor(EHE$tree_ID))[5*(i-1)+2],
                     levels(as.factor(EHE$tree_ID))[5*(i-1)+3],
                     levels(as.factor(EHE$tree_ID))[5*(i-1)+4],
                     levels(as.factor(EHE$tree_ID))[5*(i-1)+5]),
                   3,replace=F))
    ifelse(length(EHE[EHE$tree_ID==ssArb[1],]$ind_ID)>1,
           gr[3*(j-1)+1,]<-sample(EHE[EHE$tree_ID==ssArb[1],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+1,]<-rep(EHE[EHE$tree_ID==ssArb[1],]$ind_ID,
                               1000))
    ifelse(length(EHE[EHE$tree_ID==ssArb[2],]$ind_ID)>1,
           gr[3*(j-1)+2,]<-sample(EHE[EHE$tree_ID==ssArb[2],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+2,]<-rep(EHE[EHE$tree_ID==ssArb[2],]$ind_ID,
                               1000))
    ifelse(length(EHE[EHE$tree_ID==ssArb[3],]$ind_ID)>1,
           gr[3*(j-1)+3,]<-sample(EHE[EHE$tree_ID==ssArb[3],]$ind_ID,
                                  1000,replace=T),
           gr[3*(j-1)+3,]<-rep(EHE[EHE$tree_ID==ssArb[3],]$ind_ID,
                               1000))
  }
  grSI<-rbind(grSI,gr)
}

#We define a vector of population's appartenance in resampled populations for 
#the twig level
poptwr<-c()
for (i in 1:(dim(grTW)[1]/3)){
  poptwr<-c(poptwr,rep(i,3))
}
#the tree level
poptrr<-c()
for (i in 1:(dim(grTR)[1]/3)){
  poptrr<-c(poptrr,rep(i,3))
}
#and the site level
popsir<-c()
for (i in 1:(dim(grSI)[1]/3)){
  popsir<-c(popsir,rep(i,3))
}


#transform dataset in adegenet file format
datalv<-df2genind(EHElv[,14:22], sep=NULL, ncode=6, ind.names=NULL, 
                  loc.names=NULL, pop=as.factor(EHElv$leave_ID), ploidy=2, 
                  type=c("codom"))
datalv@other$genotype<-as.factor(EHElv$geno_ID)

#here is an example for one resampled dataset at the twig level
EHEtwr<-EHE[as.vector(grTW[,1]),]
datatw<-df2genind(EHEtwr[,14:22], sep=NULL, ncode=6, ind.names=NULL, 
                  loc.names=NULL, pop=as.factor(poptwr), ploidy=2, 
                  type=c("codom"))
datatw@other$genotype<-as.factor(EHEtwr$geno_ID)


###############################################################################
#Computation of the different indices: at the leaf level
###############################################################################

#loading the dataset
EHE<-read.table(file="EHEdata.txt", header=T)
#select the usefull data
EHElv<-EHE[EHE$eff_lv==3,]
#transform dataset in adegenet file format
datalv<-df2genind(EHElv[,14:22], sep=NULL, ncode=6, ind.names=NULL, 
                  loc.names=NULL, pop=as.factor(EHElv$leave_ID), ploidy=2, 
                  type=c("codom"))
datalv@other$genotype<-as.factor(EHElv$geno_ID)
#computation of the "matsite" matrix
matsitelv<-locmat(EHElv,datalv)
#computation of allelic richness at the leaf level
Arf<-AllRich(datalv,matsitelv)
#computation of Heterozygosity at the leaf level
Hsf<-HeterNei(datalv,matsitelv)
#computation of Unbiased Simpson Complement (D*) at the leaf level
Df<-DSimpson(datalv,matsitelv)
#computation of Evenness at the leaf level
Vf<-VSimpson(datalv,matsitelv)

Arf_mean<-tapply(Arf,matrix(matsitelv[,1],nrow=1,
                            ncol=length(matsitelv[,1]),byrow=T),mean)
Hsf_mean<-tapply(Hsf,matrix(matsitelv[,1],nrow=1,
                            ncol=length(matsitelv[,1]),byrow=T),mean)
Df_mean<-tapply(Df,matrix(matsitelv[,1],nrow=1,
                          ncol=length(matsitelv[,1]),byrow=T),mean)
Vf_mean<-tapply(Vf,matrix(matsitelv[,1],nrow=1,
                          ncol=length(matsitelv[,1]),byrow=T),mean)
Arf_var<-tapply(Arf,matrix(matsitelv[,1],nrow=1,
                           ncol=length(matsitelv[,1]),byrow=T),var)
Hsf_var<-tapply(Hsf,matrix(matsitelv[,1],nrow=1,
                           ncol=length(matsitelv[,1]),byrow=T),var)
Df_var<-tapply(Df,matrix(matsitelv[,1],nrow=1,
                         ncol=length(matsitelv[,1]),byrow=T),var)
Vf_var<-tapply(Vf,matrix(matsitelv[,1],nrow=1,
                         ncol=length(matsitelv[,1]),byrow=T),var)


###############################################################################
#Computation of the different indices: at the twig level
###############################################################################

Artw<-c()
Hstw<-c()
Dtw<-c()
Vtw<-c()
for (i in (1:1000)) {
  EHEtwr<-EHE[as.vector(grTW[,i]),]
  datatw<-df2genind(EHEtwr[,14:22], sep=NULL, ncode=6, ind.names=NULL, 
                    loc.names=NULL, pop=as.factor(poptwr), ploidy=2, 
                    type=c("codom"))
  datatw@other$genotype<-as.factor(EHEtwr$geno_ID)
  matsitetw<-locmat(EHEtwr,datatw)
  Art<-AllRich(datatw,matsitetw)
  Artw<-rbind(Artw,Art)
  Hst<-HeterNei(datatw,matsitetw)
  Hstw<-rbind(Hstw,Hst)
  Dt<-DSimpson(datatw,matsitetw)
  Dtw<-rbind(Dtw,Dt)
  Vt<-VSimpson(datatw,matsitetw)
  Vtw<-rbind(Vtw,Vt)
}


#to compute mean and variance by site for every resampled dataset
Artw_mean<-tapply(Artw,matrix(matsitetw[,1],nrow=dim(Artw)[1],
                              ncol=length(matsitetw[,1]),byrow=T),mean)
Hstw_mean<-tapply(Hstw,matrix(matsitetw[,1],nrow=dim(Hstw)[1],
                              ncol=length(matsitetw[,1]),byrow=T),mean)
Dtw_mean<-tapply(Dtw,matrix(matsitetw[,1],nrow=dim(Dtw)[1],
                            ncol=length(matsitetw[,1]),byrow=T),mean)
Vtw_mean<-tapply(Vtw,matrix(matsitetw[,1],nrow=dim(Vtw)[1],
                            ncol=length(matsitetw[,1]),byrow=T),mean)
Artw_var<-tapply(Artw,matrix(matsitetw[,1],nrow=dim(Artw)[1],
                             ncol=length(matsitetw[,1]),byrow=T),var)
Hstw_var<-tapply(Hstw,matrix(matsitetw[,1],nrow=dim(Hstw)[1],
                             ncol=length(matsitetw[,1]),byrow=T),var)
Dtw_var<-tapply(Dtw,matrix(matsitetw[,1],nrow=dim(Dtw)[1],
                           ncol=length(matsitetw[,1]),byrow=T),var)
Vtw_var<-tapply(Vtw,matrix(matsitetw[,1],nrow=dim(Vtw)[1],
                           ncol=length(matsitetw[,1]),byrow=T),var)


###############################################################################
#Computation of the different indices: at the tree level
###############################################################################

Artr<-c()
Hstr<-c()
Dtr<-c()
Vtr<-c()
for (i in (1:1000)) {
  EHEtrr<-EHE[as.vector(grTR[,i]),]
  datatr<-df2genind(EHEtrr[,14:22], sep=NULL, ncode=6, ind.names=NULL, 
                    loc.names=NULL, pop=as.factor(poptrr), ploidy=2, 
                    type=c("codom"))
  datatr@other$genotype<-as.factor(EHEtrr$geno_ID)
  matsitetr<-locmat(EHEtrr,datatr)
  Art<-AllRich(datatr,matsitetr)
  Artr<-rbind(Artr,Art)
  Hst<-HeterNei(datatr,matsitetr)
  Hstr<-rbind(Hstr,Hst)
  Dt<-DSimpson(datatr,matsitetr)
  Dtr<-rbind(Dtr,Dt)
  Vt<-VSimpson(datatr,matsitetr)
  Vtr<-rbind(Vtr,Vt)
}


#to compute mean and variance by site for every resampled dataset
Artr_mean<-tapply(Artr,matrix(matsitetr[,1],nrow=dim(Artr)[1],
                              ncol=length(matsitetr[,1]),byrow=T),mean)
Hstr_mean<-tapply(Hstr,matrix(matsitetr[,1],nrow=dim(Hstr)[1],
                              ncol=length(matsitetr[,1]),byrow=T),mean)
Dtr_mean<-tapply(Dtr,matrix(matsitetr[,1],nrow=dim(Dtr)[1],
                            ncol=length(matsitetr[,1]),byrow=T),mean)
Vtr_mean<-tapply(Vtr,matrix(matsitetr[,1],nrow=dim(Vtr)[1],
                            ncol=length(matsitetr[,1]),byrow=T),mean)
Artr_var<-tapply(Artr,matrix(matsitetr[,1],nrow=dim(Artr)[1],
                             ncol=length(matsitetr[,1]),byrow=T),var)
Hstr_var<-tapply(Hstr,matrix(matsitetr[,1],nrow=dim(Hstr)[1],
                             ncol=length(matsitetr[,1]),byrow=T),var)
Dtr_var<-tapply(Dtr,matrix(matsitetr[,1],nrow=dim(Dtr)[1],
                           ncol=length(matsitetr[,1]),byrow=T),var)
Vtr_var<-tapply(Vtr,matrix(matsitetr[,1],nrow=dim(Vtr)[1],
                           ncol=length(matsitetr[,1]),byrow=T),var)


###############################################################################
#Computation of the different indices: at the site level
###############################################################################

Arsi<-c()
Hssi<-c()
Dsi<-c()
Vsi<-c()
for (i in (1:1000)) {
  EHEsir<-EHE[as.vector(grSI[,i]),]
  datasi<-df2genind(EHEsir[,14:22], sep=NULL, ncode=6, ind.names=NULL, 
                    loc.names=NULL, pop=as.factor(popsir), ploidy=2, 
                    type=c("codom"))
  datasi@other$genotype<-as.factor(EHEsir$geno_ID)
  matsitesi<-locmat(EHEsir,datasi)
  Art<-AllRich(datasi,matsitesi)
  Arsi<-rbind(Arsi,Art)
  Hst<-HeterNei(datasi,matsitesi)
  Hssi<-rbind(Hssi,Hst)
  Dt<-DSimpson(datasi,matsitesi)
  Dsi<-rbind(Dsi,Dt)
  Vt<-VSimpson(datasi,matsitesi)
  Vsi<-rbind(Vsi,Vt)
}

#to compute mean and variance by site for every resampled dataset
Arsi_mean<-tapply(Arsi,matrix(matsitesi[,1],nrow=dim(Arsi)[1],
                              ncol=length(matsitesi[,1]),byrow=T),mean)
Hssi_mean<-tapply(Hssi,matrix(matsitesi[,1],nrow=dim(Hssi)[1],
                              ncol=length(matsitesi[,1]),byrow=T),mean)
Dsi_mean<-tapply(Dsi,matrix(matsitesi[,1],nrow=dim(Dsi)[1],
                            ncol=length(matsitesi[,1]),byrow=T),mean)
Vsi_mean<-tapply(Vsi,matrix(matsitesi[,1],nrow=dim(Vsi)[1],
                            ncol=length(matsitesi[,1]),byrow=T),mean)
Arsi_var<-tapply(Arsi,matrix(matsitesi[,1],nrow=dim(Arsi)[1],
                             ncol=length(matsitesi[,1]),byrow=T),var)
Hssi_var<-tapply(Hssi,matrix(matsitesi[,1],nrow=dim(Hssi)[1],
                             ncol=length(matsitesi[,1]),byrow=T),var)
Dsi_var<-tapply(Dsi,matrix(matsitesi[,1],nrow=dim(Dsi)[1],
                           ncol=length(matsitesi[,1]),byrow=T),var)
Vsi_var<-tapply(Vsi,matrix(matsitesi[,1],nrow=dim(Vsi)[1],
                           ncol=length(matsitesi[,1]),byrow=T),var)


citation('adegenet')
citation('combinat')
citation('pegas')


###############################################################################
#END
###############################################################################