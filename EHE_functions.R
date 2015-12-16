###############################################################################
###############################################################################
#Definition of the functions needed
###############################################################################
###############################################################################


###############################################################################
#matrix of site belonging
###############################################################################

#a function to determine a matrix of dimensions 
#(number of population)X(number of loci), allowing to identify the pop_ID 
#of each population X loci

#DatFrame: dataset at the basic format, after importation with 'read.table' 
#function
#data: the same dataset at the 'genind' format
locmat<-function(DatFrame,data)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean)
  #Building the matrix
  matsite<-matrix(data=c(rep("AC",dim(DatFrame[DatFrame$pop_ID=="AC",])[1]/3),
                         rep("AD",dim(DatFrame[DatFrame$pop_ID=="AD",])[1]/3),
                         rep("AN",dim(DatFrame[DatFrame$pop_ID=="AN",])[1]/3),
                         rep("AX",dim(DatFrame[DatFrame$pop_ID=="AX",])[1]/3),
                         rep("CB",dim(DatFrame[DatFrame$pop_ID=="CB",])[1]/3)),
                  nrow=(dim(datapop@tab)[1]), 
                  ncol=(length(data@loc.names)),byrow=F)
}


###############################################################################
#Allelic Richness computation
###############################################################################

#data: a dataset at the 'genind' format
#matsite: output of 'locmat' function for the same dataset
AllRich<-function(data,matsite)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean)
  #First, determining the smaller number of allele across sampled population
  matloc<-t(matrix(data=datapop@loc.fac,nrow=(dim(datapop@tab)[2]), 
                   ncol=(dim(datapop@tab)[1])))
  matpop<-matrix(data=datapop@pop.names, nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]))
  conf<-list(matpop, matloc)
  effN<-(tapply(datapop@tab, conf, sum))
  echMin<-min(effN)
  
  #Second, build of the matrix of total number of sampled allele (in fact 
  #it's always 6 since the samples are normed to 3 diploid individuals)
  truc<-t(as.matrix(table(datapop@loc.fac)))
  x<-matrix(nrow=(dim(effN)[1]), ncol=(dim(effN)[2]), data=truc,byrow=TRUE)
  effTot<-matrix(rep(t(effN),t(x)), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]), byrow=TRUE)
  
  #Third, compute the matrix of Ar for each population/loci combination
  #(see El Mousadik and Petit 1996 for details)
  CoMat<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat[i,j]<-(1-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
                        nCm(effTot[i,j],echMin)))
    }
  }
  
  #Allelic richness in each population, for each LOCUS
  ArLOC<-(tapply(CoMat, conf, sum))
  ArLOC<-ArLOC[order(as.numeric(dimnames(ArLOC)[[1]])),]
  ##determining mean Allelic Richness across site and loci
  #Ar<-(tapply(ArLOC,matsite,mean))
  #determining mean Allelic Richness across loci
  Ar<-(apply(ArLOC,1,mean))
}


###############################################################################
#Heterozygosity computation
###############################################################################

#data: a dataset at the 'genind' format
#matsite: output of 'locmat' function for the same dataset
HeterNei<-function(data,matsite)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean)
  #Heterozygosity (Nei 1987) in each population, for each LOCUS
  HsLOC<-matrix(nrow=(dim(datapop@tab)[1]),
                ncol=(length(data@loc.names)), byrow=TRUE)
  for (i in (1:(dim(datapop@tab)[1]))) {
    dataLOC<-genind2loci(data[data$pop==data$pop[3*(i-1)+1],])
    ss<-summary(dataLOC)
    HsLOC[i,]<-sapply(ss, function(x) H(x$allele))
  }
  ##determining mean Heterozygosity across site and loci
  #Hs<-(tapply(HsLOC,matsite,mean))
  #determining mean Heterozygosity across loci
  Hs<-(apply(HsLOC,1,mean))
}


###############################################################################
#Unbiased Simpson Complement (D*) computation
###############################################################################

#data: a dataset at the 'genind' format
#matsite: output of 'locmat' function for the same dataset
DSimpson<-function(data,matsite)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean)
  #First we count the number of different genotype in each population of 3 
  #individuals
  Gr<-matrix(nrow=(dim(datapop@tab)[1]), 
             ncol=1, byrow=TRUE)
  for (i in (1:(dim(datapop@tab)[1]))) {
    ddd<-data[data$pop==data$pop[3*(i-1)+1],]
    Gr[i,]<-length(levels(ddd@other$genotype))
  }
  
  #because we only have three possibilities  (ie 3 different genotypes, 2 
  #different genotypes or only one genotype), we have only three possible 
  #values for Unbiased Simpson Complement and Simpson Eveness indices
  
  Dt<-Gr
  Dt[Dt==1] <-0
  Dt[Dt==2] <-2/3
  Dt[Dt==3] <-1
  ##determining mean Unbiased Simpson Complement index across site and loci
  #D<-(tapply(Dt,matsite[,1],mean))
  #determining mean Unbiased Simpson Complement index across loci
  D<-(apply(Dt,1,mean))
}


###############################################################################
#Simpson Evenness computation
###############################################################################

#data: a dataset at the 'genind' format
#matsite: output of 'locmat' function for the same dataset
VSimpson<-function(data,matsite)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean)
  #First we count the number of different genotype in each population of 3 
  #individuals
  Gr<-matrix(nrow=(dim(datapop@tab)[1]), 
             ncol=1, byrow=TRUE)
  for (i in (1:(dim(datapop@tab)[1]))) {
    ddd<-data[data$pop==data$pop[3*(i-1)+1],]
    Gr[i,]<-length(levels(ddd@other$genotype))
  }
  
  #because we only have three possibilities  (ie 3 different genotypes, 2 
  #different genotypes or only one genotype), we have only three possible 
  #values for Unbiased Simpson Complement and Simpson Eveness indices
  
  Vt<-Gr
  Vt[Vt==1] <-0
  Vt[Vt==2] <-0
  Vt[Vt==3] <-1
  ##determining mean Evenness index across site and loci
  #V<-(tapply(Vt,matsite[,1],mean))
  #determining mean Evenness index across loci
  V<-(apply(Vt,1,mean))
}


###############################################################################
#END
###############################################################################