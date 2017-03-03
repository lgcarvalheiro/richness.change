############################################################################################## 

#                                  Multilevel.RAR_EXTR

############################################################################################### 
 ###  this script includes functions that allow to calculate relative changes in species richness through time for specific geographical locations
 ###  based on list of records of species presences
 ###  Responsible: Luísa G. Carvalheiro (Naturalis Biodiversity Center/Univeristy of Brasília), lgcarvalheiro@gmail.com
 ###  Date: 1 September 2013
 ###  People who contributed for the code are acknowledge throughout the script 
 ###  Note: after using this code, to check if bias due to differences in sampling effort were indeed corrected,
 ###  generate corrected files if necessary,
 ###  and/or to extract average trends
 ###  please use Trend.extractor script
###########################################################################################
 
 
library(vegan)
library(boot)
          
#-------------------------------------------------------------------------------
	speciesaccum<-function(vec,n,conservative=TRUE, exponential=TRUE, repl=100){ 
#-------------------------------------------------------------------------------     
# the script for this function was written by Luisa Carvalheiro 
# Description: function which decides when to use extrapolation and interpolation to calculate tau_n and SD_n according to Colwell and Mao 
# returns richness moment estimators, in Colwell et al. 2012 this is called "individual-based" multinomial rarefaction      
# returns standard deviation of species accumulation curve at number of individuals sampled n 
# returns also bootstrapped tau_n and SD_n

# n can either be larger or smaller than the actual number of individuals in the sample.
# when n > sum(vec), extrapolation is used according the equations in Colwell et al. 2012
# otherwise it uses interpolation

# Arguments: 	vec	vector of number of records per species 
#		n number of records at which richness has to be assessed
#		conservative: 	for cases of interpolation, TRUE if the maximum number of species is Inf (default)
#					          FALSE uses Chao's estimator to estimate the total number of species
#		exponential: 	for cases of extrapolation, TRUE if exponential approximation of Colwell et al. 2012 is used
#   repl: number of replicates used to calculate bootstrapped tau_n, and SD_n, number of singletons
#
#   Values: 	list
#		tau_n richness estimate at n
#		SD_n	standard deviation of tau_n
#   singletons.mean  bootstrapped number of singletons
#   tau.mean      bootstrapped tau
#   tauSD.mean    bootrapped tauSD
 
#-------------------------------------------------------------------------------

nint<-n
dataset<-rmultinom( repl,sum(vec),vec) # draws 'repl' (default=100) samples from gridsample, each sample has the number of records of original sample
tau.vector <- 1: repl
tauSD.vector <- 1: repl
singletons.vector <- 1: repl


if(n>sum(vec)){
                     
tau_n=speciesextrapolate(vec,n,exponential=F)$tau_n[1]
SD_n=speciesextrapolate(vec,n,exponential=F)$SD_n[1]

singletons<-sum(vec==1)
   if( singletons==0)       # if there are no singletons extrapolation formulas do not work !# this step 
   {tau_n =length(vec)
   SD_n=0                   # specifies extrapolated tau and SD  if there are no singletons 
   }  else{  tau_n =tau_n
   SD_n=SD_n}

}else{
tau_n=speciesinterpolate(vec,n,conservative)$tau_n[1]
SD_n=speciesinterpolate(vec,n,conservative)$SD_n[1]
}


# resample and calculate boot tau and SD_tau 
#calculation of tau based on 'repl' resamples of the data (default is 100 repl)
R=1
for (i in 1:ncol(dataset))
{subVec=dataset[,i]
 subVec=subVec[subVec>0]
 singletonsX<-sum(subVec==1)

if(n>sum(vec)){
  
tau.vector [R] <- speciesextrapolate(subVec,n,exponential=F)$tau_n
tauSD.vector [R] <- speciesextrapolate(subVec,n,exponential=F)$SD_n
singletons.vector [R]<-singletonsX
if( singletonsX==0){tau.vector [R] =length(subVec[subVec>0])}else{tau.vector [R] =tau.vector [R]}  # if there are no singletons 
                                                                                         
singletons.mean=mean( singletons.vector )
tau.mean=mean(tau.vector)
tauSD.mean=mean(tauSD.vector)
R=R+1

 }else{  

tau.vector [R] <- speciesinterpolate(subVec,n,conservative)$tau_n
tauSD.vector [R] <- speciesinterpolate(subVec,n,conservative)$SD_n
singletons.vector [R]<-singletonsX

singletons.mean=mean( singletons.vector )
tau.mean=mean(tau.vector)
tauSD.mean=mean(tauSD.vector)

R=R+1
}



}

  return (list(tau_n=tau_n,SD_n=SD_n,tau.mean=tau.mean, tauSD.mean=tauSD.mean, singletons.mean=singletons.mean, tauVector=tau.vector))
} #end of  speciesaccum


 

#-------------------------------------------------------------------------------
     alpha_sn<-function(s,n,max){
#------------------------------------------------------------------------------- 
# the script for this function was written by Tom van Dooren    
# Description: function which help to calculate the Colwell and Mao 
# richness moment estimators      
# function which calculates combinatorial coefficients.
# Arguments: 	s ,n,max
# Value: a combinatorial coefficient
#-------------------------------------------------------------------------------
  loga<-ifelse(s>(max-n),-Inf,lfactorial(max-n)+lfactorial(max-s)-lfactorial(max-n-s)-lfactorial(max)) 
            # with log factorial, can handle large numbers
          return(exp(loga))
          
          } # end of alpha_sn




#-------------------------------------------------------------------------------
	speciesinterpolate<-function(vec,n,conservative=TRUE){ 
#-------------------------------------------------------------------------------     
# the script for this function was written by Tom van Dooren and, where indicated, eddited by Luisa Carvalheiro    
# Description: function which help to calculate the Colwell and Mao 
# richness moment estimators, in Colwell et al. 2012 this is called "individual-based" multinomial rarefaction      
# returns standard deviation of species accumulation curve at number of individuals sampled n 
# n is smaller or equal to sum(vec)

# Arguments: 	vec	vector of number of records per species 
#		n number of records at which richness has to be assessed
#		conservative: 	TRUE if the maximum number of species is Inf (default)
#				FALSE uses Chao's estimator to estimate the total number fo species	
# Values: 	list
#		tau_n richness estimate at n
#		SD_n	standard deviation of tau_n
#-------------------------------------------------------------------------------

nint<-n

	# table numbers of species per occurring number of records
sindex<-sort(unique(vec))   #! j  , list of number of records in which each species occurs (fi)
s1_to_sn<-table(vec)       #! number of species which occur in j records
s1<-ifelse(sum(sindex==1)>0,s1_to_sn[sindex==1],0)  #! s1 is the nb of singletons
s2<-ifelse(sum(sindex==2)>0,s1_to_sn[sindex==2],0)  #! s2 is the nb of doubletons



Sobs<-sum(s1_to_sn)    #! total number of sps observed
Sest<-Sobs+ifelse(s2>0,s1^2/s2,(s1-1)*s1/(s2+1))/2 # Chao1 estimator
                                                   
alphas<-numeric(length(sindex))   
for(i in 1:length(alphas)) alphas[i]<-alpha_sn(s=sindex[i],n=nint,max=sum(vec))  # vectors of alphas for all js  
                                                     #! alphas will be zero when n is near to the max sampling effort  
tau_n<-Sobs-alphas%*%s1_to_sn       # sum (1 - alpha_j)*s_j
var_n<-(1-alphas)^2%*%s1_to_sn - ifelse(conservative,0,tau_n^2/Sest)

# change made by Luisa Carvalheiro
# when var_n is nearly zero the script from Tom va Dooren can give negative values of var_n. 
# this section garantees that any value that tends to zero (i.e abs(var_n)<10^-14) is transformed to zero
if (abs(var_n[1,])<10^-14) {var_n[1,]=0}  
if(is.na(var_n)){SD_n=0}else{SD_n=sqrt(var_n)}

return (list(tau_n=tau_n,SD_n=SD_n))
} #end of  speciesinterpolate



#-------------------------------------------------------------------------------
	speciesextrapolate<-function(vec,n,exponential=TRUE){ 
#-------------------------------------------------------------------------------     
# the script for this function was written by Tom van Dooren and, where indicated, eddited by Luisa Carvalheiro    
# Description: function which help to calculate the Colwell and Mao 
# richness moment estimators, in Colwell et al. 2012 this is called "individual-based" multinomial rarefaction      
# returns standard deviation of species accumulation curve at number of individuals sampled n 
# n is smaller or equal to sum(vec)

# Arguments: 	vec	vector of number of records per species 
#		n number of records at which richness has to be assessed
#		conservative: 	TRUE if the maximum number of species is Inf (default)
#				FALSE uses Chao's estimator to estimate the total number fo species	
# Values: 	list
#		tau_n richness estimate at n
#		SD_n	standard deviation of tau_n
#-------------------------------------------------------------------------------

nobs<-sum(vec)
m<-n-nobs

# re-table numbers of species per occurring number of records

findex<-sort(unique(vec))   # list of fj occurring, list of counts per species which occur
nf<-length(findex)
f1_to_fn<-table(vec)       # number of species which occur with j records, j goes from the first element of findex, to the last
f1<-ifelse(sum(findex==1)>0,f1_to_fn[findex==1],0)  # f1 is the nb of singletons
f2<-ifelse(sum(findex==2)>0,f1_to_fn[findex==2],0)  # f2 is the nb of doubletons

Sobs<-sum(f1_to_fn)    #! total number of sps observed


f0est<-ifelse(f2>0,f1^2/f2,(f1-1)*f1)/2 # Chao1 estimator

#original code from Tom van Dooren: 
#f0f1ratio<-ifelse(f2>0,f1/f2,(f1-1)/(1+f2))/2 # f0 divided by f1
# when f2 is equal to 0 the code specifies that  f0f1ratio= ((f1-1)/(1+f2))/2
# however this leeds to f0f1ratio=0 if f1=1 , 
# this caused problems further down in the code, when we assum that SD can only be zero if there are no singletons
# change of code done by Luisa Carvalheiro to make it work with no doubletons and only 1 singleton  
f0f1ratio<-if(f2>0){(f1/f2)/2} else if (f2==0&f1==1) {((f1)/(1+f2))/2} else {((f1-1)/(1+f2))/2}      
#other possibility: f0f1ratio<-ifelse(f2>0,f1/f2,(f1)/(1+f2))/2 
# which would simply assume that whenever there is 1 singleton, at least 1 doubleton should be present


expo<-exponential # turn exponential in internal variab
# derivative of S wrt f1

devStof1<-ifelse(f2>0,ifelse(expo,
f1/f2-(m*f2+nobs*f1)*exp(-2*m*f2/nobs/f1)/f2/nobs, # checked
f1/f2 - m*((1-2*f2/f1/nobs)^(m-1))/nobs - f1*((1-2*f2/f1/nobs)^m)/f2), # checked
ifelse(expo,
(2*f1-1)/(f2+1)/2 - (2*m*f1*f2+2*m*f1-3*nobs*f1+nobs+2*nobs*f1^2)*exp(-2*m*(f2+1)/nobs/(f1-1))/(f2+1)/(f1-1)/nobs/2, # checked
-m*f1/(f1-1)*((1-2*(f2+1)/(f1-1)/nobs)^(m-1))/nobs + (1-2*f1)/(f2+1)/2*((1-2*(f2+1)/(f1-1)/nobs)^m-1) ) # checked
)  
if(devStof1=="NaN"){devStof1=0}else{devStof1=devStof1}
# derivative of S wrt f2
 
devStof2<-ifelse(f2>0,ifelse(expo,-exp(-2*m*f2/nobs/f1)*f1*(nobs*f1*(exp(2*m*f2/f1/nobs)-1)-2*m*f2)/2/nobs/f2^2, # checked
m*(f1/f2)*((1-2*f2/f1/nobs)^(m-1))/nobs + (f1^2)/(f2^2)/2*((1-2*f2/f1/nobs)^m-1) ), # checked
ifelse(expo,
-exp(-2*m*(f2+1)/nobs/(f1-1)) * (nobs*(f1-1)*(exp(2*m*(f2+1)/nobs/(f1-1))-1)-2*m*(f2+1))/(2*nobs*(f2+1)^2), # checked
m*f1/(1+f2)*((1-2*(f2+1)/(f1-1)/nobs)^(m-1))/nobs + f1*(f1-1)/((f2+1)^2)/2*((1-2*(f2+1)/(f1-1)/nobs)^m-1) ) # checked
)  
if(devStof2=="NaN"){devStof2=0}else{devStof2=devStof2}
#  derivative of S wrt n

if(f0f1ratio>0)
{
devSton<-ifelse(expo,-f0est*f0f1ratio*m*exp(-m*f0f1ratio/nobs)/nobs^2,-2*f0est*m*(1-1/f0f1ratio/nobs)^(m-1)/(nobs^2)/f0f1ratio) # checked last
} else
{
devSton<-ifelse(expo,-f0est*f0f1ratio*m*exp(-m*f0f1ratio/nobs)/nobs^2,0) # if expo is FALSE and  f0f1ratio=0, devSton will be 0
}

covfifj<- matrix(rep(0,max(2,nf)^2),nrow=max(2,nf)) # covariance matrix fi fj
for(i in 1:nf)for(j in 1:nf)covfifj[i,j]<- -f1_to_fn[i]*f1_to_fn[j]/(Sobs+f0est)
diag(covfifj)<-diag(covfifj)+f1_to_fn # term added to diagonal

devfifj<-rep(1,length(findex))+findex*devSton # dSobs/dfi + DS/dn x dn/dfi

devfifj[1]<-ifelse(f1>0,devStof1+devfifj[1],devfifj[1])
devfifj[1]<-ifelse(f1==0&f2>0,devStof2+devfifj[1],devfifj[1])
devfifj[2]<-ifelse(f1>0&f2>0,devStof2+devfifj[2],ifelse(length(findex)==1,0,devfifj[2]))

# do the ratio in this estimate first, with f1 zero or not f2 zero or not



est<-Sobs+ifelse(f1==0,0, # estimate of S at nobs + m
f0est*(1-ifelse(exponential,exp(-m/f0f1ratio/nobs),(1-1/f0f1ratio/nobs)^m))
)

var<-devfifj%*%covfifj%*%devfifj # summation of variance terms
var<-ifelse(var<0,0,var)

if(is.na(var[1])){SD_n=0}else{SD_n=sqrt(var[1])}
return (list(fvals=f1_to_fn,fnames=findex,Sobs=Sobs,f0=f0est,fratio=f0f1ratio,ex=expo,d1=devStof1,d2=devStof2,dn=devSton,ds=devfifj,cv=covfifj,tau_n=est,SD_n=SD_n))

} #end of  speciesextrapolate


 

# ------------------------------------------------------------------------------
rarefy.1.unit <- function (unit.data, effort){
# ------------------------------------------------------------------------------
# script written by Petr Keil and editted by Luisa Carvalheiro with suggestions from Tom van Dooren
# Description: function similar to Petr Keil´s rarefy.1.unit but it uses the implementation
#              of rarefaction from the vegan package (Hecks fromulas) ; 
#              as well as rarefied richness and SD based on Colwell and Mao work
# Arguments: unit.data - vector containing the species names (from database);
#                        the abundance of each species is represented by the
#                        number of times it occurs in the vector
#            effort - number of records we want to rarefy for
# Value: a vector which includes estimated richness, SD, low95, up95  [according to vegan package code (Heck's formulas) ]
# Col_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean,singletons.mean [from specacum function]
# ------------------------------------------------------------------------------
      # data preparation of the sample.rarefaction
      ones <- rep(1, times=length(unit.data))
      spec.sums <- sort(tapply(ones, unit.data, sum), decreasing=TRUE)  #  table with nb records per sps
       Row.number=1:length(unit.data)
      community.matrix <- t(as.matrix(spec.sums)) # conversion to vegan rarefy format
      community.vector = 1: ncol( community.matrix )    # vector only with sps abund
      for (a in community.vector) {community.vector[a]=community.matrix[a]} 
community.vector<-community.vector[community.vector>0&!is.na(community.vector)] # zero counts removed, probably not necessary  but just in case


       
      singletons<-sum(community.vector==1)
      nrecords<-sum(community.vector)
      effort<-ifelse(is.na(effort),nrecords,effort)
      intensity<-nrecords/length(community.vector)
      
 

    # application of the rarefaction function from VEGAN package   ,
    # SE from Heck  which assumes that SE tends to zero as effort tends to maximum number of records present in the assemblage
  if (effort>sum(community.matrix)) {
    rarefied_i <- NA 
    low95 <- NA
    up95 <- NA
  }else{
       rarefied_i <- (rarefy(community.matrix, effort, se=TRUE))[1:2]
   low95 <- rarefied_i[1]-rarefied_i[2]*1.96
   up95 <- rarefied_i[1]+rarefied_i[2]*1.96
  }
  
 # application of the rarefaction using Colwell and Mao formulas  
SPSAC=speciesaccum(vec=community.vector,n=effort,conservative=F,exponential=F)
Col_MaoSD=SPSAC$SD_n 
Col_Mao_tau=SPSAC$tau_n    
tau.mean=SPSAC$tau.mean
tauSD.mean=SPSAC$tauSD.mean

singletons.mean=speciesaccum(vec=community.vector,n=effort,conservative=F,exponential=F)$singletons.mean
 
tauVector=speciesaccum(vec=community.vector,n=effort,conservative=F,exponential=F)$tauVector

   rarefied <- round(c(rarefied_i[1], low95, up95, Col_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean,singletons.mean),4)
   return(list(rarefied=rarefied,tauVector=tauVector))

} # END of rarery.1.unit





# ------------------------------------------------------------------------------
   RAR.buster <- function(RAR.data,
                          period.1,
                          period.2,
                          output,
                          MIN.rec,
                          spat.scale,
                          DIF,
                          extrapolation){
# script written by Petr Keil and Luisa Carvalheiro
# Description: function that calculates various rarefaction statistics for the
#              defined dataset ot records. The dataset needs to be a data.frame
#              arranged as follows (each row is a record):
#                     column1: species names
#                     column2: trait codes
#                     column3: year of the record
#                     column4: national grid cell code
#                     column5: vector with scale codes for the finer scale
# Arguments: RAR.data - see above
#            period.1 - vector of length 2 (i.e. c(1900,1980))
#            period.2 - vector of length 2 (i.e. c(1981,2005))
#            output - name of the file into which the results will be saved
#                     (the default value is "results.txt")
#            MIN.rec - the criterion of minimal number of records in a
#                      grid cell. If the criterion is not met in the grid cell
#                      than the rarefaction is not performed.
#            spat.scale - the name that identifies the spatial scale at which
#                         the analysis is being done
# Value: a file with all the statistics (see the headers section below)
# ------------------------------------------------------------------------------

  
    RAR.m <- as.matrix(RAR.data)
   cell.codes <- subset(unique(RAR.m[,4]),nchar(unique(RAR.m[,4]))>0)
  trait.codes <- unique(RAR.m[,2])

  # the main body of the function:
  for (trait.group in trait.codes){
    extr0 <- RAR.m[RAR.m[,2]==trait.group,] # extracts the trait.group
  if (length(extr0)>4){   #it only rund if there is more than 1 line (1 line has 4 values) 
  
        for (cell in cell.codes){
         
 
      extr <- extr0[extr0[,4]==cell,] # extracts tge grid.cell
      # extraction of the data from desired time periods:
      # (when the "extr" has only one row then R automaticly converts it to
      # vector and hece the if()s are needed to convert it back to matrix
      if (class(extr)!="matrix") {
        short.mat <- matrix(ncol=4, nrow=1)
        short.mat[] <- extr
        extr <- short.mat
        }
      # extraction of the data in the period.1
      PRE.data <- extr[(extr[,3]>=period.1[1]) & (extr[,3]<=period.1[2]),]
        if (class(PRE.data)!="matrix"){
          short.mat <- matrix(ncol=4, nrow=1)
          short.mat[] <- PRE.data
          PRE.data <- short.mat
        }
      # extraction of the data in the period.2
      POST.data <- extr[(extr[,3]>=period.2[1]) & (extr[,3]<=period.2[2]),]
        if (class(POST.data)!="matrix"){
          short.mat <- matrix(ncol=4, nrow=1)
          short.mat[] <- POST.data
          POST.data <- short.mat
        }

      # calculation of some basic values
      pre.records <- nrow(PRE.data)
      pre.species <- length(unique(PRE.data[,1]))
      pre.ratio <- round(pre.records/pre.species,2)
        if (pre.records==0){pre.ratio <- 0}
      post.records <- nrow(POST.data)
      post.species <- length(unique(POST.data[,1]))
      post.ratio <- round(post.records/post.species, 2)
        if (post.records==0){post.ratio <- 0}
      # difference in the number of records (the 10x fold diff. condition)
        pre.post <- c(pre.records, post.records)
        

   if (pre.records==0|post.records==0){fold.diff = FALSE} else if(post.records/pre.records<DIF & pre.records/post.records<DIF){fold.diff = TRUE} else {fold.diff = FALSE }
      # the condition of the 1.5 records/spec ratio
          if (pre.ratio>=1.5 & post.ratio>=1.5){ratio.1.5 = TRUE} else ratio.1.5 = FALSE


 #### assessing singletons and max rec per sps
 unit.data=PRE.data[,1]
 ones <- rep(1, times=length(unit.data))
      spec.sums <- sort(tapply(ones, unit.data, sum), decreasing=TRUE)  #  table with nb records per sps
       Row.number=1:length(unit.data)
      community.matrix <- t(as.matrix(spec.sums)) # conversion to vegan rarefy format
      community.vector = 1: ncol( community.matrix )    # vector only with sps abund
      for (a in community.vector) {community.vector[a]=community.matrix[a]} 
community.vector<-community.vector[community.vector>0&!is.na(community.vector)] # zero counts removed, probably not necessary  but just in case
PRE_singletons<-sum(community.vector==1)
  if (length(unit.data)>0) { PRE_max_ton=max(community.vector) }    else {PRE_max_ton=0}     #maximum number of records per sps
                                                                                               # if there is no data for one of the periods max_ton=0 
  unit.data=POST.data[,1]
 ones <- rep(1, times=length(unit.data))
      spec.sums <- sort(tapply(ones, unit.data, sum), decreasing=TRUE)  #  table with nb records per sps
       Row.number=1:length(unit.data)
      community.matrix <- t(as.matrix(spec.sums)) # conversion to vegan rarefy format
      community.vector = 1: ncol( community.matrix )    # vector only with sps abund
      for (a in community.vector) {community.vector[a]=community.matrix[a]} 
community.vector<-community.vector[community.vector>0&!is.na(community.vector)] # zero counts removed, probably not necessary  but just in case
  POST_singletons<-sum(community.vector==1)
 if (length(unit.data)>0) { POST_max_ton=max(community.vector) }    else {POST_max_ton=0}  #maximum number of records per sps
                                                                                          # if there is no data for one of the periods max_ton=0 
 
 
 
 
##############

# THE RAREFACTION                       
positive<-function(x)ifelse(x>0,x,10^-6) # here zeroes or negatives are turned into 10^-6
#positive<-function(x) if(is.na(x)) {10^-6} else{ifelse(x>0,x,10^-6)} # here zeroes or negatives are turned into 10^-6

# Ratio function is used to calculate bootstraped change and the confidence intervals for each gridcell
#with the boot function we pick values of pre and post richness each from a normal distribution of mean =tau and SD = tau_sd
# and apply the ratio function.
# from the range of change values we get, we calculate low and up 95%quantile
   Ratio<-function(data){
 mean(log((data$y)/(data$x)),na.rm=TRUE) # I use the log of the ratio as statistic
} 



# negative: decrease, positive: increase

ratio.rg <- function(data,mle) # to draw new random datasets
#  Function to generate random normal variates.  mle will not be used

{    out <- data
     out$y <- positive(rnorm(dim(data)[1],data$y,data$sda)) # note that these rnorm have different sd
     out$x <- positive(rnorm(dim(data)[1],data$x,data$sdb))
     out
}

      
      if (pre.records>=MIN.rec & post.records>=MIN.rec&ratio.1.5==T){

#----------------------extrapolation  ---------------------------------------------#
    
 RRpreL <- rarefy.1.unit(unit.data = PRE.data[,1], effort=max(pre.records,post.records)) 
RRposL <- rarefy.1.unit(unit.data = POST.data[,1], effort=max(pre.records,post.records))
RRpre=RRpreL$rarefied    #rarefied_i[1], low95, up95, Col_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean,singletons.mean, tau Vector
RRpos=RRposL$rarefied
       
       
#######  calculate bootstraped change and confidence intervals
       
 if(!is.na(RRpre[8])&!is.na(RRpos[8]))    # i.e, if neither pre and post tau SD are NA
 { 
 if(RRpre[8]==0&RRpos[8]==0){             # i.e, if both pre and post tau SDare zero
EXT_Change=((RRpos[4]/RRpre[4])-1)*100                # -1  so that 0 corresponds to no.change
EXT_Change_low=EXT_Change       
EXT_Change_up=EXT_Change
EXT_Change_CIspan=0
EXT_Pvalue=  0

}else    {    # if both tau_SD mean is not zero 
 datarel<-data.frame(x=RRpre[4],y=RRpos[4],z=RRpos[4]/RRpre[4],sdb=RRpre[8]*RRpre[4]/RRpre[7],sda=RRpos[8]*RRpos[4]/RRpos[7])


relboot<-boot(datarel,Ratio,R=1000,sim="parametric",ran.gen=ratio.rg,mle=1)
BOOT=boot.ci(relboot,type=c("perc")) # confidence intervals, basic and percentile type

EXT_Change=(exp(BOOT[2]$t0)-1)*100                # 1 is subtracted so that 0 corresponds to no.change
EXT_Change_low=(exp(BOOT[4]$perc[4])-1)*100       
EXT_Change_up=(exp(BOOT[4]$perc[5])-1)*100
EXT_Change_CIspan=abs(EXT_Change_up-EXT_Change_low)
EXT_Pvalue=  t.test(rnorm(100, mean=datarel$x, sd=datarel$sdb),rnorm(100, mean=datarel$y, sd=datarel$sda))$p.value

} 






}else {
EXT_Change=((RRpos[4]/RRpre[4])-1)*100                # -1  so that 0 corresponds to no.change
EXT_Change_low="NA"     
EXT_Change_up="NA"
EXT_Change_CIspan="NA"
EXT_Pvalue=  "-"}
       
#EXT_Change=(RRpos[7]-RRpre[7]) / RRpre[7]

#rarefied_i[1], low95, up95, Col_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean, tau.95CIlow,tau.95CIup,singletons.mean, tau.95CIspan
RR_extr <- c(RRpre[1],RRpre[2],RRpre[3],RRpre[4],RRpre[5],RRpre[6],RRpre[7],RRpre[8],RRpre[9],  	                                                                                                                  # tau.median,tau.95CIlow,tau.95CIup,singletons.mean, tau.95CIspan
RRpos[1],RRpos[2],RRpos[3],RRpos[4],RRpos[5],RRpos[6],RRpos[7],RRpos[8],RRpos[9]) 

 
#----------------------interpolation --------------------------------------#
    
RRpreL <- rarefy.1.unit(unit.data = PRE.data[,1], effort=min(pre.records,post.records))
RRposL <- rarefy.1.unit(unit.data = POST.data[,1], effort=min(pre.records,post.records))
RRpre=RRpreL$rarefied
RRpos=RRposL$rarefied

########################
   #######  calculate bootstraped change
  if(!is.na(RRpre[8])&!is.na(RRpos[8]))
 { 
 if(RRpre[8]==0&RRpos[8]==0){
INT_Change=(RRpos[4]/RRpre[4]-1)*100          # -1  so that 0 corresponds to no.change
INT_Change_low=INT_Change
INT_Change_up=INT_Change
INT_Change_CIspan=0
INT_Pvalue= 0
}else
 {  
 datarel<-data.frame(x=RRpre[4],y=RRpos[4],z=RRpos[4]/RRpre[4],sdb=RRpre[8]*RRpre[4]/RRpre[7],sda=RRpos[8]*RRpos[4]/RRpos[7])


relboot<-boot(datarel,Ratio,R=1000,sim="parametric",ran.gen=ratio.rg,mle=1)
BOOT=boot.ci(relboot,type=c("perc")) # confidence intervals, basic and percentile type
INT_Change=(exp(BOOT[2]$t0)-1)*100          # -1  so that 0 corresponds to no.change
INT_Change_low=(exp(BOOT[4]$perc[4])-1)*100
INT_Change_up=(exp(BOOT[4]$perc[5])-1)*100
INT_Change_CIspan=abs(INT_Change_up-INT_Change_low)
INT_Pvalue=  t.test(rnorm(100, mean=datarel$x, sd=datarel$sdb),rnorm(100, mean=datarel$y, sd=datarel$sda))$p.value
}
  } else {
INT_Change=(RRpos[4]/RRpre[4]-1)*100          # -1  so that 0 corresponds to no.change
INT_Change_low="NA"
INT_Change_up="NA"
INT_Change_CIspan="NA"
INT_Pvalue= "-"
}


    
RR_int <- c(RRpre[1],RRpre[2],RRpre[3],RRpre[4],RRpre[5],RRpre[6],RRpre[7],RRpre[8],RRpre[9],  	                                                                                                                  # tau.median,tau.95CIlow,tau.95CIup,singletons.mean, tau.95CIspan
RRpos[1],RRpos[2],RRpos[3],RRpos[4],RRpos[5],RRpos[6],RRpos[7],RRpos[8],RRpos[9]) 

#------------------------------------------------------------------------------------------------

#----------------------interpolation +extrapolation to intermediate point
# Colwell et al 2012 - the extrapolation method gives
#reasonable results for extrapolations up to about double or triple
#the original abundance or area of the reference sample.
#--------------------------------------------------------------------------    
if (max(pre.records,post.records)>3*min(pre.records,post.records))
{
RRpreL <- rarefy.1.unit(unit.data = PRE.data[,1], effort=3*min(pre.records,post.records))
RRposL <- rarefy.1.unit(unit.data = POST.data[,1], effort=3*min(pre.records,post.records))
RRpre=RRpreL$rarefied
RRpos=RRposL$rarefied


       #######  calculate bootstraped change
       if(!is.na(RRpre[8])&!is.na(RRpos[8]))
 { 
 if(RRpre[8]==0&RRpos[8]==0){
EXT.3x_Change=(RRpos[4]/RRpre[4]-1)*100         # -1  so that 0 corresponds to no.change
EXT.3x_Change_low=EXT.3x_Change
EXT.3x_Change_up=EXT.3x_Change
EXT.3x_Change_CIspan=0
EXT.3x_Pvalue=0
}else { 
 datarel<-data.frame(x=RRpre[4],y=RRpos[4],z=RRpos[4]/RRpre[4],sdb=RRpre[8]*RRpre[4]/RRpre[7],sda=RRpos[8]*RRpos[4]/RRpos[7])

relboot<-boot(datarel,Ratio,R=1000,sim="parametric",ran.gen=ratio.rg,mle=1)

BOOT=boot.ci(relboot,type=c("perc")) # confidence intervals, basic and percentile type
EXT.3x_Change=(exp(BOOT[2]$t0)-1)*100         # -1  so that 0 corresponds to no.change
EXT.3x_Change_low=(exp(BOOT[4]$perc[4])-1)*100
EXT.3x_Change_up=(exp(BOOT[4]$perc[5])-1)*100
EXT.3x_Change_CIspan=abs(EXT.3x_Change_up-EXT.3x_Change_low)
EXT.3x_Pvalue=  t.test(rnorm(100, mean=datarel$x, sd=datarel$sdb),rnorm(100, mean=datarel$y, sd=datarel$sda))$p.value

}
             
} else { 
 RR_3x <- RR_extr
 EXT.3x_Change=EXT_Change
 EXT.3x_Change_low= "NA"
EXT.3x_Change_up="NA"
EXT.3x_Change_CIspan="NA"
EXT.3x_Pvalue="-"}
}   else
#if difference in rec nb is smaller then 3x values no calculations are needed and values will be equal to EXT
 {RR_3x <- RR_extr
 EXT.3x_Change=EXT_Change
 EXT.3x_Change_low= EXT_Change_low
EXT.3x_Change_up=EXT_Change_up
EXT.3x_Change_CIspan=EXT_Change_CIspan
EXT.3x_Pvalue=EXT_Pvalue
}

RR_3x <- c(RRpre[1],RRpre[2],RRpre[3],RRpre[4],RRpre[5],RRpre[6],RRpre[7],RRpre[8],RRpre[9],  	                                                                                                                  # tau.median,tau.95CIlow,tau.95CIup,singletons.mean, tau.95CIspan
RRpos[1],RRpos[2],RRpos[3],RRpos[4],RRpos[5],RRpos[6],RRpos[7],RRpos[8],RRpos[9]) 
###############################################################################
   
 which.rarefied <- "both"     } else {
        RR_extr <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        RR_int <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        RR_3x <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        INT_Change=0
INT_Change_low=0
INT_Change_up=0
INT_Change_CIspan=0
INT_Pvalue="-"

        EXT_Change=0
EXT_Change_low=0
EXT_Change_up=0
EXT_Change_CIspan=0
EXT_Pvalue="-"
        EXT.3x_Change=0
EXT.3x_Change_low=0
EXT.3x_Change_up=0
EXT.3x_Change_CIspan=0
EXT.3x_Pvalue="-"
      
        which.rarefied <- "none"
      }

 
              
###################### Interpolation Heck  


perc.change =    round((RR_int[10]-RR_int[1])/RR_int[1]*100, 3)
      # evaluation of status change based on the comparison of 95 conf ints
      # based on vegan SE
           status.change <- "0"
        # case1
        if (which.rarefied=="none"|is.na(RR_int[1])|is.na(RR_int[11])|is.na(RR_int[12])|is.na(RR_int[2])|is.na(RR_int[3])){
          status.change <- "0"
        } else if (post.records>pre.records & RR_int[1]<RR_int[11]){
          status.change <- "increase"
        } else if (post.records<pre.records & RR_int[3]<RR_int[10]){
           status.change <- "increase"
        }  else if (post.records>pre.records & RR_int[1]>RR_int[12]){
          status.change <- "decrease"
       } else if (post.records<pre.records & RR_int[2]>RR_int[10]){
          status.change <- "decrease"
        } else {status.change <- "no.change"}
        
              
              
      # calculation of the percentual change confidence intervals
      if (post.records<pre.records){
        change.95low <- round((RR_int[10]-RR_int[2])/RR_int[2]*100, 3)
        change.95up <- round((RR_int[10]-RR_int[3])/RR_int[3]*100, 3)
      } else {
        change.95low <- round((RR_int[11]-RR_int[1])/RR_int[1]*100, 3)
        change.95up <- round((RR_int[12]-RR_int[1])/RR_int[1]*100, 3)
      }

      change.CI.span <- abs(change.95up-change.95low)


#######################################################

#######################################################
if(is.na(INT_Change_low)|is.na(INT_Change_up))
  {    status.changeMAO_int<-"-"}   else
{

if  (INT_Change_low>0)
{status.changeMAO_int<-"increase"} else if
(INT_Change_up<0)
{status.changeMAO_int<-"decrease"}  else
{status.changeMAO_int<-"no.change"}
}

      if(is.na(EXT_Change_low)|is.na(EXT_Change_up))
  {    status.changeMAO_extr<-"-"} else
{
if  (EXT_Change_low>0)
{status.changeMAO_extr<-"increase"} else if
(EXT_Change_up<0)
{status.changeMAO_extr<-"decrease"}  else
{status.changeMAO_extr<-"no.change"}
}

      if(is.na(EXT.3x_Change_low)|is.na(EXT.3x_Change_up))
  {    status.changeMAO_3x<-"-"}   else
{
if  (EXT.3x_Change_low>0)
{status.changeMAO_3x<-"increase"} else if
(EXT.3x_Change_up<0)
{status.changeMAO_3x<-"decrease"}  else
{status.changeMAO_3x<-"no.change"}
}



      RESULTS <- data.frame(MIN.rec,
                            spat.scale,
                            cell,
                            trait.group,
                            pre.records,
                            pre.species,
                            pre.ratio,
                            post.records,
                            post.species,
                            post.ratio,
# Interpolation                           
                            RR_int[1], RR_int[2],RR_int[3],  #pre_rarefied  SR  + low and up Heck CI
                            RR_int[10],RR_int[11],RR_int[12],#pos_rarefied  SR  + low and up Heck CI
                            status.change,
                            perc.change,
                            change.95low,
                            change.95up,
                            change.CI.span,
# Interpolation using Mao tau and SD and resampled tau and SD 
                            RR_int[4],RR_int[5],RR_int[6], RR_int[7],RR_int[8],RR_int[9], #PRE Col_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean, singletons.mean
                            RR_int[13],RR_int[14],RR_int[15],RR_int[16],RR_int[17],RR_int[18],  #POSCol_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean, singletons.mean
                            INT_Change,

INT_Change_low,
INT_Change_up,
INT_Change_CIspan,
INT_Pvalue,
                            status.changeMAO_int,
# Extrapolation using Mao tau and SD and resampled tau and SD 
                            RR_extr[4],RR_extr[5],RR_extr[6],RR_extr[7],RR_extr[8],RR_extr[9],  #PRE Col_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean, singletons.mean
                            RR_extr[13],RR_extr[14],RR_extr[15],RR_extr[16],RR_extr[17],RR_extr[18], #POSCol_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean, singletons.mean
                            EXT_Change,
EXT_Change_low,
EXT_Change_up,
EXT_Change_CIspan,
EXT_Pvalue,
                            status.changeMAO_extr,
# Intermediate Extrapolation using Mao tau and SD and resampled tau and SD 
                            RR_3x[4],RR_3x[5],RR_3x[6],RR_3x[7],RR_3x[8],RR_3x[9],  #PRE Col_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean, singletons.mean
                            RR_3x[13],RR_3x[14],RR_3x[15],RR_3x[16],RR_3x[17],RR_3x[18], #POSCol_Mao_tau,Col_MaoSD,singletons, tau.mean,tauSD.mean, singletons.mean
                            EXT.3x_Change,      
EXT.3x_Change_low,
EXT.3x_Change_up,
EXT.3x_Change_CIspan,
EXT.3x_Pvalue,
                            status.changeMAO_3x

 )
                
      if (which.rarefied != "none"&RESULTS$pre.records!=RESULTS$post.records){
        write.table(RESULTS, file=output, append=T, sep="\t", col.names=F)
  }  #end of rarefaction
    
    }
    } # end of the "cell" loop
  } # end of the "trait" loop
} # END of the RAR.buster function









# ------------------------------------------------------------------------------
  multilevel.RAR.buster <- function (multi.RAR.data,
                          taxon.indx,
                          trait.indx,
                          year.indx,
                          scale.data.indcs,
                          period.1m,
                          period.2m,
                          output.file,
                           min.recs,
                           MIN.PROP,
                           FIXGRID,
                           DIF,
                           extrapolation=T){
# ------------------------------------------------------------------------------
# Descritpion: This function applies the RAR.buster function, over all of the
#              spatial scales and over all values of the minimal number of records.
#              The function automatically writes the results into file. This is
#              because it may take some time (typically overnight, or longer) to
#              finish. If the computer stuck, collapse or the process is broken
#              for other reasons, the results remain in the file.
# Arguments:   multi.RAR.data - data.frame with ALL of the data
#              taxon.indx, year.indx, scale.data.indcs - these arguments are
#              simply numbers (indexes) of columns, in which the data occur
#              trait.indx - this is useful when you want to subset your analysis
#                           only for a particular functional group. If you are
#                           not interested in this, just generate a vecor of
#                           a repeating value, for example
#                           rep("allspec", times=nrow(multi.RAR.data), then
#                           add it to the multi.RAR.data and use it as the
#                           trait index
#              period.1m, period.2m - the time periods
#              output.file - that is simple
#              min.recs - vector of the minimum number of records criteria
# Value:       the function generates a file with all of the results in it
#
# ------------------------------------------------------------------------------
                  

 

           
                          headers <- data.frame("min.rec",
               "spat.scale",
               "cell",   # grid-cell id
               "trait.g",     # id of the functional group
               "pre.rec",     # number of records in period.1
               "pre.sps",     # number species in period.1
               "pre.ratio",       # ratio of records/species in period.1
               "post.rec",    # number of records in period.2
               "post.sps",    # number species in period.2
               "post.ratio",      # ratio of records/species in period.2

# Interpolation                           
               "pre.raref",
               "pre.95.low",  #Heck
               "pre.95.up",   #Heck
                 "post.raref",
               "post.95.low", #Heck
               "post.95.up", #Heck
                "status.change",    #! this indicates if the change in a given gridcell is sign or not based on overlap of 95CI based on Heck's SE
               "perc.change",
               "perc.95low",
               "perc.95up",
               "change.CI.span",
               
# Interpolation using Mao tau and SD and resampled tau and SD 
                            
               "pre.Mao.tau_int",
               "pre.Mao.SD_int",
               "pre.singletons",
               "pre.tau.mean_int",
               "pre.tauSD.mean_int",
               "pre.singletons.mean_int",   
                                           
                "post.Mao.tau_int",
               "post.Mao.SD_int",
               "post.singletons",
                "post.tau.mean_int",
               "post.tauSD.mean_int",
               "post.singletons.mean_int", 
    
                "Change.Mao_int",              
"Change.Mao_int_CIlow",
"Change.Mao_int_CIup",
"Change.Mao_int_CIspan",
"Int_Pvalue",
               "status.change.Mao_int",

# Extrapolation using Mao tau and SD and resampled tau and SD 
                    "pre.Mao.tau_extr",
               "pre.Mao.SD_extr",
               "pre.singletons",
               "pre.tau.mean_extr",
               "pre.tauSD.mean_extr",
               "pre.singletons.mean_extr",   
                                           
                "post.Mao.tau_extr",
               "post.Mao.SD_extr",
               "post.singletons",
                "post.tau.mean_extr",
               "post.tauSD.mean_extr",
               "post.singletons.mean_extr", 
    
                "Change.Mao_extr",              
"Change.Mao_extr_CIlow",
"Change.Mao_extr_CIup",
"Change.Mao_extr_CIspan",
 "Extr_Pvalue",
               "status.change.Mao_extr",

# Intermediate Extrapolation using Mao tau and SD and resampled tau and SD  
                    "pre.Mao.tau_extr_3x",
               "pre.Mao.SD_extr_3x",
               "pre.singletons",
               "pre.tau.mean_extr_3x",
               "pre.tauSD.mean_extr_3x",
               "pre.singletons.mean_extr_3x",   
                                           
                "post.Mao.tau_extr_3x",
               "post.Mao.SD_extr_3x",
               "post.singletons",
                "post.tau.mean_extr_3x",
               "post.tauSD.mean_extr_3x",
               "post.singletons.mean_extr_3x", 
    
                "Change.Mao_extr_3x",              
"Change.Mao_extr_3x_CIlow",
"Change.Mao_extr_3x_CIup",
"Change.Mao_extr_3x_CIspan",
"Extr.3x_Pvalue",
               "status.change.Mao_extr_3x"              
               )


 write.table(headers, file=output.file, append=T, sep="\t", col.names=F)

                  
 
 
  # (1) applies RAR.buster for each combination of minimal effort and scale
  
  
  ####################################################################################################
  ###### calculating maximum nb of sps per cell, so that we can calculate nb of records needed to fullfil MinPROP criteria
          
                   
         multi.RAR.data=subset(multi.RAR.data, multi.RAR.data[,trait.indx]!="NA")
         
          trait.data <- multi.RAR.data[,trait.indx]
           taxon.data <- multi.RAR.data[,taxon.indx]
            year.data  <- multi.RAR.data[,year.indx]
           SCALE_SMALL_COLUMN=scale.data.indcs[1]
           cell.data <-  multi.RAR.data[,SCALE_SMALL_COLUMN] 
         
    subset.data <-  data.frame(taxon.data,
                                    trait.data,
                                    year.data,
                                    cell.data
                                    )
            subset.data <-  subset.data[cell.data!=" ",]    
   
               YEAR= subset.data$year.data
   GRIDsmall= subset.data$cell.data
  TAXON= subset.data$taxon.data

  
    pre_vector= seq(period.1m[1], period.1m[2], by=1)      # indicate  dates of pre period
 pos_vector = seq(period.2m[1], period.2m[2], by=1)      #indicate dates of pos period
 pre_data= subset(subset.data, YEAR %in% pre_vector)
 pre_grid=  pre_data$GRIDsmall
 pos_data= subset(subset.data, YEAR %in% pos_vector& GRIDsmall%in%pre_grid )
 pos_grid=  pos_data$GRIDsmall
 pre_data= subset(pre_data, YEAR %in% pre_vector&GRIDsmall%in%pos_grid )


 pre_GRID_SPS=  with(pre_data, aggregate(YEAR/YEAR, by=list(LIST=GRIDsmall,TAXON), FUN='sum'))    #creates list of species per grid cell
 pos_GRID_SPS=  with(pos_data, aggregate(YEAR/YEAR, by=list(LIST=GRIDsmall,TAXON), FUN='sum'))    #creates list of species per grid cell
 
  pre_speciesNB=  with(pre_GRID_SPS, aggregate(x/x, by=list(LIST=LIST), FUN='sum'))                  #counts row number (i.e. species) per grid cell
 pos_speciesNB=  with(pos_GRID_SPS, aggregate(x/x, by=list(LIST=LIST), FUN='sum'))

       MAXpre=max(   pre_speciesNB$x)
MAXpos= max( pos_speciesNB$x)
RecNBprepos=round(max( MIN.PROP*MAXpre,MIN.PROP*MAXpos,min.recs),0  )

pre_GRID_records=  with(pre_data, aggregate(YEAR/YEAR, by=list(LIST=GRIDsmall), FUN='sum'))    #creates list of records per grid cell
 pos_GRID_records=  with(pos_data, aggregate(YEAR/YEAR, by=list(LIST=GRIDsmall), FUN='sum'))    #creates list of records per grid cell
GRIDS= cbind(pre_GRID_records,pos_GRID_records[,2])
colnames(GRIDS)[2]="pre_records"
colnames(GRIDS)[3]="pos_records"
good.grids=subset(GRIDS,pre_records>=RecNBprepos&pos_records>=RecNBprepos)    # selecting grids that match criterioa MinPROP

colnames (multi.RAR.data)[SCALE_SMALL_COLUMN]="SMALL_GRID"
 
 
               
          if(FIXGRID==T){  multi.RAR.data=subset( multi.RAR.data,   SMALL_GRID%in% good.grids$LIST)  } else    { multi.RAR.data= multi.RAR.data  }
 

 
        taxon.data <- multi.RAR.data[,taxon.indx]
           trait.data <- multi.RAR.data[,trait.indx]
           year.data  <- multi.RAR.data[,year.indx]   
             
  # identifying small scale grid cells with good data     
   
  for (spatial.scale in scale.data.indcs){
      for (records.no in min.recs){                          # this line is only relevant if we impose a fixed absolute nb of records instead of a relative nb
      
           cell.data <-  multi.RAR.data[,spatial.scale]
           # extract the data from the sub-scale (the finer scale) which will
           # later be used to describe "the index of dispersion of records":
           
 
 
   RAR.data <-  data.frame(taxon.data,
                                    trait.data,
                                    year.data,
                                    cell.data)
           RAR.data <- RAR.data[cell.data!=" ",] 
        
           
                                                     
               YEAR=RAR.data$year.data
   GRIDsmall=RAR.data$cell.data
  TAXON=RAR.data$taxon.data

    pre_vector= seq(period.1m[1], period.1m[2], by=1)      # indicate  dates of pre period
 pos_vector = seq(period.2m[1], period.2m[2], by=1)      #indicate dates of pos period
 pre_data= subset(RAR.data, year.data %in% pre_vector)
 pre_grid=  pre_data$cell.data
 pos_data= subset(RAR.data, year.data %in% pos_vector& cell.data%in%pre_grid )
 pos_grid=  pos_data$cell.data
 pre_data= subset(pre_data, year.data %in% pre_vector&cell.data%in%pos_grid )

 if(nrow(pre_data)>records.no&nrow(pos_data)>records.no) {

 pre_GRID_SPS=  with(pre_data, aggregate(year.data/year.data, by=list(LIST=cell.data,taxon.data), FUN='sum'))    #creates list of species per grid cell
 pos_GRID_SPS=  with(pos_data, aggregate(year.data/year.data, by=list(LIST=cell.data,taxon.data), FUN='sum'))    #creates list of species per grid cell
 
  pre_speciesNB=  with(pre_GRID_SPS, aggregate(x/x, by=list(LIST=LIST), FUN='sum'))                  #counts row number (i.e. species) per grid cell
 pos_speciesNB=  with(pos_GRID_SPS, aggregate(x/x, by=list(LIST=LIST), FUN='sum'))

       MAXpre=max(   pre_speciesNB$x)
MAXpos= max( pos_speciesNB$x)
RecNBprepos=round(max( MIN.PROP*MAXpre,MIN.PROP*MAXpos,min.recs,records.no),0  )


#RAR.data= rbind(pre_data  , mid_data)
 
           RAR.buster (RAR.data,
                        period.1=period.1m,
                        period.2=period.2m,
                        MIN.rec=RecNBprepos,
                        spat.scale=(names(multi.RAR.data))[spatial.scale],
                        output=output.file,
                        DIF=DIF,
                        extrapolation=extrapolation)
      
           }
      }
  }

} # program end
           
     