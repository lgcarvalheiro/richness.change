

####################################################################################################################################

                                                     #TREND EXTRACTOR #

#####################################################################################################################################
###  this script includes functions that allow to check if bias due to differences in sampling effort 
###  were indeed corrected and/or to extract average trends for estimates of species richness change calculated based on Multilevel.RAR_EXTR script 
###  please use Trend.extractor script
###  based on list of records of species presences
###  Responsible: Luísa G. Carvalheiro (Naturalis Biodiversity Center/Univeristy of Brasília), lgcarvalheiro@gmail.com
###  Date: 1 September 2013
###########################################################################################


library(boot)

library(lattice)
library(metafor)

# ---------------------------------------------------------------------------------------------------------------------------------
ztest<-function(mean,sd) 2* pnorm(-abs(mean/sd))  # z test= (score-mean)/sd ## 2 * pnorm(- abs(z))   for the two-tailed P-value.
# function to test if richness change is significantly different from zero
#---------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
bootstrap.weighted <- function(boot.vector, wghts, repeats=10000){
# the function calculates the weighted bootstrapped median of the boot.vector
# the function is using the built-in bootstrap function "boot" (package "boot")

# Arguments:
# boot.vector - vector of values which should be bootstrapped
# wghts - weights for each observation
# repeats - number of bootstrap repetitions


# Value:
# mean - simple mean of the boot.vector
# boot.median - bootstrapped median
# boot.median.down, boot.median.up - lower and upper 95% conf. intervals of
# the bootsrapped median
# ------------------------------------------------------------------------------


  vector.mean <- mean(boot.vector)

  my.median <- function(values, i){
    median(values[i])
  }

  boot.results <- boot(data=boot.vector,
                       statistic=my.median,
                       R=repeats,
                       weights=wghts)
  boot.median <- mean(boot.results$t)
  boot.cis <- quantile(boot.results$t, c(0.025,0.975))
  boot.down <- boot.cis[1]
  boot.up <- boot.cis[2]


  results <- c(vector.mean, boot.median, boot.down, boot.up)
  names(results) <- c("mean",
                      "boot.median",
                      "boot.median.down",
                      "boot.median.up")
  return(results)

} # END WEIGHTED BOOTSTRAP

  

# ------------------------------------------------------------------------------

trend.extractor <- function(Data,
                            weight.SE, Min.PROP,Min.REC, Min.SPS, RATIO, DIF, TOL,output.file,
                            corrected.file="corrected_file.txt",TEST=T,CORRECTED.FILE=F,LM=F,
                            Resampling_Bias_correction=T, PhC_threshold=6,rmaTHR=10^-5,Pseff=0.06 ){

# Arguments:
# Data- table of values originated from the multilevel.RAR.buster function 
#min.rec - vector of minimum number of records (tresholds)
# spat.scale - vector of codes indicating spatial scale of an analysis
# perc.change - vector of percentual changes in species richness
# weight.col - vector of measure used to weight richness change
#weight.type - SE or REC (if SE the weight measure will be inverted)
#corrected.file - path where file with estimations per cell (after post-hoc correction, if applicable)  will be saved
# TEST - indicates if we want to run meta-analyses to get mean values per spatial scale
# CORRECTED.FILE - indicates if we want to create a file with estimations per cell (after post-hoc correction, if applicable) 
#Resampling_Bias_correction - indicates if we want to consider the Bias correction due to resampling in the calculations of bootraped SD 
#relative (scaled) bias= abs((boot.tau-tau)/tau), as indicated in Walther and Moore (2005)
#PhC_threshold=indicates a threshold value of number of cells below which the post-hoc correction of sampling effort bias is not applicable
# post-hoc corrections based on very few datapoints might be highly influenced by an outlier 
# rmaTHR indicates the  threshold value of the rma function of package metafor. Usually between 0.01 and 10^-5 results are equal
# Pseff - indicates the threshold P value that determines if the bias caused by differences in sampling effort was completely corrected with interpolation or extrapolation
# if the bias persists a post-hoc correction is applied 

#-------------------------------------------------------------------------------

Data[is.na(Data)]=0
Data2=as.data.frame(subset(Data,Data$cell!=""))
bootstrap=subset(Data2,
Data2$pre.sps>Min.SPS&
Data2$post.sps>Min.SPS&
Data2$pre.ratio>RATIO&
Data2$post.ratio>RATIO)


singletonsTOTAL= as.numeric(bootstrap$pre.singletons)+as.numeric(bootstrap$post.singletons)
subset1=cbind(bootstrap,singletonsTOTAL)


DATA_SELECT =cbind (subset1$pre.rec,subset1$post.rec)
MIN_rec=apply( DATA_SELECT[,c(1,2)],1,min)
DATA_SELECT2 =cbind (subset1$pre.rec,subset1$post.rec, subset1$pre.sps,subset1$post.sps)
least.recNb=1:nrow(subset1)
for (i in 1:nrow(subset1)) {least.recNb[i]=min(subset1$pre.rec[i],subset1$post.rec[i])}

bootstrap2=cbind(subset1,MIN_rec,least.recNb)
Country.set=unique(bootstrap2$Country)

NAMES_trend.results <- data.frame("Country",
                        "Group",
                        "Period",
                        "scale",
                        "trait",
                        "TREND.mean",
                        "estimated.change",
                 #       "median corrected.change",
                  #      "corrected.5quantile",
                   #     "corrected.95quantile",
                        "weight",
                        "Nb_cells",
                        "Min.PROP",
                        "Min.REC",
                        "Min.SPS",
                        "ratio",
                        "DIF",
                        "TOL",
                        "Nb_cells_decline",
                        "Nb_cells_SIGN_decline",
                        "Nb_cells_SIGN_increase",
                        "Weighted_CorlogRatio_95CIlow" ,
                        "Weighted_CorlogRatio_95CIup" ,
                      "Weighted_CorChange_95CIlow" ,
                        "Weighted_CorChange_95CIup" ,
                        "Change Z-test",
                        "Change P-value",
                        "Seffort Z-test",
                        "Seffort P-value",
                        "METHOD",
                        "Estimated_logRatio_SE")




 write.table(NAMES_trend.results, file=output.file, col.names=F, append=T,sep="\t",)
  

if (CORRECTED.FILE==T ){

NAMES_corrected.file <- data.frame("COUNTRY","GROUP","TRAIT","PERIOD","SCALE","METHOD","CellID","pre.Mao.tau","post.Mao.tau",  "Ratio","pre.tau.mean" ,"post.tau.mean" , "Ratio.boot",
                      "pre.sd_rel","post.sd_rel",
                        "pre.tauSD.mean","post.tauSD.mean","least.recNb", "(max.rec-min.rec)/min.rec",
                       "MIN.recNb", "MAX.recNb", "pre.rec","post.rec",
                      "Change.CIspan","abs.diff.tau","diffvar","diffvar.boot",
                 "weight.null","weight1-Var.inv", "weight2-HedgesVar.inv", "weight3-HedgesBootVar.inv", "weight4-BootVar.inv", "pre.sps","post.sps",
                 "CORRECTED_logRatio","Corrected_logratioSD","PValue_cell","CORRECTED_RATIO","PERC","ZTEST_noBootSD","SElogratio","SElogratio_notBoot","CorrectedSD_notboot",
                 "sdb_boot","sdb","sda_boot","sda","Method_for_weightedAverage" ,"weight")
write.table(NAMES_corrected.file,file=corrected.file,col.names=F, append=T,sep="\t")
}

       

for (C in Country.set) {
DATAsubset1=subset(bootstrap2,bootstrap2$Country==C)
Group.set=unique(DATAsubset1$Group)
Ct=C

for (G in Group.set){
DATAsubset2=subset(DATAsubset1,DATAsubset1$Group==G&DATAsubset1$Country==C)
Period.set=unique(DATAsubset2$Period)
Gr=G

for(P in Period.set) {
DATAsubset3=subset(DATAsubset2,DATAsubset2$Period==P)
Pr=P
  Trait.set <- unique(DATAsubset3$trait.g)
  for (M in Trait.set){
    DATAsubset4=subset(DATAsubset3,DATAsubset3$trait.g==M)
Scale.set <- as.matrix(unique(DATAsubset4$spat.scale))

    for (S in Scale.set){
    
DATAsubset5=subset(DATAsubset4,DATAsubset4$spat.scale==S)



 MIN= max(Min.REC,Min.PROP*max(as.numeric(DATAsubset5$pre.sps),as.numeric(DATAsubset5$post.sps)))   #calculates min rec for a given scale
  SELECTION=1:nrow(DATAsubset5)

 # selecting gridcells that match selection criteria:
 MAX_SPS= max(as.numeric(DATAsubset4$pre.sps),as.numeric(DATAsubset4$post.sps))
 for (i in 1:nrow(DATAsubset5)){                 
 if(as.numeric(DATAsubset5$pre.rec[i])>max(100,TOL*MAX_SPS)  # here we indicate that if number of rec is high (i.e >1000 or bigger than max number of sps for the all country), DIF is irrelevant
 &as.numeric(DATAsubset5$post.rec[i])>max(100,TOL*MAX_SPS)) {SELECTION[i]="T"} else if (
as.numeric(DATAsubset5$pre.rec[i])>MIN&as.numeric(DATAsubset5$post.rec[i])>MIN&
as.numeric(DATAsubset5$pre.rec[i])/as.numeric(DATAsubset5$post.rec[i])<DIF&
round(as.numeric(DATAsubset5$post.rec[i])/as.numeric(DATAsubset5$pre.rec[i]),0)<=DIF){SELECTION[i]="T"} else {SELECTION[i]="F"}  }

             
 DATAsubset=subset(DATAsubset5, SELECTION=="T")


if(nrow( DATAsubset)>0){
TR=M

Nb_cells=nrow(DATAsubset)


############################# minimum SD is set to 0.001 to avoid cells with Inf weight 
  DATAsubset$pre.tauSD.mean_int[DATAsubset$pre.tauSD.mean_int==0]=0.001
  DATAsubset$post.tauSD.mean_int[DATAsubset$post.tauSD.mean_int==0]=0.001
  DATAsubset$pre.tauSD.mean_extr[DATAsubset$pre.tauSD.mean_extr==0]=0.001
  DATAsubset$post.tauSD.mean_extr[DATAsubset$post.tauSD.mean_extr==0]=0.001
  DATAsubset$pre.tauSD.mean_extr_3x[DATAsubset$pre.tauSD.mean_extr_3x==0]=0.001
  DATAsubset$post.tauSD.mean_extr_3x[DATAsubset$post.tauSD.mean_extr_3x==0]=0.001
  
  DATAsubset$pre.Mao.SD_int[DATAsubset$pre.Mao.SD_int==0]=0.001
  DATAsubset$post.Mao.SD_int[DATAsubset$post.Mao.SD_int==0]=0.001
  DATAsubset$pre.Mao.SD_extr[DATAsubset$pre.Mao.SD_extr==0]=0.001
  DATAsubset$post.Mao.SD_extr[DATAsubset$post.Mao.SD_extr==0]=0.001
  DATAsubset$pre.Mao.SD_extr_3x[DATAsubset$pre.Mao.SD_extr_3x==0]=0.001
  DATAsubset$post.Mao.SD_extr_3x[DATAsubset$post.Mao.SD_extr_3x==0]=0.001
  

  
# The error (SD) of an estimation at a given point of the accumulation curve (tau) is given by Colwell formulas
# However, such calculations assume random sampling and in historical datasets under/over-representation of singletons or doubletons may occur 
 #(either due to extra effort on sampling rare species or due to extra effort on getting male/female, worker/queen, juvenile/adult, etc)
 # this may lead to very small (under-representation of singletons) or very large (over-representation of singletons) SD estimations
 # bootstrapped SD (tauSD.mean) was calculated to correct for  under/over-representation of singletons or doubletons.
 # SD of tau is relative to tau, and bootSD is relative to boot_tau, which might be higher or lower than the original tau.
 
  # the difference between  the original value (tau)  and mean value obtained with re-sampling ( tau.mean )
   #can be used as an indicator of magnitude of the mean bias caused by re-sampling (Walther and Moore 2005, Ecography 28:815-829),
  # this measured should be scaled by dividing by the original value (Walther and Moore 2005)
  
 #  to take into account such bias bootstrapped SD increases proportionally to such bias (so if an estimation is 10% biased, bootrapped SD will be 10% larger)
 # so that a cell is highly affected by under/over-representation of singletons or doubletons,  will be penalized in subsequent meta-analyses
# the code provides a choice of using the bias correction or not
 
#Note that if there is a very high over-representation of singletons, SD boot will be much smaller than sd_nonboot. 
#so, although a penalization is added for cells that have a large difference between boot and non-boot,
#such penalization might still keep sd-boot smaller than sd non-boot
 
 INT_pre.bias= abs( (as.numeric(DATAsubset$pre.Mao.tau_int)-as.numeric(DATAsubset$pre.tau.mean_int)) /as.numeric(DATAsubset$pre.Mao.tau_int))
INT_post.bias= abs( (as.numeric(DATAsubset$post.Mao.tau_int)-as.numeric(DATAsubset$post.tau.mean_int)) /as.numeric(DATAsubset$post.Mao.tau_int))
EXT_pre.bias= abs( (as.numeric(DATAsubset$pre.Mao.tau_extr)-as.numeric(DATAsubset$pre.tau.mean_extr))  /as.numeric(DATAsubset$pre.Mao.tau_extr))
EXT_post.bias= abs( (as.numeric(DATAsubset$post.Mao.tau_extr)-as.numeric(DATAsubset$post.tau.mean_extr))  /as.numeric(DATAsubset$post.Mao.tau_extr))
EXT3x_pre.bias= abs( (as.numeric(DATAsubset$pre.Mao.tau_extr_3x)-as.numeric(DATAsubset$pre.tau.mean_extr_3x)) /as.numeric(DATAsubset$pre.Mao.tau_extr_3x))
EXT3x_post.bias= abs( (as.numeric(DATAsubset$post.Mao.tau_extr_3x)-as.numeric(DATAsubset$post.tau.mean_extr_3x)) /as.numeric(DATAsubset$post.Mao.tau_extr_3x))



if (Resampling_Bias_correction==T){
INT_pre.sd_rel= as.numeric(DATAsubset$pre.tauSD.mean_int)*(1+ INT_pre.bias)
INT_post.sd_rel=  as.numeric(DATAsubset$post.tauSD.mean_int)*(1+ INT_post.bias)
EXT_pre.sd_rel= as.numeric(DATAsubset$pre.tauSD.mean_extr)*(1+ EXT_pre.bias)
EXT_post.sd_rel= as.numeric(DATAsubset$post.tauSD.mean_extr)*(1+ EXT_post.bias)
EXT3x_pre.sd_rel= as.numeric(DATAsubset$pre.tauSD.mean_extr_3x)*(1+EXT3x_pre.bias)
EXT3x_post.sd_rel= as.numeric(DATAsubset$post.tauSD.mean_extr_3x)*(1+ EXT3x_post.bias)} else   {

INT_pre.sd_rel= as.numeric(DATAsubset$pre.tauSD.mean_int)
INT_post.sd_rel=  as.numeric(DATAsubset$post.tauSD.mean_int)
EXT_pre.sd_rel= as.numeric(DATAsubset$pre.tauSD.mean_extr)
EXT_post.sd_rel= as.numeric(DATAsubset$post.tauSD.mean_extr)
EXT3x_pre.sd_rel= as.numeric(DATAsubset$pre.tauSD.mean_extr_3x)
EXT3x_post.sd_rel= as.numeric(DATAsubset$post.tauSD.mean_extr_3x) }




# create a column of minimum records

MIN.col = 1: nrow(DATAsubset)
for (i in 1: nrow(DATAsubset)) {MIN.col[i] = min(DATAsubset$pre.rec[i],DATAsubset$post.rec[i])}
MAX.col = 1: nrow(DATAsubset)
for (i in 1: nrow(DATAsubset)) {MAX.col[i] = max(DATAsubset$pre.rec[i],DATAsubset$post.rec[i])}

             
  
 



  INTdatarel<-data.frame(CellID=DATAsubset$cell,x=DATAsubset$pre.Mao.tau_int,y=DATAsubset$post.Mao.tau_int,  z=as.numeric(DATAsubset$post.Mao.tau_int)/as.numeric(DATAsubset$pre.Mao.tau_int),
                         x.boot=DATAsubset$pre.tau.mean_int ,y.boot=DATAsubset$post.tau.mean_int , z.boot=as.numeric(DATAsubset$post.tau.mean_int)/as.numeric(DATAsubset$pre.tau.mean_int),
                      sdb= DATAsubset$pre.Mao.SD_int,sda= DATAsubset$post.Mao.SD_int,
                        sdb_boot= INT_pre.sd_rel,sda_boot= INT_post.sd_rel,
                      eff=DATAsubset$least.recNb, dif=abs(as.numeric(DATAsubset$post.rec)-as.numeric(DATAsubset$pre.rec))/as.numeric(MIN.col),
                       MIN.recNb= MIN.col, MAX.recNb =MAX.col, pre.rec=DATAsubset$pre.rec,post.rec=DATAsubset$post.rec,
                      Change.CIspan=DATAsubset$Change.Mao_int_CIspan,
                 abs.dif=(DATAsubset$post.Mao.tau_int)-(DATAsubset$pre.Mao.tau_int),  # absolute difference
                 diffvar= DATAsubset$post.Mao.SD_int^2+DATAsubset$pre.Mao.SD_int^2,  # non-bootstrap variance estimates
                 diffvarb=INT_post.sd_rel^2+INT_pre.sd_rel^2)

  
if  (nrow(  INTdatarel)==1) {
weight0=1
weight1o=1
weight2o=1
weight3o=1
weight4o=1
} else {


weight1o<-1/(  INTdatarel$sdb^2/(INTdatarel$x)^2+  INTdatarel$sda^2/(INTdatarel$y)^2) #inverse of VAR of the log of change
weight2o<-1/(  INTdatarel$sdb^2/(INTdatarel$x^2*INTdatarel$pre.rec)+  INTdatarel$sda^2/(INTdatarel$y^2*INTdatarel$post.rec)) #inverse of VAR of the log of change , Var calculated as described in Hedges et al 1999
weight3o<-1/(  INTdatarel$sdb_boot^2/(INTdatarel$x^2*INTdatarel$pre.rec)+  INTdatarel$sda_boot^2/(INTdatarel$y^2*INTdatarel$post.rec) )      # inverse of bootstrapped var of log of change      , Var calculated as described in Hedges et al 1999
weight4o<- 1/(  INTdatarel$sdb_boot^2/(INTdatarel$x)^2+  INTdatarel$sda_boot^2/(INTdatarel$y)^2 )          # inverse of bootstrapped var)
weight0<-rep(mean(weight4o),length(MIN.col))      # for each cell the weight nule is mean boot_variance of all cells
       
 }   

                
    INTdatarel<-data.frame(  INTdatarel,weight0=weight0,weight1=weight1o,  weight2=weight2o, weight3=weight3o, weight4=weight4o,pre.sps=DATAsubset$pre.sps,post.sps=DATAsubset$post.sps)

   EXTdatarel<-data.frame(CellID=DATAsubset$cell,x=DATAsubset$pre.Mao.tau_extr,y=DATAsubset$post.Mao.tau_extr,  z=as.numeric(DATAsubset$post.Mao.tau_extr)/as.numeric(DATAsubset$pre.Mao.tau_extr),
   x.boot=DATAsubset$pre.tau.mean_extr ,y.boot=DATAsubset$post.tau.mean_extr , z.boot=as.numeric(DATAsubset$post.tau.mean_extr)/as.numeric(DATAsubset$pre.tau.mean_extr),
                   sdb= DATAsubset$pre.Mao.SD_extr,sda= DATAsubset$post.Mao.SD_extr,
                        sdb_boot= EXT_pre.sd_rel,sda_boot= EXT_post.sd_rel,
                      eff=DATAsubset$least.recNb, dif=abs(as.numeric(DATAsubset$post.rec)-as.numeric(DATAsubset$pre.rec))/as.numeric(MIN.col),
                       MIN.recNb= MIN.col, MAX.recNb =MAX.col, pre.rec=DATAsubset$pre.rec,post.rec=DATAsubset$post.rec,
                      Change.CIspan=DATAsubset$Change.Mao_extr_CIspan,
                 abs.dif=(DATAsubset$post.Mao.tau_extr)-(DATAsubset$pre.Mao.tau_extr),  # absolute difference
                 diffvar= DATAsubset$post.Mao.SD_extr^2+DATAsubset$pre.Mao.SD_extr^2,  # non-bootstrap variance estimates
                 diffvarb=EXT_post.sd_rel^2+EXT_pre.sd_rel^2)

  
if  (nrow(  EXTdatarel)==1) {
weight0=1
weight1o=1
weight2o=1
weight3o=1
weight4o=1
} else {

  
weight1o<-1/(  EXTdatarel$sdb^2/(EXTdatarel$x)^2+  EXTdatarel$sda^2/(EXTdatarel$y)^2 ) #inverse of VAR of the log of change
weight2o<-1/(  EXTdatarel$sdb^2/(EXTdatarel$x^2*EXTdatarel$pre.rec)+  EXTdatarel$sda^2/(EXTdatarel$y^2*EXTdatarel$post.rec) ) #inverse of VAR of the log of change , Var calculated as described in Hedges et al 1999
weight3o<-1/(  EXTdatarel$sdb_boot^2/(EXTdatarel$x^2*EXTdatarel$pre.rec)+  EXTdatarel$sda_boot^2/(EXTdatarel$y^2*EXTdatarel$post.rec) )      # inverse of bootstrapped var of log of change      , Var calculated as described in Hedges et al 1999
weight4o<- 1/(  EXTdatarel$sdb_boot^2/(EXTdatarel$x)^2+  EXTdatarel$sda_boot^2/(EXTdatarel$y)^2)          # inverse of bootstrapped var)
weight0<-rep(mean(weight4o),length(MIN.col))      # weight nule is mean boot_variance for every cell
 }   

                 
  EXTdatarel<-data.frame(  EXTdatarel,weight0=weight0,weight1=weight1o,  weight2=weight2o, weight3=weight3o, weight4=weight4o,pre.sps=DATAsubset$pre.sps,post.sps=DATAsubset$post.sps)


   EXT3xdatarel<-data.frame(CellID=DATAsubset$cell,x=DATAsubset$pre.Mao.tau_extr_3x,y=DATAsubset$post.Mao.tau_extr_3x,
                      z=as.numeric(DATAsubset$post.Mao.tau_extr_3x)/as.numeric(DATAsubset$pre.Mao.tau_extr_3x),
                  x.boot=DATAsubset$pre.tau.mean_extr_3x,y.boot=DATAsubset$post.tau.mean_extr_3x , z.boot=as.numeric(DATAsubset$post.tau.mean_extr_3x)/as.numeric(DATAsubset$pre.tau.mean_extr_3x),
                  sdb= DATAsubset$pre.Mao.SD_extr_3x,sda= DATAsubset$post.Mao.SD_extr_3x,
                        sdb_boot= EXT3x_pre.sd_rel,sda_boot= EXT3x_post.sd_rel,
                      eff=DATAsubset$least.recNb, dif=abs(as.numeric(DATAsubset$post.rec)-as.numeric(DATAsubset$pre.rec))/as.numeric(MIN.col),
                       MIN.recNb= MIN.col, MAX.recNb =MAX.col, pre.rec=DATAsubset$pre.rec,post.rec=DATAsubset$post.rec,
                      Change.CIspan=DATAsubset$Change.Mao_extr_3x_CIspan,
                 abs.dif=(DATAsubset$post.Mao.tau_extr_3x)-(DATAsubset$pre.Mao.tau_extr_3x),  # absolute difference
                 diffvar= DATAsubset$post.Mao.SD_extr_3x^2+DATAsubset$pre.Mao.SD_extr_3x^2,  # non-bootstrap variance estimates
                 diffvarb=EXT3x_post.sd_rel^2+EXT3x_pre.sd_rel^2)

  
if  (nrow(  INTdatarel)==1) {
weight0=1
weight1o=1
weight2o=1
weight3o=1
weight4o=1
} else {

weight1o<-1/(  EXT3xdatarel$sdb^2/(EXT3xdatarel$x)^2+  EXT3xdatarel$sda^2/(EXT3xdatarel$y)^2 ) #inverse of VAR of the log of change
weight2o<-1/(  EXT3xdatarel$sdb^2/(EXT3xdatarel$x^2*EXT3xdatarel$pre.rec)+  EXT3xdatarel$sda^2/(EXT3xdatarel$y^2*EXT3xdatarel$post.rec) ) #inverse of VAR of the log of change , Var calculated as described in Hedges et al 1999
weight3o<-1/(  EXT3xdatarel$sdb_boot^2/(EXT3xdatarel$x^2*EXT3xdatarel$pre.rec)+  EXT3xdatarel$sda_boot^2/(EXT3xdatarel$y^2*EXT3xdatarel$post.rec) )      # inverse of bootstrapped var of log of change      , Var calculated as described in Hedges et al 1999
weight4o<- 1/(  EXT3xdatarel$sdb_boot^2/(EXT3xdatarel$x)^2+  EXT3xdatarel$sda_boot^2/(EXT3xdatarel$y)^2 )          # inverse of bootstrapped var)
weight0<-rep(mean(weight4o),length(MIN.col))      # weight nule is mean boot_variance for every cell
 }   

                   
  EXT3xdatarel<-data.frame(  EXT3xdatarel,weight0=weight0,weight1=weight1o,  weight2=weight2o, weight3=weight3o, weight4=weight4o,pre.sps=DATAsubset$pre.sps,post.sps=DATAsubset$post.sps)

datarel2<-na.omit(INTdatarel) # to remove missing values, observations which we can't bootstrap
Edatarel2<-na.omit(EXTdatarel) # to remove missing values, observations which we can't bootstrap
E3xdatarel2<-na.omit(EXT3xdatarel) # to remove missing values, observations which we can't bootstrap




   m="-"
WeightClass=c(0,1,4) # the analyses will be run only with no.weight, Inverse of Variance and Inverse of bootstrapped variance
                     # the other options (not used here) relate to Variance values using vegan package 

for (Weight in  WeightClass) {
if (Weight==1){ Weight.code="Var.Inv"}  else if (Weight==0){ Weight.code="no.weight"} else if (Weight==2){
 Weight.code="HedgesVar.inv"} else if (Weight==3){ Weight.code="HedgesBootVar.inv"}   else if(Weight==4){ Weight.code="BootVar.inv"} 





DATA_TABLE=datarel2

  WEIGHTS.TABLE=cbind(DATA_TABLE$weight1,DATA_TABLE$weight2,DATA_TABLE$weight3,DATA_TABLE$weight4) #,weight4

 PRE.VECTOR = as.matrix(DATA_TABLE$x)
 PRECellIDVECTOR=as.matrix(DATA_TABLE$CellID )
 PRErecIDVECTOR=as.matrix(DATA_TABLE$pre.rec )
PREWEIvector=as.matrix(DATA_TABLE[,22+Weight])

Period=rep("pre",length(PRE.VECTOR ))
PRE.TABLE=cbind( Period,CellID=PRECellIDVECTOR,TAU=as.numeric(PRE.VECTOR),pre.rec= PRErecIDVECTOR,WEIvector=PREWEIvector)
 PRE.TABLE[,3]= as.numeric(PRE.TABLE[,3])





 POST.VECTOR = as.matrix(DATA_TABLE$y)
 POSTCellIDVECTOR=as.matrix(DATA_TABLE$CellID )
 POSTrecIDVECTOR=as.matrix(DATA_TABLE$post.rec )
POSTWEIvector=as.matrix(DATA_TABLE[,22+Weight])

Period=rep("pre",length(POST.VECTOR ))
POST.TABLE=cbind( Period,CellID=POSTCellIDVECTOR,TAU=as.numeric(POST.VECTOR),post.rec= POSTrecIDVECTOR,WEIvector=POSTWEIvector)
 POST.TABLE[,3]= as.numeric(POST.TABLE[,3])


if(nrow(PRE.TABLE)==1){TEST.TABLEa=as.data.frame(PRE.TABLE ) }else{
TEST.TABLEa=as.data.frame(PRE.TABLE,row.names=F ) }

PRE.TAU=  as.numeric(PRE.TABLE[,3])
PRE.TAU[PRE.TAU<0.1]=0.1
POST.TAU=  as.numeric(POST.TABLE[,3])
POST.TAU[POST.TAU<0.1]=0.1

TEST.TABLEchangeDIF=POST.TAU-PRE.TAU
TEST.TABLEchangeRAT=(POST.TAU/PRE.TAU -1)
TEST.TABLEchangeRATlog=log(POST.TAU/PRE.TAU)

TEST.TABLE=cbind(TEST.TABLEa[,c(1:2)],WEIvector=as.numeric(PRE.TABLE[,5]),TEST.TABLEchangeDIF,TEST.TABLEchangeRAT,TEST.TABLEchangeRATlog,
post.rec=as.numeric(POST.TABLE[,4]),pre.rec=as.numeric(PRE.TABLE[,4]))

colnames(TEST.TABLE)[1]="Period"
colnames(TEST.TABLE)[2]="CellID"
colnames(TEST.TABLE)[3]="WEIvector"
colnames(TEST.TABLE)[4]="ChangeDIF"
colnames(TEST.TABLE)[5]="ChangeRAT"
colnames(TEST.TABLE)[6]="ChangeRATlog"
colnames(TEST.TABLE)[7]="post.rec"
colnames(TEST.TABLE)[8]="pre.rec"
TEST.TABLE1=rbind(PRE.TABLE,POST.TABLE)
TAU=as.numeric(TEST.TABLE1[1,3])
TEST.TABLE2=as.data.frame(TEST.TABLE1,row.names=F)




if ( nrow(TEST.TABLE)==1){


Nb_cells_decline="-"
Nb_cells_SIGN_decline="-"
Nb_cells_SIGN_increase="-"

if(nrow(TEST.TABLE)==1&DATA_TABLE$sdb==0&DATA_TABLE$sda==0)
{
EST=log(DATA_TABLE$z)
SEnoBoot=0.00001   # to avoid NaN when calculating pvalues
SE= 0.00001
LOG_RATIO_SE=0.00001

}else          {
EST=log(DATA_TABLE$z)
SEnoBoot= sqrt((DATA_TABLE$sdb^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda^2/DATA_TABLE$y^2))
SE= sqrt((DATA_TABLE$sdb_boot^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda_boot^2/DATA_TABLE$y^2))
LOG_RATIO_SE=SE

                }
EST_low95CI = (exp(EST-1.96*SE)-1)*100
EST_up95CI =  ((EST+1.96*SE)-1)*100
 EST_CHANGE=(exp(EST)-1)*100 
 
 
SEFFZtest="-"
SEFFPvalue="-"
     
Weighted_Corrected_est="-"
Weighted_Corrected_95CIup="-"
Weighted_Corrected_95CIlow="-"
   
Pvalue=ztest(EST,SE)    #               
Ztest=((EST)/SE)
Weighted_Corrected_95CIupBACK_T=EST_up95CI
Weighted_Corrected_95CIlowBACK_T=EST_low95CI

          }else {
#############################
# tests for overall effects #
#############################



# a. meta-analysis - fixed effect approach
# fit the meta-analytic fixed- and random-effects models with or without moderators via the linear (mixed-effects) model

SElogratio= sqrt((DATA_TABLE$sdb_boot^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda_boot^2/DATA_TABLE$y^2))
SElogratio_notBoot=sqrt((DATA_TABLE$sdb^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda^2/DATA_TABLE$y^2))

TEST.TABLE$Ratio=with(TEST.TABLE, (post.rec/pre.rec))



m ="weightedReg"           

 
if(nrow(TEST.TABLE)>PhC_threshold){     # here it is indicated the minimum number of cells needed to apply the regression

 
mmvar<-rma.uni(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,mods=log(TEST.TABLE$Ratio),method="ML",control=list(threshold=rmaTHR)) 
#yi	-vector of length k with the observed effect sizes or outcomes. 
#vi	- vector of length k with the corresponding sampling variances. 
 # tau^2 is estimated automatically, and it provides a measure of total amount of heterogeneity                 



#By default, the starting value
#is set equal to the value of the Hedges estimator and the algorithm terminates when the
#change in the estimated value of  2 is smaller than 10-5 from one iteration to the next. 


CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot

SEFFPvalue=mmvar$pval[2]
SEFFZtest=mmvar$zval[2]
EST_SEFF=   (exp(mmvar$b[2])-1)*100
SE_SEFF=    (exp(mmvar$se[2] )-1)*100

estimate_SEFF=mmvar$b[2]

# post-hoc correction is done when sampling effort has a significant effect, even after rarefaction/extrapolation
# we consider that such bias occurs (and should be corrected) 
# whenever the positive relation between sampling effort change and richness is significant: 
#P value used as threshold to determine that the relation is significant is given by Pseff
# We alert that based on several simulation and analyses with rela data we noticed that when the number of cells is very small,
#the slope could be highly influenced by an outlier, potentially leading to a misleading slope estimate 
if (estimate_SEFF<0|SEFFPvalue>Pseff){   

CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot
 
EST_SEFF=   0
SE_SEFF=   0
estimate_SEFF=0

mmvar<-rma.uni(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,method="ML",control=list(threshold=rmaTHR)) 

} }else {  
mmvar<-rma.uni(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,method="ML",control=list(threshold=rmaTHR)) 
#yi	-vector of length k with the observed effect sizes or outcomes. See Details.
#vi	- vector of length k with the corresponding sampling variances. See Details.

CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot
SEFFPvalue="-"
SEFFZtest="-"

EST_SEFF=   0
SE_SEFF=   0
estimate_SEFF=0}



CORRECTED_VALUES=TEST.TABLEchangeRATlog-estimate_SEFF*log(TEST.TABLE$Ratio)



########################
#CORRECTED VALUES

CORRECTED_RATIO=exp(CORRECTED_VALUES)
 if (S=="G10"&Weight==1){hist(CORRECTED_VALUES, main=paste(Ct,"_GR",Gr,"_P",P))}






if (TEST==T) {

Ztest=mmvar$zval[1]       #intercept
Pvalue=mmvar$pval[1]
EST_CHANGE= (exp(mmvar$b[1])-1)*100 
LOG_RATIO_SE= mmvar$se[1]
}  else {
Ztest="-"       #intercept
Pvalue="-"
EST_CHANGE= (exp(mmvar$b[1])-1)*100        

LOG_RATIO_SE= mmvar$se[1]
            }

Corrected_SDup=CORRECTED_VALUES+CorrectedSD
Corrected_SDlow=CORRECTED_VALUES-CorrectedSD



PERC= (exp(CORRECTED_VALUES)-1)*100

COUNTRY=rep(C,nrow(DATA_TABLE))
GROUP=rep(G,nrow(DATA_TABLE))
PERIOD=rep(P,nrow(DATA_TABLE))
SCALE=rep(S,nrow(DATA_TABLE))
TRAIT=rep(M,nrow(DATA_TABLE))

METHOD=rep("Interp",nrow(DATA_TABLE))

ZTEST= rep(0,nrow(DATA_TABLE))
for (i in 1:nrow(DATA_TABLE)) { ZTEST[i]=ztest(CORRECTED_VALUES[i], (CorrectedSD[i]))}

ZTEST_noBootSD= rep(0,nrow(DATA_TABLE))
for (i in 1:nrow(DATA_TABLE)) { ZTEST_noBootSD[i]=ztest(CORRECTED_VALUES[i], (CorrectedSD_notboot[i]))}

NEW_TABLE=cbind(COUNTRY,GROUP,TRAIT,PERIOD,SCALE,METHOD,DATA_TABLE,CORRECTED_VALUES,CorrectedSD,ZTEST,CORRECTED_RATIO,
PERC,ZTEST_noBootSD, SElogratio,SElogratio_notBoot,CorrectedSD_notboot,
DATA_TABLE$sdb_boot,DATA_TABLE$sdb,DATA_TABLE$sda_boot,DATA_TABLE$sda,m,Weight.code)

 
if (CORRECTED.FILE==T ){
write.table(NEW_TABLE,file=corrected.file,col.names=F, append=T,sep="\t")

                        }
                        
                         

#calculating overall 95CI based on sd or VAR
WeightZERO_ONE=(NEW_TABLE[,28+Weight])/sum(NEW_TABLE[,28+Weight]) # it ensures that weight values vary from 0 to 1
dados<-data.frame(m=CORRECTED_VALUES,w=WeightZERO_ONE)
bootfun<-function(d,i){d2<-d[i,];sum(d2$m*d2$w)/sum(d2$w)}
if( nrow(dados)>2000){bootobj<-boot(dados,bootfun,R=nrow(dados))}else {bootobj<-boot(dados,bootfun,R=2000)}

# these lines account for the fact that if std. error of bootobj is Zero 
#all generated values are equal and boot.ci function doesn't work  
if(length(unique(bootobj$t))==1){            
Weighted_Corrected_95CIup=Weighted_Corrected_est
Weighted_Corrected_95CIlow=Weighted_Corrected_est
Weighted_Corrected_95CIupBACK_T=Weighted_Corrected_est
Weighted_Corrected_95CIlowBACK_T=Weighted_Corrected_est
}else{


Weighted_Corrected_est=bootobj$t0
Weighted_Corrected_95CIup=boot.ci(bootobj,type="bca")$bca[5]
Weighted_Corrected_95CIlow=boot.ci(bootobj,type="bca")$bca[4]

Weighted_Corrected_95CIupBACK_T=(exp(boot.ci(bootobj,type="bca")$bca[5])-1)*100
Weighted_Corrected_95CIlowBACK_T=(exp(boot.ci(bootobj,type="bca")$bca[4])-1)*100
 }
 

Nb_cells_decline=nrow(NEW_TABLE[NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_decline=nrow(NEW_TABLE[NEW_TABLE$ZTEST<0.05&NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_increase=nrow(NEW_TABLE[NEW_TABLE$ZTEST<0.05&NEW_TABLE$CORRECTED_VALUES>0,])

Nb_cells_SIGN_decline_noboot=nrow(NEW_TABLE[NEW_TABLE$ZTEST_noBootSD<0.05&NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_increase_noboot=nrow(NEW_TABLE[NEW_TABLE$ZTEST_noBootSD<0.05&NEW_TABLE$CORRECTED_VALUES>0,])


                 }
########################
 INT.result <- data.frame(Ct,Gr,Pr,S,M,"-",
 EST_CHANGE,
 paste("Int_BOOT_mean_W", Weight.code),Nb_cells,Min.PROP, Min.REC, Min.SPS, RATIO, DIF,TOL,
Nb_cells_decline,
Nb_cells_SIGN_decline,
Nb_cells_SIGN_increase,
Weighted_Corrected_95CIlow,
Weighted_Corrected_95CIup,
Weighted_Corrected_95CIlowBACK_T,
Weighted_Corrected_95CIupBACK_T,
Ztest,
Pvalue,
SEFFZtest,
SEFFPvalue,
m,LOG_RATIO_SE
)
write.table(INT.result, file=output.file, col.names=F, append=T,sep="\t")


       
   


DATA_TABLE=Edatarel2

   WEIGHTS.TABLE=cbind(DATA_TABLE$weight1,DATA_TABLE$weight2,DATA_TABLE$weight3,DATA_TABLE$weight4) 


 PRE.VECTOR = as.matrix(DATA_TABLE$x)
 PRECellIDVECTOR=as.matrix(DATA_TABLE$CellID )
 PRErecIDVECTOR=as.matrix(DATA_TABLE$pre.rec )
PREWEIvector=as.matrix(DATA_TABLE[,22+Weight])

Period=rep("pre",length(PRE.VECTOR ))
PRE.TABLE=cbind( Period,CellID=PRECellIDVECTOR,TAU=as.numeric(PRE.VECTOR),pre.rec= PRErecIDVECTOR,WEIvector=PREWEIvector)
 PRE.TABLE[,3]= as.numeric(PRE.TABLE[,3])



 POST.VECTOR = as.matrix(DATA_TABLE$y)
 POSTCellIDVECTOR=as.matrix(DATA_TABLE$CellID )
 POSTrecIDVECTOR=as.matrix(DATA_TABLE$post.rec )
POSTWEIvector=as.matrix(DATA_TABLE[,22+Weight])

Period=rep("pre",length(POST.VECTOR ))
POST.TABLE=cbind( Period,CellID=POSTCellIDVECTOR,TAU=as.numeric(POST.VECTOR),post.rec= POSTrecIDVECTOR,WEIvector=POSTWEIvector)
 POST.TABLE[,3]= as.numeric(POST.TABLE[,3])



if(nrow(PRE.TABLE)==1){TEST.TABLEa=as.data.frame(PRE.TABLE ) }else{
TEST.TABLEa=as.data.frame(PRE.TABLE,row.names=F ) }

PRE.TAU=  as.numeric(PRE.TABLE[,3])
PRE.TAU[PRE.TAU<0.1]=0.1
POST.TAU=  as.numeric(POST.TABLE[,3])
POST.TAU[POST.TAU<0.1]=0.1

TEST.TABLEchangeDIF=POST.TAU-PRE.TAU
TEST.TABLEchangeRAT=(POST.TAU/PRE.TAU -1)
TEST.TABLEchangeRATlog=log(POST.TAU/PRE.TAU)

TEST.TABLE=cbind(TEST.TABLEa[,c(1:2)],WEIvector=as.numeric(PRE.TABLE[,5]),TEST.TABLEchangeDIF,TEST.TABLEchangeRAT,TEST.TABLEchangeRATlog,post.rec=as.numeric(POST.TABLE[,4]),pre.rec=as.numeric(PRE.TABLE[,4]))


TEST.TABLE1=rbind(PRE.TABLE,POST.TABLE)
TAU=as.numeric(TEST.TABLE1[1,3])
TEST.TABLE2=as.data.frame(TEST.TABLE1,row.names=F)



if ( nrow(TEST.TABLE)==1){


Nb_cells_decline="-"
Nb_cells_SIGN_decline="-"
Nb_cells_SIGN_increase="-"

if(nrow(TEST.TABLE)==1&DATA_TABLE$sdb==0&DATA_TABLE$sda==0)
{
EST=log(DATA_TABLE$z)
SEnoBoot=0.00001   # to avoid NaN when calculating pvalues
SE= 0.00001
LOG_RATIO_SE=0.00001

}else          {
EST=log(DATA_TABLE$z)
SEnoBoot= sqrt((DATA_TABLE$sdb^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda^2/DATA_TABLE$y^2))
SE= sqrt((DATA_TABLE$sdb_boot^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda_boot^2/DATA_TABLE$y^2))
LOG_RATIO_SE=SE

                }
EST_low95CI = (exp(EST-1.96*SE)-1)*100
EST_up95CI =  ((EST+1.96*SE)-1)*100
 EST_CHANGE=(exp(EST)-1)*100 
 
 
SEFFZtest="-"
SEFFPvalue="-"
     
Weighted_Corrected_est="-"
Weighted_Corrected_95CIup="-"
Weighted_Corrected_95CIlow="-"
   
Pvalue=ztest(EST,SE)    #               
Ztest=((EST)/SE)
Weighted_Corrected_95CIupBACK_T=EST_up95CI
Weighted_Corrected_95CIlowBACK_T=EST_low95CI

          }else {
#############################
# tests for overall effects #
#############################



# a. meta-analysis - fixed effect approach
# fit the meta-analytic fixed- and random-effects models with or without moderators via the linear (mixed-effects) model

SElogratio= sqrt((DATA_TABLE$sdb_boot^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda_boot^2/DATA_TABLE$y^2))
SElogratio_notBoot=sqrt((DATA_TABLE$sdb^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda^2/DATA_TABLE$y^2))

TEST.TABLE$Ratio=with(TEST.TABLE, (post.rec/pre.rec))



m="weightedReg"        
  if(nrow(TEST.TABLE)>PhC_threshold){     # here I indicate the minimum number of cells needed to apply the regression


mmvar<-rma.uni(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,mods=log(TEST.TABLE$Ratio),method="ML",control=list(threshold=rmaTHR)) 
#yi	-vector of length k with the observed effect sizes or outcomes. See Details.
#vi	- vector of length k with the corresponding sampling variances. See Details.

SEFFPvalue=mmvar$pval[2]
SEFFZtest=mmvar$zval[2]
EST_SEFF=   (exp(mmvar$b[2])-1)*100
SE_SEFF=    (exp(mmvar$se[2] )-1)*100

estimate_SEFF=mmvar$b[2]


CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot


# post-hoc correction is done when sampling effort has a significant effect, even after rarefaction/extrapolation
# to be conservative we consider that such bias occurs (and should be corrected) 
# whenever the positive relation between sampling effort change and richness change as more than 94% chance to be real (i.e. P-value<0.06&estimate_SEFF>0)  
# We alert that based on several simulation and analyses with rela data we noticed that when the number of cells is very small,
# the slope could be highly influenced by an outlier, potentially leading to a misleading slope estimate 

if (estimate_SEFF<0|SEFFPvalue>0.06){
CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot
 

EST_SEFF=   0
SE_SEFF=   0
estimate_SEFF=0

mmvar<-rma.uni(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,method="ML",control=list(threshold=rmaTHR)) 
  

} }else {  
mmvar<-rma(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,method="ML",control=list(threshold=rmaTHR)) 
#yi	-vector of length k with the observed effect sizes or outcomes. See Details.
#vi	- vector of length k with the corresponding sampling variances. See Details.

CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot
  SEFFPvalue="-"
SEFFZtest="-"

EST_SEFF=   0
SE_SEFF=   0
estimate_SEFF=0
}

    ########################
#CORRECTED VALUES


CORRECTED_VALUES=TEST.TABLEchangeRATlog-estimate_SEFF*log(TEST.TABLE$Ratio)
CORRECTED_RATIO=exp(CORRECTED_VALUES)
 if (S=="G10"&Weight==1){hist(CORRECTED_VALUES, main=paste(Ct,"_GR",Gr,"_P",P))}





if (TEST==T) {

Ztest=mmvar$zval[1]       #intercept
Pvalue=mmvar$pval[1]
EST_CHANGE= (exp(mmvar$b[1])-1)*100 
LOG_RATIO_SE= mmvar$se[1]
}  else {
Ztest="-"       #intercept
Pvalue="-"
EST_CHANGE= (exp(mmvar$b[1])-1)*100        
LOG_RATIO_SE= mmvar$se[1]}  
         



Corrected_SDup=CORRECTED_VALUES+CorrectedSD
Corrected_SDlow=CORRECTED_VALUES-CorrectedSD



PERC= (exp(CORRECTED_VALUES)-1)*100

COUNTRY=rep(C,nrow(DATA_TABLE))
GROUP=rep(G,nrow(DATA_TABLE))
PERIOD=rep(P,nrow(DATA_TABLE))
SCALE=rep(S,nrow(DATA_TABLE))
TRAIT=rep(M,nrow(DATA_TABLE))


METHOD=rep("Extrap",nrow(DATA_TABLE))

ZTEST= rep(0,nrow(DATA_TABLE))
for (i in 1:nrow(DATA_TABLE)) { ZTEST[i]=ztest(CORRECTED_VALUES[i], (CorrectedSD[i]))}

ZTEST_noBootSD= rep(0,nrow(DATA_TABLE))
for (i in 1:nrow(DATA_TABLE)) { ZTEST_noBootSD[i]=ztest(CORRECTED_VALUES[i], (CorrectedSD_notboot[i]))}

NEW_TABLE=cbind(COUNTRY,GROUP,TRAIT,PERIOD,SCALE,METHOD,DATA_TABLE,CORRECTED_VALUES,CorrectedSD,ZTEST,CORRECTED_RATIO,
PERC,ZTEST_noBootSD, SElogratio,SElogratio_notBoot,CorrectedSD_notboot,
DATA_TABLE$sdb_boot,DATA_TABLE$sdb,DATA_TABLE$sda_boot,DATA_TABLE$sda,m,Weight.code)

if (CORRECTED.FILE==T ){
write.table(NEW_TABLE,file=corrected.file,col.names=F, append=T,sep="\t")
             
          }                 

       




mmvar2<-rma(yi=NEW_TABLE$CORRECTED_VALUES,vi=1/(NEW_TABLE[,28+Weight]),method="ML") 
if (TEST==T) {

Ztest=mmvar$zval[1]       #intercept  test values
Pvalue=mmvar$pval[1]
EST_CHANGE= (exp(mmvar$b[1])-1)*100 
LOG_RATIO_SE= mmvar$se[1]
}  else {
Ztest="-"       
Pvalue="-"
EST_CHANGE= (exp(mmvar$b[1])-1)*100
LOG_RATIO_SE= mmvar$se[1]}        


#calculating overall 95CI based on sd or VAR
WeightZERO_ONE=(NEW_TABLE[,28+Weight])/sum(NEW_TABLE[,28+Weight])
dados<-data.frame(m=CORRECTED_VALUES,w=WeightZERO_ONE)
bootfun<-function(d,i){d2<-d[i,];sum(d2$m*d2$w)/sum(d2$w)}     # function of Tom van dooren
if( nrow(dados)>2000){bootobj<-boot(dados,bootfun,R=nrow(dados))}else {bootobj<-boot(dados,bootfun,R=2000)}    # R indicates the number of resamples
Weighted_Corrected_est=bootobj$t0

# these lines account for the fact that if std. error of bootobj is Zero 
#all generated values are equal and boot.ci function doesn't work    
if(length(unique(bootobj$t))==1){            
Weighted_Corrected_95CIup=Weighted_Corrected_est
Weighted_Corrected_95CIlow=Weighted_Corrected_est
Weighted_Corrected_95CIupBACK_T=Weighted_Corrected_est
Weighted_Corrected_95CIlowBACK_T=Weighted_Corrected_est
}else{

Weighted_Corrected_95CIup=boot.ci(bootobj,type="bca")$bca[5]
Weighted_Corrected_95CIlow=boot.ci(bootobj,type="bca")$bca[4]

Weighted_Corrected_95CIupBACK_T=(exp(boot.ci(bootobj,type="bca")$bca[5])-1)*100
Weighted_Corrected_95CIlowBACK_T=(exp(boot.ci(bootobj,type="bca")$bca[4])-1)*100
}


Nb_cells_decline=nrow(NEW_TABLE[NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_decline=nrow(NEW_TABLE[NEW_TABLE$ZTEST<0.05&NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_increase=nrow(NEW_TABLE[NEW_TABLE$ZTEST<0.05&NEW_TABLE$CORRECTED_VALUES>0,])

Nb_cells_SIGN_decline_noboot=nrow(NEW_TABLE[NEW_TABLE$ZTEST_noBootSD<0.05&NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_increase_noboot=nrow(NEW_TABLE[NEW_TABLE$ZTEST_noBootSD<0.05&NEW_TABLE$CORRECTED_VALUES>0,])
    
                                 }


########################
     
 EXT.result <- data.frame(Ct,Gr,Pr,S,M,"-",
 EST_CHANGE,   
 paste("Ext_BOOT_mean_W", Weight.code),Nb_cells,Min.PROP,Min.REC, Min.SPS, RATIO, DIF,TOL,
Nb_cells_decline,
Nb_cells_SIGN_decline,
Nb_cells_SIGN_increase,
Weighted_Corrected_95CIlow,
Weighted_Corrected_95CIup,
Weighted_Corrected_95CIlowBACK_T,
Weighted_Corrected_95CIupBACK_T,
Ztest,
Pvalue,
SEFFZtest,
SEFFPvalue ,m,LOG_RATIO_SE

)
write.table(EXT.result, file=output.file, col.names=F, append=T,sep="\t")
   



DATA_TABLE=E3xdatarel2
  WEIGHTS.TABLE=cbind(DATA_TABLE$weight1,DATA_TABLE$weight2,DATA_TABLE$weight3,DATA_TABLE$weight4) 


 PRE.VECTOR = as.matrix(DATA_TABLE$x)
 PRECellIDVECTOR=as.matrix(DATA_TABLE$CellID )
 PRErecIDVECTOR=as.matrix(DATA_TABLE$pre.rec )
PREWEIvector=as.matrix(DATA_TABLE[,22+Weight])

Period=rep("pre",length(PRE.VECTOR ))
PRE.TABLE=cbind( Period,CellID=PRECellIDVECTOR,TAU=as.numeric(PRE.VECTOR),pre.rec= PRErecIDVECTOR,WEIvector=PREWEIvector)
 PRE.TABLE[,3]= as.numeric(PRE.TABLE[,3])




# create a POST-period vector (POSTV) of  POSTEi, 
#composed by Ri values taken from a normal distribution of log (change) for each i cell, Ri being proportional to the weight of that cell i

 POST.VECTOR = as.matrix(DATA_TABLE$y)
 POSTCellIDVECTOR=as.matrix(DATA_TABLE$CellID )
 POSTrecIDVECTOR=as.matrix(DATA_TABLE$post.rec )
POSTWEIvector=as.matrix(DATA_TABLE[,22+Weight])

Period=rep("pre",length(POST.VECTOR ))
POST.TABLE=cbind( Period,CellID=POSTCellIDVECTOR,TAU=as.numeric(POST.VECTOR),post.rec= POSTrecIDVECTOR,WEIvector=POSTWEIvector)
 POST.TABLE[,3]= as.numeric(POST.TABLE[,3])



if(nrow(PRE.TABLE)==1){TEST.TABLEa=as.data.frame(PRE.TABLE ) }else{
TEST.TABLEa=as.data.frame(PRE.TABLE,row.names=F ) }

PRE.TAU=  as.numeric(PRE.TABLE[,3])
PRE.TAU[PRE.TAU<0.1]=0.1
POST.TAU=  as.numeric(POST.TABLE[,3])
POST.TAU[POST.TAU<0.1]=0.1

TEST.TABLEchangeDIF=POST.TAU-PRE.TAU
TEST.TABLEchangeRAT=(POST.TAU/PRE.TAU -1)
TEST.TABLEchangeRATlog=log(POST.TAU/PRE.TAU)

TEST.TABLE=cbind(TEST.TABLEa[,c(1:2)],WEIvector=as.numeric(PRE.TABLE[,5]),TEST.TABLEchangeDIF,TEST.TABLEchangeRAT,TEST.TABLEchangeRATlog,
post.rec=as.numeric(POST.TABLE[,4]),pre.rec=as.numeric(PRE.TABLE[,4]))


TEST.TABLE1=rbind(PRE.TABLE,POST.TABLE)
TAU=as.numeric(TEST.TABLE1[1,3])
TEST.TABLE2=as.data.frame(TEST.TABLE1,row.names=F)





if ( nrow(TEST.TABLE)==1){


Nb_cells_decline="-"
Nb_cells_SIGN_decline="-"
Nb_cells_SIGN_increase="-"

if(nrow(TEST.TABLE)==1&DATA_TABLE$sdb==0&DATA_TABLE$sda==0)
{
EST=log(DATA_TABLE$z)
SEnoBoot=0.00001   # to avoid NaN when calculating pvalues
SE= 0.00001
LOG_RATIO_SE=0.00001

}else          {
EST=log(DATA_TABLE$z)
SEnoBoot= sqrt((DATA_TABLE$sdb^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda^2/DATA_TABLE$y^2))
SE= sqrt((DATA_TABLE$sdb_boot^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda_boot^2/DATA_TABLE$y^2))
LOG_RATIO_SE=SE

                }
EST_low95CI = (exp(EST-1.96*SE)-1)*100
EST_up95CI =  ((EST+1.96*SE)-1)*100
 EST_CHANGE=(exp(EST)-1)*100 
 
 
SEFFZtest="-"
SEFFPvalue="-"
     
Weighted_Corrected_est="-"
Weighted_Corrected_95CIup="-"
Weighted_Corrected_95CIlow="-"
   
Pvalue=ztest(EST,SE)    #               
Ztest=((EST)/SE)
Weighted_Corrected_95CIupBACK_T=EST_up95CI
Weighted_Corrected_95CIlowBACK_T=EST_low95CI

          }else {
      
#############################
# tests for overall effects #
#############################



# a. meta-analysis 
# fit the meta-analytic fixed- and random-effects models with or without moderators via the linear (mixed-effects) model

SElogratio= sqrt((DATA_TABLE$sdb_boot^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda_boot^2/DATA_TABLE$y^2))
SElogratio_notBoot=sqrt((DATA_TABLE$sdb^2/(DATA_TABLE$x^2))+(DATA_TABLE$sda^2/DATA_TABLE$y^2))

TEST.TABLE$Ratio=with(TEST.TABLE, (post.rec/pre.rec))




 m= "weightedReg"   
 if(nrow(TEST.TABLE)>PhC_threshold){     # here I indicate the minimum number of cells needed to apply the regression


mmvar<-rma.uni(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,mods=log(TEST.TABLE$Ratio),method="ML",control=list(threshold=rmaTHR)) 
#yi	-vector of length k with the observed effect sizes or outcomes. See Details.
#vi	- vector of length k with the corresponding sampling variances. See Details.

SEFFPvalue=mmvar$pval[2]
SEFFZtest=mmvar$zval[2]
EST_SEFF=   (exp(mmvar$b[2])-1)*100
SE_SEFF=    (exp(mmvar$se[2] )-1)*100

estimate_SEFF=mmvar$b[2]


CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot


# post-hoc correction is done when sampling effort has a significant effect, even after rarefaction/extrapolation
# to be conservative we consider that such bias occurs (and should be corrected) 
#whenever the positive relation between sampling effort change and richness change as more than 94% chance to be real (i.e. P-value<0.06&estimate_SEFF>0)  
# We alert that based on several simulation and analyses with rela data we noticed that when the number of cells is very small,
#the slope could be highly influenced by an outlier, potentially leading to a misleading slope estimate 

if (estimate_SEFF<0|SEFFPvalue>0.06){
CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot

EST_SEFF=   0
SE_SEFF=   0
estimate_SEFF=0

mmvar<-rma.uni(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,method="ML",control=list(threshold=rmaTHR))
  

} }else {  
mmvar<-rma(yi=TEST.TABLEchangeRATlog,vi=1/TEST.TABLE$WEIvector,method="ML",control=list(threshold=rmaTHR)) 
#yi	-vector of length k with the observed effect sizes or outcomes. See Details.
#vi	- vector of length k with the corresponding sampling variances. See Details.

CorrectedSD=(SElogratio) 
  CorrectedSD_notboot=  SElogratio_notBoot
  SEFFPvalue="-"
SEFFZtest="-"

EST_SEFF=   0
SE_SEFF=   0
estimate_SEFF=0
}

    ########################
#CORRECTED VALUES


CORRECTED_VALUES=TEST.TABLEchangeRATlog-estimate_SEFF*log(TEST.TABLE$Ratio)
CORRECTED_RATIO=exp(CORRECTED_VALUES)
 if (S=="G10"&Weight==1){hist(CORRECTED_VALUES, main=paste(Ct,"_GR",Gr,"_P",P))}





if (TEST==T) {

Ztest=mmvar$zval[1]       #intercept
Pvalue=mmvar$pval[1]
EST_CHANGE= (exp(mmvar$b[1])-1)*100 
LOG_RATIO_SE= mmvar$se[1]
}  else {"-"    # variance from Colwell formulas
Ztest="-"       #intercept
Pvalue="-"
EST_CHANGE= (exp(mmvar$b[1])-1)*100        
LOG_RATIO_SE= mmvar$se[1]}  
         



Corrected_SDup=CORRECTED_VALUES+CorrectedSD
Corrected_SDlow=CORRECTED_VALUES-CorrectedSD




PERC= (exp(CORRECTED_VALUES)-1)*100

COUNTRY=rep(C,nrow(DATA_TABLE))
GROUP=rep(G,nrow(DATA_TABLE))
PERIOD=rep(P,nrow(DATA_TABLE))
SCALE=rep(S,nrow(DATA_TABLE))
TRAIT=rep(M,nrow(DATA_TABLE))


METHOD=rep("Extrap3x",nrow(DATA_TABLE))

ZTEST= rep(0,nrow(DATA_TABLE))
for (i in 1:nrow(DATA_TABLE)) { ZTEST[i]=ztest(CORRECTED_VALUES[i], (CorrectedSD[i]))}

ZTEST_noBootSD= rep(0,nrow(DATA_TABLE))
for (i in 1:nrow(DATA_TABLE)) { ZTEST_noBootSD[i]=ztest(CORRECTED_VALUES[i], (CorrectedSD_notboot[i]))}

NEW_TABLE=cbind(COUNTRY,GROUP,TRAIT,PERIOD,SCALE,METHOD,DATA_TABLE,CORRECTED_VALUES,CorrectedSD,ZTEST,CORRECTED_RATIO,
PERC,ZTEST_noBootSD, SElogratio,SElogratio_notBoot,CorrectedSD_notboot,
DATA_TABLE$sdb_boot,DATA_TABLE$sdb,DATA_TABLE$sda_boot,DATA_TABLE$sda,m,Weight.code)

 
if (CORRECTED.FILE==T){
write.table(NEW_TABLE,file=corrected.file,col.names=F, append=T,sep="\t")
             
  }
                             
 
#calculating overall 95CI based on sd or VAR
WeightZERO_ONE=(NEW_TABLE[,28+Weight])/sum(NEW_TABLE[,28+Weight])
dados<-data.frame(m=CORRECTED_VALUES,w=WeightZERO_ONE)
bootfun<-function(d,i){d2<-d[i,];sum(d2$m*d2$w)/sum(d2$w)}
if( nrow(dados)>2000){bootobj<-boot(dados,bootfun,R=nrow(dados))}else {bootobj<-boot(dados,bootfun,R=2000)}
Weighted_Corrected_est=bootobj$t0

# these lines account for the fact that if std. error of bootobj is Zero 
#all generated values are equal and boot.ci function doesn't work  
if(length(unique(bootobj$t))==1){         
Weighted_Corrected_95CIup=Weighted_Corrected_est
Weighted_Corrected_95CIlow=Weighted_Corrected_est
Weighted_Corrected_95CIupBACK_T=Weighted_Corrected_est
Weighted_Corrected_95CIlowBACK_T=Weighted_Corrected_est
}else{

Weighted_Corrected_95CIup=boot.ci(bootobj,type="bca")$bca[5]
Weighted_Corrected_95CIlow=boot.ci(bootobj,type="bca")$bca[4]

#Weighted_Corrected_estBACK_T=(exp(bootobj$t0)-1)*100
Weighted_Corrected_95CIupBACK_T=(exp(boot.ci(bootobj,type="bca")$bca[5])-1)*100
Weighted_Corrected_95CIlowBACK_T=(exp(boot.ci(bootobj,type="bca")$bca[4])-1)*100
 }
 
##########################

Nb_cells_decline=nrow(NEW_TABLE[NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_decline=nrow(NEW_TABLE[NEW_TABLE$ZTEST<0.05&NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_increase=nrow(NEW_TABLE[NEW_TABLE$ZTEST<0.05&NEW_TABLE$CORRECTED_VALUES>0,])

Nb_cells_SIGN_decline_noboot=nrow(NEW_TABLE[NEW_TABLE$ZTEST_noBootSD<0.05&NEW_TABLE$CORRECTED_VALUES<0,])
Nb_cells_SIGN_increase_noboot=nrow(NEW_TABLE[NEW_TABLE$ZTEST_noBootSD<0.05&NEW_TABLE$CORRECTED_VALUES>0,])
  


Graph=F
if (Weight==4&length(Scale.set)==1){Graph=T} else if (Weight==4&G<10&S==Scale.set[1]|(Weight==4&C>1&G>9&S==Scale.set[2])) {Graph=T}    else if (Weight==4&C==1&S==Scale.set[1])  {Graph=T}

if (Graph==T){
plot(NEW_TABLE$CORRECTED_VALUES~log(NEW_TABLE$post.rec/NEW_TABLE$pre.rec),data=NEW_TABLE,main=paste(Ct,"-",Gr,"-",TR,"_",Pr,"_",S),xlab="N records", ylab="log(Richness change ratio)")
lines(c(-50,50),c(0,0))
lines(c(0,0),c(-50,50))


}
                }
########################

 EXT3x.result <- data.frame(Ct,Gr,Pr,S,M,"-",
 EST_CHANGE,  
 paste("Ext3x_BOOT_mean_W", Weight.code),Nb_cells,Min.PROP,Min.REC, Min.SPS, RATIO, DIF,TOL,
Nb_cells_decline,
Nb_cells_SIGN_decline,
Nb_cells_SIGN_increase,
Weighted_Corrected_95CIlow,
Weighted_Corrected_95CIup,
Weighted_Corrected_95CIlowBACK_T,
Weighted_Corrected_95CIupBACK_T,
Ztest,
Pvalue,
SEFFZtest,
SEFFPvalue,m,LOG_RATIO_SE

)
write.table(EXT3x.result, file=output.file, col.names=F, append=T,sep="\t")


     
  
 

       }
  } 
  }}}
   }
}     
}# END OF TREND EXTRACTOR




  library (compiler)
 Trend<- cmpfun (trend.extractor)
 
 
