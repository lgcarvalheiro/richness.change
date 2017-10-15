# richness.change
 scripts to estimate richness change  based on historical data

Luísa G. Carvalheiro


Estimates of species richness change based on historical records require a methodology that corrects for the effects of sampling effort (e.g. number of records per grid cell and time period), as well as bias associated with collectors (e.g. preference for rare species, under-representation of singletons due to efforts to capture differences due to sexual dimorphism). The scripts ‘Multilevel.RAR_EXTR’ and  ‘Trend.extractor’ allow to deal with such particularities of historical datasets when estimating richness changes. See example in the end of this file.
	‘Multilevel.RAR_EXTR’ allows to obtain estimates of richness based on species accumulation curves (providing estimates of expected richness with increasing number of records) using 3 approaches: (1) only interpolation, (2) only  extrapolation, and (3) combination of extrapolation and interpolation. All these methods allow comparison of species richness estimates among regions or between time periods with unequal sampling effort. However, the shape of the accumulation curve (hence estimates in earlier stages of the curve) can be very sensitive to certain collector’s sampling bias (e.g. overrepresentation of rare species). Therefore, extrapolation is preferable when dealing with historical datasets, but only if extrapolation is up to less than threefold of the real sampling effort (Colwell et al 2012). It is hence advisable to use a combination of extrapolation and interpolation (see Carvalheiro et al 2013 for an example). Where the more sampled period had more than threefold the number of records of the less sampled period, extrapolation and interpolation are combined  comparisons: extrapolating the smaller sample’s accumulation curve up to three times its sample size, and rarefying the larger sample down to this same ‘comparison point’. Such estimates are less sensitive to differences in the shape of the accumulation curve, than estimates based solely on interpolation. The uncertainty (i.e. standard deviation) associated with such estimates is calculated also based on Colwell et al. (2012). However, as in the absence of singletons, the assemblage is considered to be well sampled, and the standard deviation  will be zero (Colwell et al. 2012), to account for any under- or overestimation of singletons and doubletons (collectors may put effort on registering the rarest species, or try to obtain two specimens to capture morphological diversity of the species, e.g. males and females; flowering and fruiting plant specimens), bootstrapping is used to calculate standard deviation of our richness estimates, and consequently of the richness change estimate (i.e. 100 replicates, with equal total number of records, of the original community samples are generated, and the frequency of each species is computed assuming multinomial probabilities). This procedure creates a corrected estimator of the number of unseen species, where the effect of under/oversampling of singletons is reduced. 

----------------------------------------------------------------------------------------------------------
Description of the main columns of the output table of ‘Multilevel.RAR_EXTR’ :

spat.scale - id of spatial scale at which calculations were done, as indicated by the user (e.g. '10x10km cell', 'county', 'state')	

cell- id of the spatial location at which calculations were done, as indicated by the user (e.g. 'TQ28', 'Sao Paulo state')	

trait.g - id of the trait group for which the calculations were done, as indicated by the user (e.g. 'oligolectic', 'monolectic')

pre.rec- original sampling efforts in the pre-period 

post.rec - original sampling efforts in the post periods	

pre.tau.mean_extr_3x - measures of estimated richness in pre time period when using combination of extrapolation and interpolation 

post.tau.mean_extr_3x - measures of estimated richness in post time period when using combination of extrapolation and interpolation

Change.Mao_extr_3x - relative richness change estimate when using combination of extrapolation and interpolation

Change.Mao_extr_3x_CIlow - 95% confidence interval (lowest limit) of relative richness change estimate when using combination of extrapolation and interpolation

Change.Mao_extr_3x_CIup - 95% confidence interval (lowest limit) of relative richness change estimate when using combination of extrapolation and interpolation

---------------------------------------------------------------------------------------------------------


For datasets where richness change estimates are done for multiple geographical regions,‘Trend.extractor’ allows to check if accumulation curve estimations did completely remove the bias due to sampling effort. Whenever the difference in records between the two time periods has a significant effect on estimated richness change across sites, partial residuals are calculated after removing the effect of sampling effort for each cell to obtain unbiased estimates of richness change for each geographical location. This was done using the rma.uni function of the R package metaphor (Viechtbauer 2010), each grid cell was weighted based on the inverse of the variance, so that cells with more reliable estimations have a higher weight in the analyses (Hartung et al. 2008). 
This analyses also allowed to obtain an overall weighted value of richness change and assess if the mean value of change at a given scale was significantly different from zero.


If you have multiple geographical locations with variable data quality, it is advisable to set a priori  selection criteria. As an example, in Carvalheiro et al (2013) we only estimated richness change for cells that: (1) had more than 15 records per time period, (2) the ratio of records/number of species was higher than 1.5 in each of the two compared time periods, and (3) if number of records was less than five-fold the total number of species in the country, we only considered cells with less than a 10-fold difference in numbers of records between periods. Applying less strict selection criteria would have led to the selection of more geographical locations, but also to results that are highly influenced by individual cells. 

Please note that the 95% confidence interval represented in Fig 1 of Carvalheiro et al (2013) was calculated based on standard deviation of the log ratio obtained when using the package metafor (In the tables produced by Trend.extractor the column used was ‘Estimated_logRatio_SE’), and after the values (logratio and 95CI) were back-transformed. 

# References

Carvalheiro L.G., et al. 2013. Species richness declines and biotic homogenization have slowed down for NW-European pollinators and plants. Ecology Letters, 16, 870-878. 

Colwell, R.K., Chao, A., Gotelli, N.J., Lin, S., Mao, C.X., Chazdon, R.L., et al. (2012). Models and estimators linking individual-based and sample-based rarefaction, extrapolation and comparison of assemblages. J. Plant Ecol., 5, 3–21.

Hartung, J., Knapp, G. & Sinha, B.K. (2008). Statistical Meta-Analysis with Applications. John Wiley & Sons, Inc., Hoboken, New Jersey, p. 248.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. J. Stat. Softw., 36, 1–48.

# ---------------------------------------------------------
# EXAMPLE - HOW TO RUN THE ANALYSIS 
# 1. copy paste all Multilevel.RAR_EXTR and Trend.extractor functions into R
# 2. upload data

# 3.   Run MultilevelRAR function

	multilevel.RAR.buster (multi.RAR.data=Data, #inform name of original dataset
	taxon.indx=7,   # column nb with info on species name
	trait.indx=16,  # column nb that has info on traits, if no trait analyses will be done the column should cointain the same value throughout
	year.indx=3,    # column nb with info on year
	scale.data.indcs=c(1,12), # columns with information on cell ID. If analyses are to be repeated at several scales, indicate several columns
	period.1m=c(1950,1969),  # lower and upper limit of pre.period 
	period.2m=c(1970,1989),  # lower and upper limit of post.period
	output.file="C:/R/RichnessChange_per_cell.txt",  #path and name of the file that will be created
	min.recs=c(20),   #  minimum number of records per gridcell/period 
	MIN.PROP=0.2,      #proportion of the max number of sps per cell that will be used as minimum number of records
	FIXGRID=F,   # specifies if the data included in the analyses should only come from cells selected in smallest scale analyses or not 
	DIF=10)  # specifies the maximum diference between rec number in pre and pos. Ratio maxRecords/minRecords < DIF
                
 
# 4. upload table generated by MultilevelRARdata
# 5. Run Trend
 # 
 
	par(mfrow=c(6,4),bty="l", mar=c(2,5,2,1)+.2, oma=c(4,2,3,0), cex.axis=1.3,lwd=2) # creates pannels to visusalize how richness change values estimates are affected by sampling effort  
  
	Trend(Data= multilevel.RAR.DATA,   
	weight.SE=T,  
	Min.PROP=0.2, #proportion of the max number of sps per cell that will be used as minimum number of reocrds
	Min.REC=15,  # minimum number of records within a cell required for each time period
	Min.SPS=2, # minimum number of species within a cell required   for a cell to be selected
	RATIO=1.5, # specifies the minimum value of the ratio between nb.records/nb.species  required for a cell to be selected
	DIF=10,   # specifies the maximum diference between rec number in pre and pos. Ratio maxRecords/minRecords < DIF
	TOL=5,    # specifies that if the total number of records is higher than TOL * maximum number of species of that group, DIF will be ignored
	PhC_threshold=6,  #specifies the minimum number of cells required to check if multilevel.RAR corrected bias do to sampling effort, and apply a pos-hoc correction if needed
	output.file="C:/R/Average_RichnessChange_perScale.txt",    #AVERAGE CHANGE PER CELL FOR A GIVEN SCALE
	corrected.file="C:/R/RichnessChange_per_cell_Corrected_file.txt",       #CHANGE PER CELL after  testing (and correcting if needed) for additional effect of sampling effort
 	TEST=T,  # specifies that    output.file will be generated
 	CORRECTED.FILE=T)   # specifies that  corrected.file will be generated
