options(scipen=999)
library(randtests)
############################################## read in output from ADMB
source("r_ReadADMB.r")

filename="sssTrout"
report=read.admb(filename)
std_est=read.std(filename) #contain est and std for estimated paras
par_est=read.par(filename) #contain est for all paras inlcuding estimated and fixed ones

############################################## model information
fleet_name=c('ComLanding') #'RecLanding','ComDiscard','RecDiscard')
index_name=c('P915NorthSpring','P915NorthFall')
Nfleet=length(fleet_name)
Nindex=length(index_name)

Nyear=report$N_Year
Nseason=report$N_Season
Nsizebin=report$N_Sizebin
Nmonth_season=report$Month_In_Season #number of months in each season

year=seq(report$Begin_Year,report$End_Year,by=1) #1991,1992,...
year_axis=seq(report$Begin_Year,report$End_Year,by=2)
sizebin=(report$Size_Bins[1:Nsizebin]+report$Size_Bins[2:(Nsizebin+1)])*0.5
size_axis=seq(0,250,by=10)

# cold_year=c(1995,1999,2000,2002,2004,2009,2010,2013,2014,2017)  #cold stun years

############################################### calculate BRP
#### user input
ThreshSPR=0.2  #SPRatio for threshold
TargetSPR=0.3  #SPRatio for target
BRPSelectivityFlag=0  #-1-input; 0-avg fleet selectivity; #?#-fleet selectivity
BRPSel_input=matrix(1,nrow=Nseason,ncol=Nsizebin)   #input selectivity if flag=-1
BRPPeriod=50  #number of years to run to reach equalibrium for calculating BRP
BRPRefYear_Ini=year[length(year)]  #which year is used for biological params (growth,natual mortality,maturity,fecundity)
BRPRefYearSele_Ini=year[length(year)]  #which year is used for selectivity
FRInSeason=c(0.4,0.6)  #F proportions for each season
BRPModelFlag=0 #Type of model projection to get BRP; These two types would give the same answer:Type 0-follow one cohort through years to accumulate the ssb and yield it can contribute;Type 1-start with one cohort and project popn until equilibrium;then get ssb and yield the equi popn can contriubte
RSFlag=1  #1 Recruitment independent of spawner; 2: Beverton-Holt Model; 3: Ricker Model

DigitalPrec=0.000001  #precision of the calculation
Finterval=0.005   #the deltaF that is considered small enough to stop searching in optimizing function
Fmin=0     #min F to start searching
Fmax=10    #max F to stop searching
NumStep=20  #number of iterations in the searching

#### calculate BRP
### SPR and YPR are mtons per 1000 recruits for numerical convienence because
### The model units are 1000s fish or mton, the SR function R unit is 1000s fish
### If SPR or YPR values are too small, much smaller than the pre-set Finterval, the optimize (bisection) may not be able to find a solution
BRPRefYear=BRPRefYear_Ini-report$Begin_Year+1      #convert year to number 1,2,3
BRPRefYearSele=BRPRefYearSele_Ini-report$Begin_Year+1      #convert year to number 1,2,3
SSBSeason=report$Spawning_Season
FracSeasonSSB=report$Spawning_Month_Fraction_In_Season
RecruitmentSeaProp=report$Recruit_Season_Ratio
RecruitSizeProp=report$Recruit_Sizebin_Ratio
GrowthMatrixBlockFlag=report$Growth_Matrix_Block_Flag
GrowthMatrixBlockNum=report$Growth_Matrix_Block_Num   #number of blocks

source("r_FuncsBRP.r")

##################BRP summary
# BRPsummary=matrix(NA,nrow=Nyear+1,ncol=21)
# colnames(BRPsummary)=c('year','F','sdF','SSB','sdSSB','spr','sprF0','spratio','SSBF0','ypr','FThreshSpr','SSBThreshSpr','FTargetSpr','SSBTargetSpr','MSY','Fmsy','SSBmsy','Fmax','SSBfmax','F01','SSBf01')  #the unit for spr and ypr is mt per 1000 recruits, the unit for ssb is mt, -999 indicates the root does not exist; last row is for BRP calculated using avg M and growth matrix
# 
# for(yr_id in 1:(Nyear+1)){
# 	if(yr_id<=Nyear){
# 		BioAvgFlag=0
# 		yr_id_bio=yr_id
# 		yr_id_sel=yr_id		
# 		BRPsummary[yr_id,'year']=year[yr_id]
# 	    BRPsummary[yr_id,'F']=report$Fishing_Mortality_Annual[yr_id]
# 	    BRPsummary[yr_id,'sdF']=std_est$std$sdr_vAnnualF[yr_id]
# 	    BRPsummary[yr_id,'SSB']=report$Spawning_Stock_Pred[yr_id]
# 	    BRPsummary[yr_id,'sdSSB']=std_est$std$sdr_vSSB[yr_id]
# 	    BRPsummary[yr_id,'spr']=SPRFunc(BRPsummary[yr_id,'F'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 	    BRPsummary[yr_id,'sprF0']=SPRFunc(0,yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 	    BRPsummary[yr_id,'spratio']=BRPsummary[yr_id,'spr']/BRPsummary[yr_id,'sprF0']
# 	    BRPsummary[yr_id,'SSBF0']=EquiSSBFun(BRPsummary[yr_id,'sprF0'],RSFlag)	    
# 	    BRPsummary[yr_id,'ypr']=YPRFunc(BRPsummary[yr_id,'F'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 	}else{
# 		BioAvgFlag=1
# 		yr_id_bio=BRPRefYear
# 		yr_id_sel=BRPRefYearSele
# 		BRPsummary[yr_id,'year']=999   #the last row is for BRP calculated using avg M and growth matrix
# 	}
# 	
# 	####calculate SPRatio based BRP
# 	BRPsummary[yr_id,'FThreshSpr']=SPRBisection(Fmin,Fmax,Finterval,ThreshSPR,yr_id_bio,yr_id_sel,NumStep,BRPModelFlag,BioAvgFlag)
# 	BRPsummary[yr_id,'FTargetSpr']=SPRBisection(Fmin,Fmax,Finterval,TargetSPR,yr_id_bio,yr_id_sel,NumStep,BRPModelFlag,BioAvgFlag)
# 	
# 	if(BRPsummary[yr_id,'FThreshSpr']==-999)
# 	{
# 		BRPsummary[yr_id,'SSBThreshSpr']=-999
# 	}else{
# 		temp=SPRFunc(BRPsummary[yr_id,'FThreshSpr'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 		BRPsummary[yr_id,'SSBThreshSpr']=EquiSSBFun(temp,RSFlag)
# 	}
# 	if(BRPsummary[yr_id,'FTargetSpr']==-999)
# 	{
# 		BRPsummary[yr_id,'SSBTargetSpr']=-999
# 	}else{
# 		temp=SPRFunc(BRPsummary[yr_id,'FTargetSpr'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 		BRPsummary[yr_id,'SSBTargetSpr']=EquiSSBFun(temp,RSFlag)
# 	}
# 	
# 	####calculate YPR based BRP
# 	BRPsummary[yr_id,'Fmax']=YPRBisection(Fmin,Fmax,Finterval,0,yr_id_bio,yr_id_sel,NumStep,BRPModelFlag,BioAvgFlag)
# 	BRPsummary[yr_id,'F01']=YPRBisection(Fmin,Fmax,Finterval,0.1,yr_id_bio,yr_id_sel,NumStep,BRPModelFlag,BioAvgFlag)
# 	
# 	if(BRPsummary[yr_id,'Fmax']==-999)
# 	{
# 		BRPsummary[yr_id,'SSBfmax']=-999
# 	}else{
# 		temp=SPRFunc(BRPsummary[yr_id,'Fmax'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 		BRPsummary[yr_id,'SSBfmax']=EquiSSBFun(temp,RSFlag)
# 	}
# 	if(BRPsummary[yr_id,'F01']==-999)
# 	{
# 		BRPsummary[yr_id,'SSBf01']=-999
# 	}else{
# 		temp=SPRFunc(BRPsummary[yr_id,'F01'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 		BRPsummary[yr_id,'SSBf01']=EquiSSBFun(temp,RSFlag)
# 	}
# 	
# 	####calculate MSY based BRP
# 	BRPsummary[yr_id,'Fmsy']=MSYBisection(Fmin,Fmax,Finterval,0,yr_id_bio,yr_id_sel,NumStep,BRPModelFlag,RSFlag,BioAvgFlag)
# 	
# 	if(BRPsummary[yr_id,'Fmsy']==-999)
# 	{
# 		BRPsummary[yr_id,'SSBmsy']=-999
# 		BRPsummary[yr_id,'MSY']=-999
# 	}else{
# 		temp=SPRFunc(BRPsummary[yr_id,'Fmsy'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 		BRPsummary[yr_id,'SSBmsy']=EquiSSBFun(temp,RSFlag)
# 		
# 		temp0=YPRFunc(BRPsummary[yr_id,'Fmsy'],yr_id_bio,yr_id_sel,BRPModelFlag,BioAvgFlag)
# 		BRPsummary[yr_id,'MSY']=temp0*BRPsummary[yr_id,'SSBmsy']/temp
# 	}	
# }
# 
# write.table(BRPsummary,"BRP_summary.txt",row.names=F,sep='\t',quote=F)
# #BRPsummary=as.matrix(read.table('summary_BRP.txt',header=T,row.names=NULL,sep='\t'))

##############################BRP for F series for reference year
Fseries=seq(Fmin,Fmax,0.1)
BRPseries=matrix(NA,nrow=length(Fseries),ncol=4)
colnames(BRPseries)=c('F','SPRatio','YPR','Yield')

BioAvgFlag=0

for(i in 1:length(Fseries))
{
	BRPseries[i,'F']=Fseries[i]
	BRPseries[i,'SPRatio']=SPRFunc(BRPseries[i,'F'],BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)/SPRFunc(0,BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)
	BRPseries[i,'YPR']=YPRFunc(BRPseries[i,'F'],BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)
	temp=SPRFunc(BRPseries[i,'F'],BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)
	temp0=EquiSSBFun(temp,RSFlag)
	temp0=temp0/temp
	BRPseries[i,'Yield']=BRPseries[i,'YPR']*temp0
}

write.table(BRPseries,"BRP_series_ref.txt",row.names=F,sep='\t',quote=F)
#BRPseries=as.matrix(read.table('BRP_series_ref.txt',header=T,row.names=NULL,sep='\t'))

##############################BRP for F series for avg M and growth matrix
Fseries=seq(Fmin,Fmax,0.1)
BRPseries=matrix(NA,nrow=length(Fseries),ncol=4)
colnames(BRPseries)=c('F','SPRatio','YPR','Yield')

BioAvgFlag=1

for(i in 1:length(Fseries))
{
	BRPseries[i,'F']=Fseries[i]
	BRPseries[i,'SPRatio']=SPRFunc(BRPseries[i,'F'],BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)/SPRFunc(0,BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)
	BRPseries[i,'YPR']=YPRFunc(BRPseries[i,'F'],BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)
	temp=SPRFunc(BRPseries[i,'F'],BRPRefYear,BRPRefYearSele,BRPModelFlag,BioAvgFlag)
	temp0=EquiSSBFun(temp,RSFlag)
	temp0=temp0/temp
	BRPseries[i,'Yield']=BRPseries[i,'YPR']*temp0
}

write.table(BRPseries,"BRP_series_avg.txt",row.names=F,sep='\t',quote=F)
#BRPseries=as.matrix(read.table('BRP_series_avg.txt',header=T,row.names=NULL,sep='\t'))

################################################ plot output
source("r_FuncsPlot.r")
#BRPsummary=as.matrix(read.table('summary_BRP.txt',header=T,row.names=NULL,sep='\t'))

# pdf("plots.pdf",width=8,height=6)

###### plot obs and pred catch; specify fleet
PlotC(1)
# PlotC(2)
# PlotC(3)
# PlotC(4)

###### plot obs and pred indices; specify index and likelihood flag:1-robust normal;2-student t;3-log-normal distribution with obs term;4-og-normal distribution without obs term; 5-normal distribution;6-Cauchy distribution;7-uniform
PlotSI(1,3)
PlotSI(2,3)

###### plot obs and pred catch composition; sepcify fleet and season
PlotCC(1,1)
PlotCC(1,2)
# PlotCC(2,1)
# PlotCC(2,2)

###### plot obs and pred index composition; specify index
PlotSIC(1)
PlotSIC(2)

###### plot fleet selectivity; specify fleet
# #########################For time blocks
# Nblock_selF=c(4,4,4,4)  #num of blocks for each fleet

# ###begin year, begin season, block id; this will represent the selectivity of this block because the sele is the same within a block
# block_selF=list()
# block_selF[[1]]=matrix(c(year[1],1,1,2009,1,2,year[1],2,3,2009,2,4),nrow=Nblock_selF[1],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id
# block_selF[[2]]=matrix(c(year[1],1,5,2009,1,6,year[1],2,7,2009,2,8),nrow=Nblock_selF[2],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id
# block_selF[[3]]=matrix(c(year[1],1,9,2009,1,10,year[1],2,11,2009,2,12),nrow=Nblock_selF[3],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id
# block_selF[[4]]=matrix(c(year[1],1,13,2009,1,14,year[1],2,15,2009,2,16),nrow=Nblock_selF[4],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id

########################for no time blocks
# Nblock_selF=c(2,2,2,2)  #num of blocks for each fleet
# 
# ###begin year, begin season, block id; this will represent the selectivity of this block because the sele is the same within a block
# block_selF=list()
# block_selF[[1]]=matrix(c(year[1],1,1,year[1],2,2),nrow=Nblock_selF[1],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id
# block_selF[[2]]=matrix(c(year[1],1,3,year[1],2,4),nrow=Nblock_selF[2],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id
# block_selF[[3]]=matrix(c(year[1],1,5,year[1],2,6),nrow=Nblock_selF[3],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id
# block_selF[[4]]=matrix(c(year[1],1,7,year[1],2,8),nrow=Nblock_selF[4],ncol=3,byrow=T) #row is block,coln is starting year, plot season, block_id

PlotSelF(1)
# PlotSelF(2)
# PlotSelF(3)
# PlotSelF(4)

##### plot survey selectivity; specify index
Nblock_selSI=c(1,1)  #num of blocks for each survey

###begin year, index of the block; this will represent the selectivity of this block because the sele is the same within a block
# block_selSI=list()
# block_selSI[[1]]=matrix(c(2004,1),nrow=Nblock_selSI[1],ncol=2,byrow=T)  #row is block,coln is starting year,block_id
# block_selSI[[2]]=matrix(c(2003,2),nrow=Nblock_selSI[1],ncol=2,byrow=T)  #row is block,coln is starting year,block_id

PlotSelSI(1) 
PlotSelSI(2)

###### plot fishing mortality; specify fleet
PlotF(1)
# PlotF(2)
# PlotF(3)
# PlotF(4)
PlotF_season()
PlotF_all()
# PlotU_all() #exploitation rate U=1-exp(-F)

###### plot natural mortality
PlotM_Year()
# PlotM_base(3,1)  #specify flag and season; flag:0-input M and not plot; >0 estimate M and plot
# PlotM_base(3,2)
# PlotM_Input()  #plot M when flag=0 input

###### plot env variables fit to Mdev if any; specify flag and season; flag: 0-no env used to fit Mdev;1-yes, env is used to fit Mdev
# PlotEnvMdev(0,1)
# PlotEnvMdev(0,2)

######## plot N at size
PlotN_atSize()
PlotN_all()
PlotN_all_noR()

######## plot Recruitment and recruitment deviation
PlotR()

###### plot env variables fit to Rdev if any; specify flag: 0-no env used to fit Rdev;1-yes, env is used to fit Rdev
# PlotEnvRdev(0)

######## plot SSB
PlotSSB()

######## plot SR (spawner-recruit) pattern
PlotSR()

######## plot SPRatio
# PlotSPRatio()

################### plot stock status
######### thresh and target using biological parameters in a specified reference year
# fthresh_BioRefYr=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'FThreshSpr']     #report$F20SPR
# if(fthresh_BioRefYr==-999){
#   fthresh_BioRefYr=NA
# }
# ftarget_BioRefYr=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'FTargetSpr']     #report$F30SPR
# if(ftarget_BioRefYr==-999){
#   ftarget_BioRefYr=NA
# }  
# sthresh_BioRefYr=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'SSBThreshSpr']     #report$SSB20SPR
# if(sthresh_BioRefYr==-999){
#   sthresh_BioRefYr=NA
# }
# starget_BioRefYr=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'SSBTargetSpr']     #report$SSB30SPR
# if(starget_BioRefYr==-999){
#   starget_BioRefYr=NA
# }
# 
# fthresh_BioRefYr=rep(fthresh_BioRefYr,Nyear)
# ftarget_BioRefYr=rep(ftarget_BioRefYr,Nyear)
# sthresh_BioRefYr=rep(sthresh_BioRefYr,Nyear)
# starget_BioRefYr=rep(starget_BioRefYr,Nyear)
# 
# PlotStockStatus(fthresh_BioRefYr,ftarget_BioRefYr,sthresh_BioRefYr,starget_BioRefYr)
# 
# ######### thresh and target using year specific biological parameters
# fthresh_Yr=BRPsummary[1:Nyear,'FThreshSpr']
# fthresh_Yr[which(fthresh_Yr==-999)]=NA
# ftarget_Yr=BRPsummary[1:Nyear,'FTargetSpr'] 
# ftarget_Yr[which(ftarget_Yr==-999)]=NA
# sthresh_Yr=BRPsummary[1:Nyear,'SSBThreshSpr']
# sthresh_Yr[which(sthresh_Yr==-999)]=NA
# starget_Yr=BRPsummary[1:Nyear,'SSBTargetSpr']
# starget_Yr[which(starget_Yr==-999)]=NA
# 
# PlotStockStatus(fthresh_Yr,ftarget_Yr,sthresh_Yr,starget_Yr)
# PlotStockStatus_ratio(fthresh_Yr,sthresh_Yr,ThreshSPR*100)
# PlotStockStatus_ratio(ftarget_Yr,starget_Yr,TargetSPR*100)
# 
# ######### thresh and target using year-average biological parameters
# fthresh_BioAvg=BRPsummary[which(BRPsummary[,'year']==999),'FThreshSpr']     
# if(fthresh_BioAvg==-999){
#   fthresh_BioAvg=NA
# }
# ftarget_BioAvg=BRPsummary[which(BRPsummary[,'year']==999),'FTargetSpr']     
# if(ftarget_BioAvg==-999){
#   ftarget_BioAvg=NA
# }  
# sthresh_BioAvg=BRPsummary[which(BRPsummary[,'year']==999),'SSBThreshSpr']     
# if(sthresh_BioAvg==-999){
#   sthresh_BioAvg=NA
# }
# starget_BioAvg=BRPsummary[which(BRPsummary[,'year']==999),'SSBTargetSpr']     
# if(starget_BioAvg==-999){
#   starget_BioAvg=NA
# }
# 
# fthresh_BioAvg=rep(fthresh_BioAvg,Nyear)
# ftarget_BioAvg=rep(ftarget_BioAvg,Nyear)
# sthresh_BioAvg=rep(sthresh_BioAvg,Nyear)
# starget_BioAvg=rep(starget_BioAvg,Nyear)
# 
# PlotStockStatus(fthresh_BioAvg,ftarget_BioAvg,sthresh_BioAvg,starget_BioAvg)
# 
# ############ depletion level SSB/SSBF0
# ssbF0=BRPsummary[1:Nyear,'SSBF0']
# PlotDeplete(ssbF0,0.2,0.3)
# 
# ############ other curves, specify AvgFlag: 0-use M and growth for reference year; 1-use avg M adn growth matrix
# PlotSPRatio_FSeries(0)
# PlotYPR_FSeries(0)
# PlotYield_FSeries(0)
# 
# PlotSPRatio_FSeries(1)
# PlotYPR_FSeries(1)
# PlotYield_FSeries(1)

###################################### other plots
#### Plot growth of a cohort and visualize growth matrices
PlotGM(year[1],year[1],10) #specify start year, year of biological parameters to use in projection and how many years to grow
# PlotGM_blk(1,year[1],10) #speicfy block id, begin year for each block, number of years to grow; use begin year biological parameters in projection

### plot VB growth paras and curve; specify total number of blocks,flag,t0
# PlotVB_Linf(GrowthMatrixBlockNum,3) #VB paras estimates; specify num of blocks and flag; flag: 0-input GM,1-VB with blk,2-VB random walk, 3-avg with random deviation
# PlotVB_K(GrowthMatrixBlockNum,3) #VB paras estimates; specify num of blocks and flag; flag as above
# PlotVB(GrowthMatrixBlockNum,3,-0.98)  #VB for length-age relationship; specify num of blocks,flag and t0; here use t0 estimated externally; flag as above
# PlotVB_deltaL(GrowthMatrixBlockNum,3) #VB for deltaL-length relationship;specify num of blocks, flag as above

### plot length-weight relationships; specify year
# PlotWL()
# 
# ### plot maturity; specify year
# PlotML() # plot maturity
# 
# dev.off()

################################################### projection
# source("r_RunProjection.r")
#run_projection(report = report, par_est = par_est, select_id = c(5,3,4), M = 0.5,
#			    n_proj_year = 50, N_rep = 10000, F_trawl_proj = rep(1.0,50), F_trap_proj = rep(0.02,50), R_proj = rep(500,50), growth_id = c(1,2), R_rand = T, R_logsd = 0.5)













