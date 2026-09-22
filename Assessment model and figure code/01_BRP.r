# Calculate Spawning Potential Ratio (SPR) biological reference points (BRPs)
# Based on BRP code developed by Yan Li (former NCDMF senior stock assessment scientist)
# Last modified by Matt Damiano (NC DMF Stock Assessment Program Manager) 09-22-2026

Nyear=report$N_year
Nseason=report$Time_step
Nsizebin=length(report$Size_bins)-1
Nmonth_season=6 #number of months in each season

year=seq(report$BeginYear,report$EndYear,by=1) #1991,1992,...
# year_axis=seq(report$Begin_Year,report$End_Year,by=2)
sizebin=size_bm

############################################### calculate BRP
#### user input
ThreshSPR=0.2  #SPRatio for threshold
TargetSPR=0.4  #SPRatio for target
# BRPSelectivityFlag=1  #-1-input; 0-avg fleet selectivity; #?#-fleet selectivity
# BRPSel_input=matrix(1,nrow=Nseason,ncol=Nsizebin)   #input selectivity if flag=-1
BRPPeriod= 50  #number of years to run to reach equalibrium for calculating BRP
BRPRefYear_Ini=year[length(year)]  #which year is used for biological params (growth,natual mortality,maturity,fecundity)
BRPRefYearSele_Ini=year[length(year)]  #which year is used for selectivity
# Ftot <- report$Fishing_Mortality # if annual
# FRInSeason = 1 # if annual
Ftot <- sum(report$Fishing_Mortality[1,n_t],report$Fishing_Mortality[2,n_t]) # if seasonal
FRInSeason=c((report$Fishing_Mortality[1,n_t])/Ftot,(report$Fishing_Mortality[2,n_t])/Ftot)  #F proportions for each season
BRPModelFlag=1 #Type of model projection to get BRP; These two types would give the same answer:Type 0-follow one cohort through years to accumulate the ssb and yield it can contribute;Type 1-start with one cohort and project popn until equilibrium;then get ssb and yield the equi popn can contriubte
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
BRPRefYear=BRPRefYear_Ini-report$BeginYear+1      #convert year to number 1,2,3
BRPRefYearSele=BRPRefYearSele_Ini-report$BeginYear+1      #convert year to number 1,2,3
SSBSeason= n_s 
FracSeasonSSB=1.0
RecruitmentSeaProp=report$Recrut_Season_ratio
RecruitSizeProp=report$Recrut_Sizebin_ratio
GrowthMatrixBlockFlag=report$GrowthMatrixBlockNum
GrowthMatrixBlockNum=n_growblock  #number of blocks
BRPSel <- array(NA,dim=c(Nseason,Nsizebin))
for(s in 1:Nseason){
  BRPSel[s,] <- report$Fleet_selectivity[1,]
}

# source("r_FuncsBRP.r")

#######################################################################
#SPR, mtons per 1000 recruits
#######################################################################
SPRFunc=function(Ftot,year_id,iModel,AvgFlag) # set defaults so you don't use time-varying code
{
  # for running components only
  # year_id = 43 # terminal year
  # iModel = 0
  # AvgFlag = 0
  #
  
	BRPSel= BRPSel # GetBRPSel(year_id_sele,BRPSelectivityFlag)   # set to commercial selectivity, which is static over seasons
	dR=1.0 #unit recruitment because we assume the balance was reached!
	# dR=exp(par_est$est$`PbnvLogRSPara[1]`)
	# dR=geomean(report$recruitment_Pred[39:43])
	BRPNAS=rep(0,Nsizebin)
	dvR=0.0
	dvF=Ftot*FRInSeason
	BRPF=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	BRPZ=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	BRPS=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	
	for(season_id in 1:Nseason)
	{
		BRPF[season_id,]=dvF[season_id]*BRPSel[season_id,] # flagging
		if(AvgFlag==0){
			itemp=Nyear*(season_id-1)+year_id
  	        mm=report$Natural_Mortality[season_id,]  # (season,year) by size - in Jie's model, it's just season x size
		}
		if(AvgFlag==1){
			mm=report$Mean_Natural_Mortality_Para_Size[season_id,]  #season by size
		}
		
		BRPZ[season_id,]=BRPF[season_id,]+mm
		BRPS[season_id,]=exp(-BRPZ[season_id,])
		if(season_id==SSBSeason)
		{
			BRPSSB_S=exp(-BRPZ[season_id,]*FracSeasonSSB)
		}
	}
	
	for(i in 1:BRPPeriod)
	{
		for(k in 1:Nseason)
		{
			if(iModel==1) #start with one cohort and project popn until equilibrium;then get ssb and yield the equi popn can contriubte
			{
				BRPNAS=BRPNAS+dR*RecruitmentSeaProp[k]*RecruitSizeProp
			}			
			if(iModel==0) #follow one cohort through years to get the ssb and yield it can contribute; these two methods should give the same answer
			{
				if(i==1)
				{
					BRPNAS=BRPNAS+dR*RecruitmentSeaProp[k]*RecruitSizeProp
				}
				if(k==SSBSeason)
				{
					dvTemp=BRPNAS*BRPSSB_S*report$Maturity_length[year_id,]*report$Weight_length[year_id,]*report$FemaleProp_at_Size[year_id,]
					dvTemp=sum(dvTemp)
					dvR=dvR+dvTemp
				}
			}
			
			if(AvgFlag==0){
			  itemp=(k-1)*Nsizebin
			    BRPNAS=(BRPNAS*BRPS[k,])%*%report$Growth[itemp+1:Nsizebin,]
			}
			if(AvgFlag==1){
				itemp=(k-1)*Nsizebin
				BRPNAS=(BRPNAS*BRPS[k,])%*%report$Mean_Growth_Matrix_Cal[itemp+1:Nsizebin,]
			}
			
		}#end k
	}#end i
	
	if(iModel==1)
	{
		for(k in 1:Nseason)
		{
			BRPNAS=BRPNAS+dR*RecruitmentSeaProp[k]*RecruitSizeProp
			
			if(k==SSBSeason)
			{
				dvTemp=BRPNAS*BRPSSB_S*report$Maturity_length[year_id,]*report$Weight_length[year_id,]*report$FemaleProp_at_Size[year_id,]
			    dvTemp=sum(dvTemp)
			    dvR=dvR+dvTemp
			}
			
			if(AvgFlag==0){
			  itemp=(k-1)*Nsizebin
			    BRPNAS=(BRPNAS*BRPS[k,])%*%report$Growth[itemp+1:Nsizebin,]
			 }
			 if(AvgFlag==1){
				itemp=(k-1)*Nsizebin
				BRPNAS=(BRPNAS*BRPS[k,])%*%report$Growth[itemp+1:Nsizebin,]
			}
		}#end k
	}
	dvR	
}
spr.f0 <- SPRFunc(0,n_t,1,0)
spr.term <- SPRFunc(Ftot,n_t,1,0)
spr.fterm <- spr.term/spr.f0
spr.fterm

#######################################################################
#equilibrium SSB, mtons
#######################################################################
EquiSSBFun=function(dvSPR,iRS=1,scen=1) # default is mean recruitment, default scenario is equilibrium recruitment
{
	if(iRS==1)    #annual average R independent of spawner
	{
	  if(scen==1){
	    alpha=par_est$est$`PbnvLogRSPara[1]`  # default #std_est$est$PbnvRecruitPara[1]
	    dvR=exp(alpha)*dvSPR
	  }
	  if(scen==2){
	    alpha = median(report$recruitment_Pred) # produce median equilibrium SSB e.g, Chesbay 2026 assessment
	    dvR = alpha*dvSPR
	  }
	  if(scen==3){
	    alpha = geomean(report$recruitment_Pred[39:43]) # assume productivity of last 5 years of the assessment
	    dvR = alpha*dvSPR
	  }
	}
	if(iRS==2)    #Beverton-Holt
	{
		alpha=par_est$PbnvRecruitPara_1   #std_est$est$PbnvRecruitPara[1]
		beta=par_est$PbnvRecruitPara_2    #std_est$est$PbnvRecruitPara[2]
		dvR=(alpha*dvSPR-1)/beta
	}
	if(iRS==3)    #Ricker
	{
		alpha=par_est$PbnvRecruitPara_1   #std_est$est$PbnvRecruitPara[1]
		beta=par_est$PbnvRecruitPara_2    #std_est$est$PbnvRecruitPara[2]
		dvR=(log(alpha)+log(dvSPR))/beta
	}
	
	dvR
}

ssb.f0 <- EquiSSBFun(spr.f0, 1)
ssb.fterm <- EquiSSBFun(spr.fterm, 1)

#######################################################################
#SPRatio Bisection
#######################################################################
SPRBisection=function(dF1,dF2,dInterval,dTarget,year_id,iStep,iModel,AvgFlag)
{
	x1=SPRFunc(dF1,year_id,iModel,AvgFlag)/SPRFunc(0,year_id,iModel,AvgFlag)
	x1=x1-dTarget
	x2=SPRFunc(dF2,year_id,iModel,AvgFlag)/SPRFunc(0,year_id,iModel,AvgFlag)
	x2=x2-dTarget
	
	if(abs(x1)<DigitalPrec)
	{
		return(dF1)
	}
	if(abs(x2)<DigitalPrec)
	{
		return(dF2)
	}
	if(x1*x2>0)
	{
		print('Root does not exist within the given range!')
		return(-999)
	}else{
		f1=dF1
		f2=dF2
		for(j in 1:iStep)
		{
			fmid=(f1+f2)/2
			xmid=SPRFunc(fmid,year_id,iModel,AvgFlag)/SPRFunc(0,year_id,iModel,AvgFlag)
			xmid=xmid-dTarget
			
			if(abs(xmid)<DigitalPrec)
			{
				rtb=fmid
				return(rtb)
			}
			if(abs(f1-f2)<dInterval)
			{
				rtb=fmid
				return(rtb)
			}
			if(xmid*x2<0)
			{
				f1=fmid
				x1=xmid
			}else{
				f2=fmid
				x2=xmid
			}
		}#end j
	}	
}

BRPsummary=matrix(NA,nrow=Nyear,ncol=2)
colnames(BRPsummary)=c('FThreshSpr','FTargetSpr')  #the unit for spr and ypr is mt per 1000 recruits, the unit for ssb is mt, -999 indicates the root does not exist; last row is for BRP calculated using avg M and growth matrix

for(yr_id in 1:(Nyear)){
  ####calculate SPRatio based BRP
  BRPsummary[yr_id,'FThreshSpr']=SPRBisection(Fmin,Fmax,Finterval,ThreshSPR,yr_id,NumStep,0,0)
  BRPsummary[yr_id,'FTargetSpr']=SPRBisection(Fmin,Fmax,Finterval,TargetSPR,yr_id,NumStep,0,0)
}


