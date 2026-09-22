#######################################################################
#BRP selectivity
#######################################################################
GetBRPSel=function(year_id,iFlag)
{
	BRPSel_out=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	if(iFlag==-1)
	{
		BRPSel_out=BRPSel_input
	}
	if(iFlag==0)
	{
		for(season_id in 1:Nseason)
		{
			temp=matrix(NA,nrow=Nfleet,ncol=Nsizebin)
			for(fleet_id in 1:Nfleet)
			{
				itemp=Nyear*Nseason*(fleet_id-1)+Nyear*(season_id-1)
				temp[fleet_id,]=report$Fleet_Selectivity[itemp+year_id,]      #(fleet,season,year) by size
			}#end fleet_id
			BRPSel_out[season_id,]=apply(temp,2,function(x){mean(x,na.rm=T)})
			BRPSel_out[season_id,]=BRPSel_out[season_id,]/max(BRPSel_out[season_id,],na.rm=T)  #scale by the max sel
		}#end season_id
	}
	if(iFlag>=1 && iFlag<=Nfleet)
	{
		fleet_id=BRPSelectivityFlag
		for(season_id in 1:Nseason)
		{
			itemp=Nyear*Nseason*(fleet_id-1)+Nyear*(season_id-1)
			BRPSel_out[season_id,]=report$Fleet_Selectivity[itemp+year_id,]      #(fleet,season,year) by size
		}
	}
	BRPSel_out	
}

#######################################################################
#SPR, mtons per 1000 recruits
#######################################################################
SPRFunc=function(Ftot,year_id,year_id_sele,iModel,AvgFlag)
{
	BRPSel=GetBRPSel(year_id_sele,BRPSelectivityFlag)   #calculate BRPSel=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	dR=1.0  #unit recruitment because we assume the balance was reached!
	BRPNAS=rep(0,Nsizebin)
	dvR=0.0
	dvF=Ftot*FRInSeason
	BRPF=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	BRPZ=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	BRPS=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	
	for(season_id in 1:Nseason)
	{
		BRPF[season_id,]=dvF[season_id]*BRPSel[season_id,]
		if(AvgFlag==0){
			itemp=Nyear*(season_id-1)+year_id
  	        mm=report$Natural_Mortality[itemp,]  # (season,year) by size
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
					dvTemp=BRPNAS*BRPSSB_S*report$Maturity_Length[year_id,]*report$Weight_Length_mtPer1000Fish[year_id,]*report$FemaleProp_At_Size[year_id,]
					dvTemp=sum(dvTemp)
					dvR=dvR+dvTemp
				}
			}
			
			if(AvgFlag==0){
				iL2=GrowthMatrixBlockFlag[year_id,k+1]
			    iL1=(iL2-1)*Nseason+k
			    itemp=(iL1-1)*Nsizebin
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
				dvTemp=BRPNAS*BRPSSB_S*report$Maturity_Length[year_id,]*report$Weight_Length_mtPer1000Fish[year_id,]*report$FemaleProp_At_Size[year_id,]
			    dvTemp=sum(dvTemp)
			    dvR=dvR+dvTemp
			}
			
			if(AvgFlag==0){
			    iL2=GrowthMatrixBlockFlag[year_id,k+1]
			    iL1=(iL2-1)*Nseason+k
			    itemp=(iL1-1)*Nsizebin
			    BRPNAS=(BRPNAS*BRPS[k,])%*%report$Growth[itemp+1:Nsizebin,]
			 }
			 if(AvgFlag==1){
				itemp=(k-1)*Nsizebin
				BRPNAS=(BRPNAS*BRPS[k,])%*%report$Mean_Growth_Matrix_Cal[itemp+1:Nsizebin,]
			}
		}#end k
	}
	
	dvR	
}

#######################################################################
#YPR, mtons per 1000 recruits
#######################################################################
YPRFunc=function(Ftot,year_id,year_id_sele,iModel,AvgFlag)
{
	BRPSel=GetBRPSel(year_id_sele,BRPSelectivityFlag)   #calculate BRPSel=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	dR=1.0  #unit recruitment because we assume the balance was reached!
	BRPNAS=rep(0,Nsizebin)
	dvR=0.0
	dvF=Ftot*FRInSeason
	BRPF=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	BRPZ=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	BRPS=matrix(NA,nrow=Nseason,ncol=Nsizebin)
	
	for(season_id in 1:Nseason)
	{
		BRPF[season_id,]=dvF[season_id]*BRPSel[season_id,]
		if(AvgFlag==0){
		    itemp=Nyear*(season_id-1)+year_id
  	        mm=report$Natural_Mortality[itemp,]  # (season,year) by size
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
				
				dvTemp=(BRPF[k,]/BRPZ[k,])*BRPNAS*(1-BRPS[k,])*report$Weight_Length_mtPer1000Fish[year_id,]
				dvTemp=sum(dvTemp)
				dvR=dvR+dvTemp				
			}
			
			if(AvgFlag==0){
			    iL2=GrowthMatrixBlockFlag[year_id,k+1]
			    iL1=(iL2-1)*Nseason+k
			    itemp=(iL1-1)*Nsizebin
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
						
			dvTemp=(BRPF[k,]/BRPZ[k,])*BRPNAS*(1-BRPS[k,])*report$Weight_Length_mtPer1000Fish[year_id,]
			dvTemp=sum(dvTemp)
			dvR=dvR+dvTemp
						
			if(AvgFlag==0){
			    iL2=GrowthMatrixBlockFlag[year_id,k+1]
			    iL1=(iL2-1)*Nseason+k
			    itemp=(iL1-1)*Nsizebin
			    BRPNAS=(BRPNAS*BRPS[k,])%*%report$Growth[itemp+1:Nsizebin,]
			}
			if(AvgFlag==1){
				itemp=(k-1)*Nsizebin
				BRPNAS=(BRPNAS*BRPS[k,])%*%report$Mean_Growth_Matrix_Cal[itemp+1:Nsizebin,]
			}
		}#end k
	}
	
	dvR	
}

#######################################################################
#equilibrium SSB, mtons
#######################################################################
EquiSSBFun=function(dvSPR,iRS)
{
	if(iRS==1)    #annual average R independent of spawner
	{
		alpha=par_est$PbnvRecruitPara_1   #std_est$est$PbnvRecruitPara[1]
		dvR=alpha*dvSPR
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

#######################################################################
#SPRatio Bisection
#######################################################################
SPRBisection=function(dF1,dF2,dInterval,dTarget,year_id,year_id_sele,iStep,iModel,AvgFlag)
{
	x1=SPRFunc(dF1,year_id,year_id_sele,iModel,AvgFlag)/SPRFunc(0,year_id,year_id_sele,iModel,AvgFlag)
	x1=x1-dTarget
	x2=SPRFunc(dF2,year_id,year_id_sele,iModel,AvgFlag)/SPRFunc(0,year_id,year_id_sele,iModel,AvgFlag)
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
			xmid=SPRFunc(fmid,year_id,year_id_sele,iModel,AvgFlag)/SPRFunc(0,year_id,year_id_sele,iModel,AvgFlag)
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

#######################################################################
#YPR Bisection
#######################################################################
YPRBisection=function(dF1,dF2,dInterval,dTarget,year_id,year_id_sele,iStep,iModel,AvgFlag)
{
	delta=0.005
	
	ds1=YPRFunc(0,year_id,year_id_sele,iModel,AvgFlag)
	ds2=YPRFunc(0+delta,year_id,year_id_sele,iModel,AvgFlag)
	slopeOrigin=(ds2-ds1)/delta
	slopeTarget=slopeOrigin*dTarget
	
	ds1=YPRFunc(dF1,year_id,year_id_sele,iModel,AvgFlag)
	ds2=YPRFunc(dF1+delta,year_id,year_id_sele,iModel,AvgFlag)
	x1=(ds2-ds1)/delta
	x1=x1-slopeTarget
	
	ds1=YPRFunc(dF2,year_id,year_id_sele,iModel,AvgFlag)
	ds2=YPRFunc(dF2+delta,year_id,year_id_sele,iModel,AvgFlag)
	x2=(ds2-ds1)/delta
	x2=x2-slopeTarget
		
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
			ds1=YPRFunc(fmid,year_id,year_id_sele,iModel,AvgFlag)
			ds2=YPRFunc(fmid+delta,year_id,year_id_sele,iModel,AvgFlag)
			xmid=(ds2-ds1)/delta
			xmid=xmid-slopeTarget
			
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

#######################################################################
#MSY Bisection
#######################################################################
MSYBisection=function(dF1,dF2,dInterval,dTarget,year_id,year_id_sele,iStep,iModel,iRS,AvgFlag)
{
	delta=0.005
	
	temp1=SPRFunc(dF1,year_id,year_id_sele,iModel,AvgFlag)
	temp2=EquiSSBFun(temp1,iRS)
	temp2=temp2/temp1
	ds1=YPRFunc(dF1,year_id,year_id_sele,iModel,AvgFlag)
	ds1=ds1*temp2        #total yield
	
	temp1=SPRFunc(dF1+delta,year_id,year_id_sele,iModel,AvgFlag)
	temp2=EquiSSBFun(temp1,iRS)
	temp2=temp2/temp1
	ds2=YPRFunc(dF1+delta,year_id,year_id_sele,iModel,AvgFlag)
	ds2=ds2*temp2        #total yield
	
	x1=(ds2-ds1)/delta
	x1=x1-dTarget
	
	temp1=SPRFunc(dF2,year_id,year_id_sele,iModel,AvgFlag)
	temp2=EquiSSBFun(temp1,iRS)
	temp2=temp2/temp1
	ds1=YPRFunc(dF2,year_id,year_id_sele,iModel,AvgFlag)
	ds1=ds1*temp2        #total yield
	
	temp1=SPRFunc(dF2+delta,year_id,year_id_sele,iModel,AvgFlag)
	temp2=EquiSSBFun(temp1,iRS)
	temp2=temp2/temp1
	ds2=YPRFunc(dF2+delta,year_id,year_id_sele,iModel,AvgFlag)
	ds2=ds2*temp2        #total yield
	
	x2=(ds2-ds1)/delta
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
			temp1=SPRFunc(fmid,year_id,year_id_sele,iModel,AvgFlag)
	        temp2=EquiSSBFun(temp1,iRS)
	        temp2=temp2/temp1
	        ds1=YPRFunc(fmid,year_id,year_id_sele,iModel,AvgFlag)
	        ds1=ds1*temp2        #total yield
	
	        temp1=SPRFunc(fmid+delta,year_id,year_id_sele,iModel,AvgFlag)
	        temp2=EquiSSBFun(temp1,iRS)
	        temp2=temp2/temp1
	        ds2=YPRFunc(fmid+delta,year_id,year_id_sele,iModel,AvgFlag)
	        ds2=ds2*temp2        #total yield
	
	        xmid=(ds2-ds1)/delta
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











