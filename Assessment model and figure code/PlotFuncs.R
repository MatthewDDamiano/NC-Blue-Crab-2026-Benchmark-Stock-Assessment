
sizebins=report$Size_bins[-1]-0.5
Nyear=report$N_year

#windows(record=T)
#X11(,type="cairo")

#######################################################################
#weight-length
#######################################################################
PlotWL=function()
{
  Year_end=report$EndYear; Year_begin=report$BeginYear
  year = seq(Year_begin,Year_end,1)
  
  plot_nrow = ceiling(Nyear/6); plot_ncol=6
  yaxis = rep(1,plot_nrow)
  
  for(pl in 1:plot_nrow)
  {
    yaxis[pl] = yaxis[1]+(pl-1)*plot_ncol
  }
  
    xaxis = seq(Nyear-plot_ncol+1,Nyear,1)
  
  par(mfrow=c(plot_nrow,plot_ncol))
  
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  for (yi in 1:Nyear)
  {
    y=report$Weight_length[yi,]
    plot(sizebins,y,col='Black',type='l',lwd=1,pch=16,axes=FALSE,ylim=c(0,35))
    mtext(paste(year[yi]),side=3,line=-1,adj=0.1,cex=0.6)
    if(yi %in% yaxis)
      axis(2,col='gray',at=seq(0,30,10))
    if(yi %in% xaxis)
      axis(1,col='gray',at=sizebins[-1])
    box(col='gray')
  }
  mtext('Dorsal Carapace Length (mm)',side=1,outer=T,line=2.2)
  mtext('Mean Weight (g)',side=2,outer=T,line=2.2)
  mtext(paste('Length-Weight Relationships'),side=3,outer=T,line=2.2)
}
######################################################################

#######################################################################
#maturity-length
#######################################################################
PlotML=function()
{
  Year_end=report$EndYear; Year_begin=report$BeginYear
  year = seq(Year_begin,Year_end,1)
  
  plot_nrow = ceiling(Nyear/6); plot_ncol=6
  yaxis = rep(1,plot_nrow)
  
  for(pl in 1:plot_nrow)
  {
    yaxis[pl] = yaxis[1]+(pl-1)*plot_ncol
  }
  
  xaxis = seq(Nyear-plot_ncol+1,Nyear,1)
  
  par(mfrow=c(plot_nrow,plot_ncol))
  
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  for (yi in 1:Nyear)
  {
    y=report$Maturity_length[yi,]
    plot(sizebins,y,col='Black',type='o',lwd=1,pch=16,axes=FALSE,ylim=c(0,1),cex=0.5)
    mtext(paste(year[yi]),side=3,line=-1,adj=0.1,cex=0.6)
    if(yi %in% yaxis)
      axis(2,col='gray',at=seq(0,0.8,0.2))
    if(yi %in% xaxis)
      axis(1,col='gray',at=sizebins[-1])
    box(col='gray')
  }
  mtext('Dorsal Carapace Length (mm)',side=1,outer=T,line=2.2)
  mtext('Frequency',side=2,outer=T,line=2.2)
  mtext(paste('Maturity'),side=3,outer=T,line=2.2)
}
######################################################################

#######################################################################
#Growth matrix,specify the year, season and how many years for growth
#######################################################################
Plot_GM_tb = function(n_tb,iYearNum)
{
  plot_nrow = 1; plot_ncol=n_tb
  par(mfrow=c(plot_nrow,plot_ncol))
  
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  for (i in 1:n_tb)
  {
    PlotGM(i,iYearNum)
  }
  mtext('Year',side=1,outer=T,line=2.2)
  mtext('Dorsal Carapace Length (mm)',side=2,outer=T,line=2.2)
  mtext(paste('Growth of a cohort with the growth matrices of year'),side=3,outer=T,line=2.2)
}

PlotGM=function(yearblock,iYearNum)
{
  season=report$Time_step
  nsizebin=length(sizebins)
  Year=yearblock
  if(season>1)
  {
    itemp=season*nsizebin*(Year-1)
    GM=report$Growth[(itemp+1):(itemp+nsizebin*season),]
    Maintemp=paste("Growth of a cohort with the Growth matrices of year",Year)
    
    BRPNAS<-rep(0,nsizebin)
    BRPNASData<-matrix(NA,ncol=nsizebin,nrow=iYearNum)
    Recruitment<-1000
    RecPrjRatio<-report$Recrut_Sizebin_ratio
    BRPS<-exp(-report$Natural_Mortality[(season*(Year-1)+1):(season*(Year-1)+season),])
    
    for (i in 1:(iYearNum-1))
    {
      for (j in 1:season)
      {
        iL=nsizebin*(j-1)
        if(i==1&&j==1)
          BRPNASData[1,]<- BRPNAS <- BRPNAS+Recruitment*RecPrjRatio;
        NASTemp<-BRPNAS*BRPS[j,];
        BRPNAS<-NASTemp%*%GM[(iL+1):(iL+nsizebin),]
        BRPNASData[i+1,]<-BRPNAS;
      }
    }
    
    iYear=seq(1,iYearNum,1);RRRData=NULL;RMainData=c();
    
    for( i in 1:iYearNum)
    {
      
      RRData<-BRPNASData[i,]
      RRData<-RRData/sum(RRData)
      RMainData[i]<-sum(RRData*sizebins)
      for(j in 1:nsizebin)
      {
        SizeNum<-RRData[j]*1000  
        if( SizeNum<1)
          next;
        XXX1<-rep(i,SizeNum)
        XXX2<-rep(j,SizeNum)
        RRRDataTemp<-cbind(iYear[XXX1],sizebins[XXX2])
        RRRData<-rbind(RRRData,RRRDataTemp)
      }
    }
    
  }
  else 
  {
    itemp=nsizebin*(Year-1)
    GM=report$Growth[(itemp+1):(itemp+nsizebin),]
    Maintemp=paste("Growth of a cohort with the Growth matrix of year",Year)
    
    BRPNAS<-rep(0,nsizebin)
    BRPNASData<-matrix(NA,ncol=nsizebin,nrow=iYearNum)
    Recruitment<-1e5
    RecPrjRatio<-report$Recrut_Sizebin_ratio
    BRPS<-exp(-report$Natural_Mortality[Year,])
    
    for( i in 1:(iYearNum-1))
    {
      if(i==1)
        BRPNASData[1,]<- BRPNAS <- BRPNAS+Recruitment * RecPrjRatio;
      #BRPNASData<-rbind(BRPNASData,BRPNAS);
      NASTemp<-BRPNAS*BRPS;
      BRPNAS<-NASTemp%*%GM
      BRPNASData[i+1,]<-BRPNAS;
    } 
    iYear=seq(1,iYearNum,1);RRRData=NULL;RMainData=c();
    
    for( i in 1:iYearNum)
    {
      
      RRData<-BRPNASData[i,]
      RRData<-RRData/sum(RRData)
      RMainData[i]<-sum(RRData*sizebins)
      for(j in 1:nsizebin)
      {
        SizeNum<-RRData[j]*1000  
        if( SizeNum<1)
          next;
        XXX1<-rep(i,SizeNum)
        XXX2<-rep(j,SizeNum)
        RRRDataTemp<-cbind(iYear[XXX1],sizebins[XXX2])
        RRRData<-rbind(RRRData,RRRDataTemp)
      }
    }
  }
  boxplot(RRRData[,2]~factor(RRRData[,1]),xlab="Year",ylab="CL(mm)",
          cex.lab=1.0, cex.axis=1.0, cex.main=1.0, cex.sub=1.0)
  #stripchart(RRRData[,2]~factor(RRRData[,1]),method='jitter',add=T,vertical=T,pch=1)
  lines(1:iYearNum,RMainData,type="l",col="blue",lwd=3)
}
######################################################################

#######################################################################
#Selectivity of fisheries
#######################################################################
Plot_sel_f = function()
{
  par(mfrow=c(2,3))
  
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  PlotSelF(1987,1,1)
  box(col='black')
  axis(2,col='black',at=seq(0,0.8,0.2))
  
  PlotSelF(2010,1,2)
  box(col='black')
  
  PlotSelF(2010,1,3)
  box(col='black')
  
  PlotSelF(2014,1,2)
  box(col='black')
  axis(1,col='black',at=sizebins[-1])
  axis(2,col='black',at=seq(0,0.8,0.2))
  
  PlotSelF(2015,1,3)
  box(col='black')
  axis(1,col='black',at=sizebins[-1])

  mtext('Carapace Length (mm)',side=1,outer=T,line=2.2)
  mtext('Selectivity',side=2,outer=T,line=2.2)
  mtext(paste('Estimated fleet selectivity'),side=3,outer=T,line=2.2)
}

PlotSelF=function(year,season,fishery)
{
  Season=report$Time_step
  nsizebin=length(sizebins)
  Year=year-report$BeginYear+1
  
  if(Season>1)
  {
    itemp=Nyear*Season*(fishery-1)+Nyear*(season-1)
    Maintemp=paste(year,";Season",season,";Fishery",fishery,sep='')
  }else
  {
    itemp=Nyear*(fishery-1)
    if (year==1987)
      fishery=2
    if (year>2000)
      fishery=fishery+1
    Maintemp=paste("Selectivity time block",fishery)
  }
  sel=report$Fleet_selectivity[itemp+Year,]
  plot(sizebins,sel,xlab="Length (mm)",ylab="Selectivity",
       cex=1.2,col='black',type='o',lwd=2,pch=16,
       cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3,axes=FALSE)
  mtext(Maintemp,side=3,line=-2,adj=0.1,cex=0.6)
}
######################################################################

#######################################################################
#Selectivity of surveys
#######################################################################
PlotSelS=function(year,survey)
{
  
  par(mfrow=c(1,1))
  
  par(mar =c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  nsizebin=length(sizebins)
  Year=year-report$BeginYear+1
  
  itemp=Nyear*(survey-1)+Nyear
  Maintemp=paste("Survey",survey)
  
  sel=report$Index_selectivity[itemp,]
  plot(sizebins,sel,xlab="Length (mm)",ylab="Selectivity",
       main=Maintemp,cex=1.2,col='black',type='o',lwd=2,pch=16,
       cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
}
######################################################################

#######################################################################
#Fishing mortality
#######################################################################
Plot_f=function()
{
  par(mfrow=c(3,2))
  
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  PlotF(1,1)
  
  PlotF(2,1)
  
  PlotF(3,1)
  
  mtext('Year',side=1,outer=T,line=2.2)
  mtext('Fishing mortality',side=2,outer=T,line=2.2)
  mtext(paste('Estimated fishing mortality'),side=3,outer=T,line=2.2)
}

PlotF=function(fishery,season)
{
  Year=seq(report$BeginYear,report$EndYear,1)
  Season=report$Time_step
  f=report$Fishing_Mortality
  if (Season>1)
  {
    #par(mfrow=c(2,2))
    for (i in 1:Season)
    {
      Maintemp=paste('Fishery',fishery,";Season",i)
      itemp=Season*(fishery-1)
      ftemp=f[itemp+i,]
      plot(Year,ftemp,xlab="Year",ylab="Fishing mortality",type='b',pch=16,
           cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3,ylim=c(0,1.2*max(ftemp)),axes=FALSE)
      box(col='black')
      mtext(Maintemp,side=3,line=-2,adj=0.1,cex=1.0)
      if (i==1)
        axis(2,col='black',at=seq(0,10,1))
      if (fishery==3)
        axis(1,col='black',at=Year)
    }
  }else
  {
    Maintemp=paste('Fishery',fishery)
    plot(Year,f[fishery,],xlab="Year",ylab="F",type='b',
         pch=16,cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
  }
  
}
######################################################################

#######################################################################
#Natural mortality
#######################################################################
PlotM=function(year)
{
  par(mfrow=c(1,1))
  par(mar = c(4, 4, 1, 0.5))
  
  Season=report$Time_step
  M=report$Natural_Mortality
  Year=year-report$BeginYear+1
  par(mfrow=c(1,1))
  if (Season>1)
  {
    itemp=Season*(Year-1)
    m=M[(itemp+1),]
    Maintemp=paste("Seasonal M; Year",year)
    plot(sizebins,m,xlab="Length (mm)",ylab="M",main=Maintemp,type='l',lwd=3,
         col='blue',cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
  }else
  {
    Maintemp=paste("Year",year)
    m=M[Year,]
    plot(sizebins,m,xlab="Length (mm)",ylab="M",main=Maintemp,type='l',lwd=3,
         col='blue',cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
  }  
}
######################################################################

#######################################################################
#Recruitment
#######################################################################
PlotR=function()
{
  Year=seq(report$BeginYear,report$EndYear,1)
  R=report$recruitment_Pred
  Rbar=report$Mean_Recruitment
  Rdev=report$recruitment_log_Dev
  RdevSD=report$recruitment_log_DevSD2
  if(report$SexAtSizeLamda>0)
    {SSB=report$Spawning_stock_Biomass
  }else
    SSB=report$Spawning_stock_Biomass_input
  par(mfrow=c(3,1))
  par(mar = c(4, 4.5, 2, 0.5))
  
  plot(Year,R,ylab="Recruits (million)",xlab="Year",type='o',
       col='black',cex=1.5,pch=20,lwd=1.5,lty='dashed',cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  plot(Year,Rdev,xlab="Year",ylab="Log recruitment deviation",type='b',
       col='black',cex=1.5,pch=20,lty='dotted',cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h=0,col="gray",lwd=1.5)
  plot(SSB,R,xlab="Spawning biomass (mt)",ylab="Recruits (million)",
       col='black',cex=1.5,pch=20,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  par(mfrow=c(1,1))
}
######################################################################

#######################################################################
#SSB
#######################################################################
PlotSSB=function()
{
  par(mfrow=c(1,1))
  Year=seq(report$BeginYear,report$EndYear,1)
  if(report$SexAtSizeLamda>0)
    {SSB=report$Spawning_stock_Biomass
  }else
    SSB=report$Spawning_stock_Biomass_input
  
  #par(mfrow=c(2,1))
  par(mar = c(4, 4.5, 2, 0.5))
  plot(Year,SSB,xlab="Year",ylab="Spawning biomass (mt)",type='o',
       col='red',cex=2,pch=20,lty='dashed',cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
}
######################################################################

#######################################################################
#Abundance at the beginning of the year
#######################################################################
PlotAbun=function()
{
  library(lattice)
  Abun=report$Abundance_at_Size # ;Abun[,1:2]=0
  Season=report$Time_step
  Year=seq(report$BeginYear,report$EndYear,1)
  
  AbunB=matrix(NA,ncol=3,nrow=length(Year)*length(sizebins))
  itemp=1
  
  for (i in Year)
  {
    iL=i-report$BeginYear+1
    
    for (j in 1:length(sizebins))
    {
      if (Season>1)
      {
        AbunB[itemp,]=c(i,sizebins[j],Abun[(iL-1)*Season+1,j])
        itemp=itemp+1
      }else
      {
        AbunB[itemp,]=c(i,sizebins[j],Abun[iL,j])
        itemp=itemp+1
      }
      
    }
  }
  AbunB<-as.data.frame(AbunB)    
  names(AbunB)=c("Year","CL","Abun")
  
  plot<-xyplot(CL ~ Year,data=AbunB,layout=c(1,1),par.strip.text=list(cex=1.3),
               as.table=T,subscript=T,z=AbunB[,3],
               main=paste("Beginning of year expected numbers at length"),
               scales=list(y=list(cex=1.3),x=list(cex=1.3)),panel=panel.bubble.for.R,scalem=4,
               ylab=list(label="Length (mm)",cex=1.3),xlab=list(label="Year", cex=1.3))
  print(plot)  
#  par(mfrow=c(1,1)
}
######################################################################

#######################################################################
#Bubble function
#######################################################################
panel.bubble.for.R<-function(x,y,z,subscripts,scalem=4,...){
  #panel.grid(col="lightgrey",lwd=1,h=-1,v=-1)
  cex<-abs(z) # size of points depends on the absolute value of z and is consistent across panels
  cex<-scalem*sqrt(cex/max(cex))
  flagplus<- z>=0
  panel.xyplot(x[flagplus],y[flagplus],pch=16,col="gray",cex=cex[flagplus])
  panel.xyplot(x[!flagplus],y[!flagplus],pch=1,col="gray",cex=cex[!flagplus])
}

######################################################################

#######################################################################
#Survey Length Comp
#######################################################################
PlotSLC<-function()
{
  N_syear=Nyear
  DimIndexIndex=report$DimIndexIndex
  if(is.vector(DimIndexIndex))
  {
    N_index=1
    DimIndexIndex=matrix(DimIndexIndex,nrow = 1,ncol = length(DimIndexIndex))
  }else{
    N_index=nrow(DimIndexIndex)
  }
  Year=seq(report$BeginYear,report$EndYear,1)
  s_length=c();ind=c()
  
  ind_year = matrix(NA,nrow=N_index,ncol=2)
  
  ind[1] = 1
  
  for (s in 1:N_index)
  {
    s_length[s] = sum(!is.na(DimIndexIndex[s,]))
    #SurveyESS_Obs[s,1:s_length[s]] = report$Survey_ESSinput[s,][1:s_length[s]]
    ind1=c(DimIndexIndex[s,][1:s_length[s]])
    ind[s+1] = ind[s]+s_length[s]
    #Scomp_Obs[c(ind1),,s] = report$Survey_Comp_Obs[c(),]
    ind_year[s,] = c(Year[c(ind1)][1],Year[c(ind1)][s_length[s]])
    
    plot_nrow = ceiling(s_length[s]/6); plot_ncol=6
    yaxis = rep(1,plot_nrow)
    
    for(pl in 1:plot_nrow)
    {
      yaxis[pl] = yaxis[1]+(pl-1)*plot_ncol
    }
    
    if (plot_nrow >1)
    {
      xaxis = seq(s_length[s]-plot_ncol+1,s_length[s],1)
    }else
    {
      xaxis = seq(1,s_length[s],1)
    }
    
    par(mfrow=c(plot_nrow,plot_ncol))
    
    par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    
    
    down_pred = (s-1)*N_syear+ind_year[s,1]-report$BeginYear+1
    up_pred = down_pred+s_length[s]-1
    
    Scomp_Obs = report$Survey_Comp_Obs[ind[s]:(ind[s+1]-1),]
    Scomp_Pred = report$Survey_Comp_Pred[down_pred:up_pred,]
    
    year_plot = seq(ind_year[s,1],ind_year[s,2],1)
    
    for (yi in 1:s_length[s])
    {
      plot(sizebins,Scomp_Obs[yi,],type='l',ylim=c(0,0.25),axes=FALSE,lwd=1,col='gray')
      x=c(sizebins,rev(sizebins))
      y=c(rep(0,length(sizebins)),rev(Scomp_Obs[yi,]))
      polygon(x,y,col="gray",border = NA)
      mtext(paste(year_plot[yi]),side=3,line=-1,adj=0.1,cex=0.6)
      lines(sizebins,Scomp_Pred[yi,],col="black",lwd=2)
      #text(x=27,y=0.15,paste("effN=",SurveyESS_Pred[i]),cex=1)
      if(yi %in% yaxis)
        axis(2,col='gray',at=seq(0,0.20,0.02))
      if(yi %in% xaxis)
        axis(1,col='gray',at=sizebins[-1])
      box(col='gray')
    }
    mtext('Length (mm)',side=1,outer=T,line=2.2)
    mtext('Proportion',side=2,outer=T,line=2.2)
    mtext(paste('Survey-',s,';black line = Pred',sep=''),side=3,outer=T,line=2.2)
  }
}
######################################################################
#######################################################################
#Survey Length Comp aggregated
#######################################################################
PlotSLCA=function()
{
  survey=nrow(report$Survey_Comp_Pred)/report$N_year
  par(mfrow=c(survey,1))
  #par(mar = c(4, 4, 1, 0.5))
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  
  Year=seq(1984,2013,1)
  
  
  for (i in 1:survey)
  {
    
    
    if(i==1)
    {
      #Year=seq(1984,2013,1)
      SurveyESS_Obs=rep(0,30)
      SurveyESS_Pred=rep(0,30)
      Scomp_Obs=matrix(0,nrow=30,ncol=25)
      Scomp_Pred=matrix(0,nrow=30,ncol=25)
      Scomp_Obs=report$Survey_Comp_Obs
      Scomp_Pred=report$Survey_Comp_Pred
      SurveyESS_Obs=report$Survey_ESSinput
      SurveyESS_Pred=round(report$Survey_ESSpred)
      #N_syear=length(Year)
      iL=0
      for(j in 1:nrow(Scomp_Obs))
      {
        if(sum(Scomp_Obs[j,]>0))
          iL=iL+1
      }
      
    }
    
    if(i==2)
    {
      #Year=seq(1984,2013,1)
      SurveyESS_Obs=rep(0,30)
      SurveyESS_Pred=rep(0,30)
      Scomp_Obs=matrix(0,nrow=30,ncol=25)
      Scomp_Pred=matrix(0,nrow=30,ncol=25)
      Scomp_Obs=report$Survey_Comp_Obs[26:55,]
      Scomp_Pred=report$Survey_Comp_Pred[31:60,]
      SurveyESS_Obs=report$Survey_ESSinput[2,][1:30]
      SurveyESS_Pred=round(report$Survey_ESSpred[2,])[1:30]
      #N_syear=length(Year)
      iL=0
      for(j in 1:nrow(Scomp_Obs))
      {
        if(sum(Scomp_Obs[j,]>0))
          iL=iL+1
      }
    }
    
    if(i==3)
    {
      #Year=seq(2005,2013,1)
      SurveyESS_Obs=rep(0,30)
      SurveyESS_Pred=rep(0,30)
      Scomp_Obs=matrix(0,nrow=30,ncol=25)
      Scomp_Pred=matrix(0,nrow=30,ncol=25)
      
	  Scomp_Obs[26:29,]=report$Survey_Comp_Obs[56:59,]
      Scomp_Pred[26:29,]=report$Survey_Comp_Pred[86:89,]
      SurveyESS_Obs[26:29]=report$Survey_ESSinput[3,][1:4]
      SurveyESS_Pred[26:29]=round(report$Survey_ESSpred[3,])[26:29]
	  
      #N_syear=length(Year)
      iL=0
      for(j in 1:nrow(Scomp_Obs))
      {
        if(sum(Scomp_Obs[j,]>0))
          iL=iL+1
      }
    }
    
    if(i==4)
    {
      #Year=seq(2005,2013,1)
      SurveyESS_Obs=rep(0,30)
      SurveyESS_Pred=rep(0,30)
      Scomp_Obs=matrix(0,nrow=30,ncol=25)
      Scomp_Pred=matrix(0,nrow=30,ncol=25)
      Scomp_Obs[26:29,]=report$Survey_Comp_Obs[67:70,]
      Scomp_Pred[26:29,]=report$Survey_Comp_Pred[116:119,]
      SurveyESS_Obs[26:29]=report$Survey_ESSinput[4,][1:4]
      SurveyESS_Pred[26:29]=round(report$Survey_ESSpred[4,])[26:29]
      #N_syear=length(Year)
      iL=0
      for(j in 1:nrow(Scomp_Obs))
      {
        if(sum(Scomp_Obs[j,]>0))
          iL=iL+1
      }
    }
    
    
    #Scomp_Obs_temp=report$Survey_Comp_Obs[(Nyear*(i-1)+1):(Nyear*(i-1)+Nyear),]
    Scomp_Obs_A=apply(Scomp_Obs,2,sum)/iL
    #Scomp_Pred_temp=report$Survey_Comp_Pred[(Nyear*(i-1)+1):(Nyear*(i-1)+Nyear),]
    Scomp_Pred_A=apply(Scomp_Pred,2,sum)/iL
    
    SurveyESS_Obs=sum(SurveyESS_Obs)
    SurveyESS_Pred=sum(SurveyESS_Pred)
    
    plot(sizebins,Scomp_Obs_A,type='l',ylim=c(0,0.16),axes=FALSE,lwd=1,col="gray")
    x=c(sizebins,rev(sizebins))
    y=c(rep(0,length(sizebins)),rev(Scomp_Obs_A))
    polygon(x,y,col="gray",border=NA)
    legend('topleft',paste("Survey",i),cex=1.2,bty='n')
    rp=vector('expression',2)
    rp[1]=substitute(expression(N==Value),list(Value=SurveyESS_Obs))[2]
    rp[2]=substitute(expression(effN==Value),list(Value=SurveyESS_Pred))[2]
    legend('topright',legend=rp,bty='n',y.intersp=0.8)
    lines(sizebins,Scomp_Pred_A,col='black',lwd=2)
    axis(2,col='gray',at=seq(0,0.15,0.02))
    if(i==survey)
      axis(1,col='gray',at=sizebins)
    box(col='gray')
    
  }
  mtext('Length (mm)',side=1,outer=T,line=2.2)
  mtext('Proportion',side=2,outer=T,line=2.2)
  mtext(paste('Length Comp,','aggregated across time by survey',',black line=Pred'),
        side=3,outer=T,line=2.2)
}

######################################################################

#######################################################################
#Survey index
#######################################################################
PlotSI=function(n)
{
  
  if(n>2){
    par(mfrow=c(2,2))
  }else{
    par(mfrow=c(1,2))
  }
  par(mar = c(4, 4, 1, 0.5))
#   survey=1
#   
#     Year=seq(1984,2008,1)
#     SI_Obs=log(report$Survey_Index_Obs[1,][1:25])
#     SI_Pred=log(report$Survey_Index_Pred[1,][1:25])
#     SIESS_input=report$Survey_ESSinput[1,][8:25]
#     SIESS_pred=report$Survey_ESSpred[1,][8:25]
# 	scaler=1.3
#   
#   
#   plot(Year,SI_Pred,ylab='Survey Index',type='b',pch=16,
#        col='black',cex=1.5,ylim=c(0,scaler*max(SI_Obs)))
# 	   legend('topright',paste('Survey',survey),bty='n')
#   points(Year[1:length(SI_Obs)],SI_Obs,pch=16,col='red',cex=1.3)
#   legend('topleft',c('Obs','Pred'),pch=c(16,16),col=c('red','black'),bty='n',y.intersp=0.8,xjust=-1)
#  
  
  survey=1
  
    Year=seq(1984,report$EndYear,1)
    SI_Obs=report$Survey_Index_Obs[1,]
    SI_Pred=report$Survey_Index_Pred[1,]
    #SIESS_input=report$Survey_ESSinput[1,]
    #SIESS_pred=report$Survey_ESSpred[1,]
	scaler=1.1
	ylim=max(max(SI_Obs),max(SI_Pred))
  plot(Year,SI_Pred,ylab=paste('Index',sep=''),type='l',lwd=2,
       col='black',cex=1.1,ylim=c(0,scaler*ylim),main='Summer survey')
  #legend(2014,15,paste('Survey',survey),bty='n')
  points(Year[1:length(SI_Obs)],SI_Obs,pch=16,col='red',cex=1.0)
  #legend('topleft',c('Obs','Pred'),pch=c(16,16),col=c('red','black'),bty='n',y.intersp=0.8,xjust=-1)
 
   survey=2
   
   Year=seq(1986,2016,1)
   SI_Obs=report$Survey_Index_Obs[2,][1:31]
   SI_Pred=report$Survey_Index_Pred[2,][1:31]
   ylim=max(max(SI_Obs),max(SI_Pred))
   plot(Year,SI_Pred,ylab=paste('Index',sep=''),type='l',lwd=2,
        col='black',cex=1.1,ylim=c(0,scaler*ylim),main='Fall survey')
   #legend('topright',paste('Survey',survey),bty='n')
   points(Year[1:length(SI_Obs)],SI_Obs,pch=16,col='red',cex=1.0)
   
   if (n>2){
     # survey=3
     # 
     # Year=seq(2009,2016,1)
     # SI_Obs=report$Survey_Index_Obs[3,][1:8]
     # SI_Pred=report$Survey_Index_Pred[3,][26:33]
     # 
     # plot(Year,SI_Pred,ylab=paste('Fall Index',survey,sep=''),type='l',lwd=2,
     #      col='black',cex=1.5,ylim=c(0,scaler*max(SI_Obs)),main='')
     # #legend('topright',paste('Survey',survey),bty='n')
     # points(Year[1:length(SI_Obs)],SI_Obs,pch=16,col='red',cex=1.0)
     # 
     survey=3
     
     Year=seq(2003,2016,1)
     SI_Obs=report$Survey_Index_Obs[3,][1:14]
     SI_Pred=report$Survey_Index_Pred[3,][20:33]
     ylim=max(max(SI_Obs),max(SI_Pred))
     plot(Year,SI_Pred,ylab=paste('Index',sep=''),type='l',lwd=2,
          col='black',cex=1.5,ylim=c(0,scaler*ylim),main='Spring survey')
     #legend('topright',paste('Survey',survey),bty='n')
     points(Year[1:length(SI_Obs)],SI_Obs,pch=16,col='red',cex=1.0)
   }
   
   mtext('black line ~ Pred',side=1,outer=T,line=2.2)
}

PlotSIESS=function()
{
  par(mfrow=c(3,1))
  par(mar = c(4, 4, 1, 0.5))
  survey=1
  
    Year=seq(1984,2008,1)
    SI_Obs=log(report$Survey_Index_Obs[1,][1:25])
    SI_Pred=log(report$Survey_Index_Pred[1,][1:25])
    SIESS_input=report$Survey_ESSinput[1,][8:25]
    SIESS_pred=report$Survey_ESSpred[1,][8:25]
	scaler=1.3
  
  plot(SIESS_input[1:length(SI_Obs)],SIESS_pred[1:length(SI_Obs)],type='p',pch=16,
       xlab='Input effective sample size',ylab='Estimated effective sample size')
    legend('topleft',paste('Survey',survey),bty='n')  
  #lines(SIESS_input,SIESS_input)
  
  survey=2
  
    Year=seq(1984,2013,1)
    SI_Obs=log(report$Survey_Index_Obs[2,])
    SI_Pred=log(report$Survey_Index_Pred[2,])
    SIESS_input=report$Survey_ESSinput[2,]
    SIESS_pred=report$Survey_ESSpred[2,]
	scaler=1.3
  
  plot(SIESS_input[1:length(SI_Obs)],SIESS_pred[1:length(SI_Obs)],type='p',pch=16,
       xlab='Input effective sample size',ylab='Estimated effective sample size')
  legend('topleft',paste('Survey',survey),bty='n') 
  #lines(SIESS_input,SIESS_input)
  
  survey=3
  
    Year=seq(2009,2012,1)
    SI_Obs=log(report$Survey_Index_Obs[3,][1:4])
    SI_Pred=log(report$Survey_Index_Pred[3,][26:29])
    SIESS_input=report$Survey_ESSinput[3,][1:4]
    SIESS_pred=report$Survey_ESSpred[3,][26:29]
	scaler=2
  
  plot(SIESS_input[1:length(SI_Obs)],SIESS_pred[1:length(SI_Obs)],type='p',pch=16,
       xlab='Input effective sample size',ylab='Estimated effective sample size')
  legend('topleft',paste('Survey',survey),bty='n')   
  #lines(SIESS_input,SIESS_input)
  
  #if(survey==4)
  
    #Year=seq(2009,2012,1)
    #SI_Obs=report$Survey_Index_Obs[4,][1:4]
    #SI_Pred=report$Survey_Index_Pred[4,][26:29]
    #SIESS_input=report$Survey_ESSinput[4,][1:4]
    #SIESS_pred=report$Survey_ESSpred[4,][26:29]
  
  
  #Year=seq(report$BeginYear,report$EndYear,1)
  #SI_Obs=report$Survey_Index_Obs#[survey,]
  #SI_Pred=report$Survey_Index_Pred[1:length(SI_Obs)]#[survey,]
  #SIESS_input=report$Survey_ESSinput#[survey,]
  #SIESS_pred=report$Survey_ESSpred[1:length(SI_Obs)]#[survey,]
  
  
  #plot(Year,SI_Pred,ylab='Survey Index',type='b',pch=16,
       #col='black',cex=1.5,ylim=c(0,scaler*max(SI_Obs)),
       #main=paste('Survey',survey))
  #points(Year[1:length(SI_Obs)],SI_Obs,pch=16,col='red',cex=1.3)
  #legend('topleft',c('Obs','Pred'),pch=c(16,16),col=c('red','black'),bty='n',y.intersp=0.8,xjust=-1)
  #plot(SIESS_input[1:length(SI_Obs)],SIESS_pred[1:length(SI_Obs)],type='p',pch=16,
       #xlab='Observed sample size',ylab='Effective sample size')
  #lines(SIESS_input,SIESS_input)
}


######################################################################

#######################################################################
#Total Catch
#######################################################################
PlotTC=function()
{
  par(mfrow=c(1,1))
  par(mar = c(4, 4, 1.5, 1))
  #par(tcl = -0.25)
  #par(mgp = c(2, 0.6, 0))

  bioORnum=report$BiomassORNum
  if (sum(bioORnum)==0)
  {
    ylab=paste('Catch (million)')
  }else
  {
    ylab='Catch (mt)'
  }
  
  total_catch_obs=apply(report$Catch_Obs,2,sum)
  total_catch_pred=apply(report$Catch_Pred,2,sum)
  ylim=max(max(total_catch_obs),max(total_catch_pred))*1.1
  Year=seq(report$BeginYear,report$EndYear,1)
  plot(Year,total_catch_obs,type='l',lwd=2,col='black',xlab='Year',ylab=ylab,main='Total catch',ylim=c(0,ylim))
  points(Year,total_catch_pred,pch=20,col='red',cex=1.5)
  legend('topleft',c('Pred','Obs'),lty=c(NA,1),pch=c(20,NA),lwd=c(NA,3),cex=1.0,
         col=c('red','black'),bty='n',y.intersp=0.8,xjust=1)
}
######################################################################

#######################################################################
#Fishery-specific Catch
#######################################################################
Plot_f_c=function()
{
  par(mfrow=c(3,2))
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  PlotC(1) # plot fishery-specific catch
  PlotC(2)
  PlotC(3)
}

PlotC=function(fishery)
{
  Season=report$Time_step
  Year=seq(report$BeginYear,report$EndYear,1)
  bioORnum=report$BiomassORNum
  if (sum(bioORnum)==0)
  {
    ylab=paste('Catch (million)')
  }else
  {
    ylab='Catch (mt)'
  }
  if (Season>1)
  {
    #par(mfrow=c(1,2))
    #par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
    #par(tcl = -0.25)
    #par(mgp = c(2, 0.6, 0))
    ylim=ceiling(max(max(report$Catch_Obs),max(report$Catch_Pred)))
    #Season=1
    for (i in 1:Season)
    {
      total_catch_obs=report$Catch_Obs[Season*(fishery-1)+i,]
      total_catch_pred=report$Catch_Pred[Season*(fishery-1)+i,]
      #ylim=ceiling(max(max(total_catch_obs),max(total_catch_pred)))
      plot(Year,total_catch_obs,type='l',lwd=3,col='black',axes=FALSE,ylim=c(0,ylim))
      points(Year,total_catch_pred,pch=20,col='red',cex=1.5)
      legend('topright',paste('Fishery',fishery,'Season',i),bty='n',cex=1.2)
      if (i==4)
        legend('topleft',c('Pred','Obs'),lty=c(NA,1),pch=c(20,NA),lwd=c(NA,2),cex=1,
               col=c('red','black'),bty='n',y.intersp=0.8,xjust=-1)
      if (i %in% c(1))
        axis(2,col='black',at=seq(0,ylim,50))
      if(fishery %in% c(3))
        axis(1,col='black',at=Year)
      box(col='black')
    }
    mtext('Year',side=1,outer=T,line=2.2)
    mtext(paste(ylab),side=2,outer=T,line=2.2)
    #mtext(paste('Seasonal catch,','fishery',fishery),side=3,outer=T,line=2.2,cex=1.5)
    
  }else
  {
    
    par(mfrow=c(fishery,1))
    par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    ylim=ceiling(max(max(report$Catch_Obs),max(report$Catch_Pred)))
    for (i in 1:fishery)
    {
      total_catch_obs=report$Catch_Obs[i,]
      total_catch_pred=report$Catch_Pred[i,]
      plot(Year,total_catch_obs,type='l',lwd=3,col='black',axes=FALSE,ylim=c(0,ylim))
      points(Year,total_catch_pred,pch=20,col='red',cex=1.5)
      legend('top',paste('Fishery',i),bty='n',cex=1.5)
      if (i==1)
        legend('topleft',c('Pred','Obs'),lty=c(NA,1),pch=c(20,NA),lwd=c(NA,3),cex=1,
               col=c('red','black'),bty='n',y.intersp=0.8,xjust=-1)
      if (i %in% seq(1,fishery,1))
        axis(2,col='gray',at=seq(0,ylim,50))
      if(i==fishery)
        axis(1,col='gray',at=Year)
      box(col='gray')
    }
    mtext('Year',side=1,outer=T,line=2.2)
    mtext(paste(ylab),side=2,outer=T,line=2.2)
    mtext(paste('Total catch by fleet'),side=3,outer=T,line=2.2,cex=1.5)
  }    
}
######################################################################

#######################################################################
#Catch comp 
#######################################################################
PlotCC=function(fishery,season)
{
  Year=seq(report$BeginYear,report$EndYear,1)
  Nyear=length(Year)
  n_year = Nyear
  Season=report$Time_step
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  if (Season>1)
  {
    temp=Nyear*Season*(fishery-1)+Nyear*(season-1)
    temp1=Season*(fishery-1)+season
    Ccomp_Obs=report$Catch_Comp_Obs[(temp+1):(temp+Nyear),]
    Ccomp_Pred=report$Catch_Comp_Pred[(temp+1):(temp+Nyear),]
    CatchESS_Obs=report$Catch_ESSinput[temp1,]
    CatchESS_Pred=report$Catch_ESSpred[temp1,]
    
    if (Nyear<25)
    {
      par(mfrow=c(5,5))
      temp=(Nyear%/%5*(5-1)+Nyear%%5):(Nyear%/%5*(5-1)+5)
      yaxis=c(1,6,11,16,21)
      xaxis=c(21,22,23,24,25,temp)
    }
    if (Nyear>25)
    {
      par(mfrow=c(5,6))
      yaxis=c(1,7,13,19,25,31)
      xaxis=c(25,26,27,28,29,30)
    }
    
    if (Nyear>25)
    {
      par(mfrow=c(6,6))
      yaxis=c(1,7,13,19,25,31)
      xaxis=c(28,29,30,31,32,33)
    }
    
    # if (fishery==1&season==1)
    # {
    #   Ccomp_Obs = Ccomp_Obs[1:16,]
    #   Ccomp_Pred = Ccomp_Pred[1:16,]
    #   n_year = 16
    #   Year = seq(1984,2000,1)
    #   xaxis=c(11,12,13,14,15,16)
    # }
    # 
    # if (fishery==2&season==1)
    # {
    #   Ccomp_Obs = Ccomp_Obs[17:Nyear,]
    #   Ccomp_Pred = Ccomp_Pred[17:Nyear,]
    #   n_year = 17
    #   Year = seq(2000,report$EndYear,1)
    #   xaxis=c(12,13,14,15,16,17)
    # }
    # 
    # if (fishery==3&season==1)
    # {
    #   Ccomp_Obs = Ccomp_Obs[17:Nyear,]
    #   Ccomp_Pred = Ccomp_Pred[17:Nyear,]
    #   n_year = 17
    #   Year = seq(2000,report$EndYear,1)
    #   xaxis=c(12,13,14,15,16,17)
    # }
    # 
    par(mfrow=c(3,6))
    yaxis=c(1,7,13)
    
    
    for (i in 1:n_year)
    {
      plot(sizebins,Ccomp_Obs[i,],type='l',ylim=c(0,0.35),axes=FALSE,lwd=1,col='gray')
      x=c(sizebins,rev(sizebins))
      y=c(rep(0,length(sizebins)),rev(Ccomp_Obs[i,]))
      polygon(x,y,col="gray",border=NA)
      mtext(paste(Year[i]),side=3,line=-1,adj=0.1,cex=0.6)
      #legend('topleft',paste(Year[i]),cex=1.1,bty='n')
      #rp=vector('expression',2)
      #rp[1]=substitute(expression(N==Value),list(Value=CatchESS_Obs[i]))[2]
      #rp[2]=substitute(expression(effN==Value),list(Value=CatchESS_Pred[i]))[2]
      #legend('topright',legend=rp,bty='n',y.intersp=0.8,yjust=1)
      lines(sizebins,Ccomp_Pred[i,],col='black',lwd=2)
      #text(x=27,y=0.15,paste("effN=",SurveyESS_Pred[i]),cex=1)
      if(i %in% yaxis)
        axis(2,col='gray',at=seq(0,0.2,0.05))
      if(i %in% xaxis)
        axis(1,col='gray',at=sizebins)
      box(col='gray')
    }
    mtext('Length (mm)',side=1,outer=T,line=2.2)
    mtext('Proportion',side=2,outer=T,line=2.2)
    mtext(paste('Length Comp,','fishery',fishery,'season',season,',black line=Pred'),
          side=3,outer=T,line=2.2)
  }else
  {
    temp=Nyear*(fishery-1)
    
    Ccomp_Obs=report$Catch_Comp_Obs[(temp+1):(temp+Nyear),]
    Ccomp_Pred=report$Catch_Comp_Pred[(temp+1):(temp+Nyear),]
    CatchESS_Obs=report$Catch_ESSinput[fishery,]
    CatchESS_Pred=report$Catch_ESSpred[fishery,]
    
    # if (Nyear<25)
    # {
    #   par(mfrow=c(5,5))
    #   temp=(Nyear%/%5*(5-1)+Nyear%%5):(Nyear%/%5*(5-1)+5)
    #   yaxis=c(1,6,11,16,21)
    #   xaxis=c(21,22,23,24,25,temp)
    # }
    # if (Nyear>25)
    # {
    #   par(mfrow=c(5,6))
    #   yaxis=c(1,7,13,19,25,31)
    #   xaxis=c(25,26,27,28,29,30)
    # }
    for (i in 1:Nyear)
    {
      plot(sizebins,Ccomp_Pred[i,],type='l',ylim=c(0,0.3),axes=FALSE)
      x=c(sizebins,rev(sizebins))
      y=c(rep(0,length(sizebins)),rev(Ccomp_Pred[i,]))
      polygon(x,y,col="gray")
      mtext(paste(Year[i]),side=3,line=-1,adj=0.1,cex=0.6)
      #legend('topleft',paste(Year[i]),cex=1.1,bty='n')
      rp=vector('expression',2)
      rp[1]=substitute(expression(N==Value),list(Value=CatchESS_Obs[i]))[2]
      rp[2]=substitute(expression(effN==Value),list(Value=CatchESS_Pred[i]))[2]
      legend('topright',legend=rp,bty='n',y.intersp=0.8,yjust=1)
      lines(sizebins,Ccomp_Obs[i,],col='red')
      #text(x=27,y=0.15,paste("effN=",SurveyESS_Pred[i]),cex=1)
      if(i %in% yaxis)
        axis(2,col='gray',at=seq(0,0.2,0.05))
      if(i %in% xaxis)
        axis(1,col='gray',at=sizebins)
      box(col='gray')
    }
    mtext('Length (mm)',side=1,outer=T,line=2.2)
    mtext('Proportion',side=2,outer=T,line=2.2)
    mtext(paste('Length Comp,','fishery',fishery,',red line=Obs'),
          side=3,outer=T,line=2.2)
  }
}
######################################################################

#######################################################################
#Catch comp aggregated
#######################################################################
PlotCCA=function(Nfishery)
{
  par(mfrow=c(Nfishery,1))
  #par(mar = c(4, 4, 1, 0.5))
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  Season=report$Time_step
  for (i in 1:Nfishery)
  {
    if (Season>1)
    {
      Ccomp_Obs_temp=report$Catch_Comp_Obs[(Nyear*Season*(i-1)+1):(Nyear*Season*(i-1)+Nyear*Season),]
      Ccomp_Obs_A=apply(Ccomp_Obs_temp,2,sum)/(report$N_year*Season)
      
      Ccomp_Pred_temp=report$Catch_Comp_Pred[(Nyear*Season*(i-1)+1):(Nyear*Season*(i-1)+Nyear*Season),]
      Ccomp_Pred_A=apply(Ccomp_Pred_temp,2,sum)/(report$N_year*Season)
      
      CatchESS_Obs=sum(report$Catch_ESSinput[(Season*(i-1)+1):(Season*(i-1)+1),])
      CatchESS_Pred=round(sum(report$Catch_ESSpred[(Season*(i-1)+1):(Season*(i-1)+1),]),1)
    }else
    {
      Ccomp_Obs_temp=report$Catch_Comp_Obs[(Nyear*(i-1)+1):(Nyear*(i-1)+Nyear),]
      Ccomp_Obs_A=apply(Ccomp_Obs_temp,2,sum)/report$N_year
      
      Ccomp_Pred_temp=report$Catch_Comp_Pred[(Nyear*(i-1)+1):(Nyear*(i-1)+Nyear),]
      Ccomp_Pred_A=apply(Ccomp_Pred_temp,2,sum)/report$N_year
      
      CatchESS_Obs=sum(report$Catch_ESSinput[i,])
      CatchESS_Pred=round(sum(report$Catch_ESSpred[i,]),1)
    }
    
    plot(sizebins,Ccomp_Obs_A,type='l',ylim=c(0,max(Ccomp_Obs_A)),axes=FALSE,lwd=1,col="gray")
    x=c(sizebins,rev(sizebins))
    y=c(rep(0,length(sizebins)),rev(Ccomp_Obs_A))
    polygon(x,y,col="gray",border=NA)
    legend('topleft',paste("Fishery",i),cex=1.2,bty='n')
    rp=vector('expression',2)
    rp[1]=substitute(expression(N==Value),list(Value=CatchESS_Obs))[2]
    rp[2]=substitute(expression(effN==Value),list(Value=CatchESS_Pred))[2]
    legend('topright',legend=rp,bty='n',y.intersp=0.8)
    lines(sizebins,Ccomp_Pred_A,col='black',lwd=2)
    axis(2,col='gray',at=seq(0,0.15,0.02))
    if(i==Nfishery)
      axis(1,col='gray',at=sizebins)
    box(col='gray')
    
  }
  mtext('Length (mm)',side=1,outer=T,line=2.2)
  mtext('Proportion',side=2,outer=T,line=2.2)
  mtext(paste('Catch Length Comp,','aggregated cross time by fleet',',black line=Pred'),
        side=3,outer=T,line=2.2)
}

######################################################################
#Sex COMP
######################################################################
PlotSexComp=function(season,ylim_scaler)
{
  Year=seq(report$BeginYear,report$EndYear,1)
  Season=season
  
  if(Nyear<25)
  {
  par(mfrow=c(5,5))
  temp=(Nyear%/%5*(5-1)+Nyear%%5):(Nyear%/%5*(5-1)+5)
  yaxis=c(1,6,11,16,21)
  xaxis=c(21,22,23,24,25,temp)
  }
  
  if (Nyear>25)
  {
    par(mfrow=c(6,6))
    yaxis=c(1,7,13,19,25,31)
    xaxis=c(28,29,30,31,32,33)
  }
  
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  Abun=report$Abundance_at_Size[(Nyear*(Season-1)+1):(Nyear*Season),]
  if(season==1){Abun[,1:2]=0}
  Non_female=(1-report$FemaleProp_at_Size)*Abun
  Female=report$FemaleProp_at_Size*Abun
  ylim=ylim_scaler*max(Abun)
  
  for (i in 1:Nyear)
  {
    plot(sizebins,Non_female[i,],col='lightblue',type='l',ylim=c(0,ylim),axes=FALSE)
    xx=c(sizebins,rev(sizebins))
    yy=c(rep(0,length(sizebins)),rev(Non_female[i,]))
    polygon(xx, yy, col='lightblue',border=NA)  
    mtext(paste(Year[i]),side=3,line=-1,adj=0.1,cex=0.6)    
    yy2 <- c(Non_female[i,], rev(Female[i,]) + rev(Non_female[i,])) 
    polygon(xx, yy2, col='lightpink',border=NA) 
    if(i %in% yaxis)
      axis(2,col='gray',at=seq(0,ylim,50))
    if(i %in% xaxis)
      axis(1,col='gray',at=sizebins)
    box(col='gray')
  }
  mtext('Length (mm)',side=1,outer=T,line=2.2)
  mtext('Abundance (million)',side=2,outer=T,line=2.2)
  mtext(paste('Numbers at stage and size'),side=3,outer=T,line=2.2)
}
######################################################################
#Sex COMP Fit
######################################################################
PlotFfit=function()
{
  Year=seq(report$BeginYear,report$EndYear,1)
  
  Non_female=report$FemaleProp_at_Size
  Non_female_obs=report$FemaleProp_at_Size_obs
  
  Scomp_Obs=Non_female_obs
  Scomp_Pred=Non_female
  
  if (Nyear<25)
  {
    par(mfrow=c(5,5))
    temp=(Nyear%/%5*(5-1)+Nyear%%5):(Nyear%/%5*(5-1)+5)
    yaxis=c(1,6,11,16,21)
    xaxis=c(21,22,23,24,25,temp)
  }
  if (Nyear>25)
  {
    par(mfrow=c(6,6))
    yaxis=c(1,7,13,19,25,31)
    xaxis=c(28,29,30,31,32,33)
  }
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  for (i in 1:Nyear)
  {
    
    plot(sizebins,Scomp_Pred[i,],type='l',ylim=c(0,1),axes=FALSE)
    x=c(sizebins,rev(sizebins))
    y=c(rep(0,length(sizebins)),rev(Scomp_Pred[i,]))
    polygon(x,y,col="gray")
    mtext(paste(Year[i]),side=3,line=-1,adj=0.1,cex=0.6)
    #legend('topleft',paste(Year[i]),cex=1.1,bty='n')
    #rp=vector('expression',2)
    #rp[1]=substitute(expression(N==Value),list(Value=SurveyESS_Obs[i]))[2]
    #rp[2]=substitute(expression(effN==Value),list(Value=SurveyESS_Pred[i]))[2]
    #legend('topright',legend=rp,bty='n',y.intersp=0.8,yjust=1)
    lines(sizebins,Scomp_Obs[i,],col='red')
    #text(x=27,y=0.15,paste("effN=",SurveyESS_Pred[i]),cex=1)
    if(i %in% yaxis)
      axis(2,col='gray',at=seq(0,1,0.1))
    if(i %in% xaxis)
      axis(1,col='gray',at=sizebins)
    box(col='gray')
    
  }
  mtext('Length (mm)',side=1,outer=T,line=2.2)
  mtext('Proportion',side=2,outer=T,line=2.2)
  mtext(paste('Proportion of change sex for a give size,','red line=Obs'),
        side=3,outer=T,line=2.2)
}

######################################################################
#L50 COMP
######################################################################
PlotLfifty=function()
{
  par(mfrow=c(1,1))
  Year=seq(report$BeginYear,report$EndYear,1)
  L50=report$Lfifty
  #par(mfrow=c(2,1))
  par(mar = c(4, 4.5, 2, 0.5))
  plot(Year,L50,xlab="Year",ylab="Size (mm)",type='o',
       col='red',cex=2,pch=20,lty='dashed',cex.lab=1.3, cex.axis=1.3, 
       cex.main=1.3, cex.sub=1.3, main='L50')
}

#######################################################################
#Stock status
#######################################################################
PlotStockStatus=function()
{
  par(mfrow=c(1,1))
  Year=seq(report$BeginYear,report$EndYear,1)
  f=report$Fishing_Mortality 
  N=report$Abundance_at_Size[1:Nyear,]
  W=report$Weight_length
  B=apply(N*W,1,sum)
  R=report$Mean_Recruitment
  
  F=apply(f,2,sum)
  Fmax=report$Fmax
  F0.1=report$F0.1
  Fmsy=report$FMSY
  F30=report$F30SPR
  Bmsy=report$Bmsy
  
  SSB=report$Spawning_stock_Biomass
  SSBmsy=report$SSBmsy
  
  xlim=1.2*max(SSB,SSBmsy)
  ylim=1.2*max(F,Fmax,F0.1,Fmsy,F30)
  plot(SSB,F,type='o',col='black',xlim=c(0,xlim),ylim=c(0,ylim),
       xlab='Spawning stock biomass',ylab='Fishing mortality')
  
  abline(h=Fmsy,col='red',lwd=2)
  abline(v=Bmsy,col='blue',lwd=2)
  text(0.9*xlim,Fmsy,'Fmsy')
  text(Bmsy,0.9*ylim,'Bmsy')
  text(1700,0.1,'2013')
}
#######################################################################
#Stage specific biomass
#######################################################################

PlotSpB=function(season)
{
  Year=seq(report$BeginYear,report$EndYear,1)
  W=report$Weight_length
  Abun=report$Abundance_at_Size[(Nyear*(season-1)+1):(Nyear*season),]
  
  Non_female=apply((1-report$FemaleProp_at_Size)*Abun*W,1,sum)
  Female=apply(report$FemaleProp_at_Size*Abun*W,1,sum)
  Data=matrix(NA,ncol=length(Year),nrow=2)
  Data[1,]=Non_female
  Data[2,]=Female
  
  par(xpd=T, mar=par()$mar+c(0,0,0,4))
  barplot(Data, main="", ylab="Biomass (mt)", 
          col=heat.colors(2), space=0.1, cex.axis=0.6, las=2,
          names.arg=Year, cex=0.8,border=NA) 
  legend('topleft', c('Non-females','Females'), cex=0.8, fill=heat.colors(2),bty='n')
  par(mar=c(5, 4, 4, 2) + 0.1)
}

#######################################################################
#Retro SSB
#######################################################################
PlotRo_scale=function(beginYear,referenceYear,endYear,timestep)
{ 
  par(mar = c(4, 4, 4, 2), oma = c(3, 3, 2, 0.5))
  #par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  Nyear=endYear-beginYear+1
  Year=seq(beginYear,endYear,1)
  Nrow=endYear-referenceYear+1
  Nthcol=Nyear-(endYear-referenceYear)
  SSB=matrix(NA,ncol=Nyear,nrow=Nrow)
  SSB_s=matrix(NA,ncol=Nyear,nrow=Nrow)
  Recruit=matrix(NA,ncol=Nyear,nrow=Nrow)
  ER=matrix(NA,ncol=Nyear,nrow=Nrow)
  Recruit_s=matrix(NA,ncol=Nyear,nrow=Nrow)
  ER_s=matrix(NA,ncol=Nyear,nrow=Nrow)
  
  par(mfrow=c(1,1))
  for (i in 1:Nrow)
  {
    filename=paste("NSLSAP01_",beginYear,"_",endYear+1-i,"_",timestep,".rep",sep='')
    
	if(report$SexAtSizeLamda>0)
    {SSB[i,1:(Nyear+1-i)]=read.rep(filename)$Spawning_stock_Biomass
    }else
	SSB[i,1:(Nyear+1-i)]=read.rep(filename)$Spawning_stock_Biomass_input
	
    SSB_s[i,1:(Nyear+1-i)]=SSB[i,1:(Nyear+1-i)]/SSB[1,1:(Nyear+1-i)]
    if (i==1)
    {
      plot(Year,SSB_s[i,],type='l',lty=1,col=1,xlab='Year',lwd=2,
           ylab='Estiamte/Terminal year estimate',ylim=c(0,2),
           main="SSB")
      #points(Year[Nyear],SSB[i,Nyear],pch=16)
    }else
    {
      lines(Year,SSB_s[i,],lty=1,col=1,lwd=2)
      #points(Year[Nyear-i+1],SSB[i,Nyear-i+1],pch=16)
    }
  }
  
  
  par(mfrow=c(1,1))
  for (i in 1:Nrow)
  {
    filename=paste("NSLSAP01_",beginYear,"_",endYear+1-i,"_",timestep,".rep",sep='')
    
    Recruit[i,1:(Nyear+1-i)]=read.rep(filename)$recruitment_Pred
    Recruit_s[i,1:(Nyear+1-i)]=Recruit[i,1:(Nyear+1-i)]/Recruit[1,1:(Nyear+1-i)]
    
    if (i==1)
    {
      plot(Year,Recruit_s[i,],type='l',lty=1,col=1,xlab='Year',lwd=2,
           ylab='Estiamte/Terminal year estimate',ylim=c(0,3),
           main="Recruitment")
      #points(Year[Nyear],Recruit[i,Nyear],pch=16)
    }else
    {
      lines(Year,Recruit_s[i,],lty=1,col=1,lwd=2)
      #points(Year[Nyear-i+1],Recruit[i,Nyear-i+1],pch=16)
    }
  }
  
  
    par(mfrow=c(1,1))
    for (i in 1:Nrow)
    {
      filename=paste("NSLSAP01_",beginYear,"_",endYear+1-i,"_",timestep,".rep",sep='')
      Prop_female=read.rep(filename)$FemaleProp_at_Size
      catch1=read.rep(filename)$Catch_Pred[1,]*read.rep(filename)$Catch_Comp_Pred[1:(Nyear+1-i),]
      catch2=read.rep(filename)$Catch_Pred[2,]*read.rep(filename)$Catch_Comp_Pred[(Nyear+2-i):(2*(Nyear+1-i)),]
      catch3=read.rep(filename)$Catch_Pred[3,]*read.rep(filename)$Catch_Comp_Pred[(2*(Nyear+1-i)+1):(3*(Nyear+1-i)),]
      TCB=apply(read.rep(filename)$Weight_length[1:(Nyear+1-i),]*(catch1+catch2+catch3)*Prop_female,1,sum)
      TBF=read.rep(filename)$Abundance_at_Size[1:(Nyear+1-i),]*Prop_female*read.rep(filename)$Weight_length[1:(Nyear-i+1),]
      Exploitation_rate_TB=TCB/apply(TBF,1,sum)
      
      ER[i,1:(Nyear+1-i)]=Exploitation_rate_TB
      ER_s[i,1:(Nyear+1-i)]=ER[i,1:(Nyear+1-i)]/ER[1,1:(Nyear+1-i)]
      
      if (i==1)
      {
        plot(Year,ER_s[i,],type='l',lty=1,col=1,xlab='Year',lwd=2,
             ylab='Estiamte/Terminal year estimate',ylim=c(0,2),
             main="Exploitation rate")
        #points(Year[Nyear],ER[i,Nyear],pch=16)
      }else
      {
        lines(Year,ER_s[i,],lty=1,col=1,lwd=2)
        #points(Year[Nyear-i+1],ER[i,Nyear-i+1],pch=16)
      }
    }
  
}

PlotRo=function(beginYear,referenceYear,endYear,timestep)
{ 
  par(mar = c(4, 4, 4, 2), oma = c(3, 3, 2, 0.5))
  #par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  Nyear=endYear-beginYear+1
  Year=seq(beginYear,endYear,1)
  Nrow=endYear-referenceYear+1
  Nthcol=Nyear-(endYear-referenceYear)
  SSB=matrix(NA,ncol=Nyear,nrow=Nrow)
  Recruit=matrix(NA,ncol=Nyear,nrow=Nrow)
  ER=matrix(NA,ncol=Nyear,nrow=Nrow)
  
  par(mfrow=c(1,1))
  for (i in 1:Nrow)
  {
    filename=paste("NSLSAP01_",beginYear,"_",endYear+1-i,"_",timestep,".rep",sep='')
    
	if(report$SexAtSizeLamda>0)
    {SSB[i,1:(Nyear+1-i)]=read.rep(filename)$Spawning_stock_Biomass
    }else
	SSB[i,1:(Nyear+1-i)]=read.rep(filename)$Spawning_stock_Biomass_input
	
    if (i==1)
    {
      plot(Year,SSB[i,],type='l',lty=1,col=1,xlab='Year',ylab='Spawning biomass (mt)')
      points(Year[Nyear],SSB[i,Nyear],pch=16)
    }else
    {
      lines(Year,SSB[i,],lty=1,col=1)
      points(Year[Nyear-i+1],SSB[i,Nyear-i+1],pch=16)
    }
  }
  
  Mohn_SSB=sum((SSB[,Nthcol]-SSB[Nrow,Nthcol])/SSB[,Nthcol])
  print(Mohn_SSB)
  
  
  par(mfrow=c(1,1))
  for (i in 1:Nrow)
  {
    filename=paste("NSLSAP01_",beginYear,"_",endYear+1-i,"_",timestep,".rep",sep='')
    
    Recruit[i,1:(Nyear+1-i)]=read.rep(filename)$recruitment_Pred
    if (i==1)
    {
      plot(Year,Recruit[i,],type='l',lty=1,col=1,xlab='Year',ylab='Recruitment')
      points(Year[Nyear],Recruit[i,Nyear],pch=16)
    }else
    {
      lines(Year,Recruit[i,],lty=1,col=1)
      points(Year[Nyear-i+1],Recruit[i,Nyear-i+1],pch=16)
    }
  }
  
  Mohn_R=sum((Recruit[,Nthcol]-Recruit[Nrow,Nthcol])/Recruit[,Nthcol])
  print(Mohn_R)
  
  
  par(mfrow=c(1,1))
  for (i in 1:Nrow)
  {
    filename=paste("NSLSAP01_",beginYear,"_",endYear+1-i,"_",timestep,".rep",sep='')
    Prop_female=read.rep(filename)$FemaleProp_at_Size
    catch1=read.rep(filename)$Catch_Pred[1,]*read.rep(filename)$Catch_Comp_Pred[1:(Nyear+1-i),]
    catch2=read.rep(filename)$Catch_Pred[2,]*read.rep(filename)$Catch_Comp_Pred[(Nyear+2-i):(2*(Nyear+1-i)),]
    catch3=read.rep(filename)$Catch_Pred[3,]*read.rep(filename)$Catch_Comp_Pred[(2*(Nyear+1-i)+1):(3*(Nyear+1-i)),]
    TCB=apply(read.rep(filename)$Weight_length[1:(Nyear+1-i),]*(catch1+catch2+catch3)*Prop_female,1,sum)
    TBF=read.rep(filename)$Abundance_at_Size[1:(Nyear+1-i),]*Prop_female*read.rep(filename)$Weight_length[1:(Nyear-i+1),]
    Exploitation_rate_TB=TCB/apply(TBF,1,sum)
    
    
    ER[i,1:(Nyear+1-i)]=Exploitation_rate_TB
    if (i==1)
    {
      plot(Year,ER[i,],type='l',lty=1,col=1,xlab='Year',ylab='Exploitation rate',ylim=c(0,1))
      points(Year[Nyear],ER[i,Nyear],pch=16)
    }else
    {
      lines(Year,ER[i,],lty=1,col=1)
      points(Year[Nyear-i+1],ER[i,Nyear-i+1],pch=16)
    }
  }
  
  Mohn_ER=sum((ER[,Nthcol]-ER[Nrow,Nthcol])/ER[,Nthcol])
  print(Mohn_ER)
  
}

#######################################################################
#Get likelihoods
#######################################################################
Getlike=function()
{
  print(c(report$likelihood_total,report$likelihood_tcatch,report$likelihood_pcatch,
          report$likelihood_index,report$likelihood_pindex,
          report$likelihood_Recruit,report$likelihood_pfemal))
}


#####################exploitation rate#######################
PlotER=function()
{
  my.settings <- list(
    par.main.text = list(font = 1, # make it bold
                         just = "left", 
                         x = grid::unit(5, "mm")))
  
  library(lattice)
  ABUN=report$Abundance_at_Size[1:Nyear,]
  
  TN=apply(ABUN,1,sum)
  Exploitation_rate_TN=apply(report$Catch_Pred,2,sum)/TN

  Prop_female=report$FemaleProp_at_Size
  
  catch1=report$Catch_Pred[1,]*report$Catch_Comp_Pred[1:33,]*Prop_female
  catch2=report$Catch_Pred[2,]*report$Catch_Comp_Pred[67:99,]*Prop_female
  catch3=report$Catch_Pred[3,]*report$Catch_Comp_Pred[133:165,]*Prop_female
  
  TCFB=apply(report$Weight_length*(catch1+catch2+catch3),1,sum)
  TCFN=apply(catch1+catch2+catch3,1,sum)
  TNF=ABUN*Prop_female
  TBF=ABUN*Prop_female*report$Weight_length
  TBNF=ABUN*(1-Prop_female)*report$Weight_length
  
  Exploitation_rate_TFB=TCFB/apply(TBF,1,sum)
  x=report$BeginYear:report$EndYear
  
  catch1=report$Catch_Pred[1,]*report$Catch_Comp_Pred[1:33,]
  catch2=report$Catch_Pred[2,]*report$Catch_Comp_Pred[67:99,]
  catch3=report$Catch_Pred[3,]*report$Catch_Comp_Pred[133:165,]
  
  TCB=apply(report$Weight_length*(catch1+catch2+catch3),1,sum)
  TB=apply(report$Weight_length*ABUN,1,sum)
  Exploitation_rate_TB=TCFB/TB
  
  Exploitation_rate_CTB_TBF=TCB/apply(TBF,1,sum)
  Exploitation_rate_CTN_TNF=apply(report$Catch_Pred,2,sum)/apply(TNF,1,sum)
  Exploitation_rate_CTFB_TBF=Exploitation_rate_TFB
  
  ER<-xyplot(Exploitation_rate_TN+Exploitation_rate_TFB+Exploitation_rate_TB~x,
         t=c("l","l","l"),lty=c(1,2,3),xlab=list(label='Year',cex=1.3),
         ylab=list(label='Exploitation rate',cex=1.3),
         col.line = c('red', 'green',"blue"),lwd=2,
         scales=list(cex=1.3),par.settings=my.settings,
         main="red~TN; green~TFB; blue~TB")
  
  #print(c(Exploitation_rate_CTN_TNF[30],Exploitation_rate_CTB_TBF[30],Exploitation_rate_TFB[30]))
  Bmsy=report$Bmsy
  ratio=TB/Bmsy
  print(ER)
}
#############################################################
  
  