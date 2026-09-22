#######################################################################
#Fishery-specific Catch
#######################################################################
PlotC=function(fleet_id)
{
  bioORnum=report$BiomassORNum
  if (sum(bioORnum)==0)
  {
    ylab.name='Number (thousands)'
  }else{
    ylab.name='Weight (metric tons)'
  }  
  
  n.col=1 #num of coln per plot
  n.row=2 #num of row per plot
  n.pannel=n.col*n.row #num of pannels per plot
  n.plot=ceiling(Nseason/n.pannel)
  x.axis=rep(0,n.col)
    
  for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		par(mfrow=c(n.row,n.col))
		par(mar=c(1,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
	   }
	
	  if(hh==n.plot){
		n.fig=Nseason-(n.plot-1)*n.pannel
		par(mfrow=c(n.row,n.col))
		par(mar=c(1,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
	  }
	  
	  x.axis=seq(n.fig-n.col+1,n.fig,1)
	  
	  for(pp in 1:n.fig){
	  	season_id=(hh-1)*n.pannel+pp
	  	obs=report$Catch_Obs[Nseason*(fleet_id-1)+season_id,]        #fleet (season by year)
        pred=report$Catch_Pred[Nseason*(fleet_id-1)+season_id,]
        ylim_plot=c(0,max(max(obs),max(pred))*1.3)
        
        plot(year,pred,type='l',xlab='',ylab='',main=NULL,lwd=2,col='black',ylim=ylim_plot,xaxt='n')
        points(year,obs,pch=1,col='red')
        legend('top',paste(fleet_name[fleet_id],'Season',season_id,by=' '),bty='n')
        
        if(pp %in% x.axis){
         	axis(1,at=year_axis,las=3)
        }
        
        if (pp==1)
        {
      	   legend('topleft',c('Obs','Pred'),lty=c(NA,1),pch=c(1,NA),lwd=c(NA,2),col=c('red','black'),bty='o',y.intersp=0.8)
        }	  	
	  } #end pp
	  
	  mtext('Year',side=1,outer=T,line=3)
      mtext(ylab.name,side=2,outer=T,line=2.5) 
	  
  } #end hh
}
######################################################################

#######################################################################
#Survey index
#######################################################################
PlotSI=function(index_id,flag)
{
  par(mfrow=c(2,1))
  par(oma = c(3, 0.5, 0.5, 0.5))
  par(mar = c(1, 4, 0, 0))
  par(tcl = -0.25)
    
  #########obs and pred
  obs=report$Survey_Index_Obs[index_id,]      #index by year
  obs[which(obs==0)]=NA                       # zero obs will not be predicted
  pred=report$Survey_Index_Pred[index_id,]
  pred[which(pred==0)]=NA
  ylim_plot=c(0,max(max(obs,na.rm=T),max(pred,na.rm=T))*1.3)
  
  plot(1,type='n',xlab='Year',ylab='Abundance Index',main=NULL,xlim=c(year[1],year[Nyear]),ylim=ylim_plot,xaxt='n')
  lines(year,pred,type='l',lwd=2,col='black')
  if(length(which(is.na(pred)))>0){      #if there are 0 obs/pred in this index, add points to show pred for obs next to zeros because line plot cannot show
  	points(year,pred,pch=20,cex=0.3,col='black')
  }
  points(year,obs,pch=1,col='red')
  axis(1,at=year_axis,las=3)
  legend('top',index_name[index_id],bty='n')
  legend('topleft',c('Obs','Pred'),lty=c(NA,1),pch=c(1,NA),lwd=c(NA,2),col=c('red','black'),bty='o',y.intersp=0.8)
  
  ########deviation   
  if(flag<=4){     #deviation on log-scale: log(pred)-log(obs)
  	res=log(pred)-log(obs)
  	ylab.name='Deviation (log-scale)'
  }else{           #deviation on regular scale: pred-obs
  	res=pred-obs
  	ylab.name='Deviation'
  }
  
  #########test
  test.result=matrix(NA,nrow=1,ncol=4+1)
  colnames(test.result)=c('Index','Runs','RunsPvalue','SW','SWPvalue')
  test.result[1,'Index']=index_name[index_id]

  ####Runs randomness test 
  r.test=runs.test(res)
  test.result[1,'Runs']=r.test$statistic
  test.result[1,'RunsPvalue']=r.test$p.value
  
  ####Shapiro-Wilk normality test
  sw.test=shapiro.test(res)
  test.result[1,'SW']=sw.test$statistic
  test.result[1,'SWPvalue']=sw.test$p.value
  write.table(test.result,paste('IndexTest_',index_id,'.txt',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
  
  ####plot 
  par(mar = c(1, 4, 4, 0))
  par(tcl = -0.25)
  res[which(is.na(res))]=0    #for plotting purpose
  ylim_plot=c(min(res),max(res))*1.3
  barplot(res,xlab='Year',ylab=ylab.name,ylim=ylim_plot,space=2,col='black',xaxt='n')
  box(col='black')
  mtext('Year',side=1,outer=T,line=0.5)
}

#######################################################################
#Catch length comp 
#######################################################################
PlotCC=function(fleet_id,season_id)
{
    n.col=5 #num of coln per plot
    n.row=6 #num of row per plot
    n.pannel=n.col*n.row #num of pannels per plot
    n.plot=ceiling(Nyear/n.pannel)
    y.axis=rep(0,n.row)
    x.axis=rep(0,n.col)

    for(i in 1:n.row){
    	y.axis[i]=(i-1)*n.col+1
    }

    temp1=Nseason*(fleet_id-1)+season_id     
    temp2=Nyear*Nseason*(fleet_id-1)+Nyear*(season_id-1)
    ylim_max=max(max(report$Catch_Comp_Obs[temp2+1:Nyear,]),max(report$Catch_Comp_Pred[temp2+1:Nyear,]))

    for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
	   }
	
	  if(hh==n.plot){
		n.fig=Nyear-(n.plot-1)*n.pannel
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
	  }
	  
      x.axis=seq(n.fig-n.col+1,n.fig,1)
      
	  for(pp in 1:n.fig){
	  	 year_id=(hh-1)*n.pannel+pp

  	     temp=temp2+year_id  
  	     comp_obs=report$Catch_Comp_Obs[temp,]  #(fleet,season,year) by size
         comp_pred=report$Catch_Comp_Pred[temp,]
         ess_obs=report$Catch_ESSinput[temp1, year_id] #(fleet,season) by year
         ess_pred=round(report$Catch_ESSpred[temp1, year_id],2)
         #ylim_plot=c(0,max(max(comp_obs),max(comp_pred))*1.1)
         ylim_plot=c(0,ylim_max*1.1)
    
  	     plot(1,type='n',xlab='',ylab='',main=NULL,lwd=1,col='black',xlim=c(sizebin[1],sizebin[Nsizebin]),ylim=ylim_plot,xaxt='n',yaxt='n')
         polygon(c(sizebin,rev(sizebin)),c(rep(0,Nsizebin),rev(comp_obs)),col='grey',border=F) 
         lines(sizebin,comp_pred,type='l',lwd=2,col='black')
         mtext(year[year_id],side=3,line=-1,adj=0.1,cex=0.6)
         
         if(pp==1){
         	legend('topright',c(paste('ESS_Obs=',ess_obs,sep=''),paste('ESS_Pred=',ess_pred,sep=''),'Obs','Pred'),lty=c(NA,NA,NA,1),pch=c(NA,NA,15,NA),lwd=c(NA,NA,NA,2),col=c('black','black','grey','black'),bty='o', y.intersp=0.8,cex=0.6)
         }else{   
            legend('topright',c(paste('ESS_Obs=',ess_obs,sep=''),paste('ESS_Pred=',ess_pred,sep='')),bty='o', y.intersp=0.8,cex=0.6)
	  	 }  		  	 
	  	 if(pp %in% y.axis){
	  	 	axis(2,col='gray',at=seq(ylim_plot[1],ylim_plot[2]-0.2,by=0.2))
	  	 }      
         if(pp %in% x.axis){
         	axis(1,col='gray',at=size_axis)
         }
         box(col='gray')
	  }

      mtext('Mid-Length (mm)',side=1,outer=T,line=2.5)
      mtext('Proportion',side=2,outer=T,line=2.5)
	  mtext(paste(fleet_name[fleet_id],'Season',season_id,by=' '),3,outer=T,line=2)
  }  
}
######################################################################

#######################################################################
#Survey Length Comp
#######################################################################
PlotSIC=function(index_id)
{
    n.col=5 #num of coln per plot
    n.row=6 #num of row per plot
    n.pannel=n.col*n.row #num of pannels per plot
    n.plot=ceiling(Nyear/n.pannel)
    y.axis=rep(0,n.row)
    x.axis=rep(0,n.col)

    for(i in 1:n.row){
    	y.axis[i]=(i-1)*n.col+1
    }

    temp1=index_id
    temp2=Nyear*(index_id-1)
    ylim_max=max(max(report$Survey_Comp_Obs[temp2+1:Nyear,]),max(report$Survey_Comp_Pred[temp2+1:Nyear,]))

    for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
	   }
	
	  if(hh==n.plot){
		n.fig=Nyear-(n.plot-1)*n.pannel
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
	  }
	  
      x.axis=seq(n.fig-n.col+1,n.fig,1)
      
	  for(pp in 1:n.fig){
	  	 year_id=(hh-1)*n.pannel+pp

  	     temp=temp2+year_id
  	     comp_obs=report$Survey_Comp_Obs[temp,]
         comp_pred=report$Survey_Comp_Pred[temp,]
         ess_obs=report$Survey_ESSinput[temp1,year_id]
         ess_pred=round(report$Survey_ESSpred[temp1,year_id],2)
    
         ylim_plot=c(0,ylim_max*1.1)
        
  	     plot(1,type='n',xlab='',ylab='',main=NULL,lwd=1,col='black',xlim=c(sizebin[1],sizebin[Nsizebin]),ylim=ylim_plot,xaxt='n',yaxt='n')
         polygon(c(sizebin,rev(sizebin)),c(rep(0,Nsizebin),rev(comp_obs)),col='grey',border=F) 
         lines(sizebin,comp_pred,type='l',lwd=2,col='black')
         mtext(year[year_id],side=3,line=-1,adj=0.1,cex=0.6)
         
         if(pp==1){
         	legend('topright',c(paste('ESS_Obs=',ess_obs,sep=''),paste('ESS_Pred=',ess_pred,sep=''),'Obs','Pred'),lty=c(NA,NA,NA,1),pch=c(NA,NA,15,NA),lwd=c(NA,NA,NA,2),col=c('black','black','grey','black'),bty='o', y.intersp=0.8,cex=0.6)
         }else{   
            legend('topright',c(paste('ESS_Obs=',ess_obs,sep=''),paste('ESS_Pred=',ess_pred,sep='')),bty='o', y.intersp=0.8,cex=0.6)
	  	 }  		  	 
	  	 if(pp %in% y.axis){
	  	 	axis(2,col='gray',at=seq(ylim_plot[1],ylim_plot[2]-0.2,by=0.2))
	  	 }      
         if(pp %in% x.axis){
         	axis(1,col='gray',at=size_axis)
         }
         box(col='gray')
	  }

      mtext('Mid-Length (mm)',side=1,outer=T,line=2.5)
      mtext('Proportion',side=2,outer=T,line=2.5)
	  mtext(index_name[index_id],3,outer=T,line=2)
  }  
}
######################################################################

#######################################################################
#Selectivity of fisheries
#######################################################################
PlotSelF=function(fleet_id)
{
  par(mfrow=c(1,1))
  par(mar = c(4, 4, 4, 0), oma = c(0.5, 0.5, 0.5, 0.5))
  par(tcl = -0.25)
  
  ylim_plot=c(0,1)
  
  plot(1,type='n',xlab='Mid-length (mm)',ylab='Selectivity',main=fleet_name[fleet_id],xlim=c(min(size_axis),max(size_axis)),ylim=ylim_plot,xaxt='n',yaxt='n')
  axis(1,at=size_axis)
  axis(2,at=seq(0,1,by=0.2))
  
  for(i in 1: Nblock_selF[fleet_id]){
  	year_id_temp=block_selF[[fleet_id]][i,1]
    season_id=block_selF[[fleet_id]][i,2]
    # block_id=block_selF[[fleet_id]][i,3]
    
    year_id=year_id_temp-report$Begin_Year+1    
    itemp=Nyear*Nseason*(fleet_id-1)+Nyear*(season_id-1)
    sel=report$Fleet_Selectivity[itemp+year_id,]      #(fleet,season,year) by size
    
    lines(sizebin,sel,type='o',col=i,lwd=2,pch=16,lty=i)
  }
  
  legend.name=rep(NA,Nblock_selF[fleet_id])
  for(j in 1:length(legend.name)){
  	legend.name[j]=paste("Season",block_selF[[fleet_id]][j,2],block_selF[[fleet_id]][j,1],sep=' ')
  }
  
  legend('right',legend=legend.name,lty=1:Nblock_selF[fleet_id],lwd=rep(2,Nblock_selF[fleet_id]),col=1:Nblock_selF[fleet_id],bty='o',y.intersp=0.8)
}
######################################################################

#######################################################################
#Selectivity of surveys
#######################################################################
PlotSelSI=function(index_id)
{
  par(mfrow=c(1,1))
  par(mar = c(4, 4, 4, 0), oma = c(0.5, 0.5, 0.5, 0.5))
  par(tcl = -0.25)
  
  ylim_plot=c(0,1)
  maintemp=index_name[index_id]
  
  plot(1,type='n',xlab='Mid-length (mm)',ylab='Selectivity',main=maintemp,xlim=c(min(size_axis),max(size_axis)),ylim=ylim_plot,xaxt='n',yaxt='n')
  axis(1,at=size_axis)
  axis(2,at=seq(0,1,by=0.2))

  for(i in 1: Nblock_selSI[index_id]){

  year_id_temp=block_selSI[[index_id]][i,1]
  # block_id=block_selSI[[index_id]][i,2]
  year_id=year_id_temp-report$Begin_Year+1
  
  itemp=Nyear*(index_id-1)
  sel=report$Index_Selectivity[itemp+year_id,]   #(index,year) by size
  
  lines(sizebin,sel,type='o',col=i,lwd=2,pch=16,lty=i)
  }
  
  legend('topleft',legend=block_selSI[[index_id]][,1],lty=1:Nblock_selSI[index_id],lwd=rep(2,Nblock_selSI[index_id]),col=1:Nblock_selSI[index_id],bty='o',y.intersp=0.8)
}
######################################################################

#######################################################################
#Fishing mortality
#######################################################################
PlotF=function(fleet_id)
{
	par(mfrow=c(2,1))
    par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
        
    itemp=Nseason*(fleet_id-1)
    
    for (i in 1:Nseason)
    {   	
    	  ff=report$Fishing_Mortality[itemp+i,]   #(fleet,season) by year
    	  ylim_plot=c(0,max(ff)*1.3)
    	
    	  plot(year,ff,type='o',xlab='',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
    	  
    	  if(i==Nseason)
      {
      	 axis(1,at=year_axis,las=3)
      }
      
      legend('top',paste(fleet_name[fleet_id],'Season',i,by=' '),bty='n')
      
      if (i==1)
      {
      	 legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
      }
    }
    
    mtext('Year',side=1,outer=T,line=2.5)
    mtext('Fishing mortality (per season)',side=2,outer=T,line=2.5)   
}

PlotF_season=function()
{
	par(mfrow=c(2,1))
    par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    
    for (i in 1:Nseason)
    {
      # ff=rep(0,Nyear)
      # for(fleet_id in 1:Nfleet)
      # {
    	   # itemp=Nseason*(fleet_id-1)
    	   # ff_temp=report$Fishing_Mortality[itemp+i,]   #(fleet,season) by year
    	   # ff=ff+ff_temp
      # }
    	   
       ff=report$Fishing_Mortality_Season[i,]   	
       ylim_plot=c(0,max(ff)*1.3)
    	
       plot(year,ff,type='o',xlab='',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
    	
       if(i==Nseason)
       {
      	  axis(1,at=year_axis,las=3)
       }
      
       legend('top',paste('Season',i,by=' '),bty='n')
      
       if (i==1)
       {
      	  legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
       }
    }
    
    mtext('Year',side=1,outer=T,line=2.5)
    mtext('Fishing mortality (per season)',side=2,outer=T,line=2.5)   
}

PlotF_all=function()
{
	par(mfrow=c(1,1))
	par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
	par(tcl = -0.25)
    
    # ff=rep(0,Nyear)
    # for (i in 1:Nseason)
    # {    	   
    	   # for(fleet_id in 1:Nfleet)
    	   # {
    	   	  # itemp=Nseason*(fleet_id-1)
    	   	  # ff_temp=report$Fishing_Mortality[itemp+i,]   #(fleet,season) by year
    	   	  # ff=ff+ff_temp
    	   # }
    # }
    
    ff=report$Fishing_Mortality_Annual
    std=std_est$std$sdr_vAnnualF
    up=ff+2*std
    lo=ff-2*std
    lo[which(lo<0)]=0
       	
    ylim_plot=c(0,max(up)*1.3)
    	
    plot(year,ff,type='o',xlab='',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
    lines(year,up,lty=3,lwd=2,col='black')
    lines(year,lo,lty=3,lwd=2,col='black')
    axis(1,at=year_axis,las=3)
    legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
    
    mtext('Year',side=1,outer=T,line=2.5)
    mtext('Fishing mortality (per year)',side=2,outer=T,line=2.5)
 
}

PlotU_all=function()
{
	par(mfrow=c(1,1))
	par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
	par(tcl = -0.25)
    
    # ff=rep(0,Nyear)
    # for (i in 1:Nseason)
    # {    	   
    	   # for(fleet_id in 1:Nfleet)
    	   # {
    	   	  # itemp=Nseason*(fleet_id-1)
    	   	  # ff_temp=report$Fishing_Mortality[itemp+i,]   #(fleet,season) by year
    	   	  # ff=ff+ff_temp
    	   # }
    # }
    
    ff=report$Fishing_Mortality_Annual
    uu=1-exp(-ff)
    	   	
    ylim_plot=c(min(uu),1)
    	
    plot(year,uu,type='o',xlab='',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
    axis(1,at=year_axis,las=3)
    legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
    
    mtext('Year',side=1,outer=T,line=2.5)
    mtext('Exploitation rate',side=2,outer=T,line=2.5)
}

######################################################################

#######################################################################
#Natural mortality
#######################################################################
PlotM_Year=function()
{
    n.col=5 #num of coln per plot
    n.row=6 #num of row per plot
    n.pannel=n.col*n.row #num of pannels per plot
    n.plot=ceiling(Nyear/n.pannel)
    y.axis=rep(0,n.row)
    x.axis=rep(0,n.col)
    
    for(i in 1:n.row){
    	y.axis[i]=(i-1)*n.col+1
    }    
    
    ylim_plot=c(0,max(report$Natural_Mortality)*1.3)
    
    name_temp=rep(NA,Nseason)
    for(j in 1:Nseason){
    	name_temp[j]=paste('Season ',j,sep='')
    }
    
    for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
	   }
	
	  if(hh==n.plot){
		n.fig=Nyear-(n.plot-1)*n.pannel
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
	  }
	  
	  x.axis=seq(n.fig-n.col+1,n.fig,1)
	  
	  for(pp in 1:n.fig){
	  	 year_id=(hh-1)*n.pannel+pp
	  	 	  	 
	  	 # plot(age,y,type='l',xlab="",ylab="",main=NULL,col='black',lwd=1,ylim=ylim_plot,axes=F)
	  	 plot(1,type='n',xlab='',ylab='',main=NULL,xlim=c(min(size_axis),max(size_axis)),ylim=ylim_plot,axes=F)	  	 
	  	 mtext(year[year_id],side=3,line=-1,adj=0.1,cex=0.6)
         
         for(i in 1:Nseason){
         	itemp=Nyear*(i-1)+year_id
  	        mm=report$Natural_Mortality[itemp,]  # (season,year) by size
    	    lines(sizebin,mm,type='l',col=i,lty=i,lwd=1)
         }
	  	 	  	 
	  	 if(pp %in% y.axis){
	  	 	axis(2,col='gray',at=size_axis)
	  	 }      
         if(pp %in% x.axis){
         	axis(1,col='gray',at=size_axis)
         }
         box(col='gray')
        
         if(pp==1){
        	legend('topright',name_temp,lty=1:Nseason,lwd=rep(1,Nseason),col=1:Nseason,bty='o',y.intersp=0.8,cex=0.8,box.col='grey')
         }        
	  }
	
	  mtext('Mid-length (mm)',1,outer=T,line=2.5)
	  mtext('Natural mortality (per season)',2,outer=T,line=2.5)
	  mtext('Natural mortality by size',3,outer=T,line=2)
    }
}

PlotM_base=function(flag,season_id)  #specify flag and season; flag: specify flag:0-input M and not plot; >0 estimate M and plot
{
  if(flag>0)
  {  	
	par(mfrow=c(2,1))
    par(mar = c(1, 4, 0, 0), oma = c(4, 0.5, 0.5, 0.5))
    par(tcl = -0.25)     
       
    yy=report$Natural_Mortality_Base_With_YrScalar[season_id,]
    ylim_plot=c(0,max(yy)*1.3)
    
    if(flag>1)
    {
       yy_dev=report$Natural_Mortality_Dev_Logscale[season_id,]
       ylim_plot_dev=c(min(yy_dev)*1.3,max(yy_dev)*1.3) 
    }else{
       yy_dev=rep(NA,Nyear)
       ylim_plot_dev=c(-1,1)
    } 	       	   	       	
    
    plot(year,yy,type='o',xlab='',ylab='Natural mortality (per season)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
    legend('top',paste('Season',season_id,by=' '),bty='n')
    
    plot(year,yy_dev,type='o',xlab='',ylab='Deviation (log-scale)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot_dev,xaxt='n')
    abline(h=0,col="black",lwd=1)
    axis(1,at=year_axis,las=3)
    
    mtext('Year',side=1,outer=T,line=2.5)
  }
}

#######################################################################
#Env variables fit to Mdev
#######################################################################
PlotEnvMdev=function(flag,season_id)
{
  if(flag>0){
  par(mfrow=c(2,1))
  par(oma = c(3, 0.5, 0.5, 0.5))
  par(mar = c(1, 4, 0, 0))
  par(tcl = -0.25)
  
  #########obs and pred
  miss_index=report$Env_Fit_Mdev_Zscore_Missing_Data_Index[season_id,]  #index for years with missing data;0-missing data;1-no missing data
  obs=report$Env_Fit_Mdev_Zscore[season_id,]
  obs[which(miss_index==0)]=NA      
  pred=report$Natural_Mortality_Dev_Logscale[season_id,]
  ylim_plot=c(min(min(obs,na.rm=T),min(pred))*1.3,max(max(obs,na.rm=T),max(pred))*1.3)
  env_name=paste('Natural mortality env fit (Season ',season_id,')',sep='')
  
  plot(1,type='n',xlab='Year',ylab='Z-score',main=NULL,xlim=c(year[1],year[Nyear]),ylim=ylim_plot,xaxt='n')
  lines(year,pred,type='l',lwd=2,col='black')
  points(year,obs,pch=1,col='red')
  # abline(h=0,col="gray",lwd=2)
  axis(1,at=year_axis,las=3)
  legend('top',env_name,bty='n')
  legend('topleft',c('Obs Z-score','Pred deviation'),lty=c(NA,1),pch=c(1,NA),lwd=c(NA,2),col=c('red','black'),bty='o',y.intersp=0.8)
  
  ########residuals
  par(mar = c(1, 4, 4, 0))
  par(tcl = -0.25)
  
  res=pred-obs
  res[which(miss_index==0)]=0       #for ploting purpose
  ylim_plot=c(min(res),max(res))*1.3
  barplot(res,xlab='Year',ylab='Residuals',ylim=ylim_plot,space=2,col='black',xaxt='n')
  box(col='black')
  mtext('Year',side=1,outer=T,line=0.5)
  }
}
######################################################################

# PlotM_Input=function()
# {
  # par(mfrow=c(1,1))
  # par(mar = c(4, 4, 0, 0), oma = c(0.5, 0.5, 0.5, 0.5))
  # par(tcl = -0.25)
  
  # year_id=c(1991,1991,cold_year[1],cold_year[1])  #years that are being ploted
  # year_id=year_id-report$Begin_Year+1   #convert years to year id
  # season_id=c(1,2,1,2) #seasons of the above years that are being ploted
  
  # ylim_plot=c(0,0.4)
  # legend_name=c('Non-winter season, non-cold-stun year','Winter season, non-cold-stun year','Non-winter season, cold-stun year','Winter season, cold-stun year')
  
  # plot(1,type='n',xlab='Mid-length (mm)',ylab='Natural mortality (per month)',main=NULL,xlim=c(min(sizebin),max(sizebin)),ylim=ylim_plot,xaxt='n',yaxt='n')
  
  # for(i in 1:length(year_id))
  # {
  	# itemp=Nyear*(season_id[i]-1)+year_id[i]
  	# mm=report$Natural_Mortality[itemp,]  # (season,year) by size
  	# mm=round(mm/Nmonth_season[season_id[i]],3)   #scale to monthly M for ploting and comparing
  	
  	# lines(sizebin,mm,type='o',lwd=2,col=i,pch=i,lty=i)
  # }
    # axis(1,at=size_axis)
    # axis(2,at=seq(0,ylim_plot[2],by=0.1))
    # legend('topright',legend=legend_name,lty=1:length(year_id),lwd=2,col=1:length(year_id),pch=1:length(year_id),bty='o',y.intersp=0.8)
     
# }
######################################################################

######################################################################
#N at size at the begining of a year
#####################################################################

##########Bubble function
panel.bubble.for.R<-function(x,y,z,subscripts,scalem=4,...){
  #panel.grid(col="lightgrey",lwd=1,h=-1,v=-1)
  cex<-abs(z) # size of points depends on the absolute value of z and is consistent across panels
  cex<-scalem*sqrt(cex/max(cex))
  flagplus<- z>=0
  panel.xyplot(x[flagplus],y[flagplus],pch=16,col="gray",cex=cex[flagplus])
  panel.xyplot(x[!flagplus],y[!flagplus],pch=1,col="gray",cex=cex[!flagplus])
}

PlotN_atSize=function()
{
  library(lattice)
  
  nn=matrix(NA,ncol=3,nrow=Nyear*Nsizebin)
  itemp=1
  
  for (i in 1:Nyear)
  {
    for (j in 1:Nsizebin)
    {      
        nn[itemp,]=c(year[i],sizebin[j],report$Abundance_At_Size[(1-1)*Nyear+i,j]) #(season, year) by size
        itemp=itemp+1            
    }
  }
  
  nn=as.data.frame(nn)    
  names(nn)=c("year","size","abun")
  
  plot<-xyplot(size ~ year,data=nn,layout=c(1,1),par.strip.text=list(cex=1.3),
               as.table=T,subscript=T,z=nn[,'abun'],
               main="Abundance at the beginning of year (thousands)",
               scales=list(y=list(cex=1.3),x=list(cex=1.3)),panel=panel.bubble.for.R,scalem=4,
               ylab=list(label="Mid-Length (mm)",cex=1.3),xlab=list(label="Year", cex=1.3))
  print(plot)  
}

PlotN_all=function()
{
  par(mfrow=c(1,1))
  par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
    
  # nn=rep(0,Nyear)  
  # for (i in 1:Nyear)
  # {
  	# nn_temp=report$Abundance_At_Size[(1-1)*Nyear+i,]  #(season, year) by size
  	# nn[i]=round(sum(nn_temp),3)
  # }

  nn=report$Annual_Abundance_At_Beginning_Of_Year
  std=std_est$std$sdr_vAnnualN
  up=nn+2*std
  lo=nn-2*std
  lo[which(lo<0)]=0

    
  ylim_plot=c(0,max(up)*1.3)
  
  plot(year,nn,type='o',xlab='',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black')
  axis(1,at=year_axis,las=3)  
  legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
  
  mtext('Year',side=1,outer=T,line=2.5)
  mtext('Abundance at the beginning of year (thousands)',side=2,outer=T,line=2.5)
}


PlotN_all_noR=function()
{
  par(mfrow=c(1,1))
  par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
    
  nn=report$Annual_Abundance_At_Beginning_Of_Year_No_Recruits
  std=std_est$std$sdr_vAnnualN_NoR
  up=nn+2*std
  lo=nn-2*std
  lo[which(lo<0)]=0

    
  ylim_plot=c(0,max(up)*1.3)
  
  plot(year,nn,type='o',xlab='Year',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black')
  axis(1,at=year_axis,las=3)  
  legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
  
  mtext('Year',side=1,outer=T,line=2.5)
  mtext('Abundance at the beginning of year without recruits (thousands)',side=2,outer=T,line=2.5)
}

######################################################################

#######################################################################
#Recruitment
#######################################################################
PlotR=function()
{
  par(mfrow=c(2,1))
  par(mar = c(1, 4, 0, 0), oma = c(4, 0.5, 0.5, 0.5))
  par(tcl = -0.25)
  
  rr=report$Recruitment_Pred
  Rdev=report$Recruitment_Dev_Logscale
  
  std.r=std_est$std$sdr_vRecruitmentP
  up.r=rr+2*std.r
  lo.r=rr-2*std.r
  lo.r[which(lo.r<0)]=0
  
  std.dev=std_est$std$PbnvRecruitDevPara
  up.dev=Rdev+2*std.dev
  lo.dev=Rdev-2*std.dev
    
  ylim_plot=c(0,max(up.r)*1.3)
  ylim_plot_dev=c(min(lo.dev)*1.3,max(up.dev)*1.3)
   
  plot(year,rr,type='o',xlab='',ylab='Recruits (thousands)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up.r,lty=3,lwd=2,col='black')
  lines(year,lo.r,lty=3,lwd=2,col='black')
  legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
  
  plot(year,Rdev,type='o',xlab='',ylab='Deviation (log-scale)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot_dev,xaxt='n')
  lines(year,up.dev,lty=3,lwd=2,col='black')
  lines(year,lo.dev,lty=3,lwd=2,col='black')  
  abline(h=0,col="black",lwd=1)  
  axis(1,at=year_axis,las=3)
  
  mtext('Year',side=1,outer=T,line=2.5)
  }
######################################################################

#######################################################################
#Env variables fit to Rdev
#######################################################################
PlotEnvRdev=function(flag)
{
  if(flag>0){
  par(mfrow=c(2,1))
  par(oma = c(3, 0.5, 0.5, 0.5))
  par(mar = c(1, 4, 0, 0))
  par(tcl = -0.25)
  
  #########obs and pred
  miss_index=report$Env_Fit_Rdev_Zscore_Missing_Data_Index  #index for years with missing data;0-missing data;1-no missing data
  obs=report$Env_Fit_Rdev_Zscore
  obs[which(miss_index==0)]=NA      
  pred=report$Recruitment_Dev_Logscale
  ylim_plot=c(min(min(obs,na.rm=T),min(pred))*1.3,max(max(obs,na.rm=T),max(pred))*1.3)
  env_name='Recruitment env fit'
  
  plot(1,type='n',xlab='Year',ylab='Z-score',main=NULL,xlim=c(year[1],year[Nyear]),ylim=ylim_plot,xaxt='n')
  lines(year,pred,type='l',lwd=2,col='black')
  points(year,obs,pch=1,col='red')
  # abline(h=0,col="gray",lwd=2)
  axis(1,at=year_axis,las=3)
  legend('top',env_name,bty='n')
  legend('topleft',c('Obs Z-score','Pred deviation'),lty=c(NA,1),pch=c(1,NA),lwd=c(NA,2),col=c('red','black'),bty='o',y.intersp=0.8)
  
  ########residuals
  par(mar = c(1, 4, 4, 0))
  par(tcl = -0.25)
  
  res=pred-obs
  res[which(miss_index==0)]=0  #for plotting purpose
  ylim_plot=c(min(res),max(res))*1.3
  barplot(res,xlab='Year',ylab='Residuals',ylim=ylim_plot,space=2,col='black',xaxt='n')
  box(col='black')
  mtext('Year',side=1,outer=T,line=0.5)
  }
}
######################################################################

#######################################################################
#SSB
#######################################################################
PlotSSB=function()
{
  par(mfrow=c(1,1))
  par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
  
  ssb=report$Spawning_Stock_Pred
  std=std_est$std$sdr_vSSB
  up=ssb+2*std
  lo=ssb-2*std
  lo[which(lo<0)]=0
    
  ylim_plot=c(0,max(up)*1.3)
  
  plot(year,ssb,type='o',xlab='',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black') 
  axis(1,at=year_axis,las=3)
  legend('topleft',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
  
  mtext('Year',side=1,outer=T,line=2.5)
  mtext('Spawning stock biomass (metric tons)',side=2,outer=T,line=2.5)
}


######################################################################

#######################################################################
#SR
#######################################################################
PlotSR=function()
{
  par(mfrow=c(1,1))
  par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
  
  rr=report$Recruitment_Pred
  ssb=report$Spawning_Stock_Pred
  
  ylim_plot=c(0,max(rr)*1.3)
  
  plot(ssb,rr,type='p',xlab='',ylab='',main=NULL,col='black',pch=1,ylim=ylim_plot)
  
  mtext('Spawning stock biomass (metric tons)',side=1,outer=T,line=2.5)
  mtext('Recruits (thousands)',side=2,outer=T,line=2.5)
}

#######################################################################
###SPRatio (spawning potential ratio)
#######################################################################
PlotSPRatio=function()
{
    par(mfrow=c(1,1))
    par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    
	thresh=ThreshSPR            #report$ThreshSPR
	target=TargetSPR            #report$TargetSPR
	yy=BRPsummary[1:Nyear,'spratio']   #report$Annual_SPRatio
	ylim_plot=c(0,0.6)
	
	plot(year,yy,type='o',xlab='',ylab='',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n',yaxt='n')
    axis(1,at=year_axis,las=3)
    axis(2,at=seq(0,ylim_plot[2],by=0.2))
    abline(h=thresh,col='black',lwd=1,lty=1)
    abline(h=target,col='black',lwd=1,lty=2)
    legend('topleft',c(paste('Threshold ',thresh*100,'%',sep=''),paste('Target ',target*100,'%',sep='')),lty=c(1,2),lwd=c(1,1),col=c('black','black'),bty='o',y.intersp=0.8) 
    
    mtext('Year',side=1,outer=T,line=2.5)
    mtext('Spawning potential ratio',side=2,outer=T,line=2.5)    
}

#######################################################################
#Stock status
#######################################################################
BRP_names=c(paste(ThreshSPR*100,'%SPR',sep=''),paste(TargetSPR*100,'%SPR',sep=''))

PlotStockStatus=function(fthresh,ftarget,sthresh,starget)
{  
  par(mfrow=c(2,1))
  par(mar = c(0.5, 4, 0, 0), oma = c(4, 0.5, 0.5, 0.5))
  par(tcl = -0.25)
  
  ff=report$Fishing_Mortality_Annual
  std=std_est$std$sdr_vAnnualF
  up=ff+2*std
  lo=ff-2*std
  lo[which(lo<0)]=0
       	
  ylim_plot=c(0,max(up,fthresh,ftarget,na.rm=T)*1.3)
    	
  plot(year,ff,type='o',xlab='',ylab='Fishing mortality (per year)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black')
  lines(year,fthresh,type='l',lty=1,lwd=1,col='black')
  lines(year,ftarget,lty=2,lwd=1,col='black')
    
  # legend('topleft',c(paste('Threshold ',thresh*100,'%',sep=''),paste('Target ',target*100,'%',sep='')),lty=c(1,2),lwd=c(1,1),col=c('blue','blue'),bty='o',y.intersp=0.8)
  legend('topleft',c(paste('Threshold ',BRP_names[1],sep=''),paste('Target ',BRP_names[2],sep='')),lty=c(1,2),lwd=c(1,1),col=c('black','black'),bty='o',y.intersp=0.8)     
  legend('topright',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)

  ssb=report$Spawning_Stock_Pred
  std=std_est$std$sdr_vSSB
  up=ssb+2*std
  lo=ssb-2*std
  lo[which(lo<0)]=0
    
  ylim_plot=c(0,max(up,sthresh,starget,na.rm=T)*1.3)
  
  plot(year,ssb,type='o',xlab='',ylab='Spawning stock biomass (metric tons)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black') 
  axis(1,at=year_axis,las=3)
  lines(year,sthresh,lty=1,lwd=1,col='black')
  lines(year,starget,lty=2,lwd=1,col='black')
  
  mtext('Year',side=1,outer=T,line=2.5)  
}

PlotStockStatus_ratio=function(fthresh,sthresh,thresh)
{  
  par(mfrow=c(2,1))
  par(mar = c(0.5, 4, 0, 0), oma = c(4, 0.5, 0.5, 0.5))
  par(tcl = -0.25)
  
  ff=report$Fishing_Mortality_Annual
  std=std_est$std$sdr_vAnnualF
  up=ff+2*std
  lo=ff-2*std
  lo[which(lo<0)]=0
  
  ff=ff/fthresh
  up=up/fthresh
  lo=lo/fthresh
       	
  ylim_plot=c(0,max(up,1,na.rm=T)*1.3)
    	
  plot(year,ff,type='o',xlab='',ylab=paste('F/F',thresh,sep=''),main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black')
  abline(h=1,lty=1,lwd=1,col='black')    
  legend('topright',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)

  ssb=report$Spawning_Stock_Pred
  std=std_est$std$sdr_vSSB
  up=ssb+2*std
  lo=ssb-2*std
  lo[which(lo<0)]=0
  
  ssb=ssb/sthresh
  up=up/sthresh
  lo=lo/sthresh
     
  ylim_plot=c(0,max(up,1,na.rm=T)*1.3)
  
  plot(year,ssb,type='o',xlab='',ylab=paste('SSB/SSB',thresh,sep=''),main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black') 
  axis(1,at=year_axis,las=3)
  abline(h=1,lty=1,lwd=1,col='black')
  
  mtext('Year',side=1,outer=T,line=2.5)  
}

PlotDeplete=function(ssbF0,sthresh,starget)
{  
  par(mfrow=c(1,1))
  par(mar = c(0.5, 4, 0, 0), oma = c(4, 0.5, 0.5, 0.5))
  par(tcl = -0.25)
  
  ssb=report$Spawning_Stock_Pred
  std=std_est$std$sdr_vSSB
  up=ssb+2*std
  lo=ssb-2*std
  lo[which(lo<0)]=0
  
  ssb=ssb/ssbF0
  up=up/ssbF0
  lo=lo/ssbF0
    
  ylim_plot=c(0,max(up,sthresh,starget,na.rm=T)*1.3)
  
  plot(year,ssb,type='o',xlab='',ylab='SSB/SSBF0',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')
  lines(year,up,lty=3,lwd=2,col='black')
  lines(year,lo,lty=3,lwd=2,col='black') 
  axis(1,at=year_axis,las=3)
  abline(h=sthresh,lty=1,lwd=1,col='black')
  abline(h=starget,lty=2,lwd=1,col='black')
  
  legend('topleft',c(paste('Threshold ',sthresh,sep=''),paste('Target ',starget,sep='')),lty=c(1,2),lwd=c(1,1),col=c('black','black'),bty='o',y.intersp=0.8)     
  legend('topright',c('Mean','±2SD'),lty=c(1,3),lwd=c(2,2),col=c('black','black'),bty='o',y.intersp=0.8)
 
  mtext('Year',side=1,outer=T,line=2.5)  
}

PlotSPRatio_FSeries=function(AvgFlag)    #specify AvgFlag: 0-use M and growth for reference year; 1-use avg M adn growth matrix
{
    if(AvgFlag==0){
    	BRPseries=as.matrix(read.table('BRP_series_ref.txt',header=T,row.names=NULL,sep='\t'))
    }
    if(AvgFlag==1){
    	BRPseries=as.matrix(read.table('BRP_series_avg.txt',header=T,row.names=NULL,sep='\t'))
    }
    
    par(mfrow=c(1,1))
    par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    
	thresh=ThreshSPR       #report$ThreshSPR
	target=TargetSPR       #report$TargetSPR
	
	fthresh=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'FThreshSpr']     #report$F20SPR
	if(fthresh==-999){
		fthresh=NA
	}
	ftarget=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'FTargetSpr']     #report$F30SPR
	if(ftarget==-999){
		ftarget=NA
	}
	xx=BRPseries[,'F']            #report$F_Series
	yy=BRPseries[,'SPRatio']      #report$SPRatio_Series
	ylim_plot=c(0,max(yy)*1.3)
	
	plot(xx,yy,type='l',xlab='',ylab='',main=NULL,lwd=2,col='black',ylim=ylim_plot,xaxt='n',yaxt='n')
    axis(1,at=seq(0,max(xx),by=1))
    axis(2,at=seq(0,1,by=0.2))
    abline(v=fthresh,col='black',lwd=1,lty=1)
    abline(v=ftarget,col='black',lwd=1,lty=2)
    legend('topright',BRP_names,lty=c(1,2),lwd=c(1,1),col=c('black','black'),bty='o',y.intersp=0.8) 
    
    mtext('Fishing mortality (per year)',side=1,outer=T,line=2.5)
    mtext('Spawning potential ratio',side=2,outer=T,line=2.5)    
}

PlotYPR_FSeries=function(AvgFlag)   #specify AvgFlag: 0-use M and growth for reference year; 1-use avg M adn growth matrix
{
    if(AvgFlag==0){
    	BRPseries=as.matrix(read.table('BRP_series_ref.txt',header=T,row.names=NULL,sep='\t'))
    }
    if(AvgFlag==1){
    	BRPseries=as.matrix(read.table('BRP_series_avg.txt',header=T,row.names=NULL,sep='\t'))
    }

    par(mfrow=c(1,1))
    par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    
	fthresh=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'F01']      #report$F01
	if(fthresh==-999){
		fthresh=NA
	}
	ftarget=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'Fmax']     #report$Fmax
	if(ftarget==-999){
		ftarget=NA
	}
	xx=BRPseries[,'F']           #report$F_Series
	yy=BRPseries[,'YPR']/1000    #report$YPR_Series/1000    #output is yield (metric tons) per 1000 recruits
	ylim_plot=c(0,max(yy)*1.3)
	
	plot(xx,yy,type='l',xlab='',ylab='',main=NULL,lwd=2,col='black',ylim=ylim_plot,xaxt='n',yaxt='n')
    axis(1,at=seq(0,max(xx),by=1))
    axis(2,at=seq(0,ylim_plot[2],by=0.0001))
    abline(v=fthresh,col='black',lwd=1,lty=1)
    abline(v=ftarget,col='black',lwd=1,lty=2)
    legend('topright',c('F01','Fmax'),lty=c(1,2),lwd=c(1,1),col=c('black','black'),bty='o',y.intersp=0.8) 
    
    mtext('Fishing mortality (per year)',side=1,outer=T,line=2.5)
    mtext('Yield (metric tons) per recruit',side=2,outer=T,line=2.5)   
}

PlotYield_FSeries=function(AvgFlag)   #specify AvgFlag: 0-use M and growth for reference year; 1-use avg M adn growth matrix
{
    if(AvgFlag==0){
    	BRPseries=as.matrix(read.table('BRP_series_ref.txt',header=T,row.names=NULL,sep='\t'))
    }
    if(AvgFlag==1){
    	BRPseries=as.matrix(read.table('BRP_series_avg.txt',header=T,row.names=NULL,sep='\t'))
    }

    par(mfrow=c(1,1))
    par(mar = c(1, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    par(tcl = -0.25)
    
	fthresh=BRPsummary[which(BRPsummary[,'year']==BRPRefYear_Ini),'Fmsy']   #report$Fmsy
	if(fthresh==-999){
		fthresh=NA
	}
	xx=BRPseries[,'F']          #report$F_Series
	yy=BRPseries[,'Yield']      #report$TotalYield_Series
	ylim_plot=c(0,max(yy)*1.3)
	
	plot(xx,yy,type='l',xlab='',ylab='',main=NULL,lwd=2,col='black',ylim=ylim_plot,xaxt='n',yaxt='n')
    axis(1,at=seq(0,max(xx),by=1))
    axis(2,at=seq(0,ylim_plot[2],by=500))
    abline(v=fthresh,col='black',lwd=1,lty=2)
    legend('topright',c('Fmsy'),lty=c(2),lwd=c(2),col=c('black'),bty='o',y.intersp=0.8) 
    
    mtext('Fishing mortality (per year)',side=1,outer=T,line=2.5)
    mtext('Total yield (metric tons)',side=2,outer=T,line=2.5)    
}

#######################################################################

#######################################################################
####Growth
#######################################################################
####plot estmated Linf and K
PlotVB_Linf=function(Nblock_vb,flag)  #specify number of blocks and flag; flag: 0-input Paras and not plot; >0 estimate paras and plot
{
  if(flag>0)
  {  	
	par(mfrow=c(2,1))
    par(mar = c(1, 4, 0, 0), oma = c(4, 0.5, 0.5, 0.5))
    par(tcl = -0.25)
           
    yy=report$Growth_Block_Linf
    ylim_plot=c(0,max(yy)*1.3)
    
    if(flag>1)
    {
    	yy_dev=report$Growth_Linf_Dev_Logscale
    	ylim_plot_dev=c(min(yy_dev)*1.3,max(yy_dev)*1.3) 
    }else{
    	yy_dev=rep(NA,Nblock_vb)
    	ylim_plot_dev=c(-1,1)
    }
    
    plot(1:Nblock_vb,yy,type='o',xlab='',ylab='Linf (mm)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')    
    plot(1:Nblock_vb,yy_dev,type='o',xlab='',ylab='Deviation (log-scale)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot_dev,xaxt='n')
    abline(h=0,col="black",lwd=1)
    axis(1,at=1:Nblock_vb)
    
    mtext('Block',side=1,outer=T,line=2.5)
  }
}

PlotVB_K=function(Nblock_vb,flag)  #specify number of blocks and flag; flag: 0-input Paras and not plot; >0 estimate paras and plot
{
  if(flag>0)
  {  	
	par(mfrow=c(2,1))
    par(mar = c(1, 4, 0, 0), oma = c(4, 0.5, 0.5, 0.5))
    par(tcl = -0.25)
           
    yy=report$Growth_Block_K
    ylim_plot=c(0,max(yy)*1.3)
    
    if(flag>1)
    {
    	yy_dev=report$Growth_K_Dev_Logscale
    	ylim_plot_dev=c(min(yy_dev)*1.3,max(yy_dev)*1.3) 
    }else{
    	yy_dev=rep(NA,Nblock_vb)
    	ylim_plot_dev=c(-1,1)
    }
    
    plot(1:Nblock_vb,yy,type='o',xlab='',ylab='K (per year)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot,xaxt='n')    
    plot(1:Nblock_vb,yy_dev,type='o',xlab='',ylab='Deviation (log-scale)',main=NULL,lwd=2,col='black',pch=16,ylim=ylim_plot_dev,xaxt='n')
    abline(h=0,col="black",lwd=1)
    axis(1,at=1:Nblock_vb)
    
    mtext('Block',side=1,outer=T,line=2.5)
  }
}

###Plot growth of a cohort,specify start year, year of biological parameters to use in projection and how many years to grow
PlotGM=function(year_start,year_bio,YearNum)
{
    par(mfrow=c(1,1)) 
    par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
    par(tcl = -0.25)
  
    year_end=year_start+YearNum-1
    
    if(year_end<=report$End_Year){
        year_id=c(year_start:year_end)-report$Begin_Year+1
    }else{
    	year_id=c(year_start:year_end,rep(year_end,year_end-report$End_Year))-report$Begin_Year+1
    }
    
    year_bio_id=year_bio-report$Begin_Year+1
        
    BRPNASData=matrix(NA,ncol=Nsizebin,nrow=YearNum)
    Recruitment=1000
    RecPrjRatio=report$Recruit_Sizebin_Ratio
    RecPrjSeaRatio=report$Recruit_Season_Ratio

    BRPNAS=rep(0,Nsizebin)
    for (i in 1:YearNum)
    {           
      for (j in 1:Nseason)
      {
        
        BRPS=exp(-report$Natural_Mortality[(j-1)*Nyear+year_bio_id,])  #natrual mortality: (season, year) by size
        
        iL2=GrowthMatrixBlockFlag[year_bio_id,j+1]
	    iL1=(iL2-1)*Nseason+j
		itemp=(iL1-1)*Nsizebin
        GM=report$Growth[itemp+1:Nsizebin,]
                
        if(i==1){
        	BRPNAS=BRPNAS+Recruitment*RecPrjRatio*RecPrjSeaRatio[j]
        }
        
        if(j==1){
        	BRPNASData[i,]=BRPNAS  #record abundance at the beginning of the year,i.e., the begning of the 1st season
        }              
        
        NASTemp=BRPNAS*BRPS
        BRPNAS=NASTemp%*%GM
      }
    }
    
    iYear=seq(1,YearNum,1)
    RRRData=NULL
    RMainData=c()
    
    for( i in 1:YearNum)
    {     
      RRData=BRPNASData[i,]
      if(sum(RRData)>0){
      	RRData=RRData/sum(RRData)
      }      
      RMainData[i]=sum(RRData*sizebin) #number-weighted mid-size
      for(j in 1:Nsizebin)
      {
        SizeNum=RRData[j]*1000  
        if(SizeNum<1)
        {
        	next;
        }
          
        XXX1=rep(i,SizeNum)
        XXX2=rep(j,SizeNum)
        RRRDataTemp=cbind(iYear[XXX1],sizebin[XXX2])
        RRRData=rbind(RRRData,RRRDataTemp)
      }
    }
    
    boxplot(RRRData[,2]~factor(RRRData[,1]),xlab="",ylab="",main=NULL,col=NULL)
    #stripchart(RRRData[,2]~factor(RRRData[,1]),method='jitter',add=T,vertical=T,pch=1)
    if(RMainData[1]==0){    #this is the case when RecPrjSeaRatio[1]=0 when N at the begining of the initial year is 0, skip plotting for year 1
    	lines(1:(YearNum-1),RMainData[-1],type="l",col="black",lwd=2)
    }else{                 #this is the case when RecPrjSeaRatio[1]=1
    	lines(1:YearNum,RMainData,type="l",col="black",lwd=2)
    }
    
  
    mtext('Year',side=1,outer=T,line=2.5)
    mtext('Mid-length (mm)',side=2,outer=T,line=2.5)
    mtext('Example of the growth of a cohort',side=3,outer=T,line=2)
}

######Plot growth of a cohort using a specified growth block
###speicfy block id, begin year for each block, number of years to grow
###use begin year biological parameters M in projection

PlotGM_blk=function(block_id,block_year,YearNum)
{
    par(mfrow=c(1,1)) 
    par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
    par(tcl = -0.25)
  
    year_id=block_year[block_id]-report$Begin_Year+1
  
    itemp=Nseason*Nsizebin*(year_id-1)
    GM=report$Growth[(itemp+1):(itemp+Nsizebin*Nseason),]
    
    BRPNASData=matrix(NA,ncol=Nsizebin,nrow=YearNum)
    Recruitment=1000
    RecPrjRatio=report$Recruit_Sizebin_Ratio
    RecPrjSeaRatio=report$Recruit_Season_Ratio
    BRPS=exp(-report$Natural_Mortality[c((1-1)*Nyear+year_id,(2-1)*Nyear+year_id),])  #natrual mortality: (season, year) by size

    BRPNAS=rep(0,Nsizebin)
    for (i in 1:YearNum)
    {
      for (j in 1:Nseason)
      {
        if(i==1){
        	BRPNAS=BRPNAS+Recruitment*RecPrjRatio*RecPrjSeaRatio[j]
        }
        
        if(j==1){
        	BRPNASData[i,]=BRPNAS  #record abundance at the beginning of the year,i.e., the begning of the 1st season
        }
              
        iL=Nsizebin*(j-1)
        NASTemp=BRPNAS*BRPS[j,]
        BRPNAS=NASTemp%*%GM[(iL+1):(iL+Nsizebin),]
      }
    }
    
    iYear=seq(1,YearNum,1)
    RRRData=NULL
    RMainData=c()
    
    for( i in 1:YearNum)
    {     
      RRData=BRPNASData[i,]
      RRData=RRData/sum(RRData)
      RMainData[i]=sum(RRData*sizebin) #number-weighted mid-size
      for(j in 1:Nsizebin)
      {
        SizeNum=RRData[j]*1000  
        if(SizeNum<1)
        {
        	next;
        }
          
        XXX1=rep(i,SizeNum)
        XXX2=rep(j,SizeNum)
        RRRDataTemp=cbind(iYear[XXX1],sizebin[XXX2])
        RRRData=rbind(RRRData,RRRDataTemp)
      }
    }
    
    boxplot(RRRData[,2]~factor(RRRData[,1]),xlab="",ylab="",main=NULL,col=NULL)
    #stripchart(RRRData[,2]~factor(RRRData[,1]),method='jitter',add=T,vertical=T,pch=1)
    lines(1:YearNum,RMainData,type="l",col="black",lwd=2)
  
    mtext('Year',side=1,outer=T,line=2.5)
    mtext('Mid-length (mm)',side=2,outer=T,line=2.5)
    mtext(paste('Example of the growth of a cohort with Block',block_id,sep=''),side=3,outer=T,line=2)
}

#####VB function
PlotVB=function(Nblock_vb,flag,t0)  #specify num of blocks,flag and t0; here use t0 estimated externally; flag: 0-input GM,1-VB with blk,2-VB random walk, 3-avg with random deviation
{
  if(flag>0){
    n.col=5 #num of coln per plot
    n.row=6 #num of row per plot
    n.pannel=n.col*n.row #num of pannels per plot
    n.plot=ceiling(Nblock_vb/n.pannel)
    y.axis=rep(0,n.row)
    x.axis=rep(0,n.col)
    
    for(i in 1:n.row){
    	y.axis[i]=(i-1)*n.col+1
    }
    
    ylim_plot=c(0,max(size_axis))
    
    for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		# layout(matrix(1:n.fig,nrow=n.row,ncol=n.col,byrow=T))
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
		# par(mgp = c(2, 0.6, 0))
	   }
	
	  if(hh==n.plot){
		n.fig=Nblock_vb-(n.plot-1)*n.pannel
		# layout(matrix(c(1:n.fig,rep(0,(n.pannel-n.fig))),nrow=n.row,ncol=n.col,byrow=T))
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
		# par(mgp = c(2, 0.6, 0))
	  }
	  
	  x.axis=seq(n.fig-n.col+1,n.fig,1)
	  
	  for(pp in 1:n.fig){
	  	 block_id=(hh-1)*n.pannel+pp
	  	 	  	 
         ##linf=par_est[[paste('PbnvGrowthPara_',(block_id-1)*5+1,sep='')]]
         ##k=par_est[[paste('PbnvGrowthPara_',(block_id-1)*5+2,sep='')]]
                  
         linf=report$Growth_Block_Linf[block_id]
         k=report$Growth_Block_K[block_id]
                 
         age=0:10
    
         y=linf*(1-exp(-k*(age-t0)))
         
	  	 plot(age,y,type='l',xlab="",ylab="",main=NULL,col='black',lwd=1,ylim=ylim_plot,axes=F)
	  	 mtext(paste('Block ',block_id,sep=''),side=3,line=-1,adj=0.1,cex=0.6)
	  	 
	  	 if(pp %in% y.axis){
	  	 	axis(2,col='gray',at=size_axis)
	  	 }      
         if(pp %in% x.axis){
         	axis(1,col='gray',at=age)
         }
         box(col='gray')
	  }
	
	  mtext('Age',1,outer=T,line=2.5)
	  mtext('Mean length (mm)',2,outer=T,line=2.5)
	  mtext('VB Growth',3,outer=T,line=2)
    }
  }
}

PlotVB_deltaL=function(Nblock_vb,flag)  #specify num of blocks, flag as above
{
  if(flag>0){
    n.col=5 #num of coln per plot
    n.row=6 #num of row per plot
    n.pannel=n.col*n.row #num of pannels per plot
    n.plot=ceiling(Nblock_vb/n.pannel)
    y.axis=rep(0,n.row)
    x.axis=rep(0,n.col)
    
    for(i in 1:n.row){
    	y.axis[i]=(i-1)*n.col+1
    }    
    
    ylim_plot=c(0,max(size_axis)/4.5)
    
    name_temp=rep(NA,Nseason)
    for(j in 1:Nseason){
    	name_temp[j]=paste('Season ',j,sep='')
    }
    
    for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
	   }
	
	  if(hh==n.plot){
		n.fig=Nblock_vb-(n.plot-1)*n.pannel
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
	  }
	  
	  x.axis=seq(n.fig-n.col+1,n.fig,1)
	  
	  for(pp in 1:n.fig){
	  	 block_id=(hh-1)*n.pannel+pp
	  	 
         ##linf=par_est[[paste('PbnvGrowthPara_',(block_id-1)*5+1,sep='')]]
         ##k=par_est[[paste('PbnvGrowthPara_',(block_id-1)*5+2,sep='')]]
         
         linf=report$Growth_Block_Linf[block_id]
         k=report$Growth_Block_K[block_id]
                  
         alpha=report$Growth_Proportion
         
	  	 # plot(age,y,type='l',xlab="",ylab="",main=NULL,col='black',lwd=1,ylim=ylim_plot,axes=F)
	  	 plot(1,type='n',xlab='',ylab='',main=NULL,xlim=c(min(size_axis),max(size_axis)),ylim=ylim_plot,axes=F)	  	 
	  	 mtext(paste('Block ',block_id,sep=''),side=3,line=-1,adj=0.1,cex=0.6)
         
         for(i in 1:Nseason){
    	    y=(linf-sizebin)*(1-exp(-alpha[i]*k))
    	    lines(sizebin,y,type='l',col=i,lty=i,lwd=1)
         }
	  	 	  	 
	  	 if(pp %in% y.axis){
	  	 	axis(2,col='gray',at=size_axis)
	  	 }      
         if(pp %in% x.axis){
         	axis(1,col='gray',at=size_axis)
         }
         box(col='gray')
        
         if(pp==1){
        	legend('topright',name_temp,lty=1:Nseason,lwd=rep(1,Nseason),col=1:Nseason,bty='o',y.intersp=0.8,cex=0.8,box.col='grey')
         }        
	  }
	
	  mtext('Mid-length (mm)',1,outer=T,line=2.5)
	  mtext('DeltaL(mm)',2,outer=T,line=2.5)
	  mtext('VB Growth DeltaL',3,outer=T,line=2)
    }
  }
}

######################################################################


#######################################################################
#weight-length
#######################################################################
PlotWL=function()
{      
    n.col=5 #num of coln per plot
    n.row=6 #num of row per plot
    n.pannel=n.col*n.row #num of pannels per plot
    n.plot=ceiling(Nyear/n.pannel)
    y.axis=rep(0,n.row)
    x.axis=rep(0,n.col)
    
    for(i in 1:n.row){
    	y.axis[i]=(i-1)*n.col+1
    }
    
    ylim_plot=c(0,max(report$Weight_Length_g[1,])*1.3)
    
    for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
	   }
	
	  if(hh==n.plot){
		n.fig=Nyear-(n.plot-1)*n.pannel
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
	  }

      x.axis=seq(n.fig-n.col+1,n.fig,1)
      
	  for(pp in 1:n.fig){
	  	 year_id=(hh-1)*n.pannel+pp
	  	 
	  	 y=report$Weight_Length_g[year_id,]
	  	 # ylim_plot=c(0,max(y)*1.3)
         
         plot(sizebin,y,type='l',xlab="",ylab="",main=NULL,col='black',lwd=1,ylim=ylim_plot,axes=F)	  	 
	  	 mtext(year[year_id],side=3,line=-1,adj=0.1,cex=0.6)
	  	 
	  	 if(pp %in% y.axis){
	  	 	axis(2,col='gray',at=seq(ylim_plot[1],ylim_plot[2],by=500))
	  	 }      
         if(pp %in% x.axis){
         	axis(1,col='gray',at=size_axis)
         }
         box(col='gray')
	  }
	
	  mtext('Mid-length (mm)',1,outer=T,line=2.5)
	  mtext('Mean weight (g)',2,outer=T,line=2.5)
	  mtext('Weight at size',3,outer=T,line=2)
  }  
}
######################################################################

#######################################################################
#maturity-length
#######################################################################
PlotML=function()
{      
    n.col=5 #num of coln per plot
    n.row=6 #num of row per plot
    n.pannel=n.col*n.row #num of pannels per plot
    n.plot=ceiling(Nyear/n.pannel)
    y.axis=rep(0,n.row)
    x.axis=rep(0,n.col)
    
    for(i in 1:n.row){
    	y.axis[i]=(i-1)*n.col+1
    }
    
    ylim_plot=c(0,1)
    
    for(hh in 1:n.plot){    
      if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))	
		par(tcl = -0.25)
	   }
	
	  if(hh==n.plot){
		n.fig=Nyear-(n.plot-1)*n.pannel
		par(mfrow=c(n.row,n.col))
		par(mar=c(0,0,0,0),oma=c(4,4,4,0.5))
		par(tcl = -0.25)
	  }
	  
      x.axis=seq(n.fig-n.col+1,n.fig,1)
      
	  for(pp in 1:n.fig){
	  	 year_id=(hh-1)*n.pannel+pp
	  	 
	  	 y=report$Maturity_Length[year_id,]
         
         plot(sizebin,y,type='l',xlab="",ylab="",main=NULL,col='black',lwd=1,ylim=ylim_plot,axes=F)	  	 
	  	 mtext(year[year_id],side=3,line=-1,adj=0.1,cex=0.6)
	  	 
	  	 if(pp %in% y.axis){
	  	 	axis(2,col='gray',at=seq(ylim_plot[1],ylim_plot[2],by=0.2))
	  	 }      
         if(pp %in% x.axis){
         	axis(1,col='gray',at=size_axis)
         }
         box(col='gray')
	  }
	
	  mtext('Mid-length (mm)',1,outer=T,line=2.5)
	  mtext('Maturity (proportion)',2,outer=T,line=2.5)
	  mtext('Maturity at size',3,outer=T,line=2)
  }  
}
######################################################################

