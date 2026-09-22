source("reptoRlist.r");source("01_ReadADMB.r")
filename="NSLSAP01"
report<-read.admb(filename);par_est <- std_read("nslsap01.std")

# Geometric mean function for calculating terminal recruitment and fishing mortality
geomean <- function(x) {
  ## compute geometric mean of a time series x
  nx <- length(x)
  y <- prod(x)^(1 / nx)
  return(y)
}
# Call functions from Cao et al.'s old code and modified version of Yan Li's BRP code
source("PlotFuncs.r")
source('01_BRP.R')


# Plot observed vs. predicted catch
windows(10,10)
par(mfrow=c(1,2))
plot(years, report$Catch_Obs[1,], type = "b", lwd = 1.5, col = "red", ylab = "Catch (1000 mt)", xlab = "Years", main = "Season 1")
lines(years, report$Catch_Pred[1,], lwd = 1.5, col = "black")
plot(years, report$Catch_Obs[2,], type = "b", lwd = 1.5, col = "red", ylab = "Catch (1000 mt)", xlab = "Years", main = "Season 2")
lines(years, report$Catch_Pred[2,], lwd = 1.5, col = "black")
legend("topright", cex = 0.75, c("Observed", "Predicted"), col = c("red", "black"), lty = 1)
# Plot observed vs. predicted indices of relative abundance
windows(10,10)
par(mfrow=c(1,2))
plot(years, report$Survey_Index_Obs[1,], type = "b", lwd = 1.5, col = "red", ylab = "Survey index", xlab = "Years", main = "Season 1")
lines(years, report$Survey_Index_Pred[1,], lwd = 1.5, col = "black")
legend("topright", cex = 0.75, c("Observed", "Predicted"), col = c("red", "black"), lty = 1)
plot(years, report$Survey_Index_Obs[2,], type = "b", lwd = 1.5, col = "red", ylab = "Survey index", xlab = "Years", main = "Season 2")
lines(years, report$Survey_Index_Pred[2,], lwd = 1.5, col = "black")

# Plot starting selectivity vs. predicted selectivity by fleet and index
windows(10,10)
par(mar=c(4,4,4,4))
par(mfrow=c(1,3))
plot(size_bm, sels, type ="b", col = "red", lwd = 1.5, ylab = "Fleet selectivity", xlab = "Size (mm)", ylim=c(0,1))
lines(size_bm, report$Fleet_selectivity[1,], col = "black", lwd = 1.5)
legend("topleft", cex = 0.8, c("Initial", "Estimated"), col=c("red","black"), lty = 1)
plot(size_bm, survey_sels, type = "b",col = "red", ylab = "Survey selectivity (season 1)", xlab = "Size (mm)",ylim=c(0,1))
lines(size_bm, report$Index_selectivity[1,], col = "black", lwd = 1.5)
plot(size_bm, survey_sels_s2, type = "b", col = "red", ylab = "Survey selectivity (season 2)", xlab = "Size (mm)",ylim=c(0,1))
lines(size_bm, report$Index_selectivity[(n_t+1),], col = "black", lwd = 1.5)

# Some calculated quantities
N.est <- rowSums(report$Abundance_at_Size)
N.est.sum <- N.est[1:n_t]+N.est[(n_t+1):(n_t*2)]
F_annual <- c((report$Fishing_Mortality[1,1:42] + report$Fishing_Mortality[2,1:42]),geomean(report$Fishing_Mortality[1,41:43])+(geomean(report$Fishing_Mortality[2,41:43])))

# Vital rates plot: recruitment, abundance, SSB, F (season one), F (season two), F (annual)
windows(10,10)
par(mar=c(4,4,4,4))
par(mfrow=c(2,3))
plot(years,c(report$recruitment_Pred[1:42]*10,geomean(report$recruitment_Pred[41:43]*10)), type = "o",  col = "black",lwd = 1.5, ylab = "Recruits (million crabs)", xlab = "Years",main = "Recruitment",ylim = c(0, max(report$recruitment_Pred*10)))
plot(years,c(N.est.sum[1:42]*10,geomean(N.est.sum[41:43]*10)), type = "o",  col = "black", lwd = 1.5, xlab = "Years", ylab = "Abundance (million crabs)", main = "Abundance", ylim = c(0, max(N.est.sum*10)))
plot(years,report$Spawning_stock_Biomass, type = "o",  col = "black", lwd = 1.5, xlab = "Years", ylab = "SSB (1000 mt)", main = "Spawning Stock Biomass", ylim = c(0, max(report$Spawning_stock_Biomass)))
plot(years,report$Fishing_Mortality[1,], type = "o", lwd = 1.5, col = "black", ylab = "F (season 1)", xlab = "Years", main = "Season 1 Fishing Mortality") #, ylim = c(0,1.0))
plot(years,report$Fishing_Mortality[2,], type = "o", lwd = 1.5, col = "black", ylab = "F (season 2)", xlab = "Years", main = "Season 2 Fishing Mortality") #, ylim = c(0,1.0))
plot(years, F_annual, type = "o", lwd = 1.5, col = "black", ylab = "F (annual)", xlab = "Years", main = "Total annual fishing mortality") #, ylim = c(0,1.0))

# Get stock status
ssb.f0 <- EquiSSBFun(spr.f0,1,1) # SSBF0
ssb.fterm <- EquiSSBFun(spr.term,1,1) #SSBFterminal year
spr.fthresh <- SPRFunc(BRPsummary[1,'FTargetSpr'],n_t,0,0) #SPR F40
ssb.fthresh <- EquiSSBFun(spr.fthresh,1,1) #SSB F40 
spr.ftarg <- SPRFunc((BRPsummary[1,'FTargetSpr']*0.75),n_t,0,0) # 75% SPR F40
ssb.ftarg <- EquiSSBFun(spr.ftarg,1,1)  # SSB 75% F40
spr.f20 <- SPRFunc((BRPsummary[1,'FThreshSpr']),n_t,0,0) #SPR F20
ssb.f20 <- EquiSSBFun(spr.f20,1,1) # SSBF20, almost exactly the same as 50% of SSBF40

# Plot stock status
windows(10,10)
par(mfrow=c(1,2))
plot(years, F_annual, type = "o", lwd = 1.5, col = "black", ylab = "Total F (annual)", xlab = "Years", ylim=c(0,1.5))
lines(years, BRPsummary[,'FTargetSpr'], col = "blue")
lines(years, BRPsummary[,'FTargetSpr']*0.75, col = "green")
# lines(years, BRPsummary[,'FThreshSpr'], col = "red")
legend("topleft", cex = 0.6, c("F40", "75%F40"), col=c("blue","green"), lty = 1)
plot(years, report$Spawning_stock_Biomass, type = "o", lwd = 1.5, col = "black", ylab = "Spawning stock biomass (1000 mt)", xlab = "Years",, ylim=c(0,max(report$Spawning_stock_Biomass)))
lines(years, rep(ssb.fthresh,n_t), col = "blue")
lines(years, rep(ssb.ftarg,n_t), col = "green")
# lines(years, rep(ssb.f20,n_t), col = "red")
lines(years, rep((0.5*ssb.fthresh),n_t), col = "black")
# lines(years, rep((0.25*ssb.f0), n_t), col = "orange") # 25% SSBF0
# lines(years, rep((0.75*ssb.fthresh), n_t), col = "pink") # 75% SSBF40
legend("topright", cex = 0.6, c("SSBF40", "SSB 75%F40", "50% SSBF40"), col=c("blue","green", "black" ), lty = 1)


# calculate dynamic B0 
# rec_pred <- c(report$recruitment_Pred[1:42],geomean(report$recruitment_Pred[41:43]))
# dSSB75F40 <- array(NA, dim = c(1,length(rec_pred)))
# dSSBF40 <- array(NA, dim = c(1,length(rec_pred)))
# dSSBF0 <- array(NA, dim = c(1,length(rec_pred)))
# # modify function to loop through recruitment vector
# DynSSBFun=function(dvSPR,iRS=1) # default is mean recruitment
# {
#   if(iRS==1)    #annual average R independent of spawner
#   {
#     alpha = rec_pred[i] # produce median equilibrium SSB e.g, Chesbay 2026 assessment
#     dvR = alpha*dvSPR
#   }
#   dvR
# }
# # loop through function
# for(i in 1:length(rec_pred)){
#   dSSBF40[i] <- DynSSBFun(spr.fthresh, 1)
#   dSSB75F40[i] <- DynSSBFun(spr.ftarg, 1)
#   dSSBF0[i] <- DynSSBFun(spr.f0, 1)
# }
# # Dynamic B0 SSB plot
# windows(10,10)
# par(mfrow=c(1,1))
# plot(years, report$Spawning_stock_Biomass, type = "o", lwd = 1.5, col = "black", ylab = "Spawning stock biomass (1000 mt)", xlab = "Years",, ylim=c(0,max(dSSBF0)))
# lines(years, dSSBF40, col = "blue")
# lines(years, dSSB75F40, col = "green")
# lines(years, 0.5*dSSBF40, col = "black")
# lines(years, dSSBF0, col = "purple")
# lines(years, rep((0.5*ssb.fthresh),n_t), col = "red")
# legend("topright", cex = 0.6, c("dSSBF40", "dSSB 75%F40", "50% dSSBF40", "dSSBF0", "50% sSSBF40"), col=c("blue","green", "black","purple"),lty = 1)

# Calculate CIs for log-scale recruitment and Fdevs
# extract SEs from std object
rec.se <- par_est$std$PbdvLogRecruitDevs
rec.upper.CI <- (par_est$est$`PbnvLogRSPara[1]`+par_est$est$PbdvLogRecruitDevs)+(1.96*rec.se)
rec.lower.CI <- (par_est$est$`PbnvLogRSPara[1]`+par_est$est$PbdvLogRecruitDevs)-(1.96*rec.se)
f.se <- as.numeric(par_est$std[1:86])
f.est <- as.numeric(par_est$est[1:86])
f.upper.CI <- f.est+1.96*f.se
f.lower.CI <- f.est-1.96*f.se

# Plot log-scale recruitment and F devs with 95% CIs
windows(10,10)
par(mar = c(4,4,4,4))
par(mfrow=c(1,2))
plot(seq(1982,2024,1),(par_est$est$`PbnvLogRSPara[1]`+par_est$est$PbdvLogRecruitDevs), lwd = 1.5, type = "l", ylim = c(min(rec.lower.CI),3), ylab = "Log-scale recruitment", xlab = "Years")
lines(seq(1982,2024,1),rec.upper.CI, lty = "dashed")
lines(seq(1982,2024,1),rec.lower.CI, lty = "dashed")
plot(seq(1982,2024.5,0.5), f.est, type = "l", lwd = 1.5, ylab = "log-scale F deviations", xlab = "Year-season", ylim=c(min(f.lower.CI),max(f.upper.CI)))
lines(seq(1982,2024.5,0.5), f.upper.CI, lty = "dashed")
lines(seq(1982,2024.5,0.5), f.lower.CI, lty = "dashed")

# Plot survey index size composition
windows(10,10)
PlotSLC() # if using the windows() function, will only show second season's comps
# Plot abundance at size (from Cao et al.'s original code)
windows(10,10)
PlotAbun()
# Modify Cao et al.'s fishery-dependent catch composition function to work with current model
PlotCC<-function()
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

    Scomp_Obs = report$Catch_Comp_Obs[ind[s]:(ind[s+1]-1),]
    Scomp_Pred = report$Catch_Comp_Pred[down_pred:up_pred,]

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
    mtext(paste('Catch Season-',s,';black line = Pred',sep=''),side=3,outer=T,line=2.2)
  }
}
windows(10,10)
PlotCC()

# Modify Cao et al.'s stock recruitment plot to give correct scale
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

  plot(Year,R,ylab="Recruits (1000 mt)",xlab="Year",type='o',
       col='black',cex=1.5,pch=20,lwd=1.5,lty='dashed',cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  plot(Year,Rdev,xlab="Year",ylab="Log recruitment deviation",type='b',
       col='black',cex=1.5,pch=20,lty='dotted',cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h=0,col="gray",lwd=1.5)
  plot(SSB,R,xlab="Spawning biomass (1000 mt)",ylab="Recruits (1000 mt)",
       col='black',cex=1.5,pch=20,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  par(mfrow=c(1,1))
}
# Stock recruitment plot from Cao et al.'s original code
windows(10,10)
PlotR()
# Plot to compare NAO index with log-scale recruitment deviations
#windows(10,10)
#par(mar=c(3,3,3,3))
#par(mfrow=c(1,2))
#plot(years,env.data[,2], type = "o", lwd = 1.5, ylab = "NAO index", xlab = "Years")
#lines(years,rep(0,length(years)), col = "black")
#plot(years,par_est$est$PbdvLogRecruitDevs, type = "o", lwd = 1.5, ylab = "log-scale recruitment deviations", xlab = "Years")
#lines(years,rep(0,length(years)), col = "black")


# Calculate index residuals
index_resids <- array(NA, dim = c(length(years),2))
for(i in 1:(length(years))){
  index_resids[i,1] <- report$Survey_Index_Obs[1,i]-report$Survey_Index_Pred[1,i]
  index_resids[i,2] <- report$Survey_Index_Obs[2,i]-report$Survey_Index_Pred[2,i]
}

# Runs Tests
library(randtests)
runs.test(index_resids[,1])
runs.test(index_resids[,2])
library(DescTools)
RunsTest(index_resids[,1])
RunsTest(index_resids[,2])

# index_resids[29,2] <- 0.0 # test to see if failure s driven by large outliers
# index_resids[37,2] <- 0.0
#
# runs.test(index_resids[,2])

sigma_res1 <- sum(index_resids[,1])/43
sd_res1 <- sqrt(sigma_res1)
sigma_res2 <- sum(index_resids[,2])/43
sd_res2 <- sqrt(sigma_res2)

index_resids[,1] <- index_resids[,1]/sd_res1
index_resids[,2] <- index_resids[,2]/sd_res2

# Index residual plots
windows(10,10)
par(mfrow=c(1,2))
plot(years, index_resids[,1], type = "o", lwd = 1.5, ylab = "Std residuals", xlab = "Years", main = "Season 1")
lines(years, rep(0,length(years)))
plot(years, index_resids[,2],  type = "o", lwd = 1.5, ylab = "Std residuals", xlab = "Years", main = "Season 2")
lines(years, rep(0,length(years)))

# The Graveyard

# Annual plots if using annual model
# years = seq(1982,2024,1)
# windows(10,10)
# par(mfrow=c(1,2))
# plot(years, report$Catch_Obs)
# lines(years, report$Catch_Pred)
# plot(years, report$Survey_Index_Obs)
# lines(years, report$Survey_Index_Pred)
# 
# windows(10,10)
# par(mfrow=c(1,2))
# plot(size_bm, sels)
# lines(size_bm, report$Fleet_selectivity[1,])
# plot(size_bm, survey_sels)
# lines(size_bm, report$Index_selectivity[1,])
# 
# #
# windows(10,10)
# par(mfrow=c(1,2))
# plot(years, c(report$Fishing_Mortality[1:42],geomean(report$Fishing_Mortality[40:42])), lwd = 1.5, type = "o", ylab = "F", xlab = "Years", main = "Annual Fishing Mortality", ylim = c(0,2.0) ) # if annual
# plot(years, report$Spawning_stock_Biomass, type = "o", lwd = 1.5, col = "black", ylab = "SSB (1000 mt)", xlab = "Years", main = "Spawning Stock Biomass", ylim=c(0,40))

# Other stuff
# View(par_est$est)
# write.csv(par_est, filename = "base_parestest.csv")
# 
# Fdevs_est <- as.numeric(par_est$est[1:86])
# Fdevs_sd <- as.numeric(par_est$std[1:86])
# Fdevs <- cbind(Fdevs_est,Fdevs_sd)
# write.csv(Fdevs, "fdevs.csv")
# Rdevs_est <- par_est$est$PbdvLogRecruitDevs
# Rdevs_sd <- par_est$std$PbdvLogRecruitDevs
# Rdevs <- cbind(Rdevs_est, Rdevs_sd)
# write.csv(Rdevs, "rdevs.csv")
# Initcomp <- par_est$est$MvNASPJD_1
# Initcomp_sd <- par_est$std$MvNASPJD_1
# initcomp <- cbind(Initcomp, Initcomp_sd)
# write.csv(initcomp, "initcomps.csv")
