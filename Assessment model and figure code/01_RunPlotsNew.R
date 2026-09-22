source("reptoRlist.r");source("r_ReadADMB.r")
filename="sssTrout"
report<-read.admb(filename) #;par_est <- std_read("filename.std")
source("r_Output.r") # will output most basic diagnostic plots

# Some quick and dirty diagnostics
# windows(10,10)
# par(mfrow=c(1,2))
# plot(years, report$Catch_Obs[1,], type = "b", col = "red", ylab = "Catch (1000 mt)", xlab = "Years", main = "Season 1")
# lines(years, report$Catch_Pred[1,], col = "black")
# plot(years, report$Catch_Obs[2,], type = "b", col = "red", ylab = "Catch (1000 mt)", xlab = "Years", main = "Season 2")
# lines(years, report$Catch_Pred[2,], col = "black")
# windows(10,10)
# par(mfrow=c(1,2))
# plot(years, report$Survey_Index_Obs[1,], type = "b", col = "red", ylab = "Survey index", xlab = "Years", main = "Season 1")
# lines(years, report$Survey_Index_Pred[1,], col = "black")
# plot(years, report$Survey_Index_Obs[2,], type = "b", col = "red", ylab = "Survey index", xlab = "Years", main = "Season 2")
# lines(years, report$Survey_Index_Pred[2,], col = "black")

# windows(10,10)

# windows(10,10)
par(mfrow=c(1,4))
# plot(report$Abundance_at_Size[,1], type = "l", col = "black", ylab = "Total abundance", xlab = "Years")
plot(years,report$Recruitment_Pred, type = "o", lwd = 1.5, col = "black", ylab = "Recruitment", xlab = "Years", ylim = c(0, max(report$Recruitment_Pred)), main = "Estimated recruitment (1000 mt)")
plot(years,report$Fishing_Mortality_Season[1,], type = "o", lwd = 1.5, col = "black", ylab = "F (season 1)", xlab = "Years", main = "Season 1 Fishing Mortality")
plot(years,report$Fishing_Mortality_Season[2,], type = "o", lwd = 1.5, col = "black", ylab = "F (season 2)", xlab = "Years", main = "Season 2 Fishing Mortality")
plot(years,report$Spawning_Stock_Pred, type = "o", lwd = 1.5, col = "black", ylab = "SSB", xlab = "Years", ylim = c(0, max(report$Spawning_Stock_Pred)), main = "Spawning stock biomass (1000 mt)")
# pdf("plots.pdf",width=8,height=6)

# # Try some of Yan's plots
# 
# PlotC(1)
# PlotSR() # so far, only this one works, otherwise you need a whole host of other variables declared for them to work. FUCK. 

# windows(10,10)
# par(mfrow=c(3,1))
# plot(years,report$Spawning_stock_Biomass, type = "o", lwd = 1.5, col = "black", ylab = "SSB", xlab = "Years", ylim = c(0, max(report$Spawning_stock_Biomass)), main = "Spawning stock biomass (1000 mt)")
# plot(years,report$Fishing_Mortality[1,], type = "o", lwd = 1.5, col = "black", ylab = "F (season 1)", xlab = "Years", main = "Season 1 Fishing Mortality")
# plot(years,report$Fishing_Mortality[2,], type = "o", lwd = 1.5, col = "black", ylab = "F (season 2)", xlab = "Years", main = "Season 2 Fishing Mortality")
# 

# PlotSLC()
# Plot_GM_tb(1,2)
# # windows(10,10)
# # PlotAbun()
# 
# # PlotWL() # plot length-weight relationships
# # PlotML() # plot maturity
#  # visualize growth matrices (user defines the number of year-specific time-blocks for growth and how many year to grow)
# # Plot_sel_f() # plot fleet selectivity
# # Plot_f() # plot fishing mortality
# # PlotTC() # plot total catch
# # Plot_f_c() # plot fishery-specific catch
# # PlotCC(3,1)
# # PlotCCA() # plot aggregated fishery catch composition
# # PlotSelS(2000,1) # plot survey selectivity; specify year and survey
# # PlotSelS(2000,2)
# # #PlotSelS(2000,3)
# # windows(10,10)
#  # plot survey size composition
# # PlotSLCA() # plot aggregated survey size composition
# # PlotSI(3) # plot survey index;n= # of indices
# PlotM(2000)
# windows(10,10)
# PlotR()
# 
# PlotSSB()
# # windows(10,10)
# PlotSpB(1)
# # PlotSexComp(1,0.5)
# # PlotFfit()
# # PlotLfifty()
# # PlotER()
# 
# # PlotStockStatus()
# #PlotRoSSB(1985,2002,2006,4)
# #PlotRo(1984,2012,2016,2)
# #PlotRo_scale(1984,2008,2013,2)
# #Getlike()
# # dev.off()
# PlotCC_seas=function(fishery,season)
# {
#   Year=seq(report$BeginYear,report$EndYear,1)
#   Season=report$Time_step
#   par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
#   par(tcl = -0.25)
#   par(mgp = c(2, 0.6, 0))
#   if (Season>1)
#   {
#     temp=Nyear*Season*(fishery-1)+Nyear*(season-1)
#     temp1=Season*(fishery-1)+season
#     Ccomp_Obs=report$Catch_Comp_Obs[(temp+1):(temp+Nyear),]
#     Ccomp_Pred=report$Catch_Comp_Pred[(temp+1):(temp+Nyear),]
#     CatchESS_Obs=report$Catch_ESSinput[temp1,]
#     CatchESS_Pred=report$Catch_ESSpred[temp1,]
# 
#     if (Nyear<25)
#     {
#       par(mfrow=c(5,5))
#       temp=(Nyear%/%5*(5-1)+Nyear%%5):(Nyear%/%5*(5-1)+5)
#       yaxis=c(1,6,11,16,21)
#       xaxis=c(21,22,23,24,25,temp)
#     }
#     if (Nyear>25)
#     {
#       par(mfrow=c(5,6))
#       yaxis=c(1,7,13,19,25,31)
#       xaxis=c(25,26,27,28,29,30)
#     }
#     
#     if (Nyear>25)
#     {
#       par(mfrow=c(6,6))
#       yaxis=c(1,7,13,19,25,31)
#       xaxis=c(28,29,30,31,32,33)
#     }
#     
#     if (fishery==1&season==1)
#     {
#       Ccomp_Obs = Ccomp_Obs[1:(length(years)),]
#       Ccomp_Pred = Ccomp_Pred[1:length(years),]
#       n_year = 51
#       Year = years
#       xaxis=c(11,12,13,14,15,16)
#     }
#     
#     # if (fishery==1&season==2)
#     # {
#     #   Ccomp_Obs = Ccomp_Obs[52:102,]
#     #   Ccomp_Pred = Ccomp_Pred[52:102,]
#     #   n_year = 51
#     #   Year = seq(1974,2024,1)
#     #   xaxis=c(12,13,14,15,16,17)
#     # }
#     # 
#     # if (fishery==3&season==1)
#     # {
#     #   Ccomp_Obs = Ccomp_Obs[17:Nyear,]
#     #   Ccomp_Pred = Ccomp_Pred[17:Nyear,]
#     #   n_year = 17
#     #   Year = seq(2000,report$EndYear,1)
#     #   xaxis=c(12,13,14,15,16,17)
#     # }
#     
#     par(mfrow=c(3,6))
#     yaxis=c(1,7,13)
#     
#     
#     for (i in 1:Nyear)
#     {
#       plot(sizebins,Ccomp_Obs[i,],type='l',ylim=c(0,0.35),axes=FALSE,lwd=1,col='gray')
#       x=c(sizebins,rev(sizebins))
#       y=c(rep(0,length(sizebins)),rev(Ccomp_Obs[i,]))
#       polygon(x,y,col="gray",border=NA)
#       mtext(paste(Year[i]),side=3,line=-1,adj=0.1,cex=0.6)
#       #legend('topleft',paste(Year[i]),cex=1.1,bty='n')
#       #rp=vector('expression',2)
#       #rp[1]=substitute(expression(N==Value),list(Value=CatchESS_Obs[i]))[2]
#       #rp[2]=substitute(expression(effN==Value),list(Value=CatchESS_Pred[i]))[2]
#       #legend('topright',legend=rp,bty='n',y.intersp=0.8,yjust=1)
#       lines(sizebins,Ccomp_Pred[i,],col='black',lwd=2)
#       #text(x=27,y=0.15,paste("effN=",SurveyESS_Pred[i]),cex=1)
#       if(i %in% yaxis)
#         axis(2,col='gray',at=seq(0,0.2,0.05))
#       if(i %in% xaxis)
#         axis(1,col='gray',at=sizebins)
#       box(col='gray')
#     }
#     mtext('Length (mm)',side=1,outer=T,line=2.2)
#     mtext('Proportion',side=2,outer=T,line=2.2)
#     mtext(paste('Length Comp,','fishery',fishery,'season',season,',black line=Pred'),
#           side=3,outer=T,line=2.2)
#   }else
#   {
#     temp=Nyear*(fishery-1)
#     
#     Ccomp_Obs=report$Catch_Comp_Obs[(temp+1):(temp+Nyear),]
#     Ccomp_Pred=report$Catch_Comp_Pred[(temp+1):(temp+Nyear),]
#     CatchESS_Obs=report$Catch_ESSinput[fishery,]
#     CatchESS_Pred=report$Catch_ESSpred[fishery,]
#     
#     if (Nyear<25)
#     {
#       par(mfrow=c(5,5))
#       temp=(Nyear%/%5*(5-1)+Nyear%%5):(Nyear%/%5*(5-1)+5)
#       yaxis=c(1,6,11,16,21)
#       xaxis=c(21,22,23,24,25,temp)
#     }
#     if (Nyear>25)
#     {
#       par(mfrow=c(5,6))
#       yaxis=c(1,7,13,19,25,31)
#       xaxis=c(25,26,27,28,29,30)
#     }
#     for (i in 1:Nyear)
#     {
#       plot(sizebins,Ccomp_Pred[i,],type='l',ylim=c(0,0.3),axes=FALSE)
#       x=c(sizebins,rev(sizebins))
#       y=c(rep(0,length(sizebins)),rev(Ccomp_Pred[i,]))
#       polygon(x,y,col="gray")
#       mtext(paste(Year[i]),side=3,line=-1,adj=0.1,cex=0.6)
#       #legend('topleft',paste(Year[i]),cex=1.1,bty='n')
#       rp=vector('expression',2)
#       rp[1]=substitute(expression(N==Value),list(Value=CatchESS_Obs[i]))[2]
#       rp[2]=substitute(expression(effN==Value),list(Value=CatchESS_Pred[i]))[2]
#       legend('topright',legend=rp,bty='n',y.intersp=0.8,yjust=1)
#       lines(sizebins,Ccomp_Obs[i,],col='red')
#       #text(x=27,y=0.15,paste("effN=",SurveyESS_Pred[i]),cex=1)
#       if(i %in% yaxis)
#         axis(2,col='gray',at=seq(0,0.2,0.05))
#       if(i %in% xaxis)
#         axis(1,col='gray',at=sizebins)
#       box(col='gray')
#     }
#     mtext('Length (mm)',side=1,outer=T,line=2.2)
#     mtext('Proportion',side=2,outer=T,line=2.2)
#     mtext(paste('Length Comp,','fishery',fishery,',red line=Obs'),
#           side=3,outer=T,line=2.2)
#   }
# }
# PlotCC_seas(1,1)
# PlotCC_seas(1,2)
# # 
# # # projection
# # run_projection(report = report, par_est = par_est, select_id = c(5,3,4), M = 0.35,
# # 			    n_proj_year = 5, N_rep = 10000, F_trawl_proj = c(0.6,0.6,0.6,0.6,0.6), F_trap_proj = c(0.01,0.01,0.01,0.01,0.01), R_proj = c(500,500,500,500,500), growth_id = c(1,2), R_rand = T, R_sd = 0.5)
# # 
# # # 
# # 
# # Only major diagnostic left are residual plots for index fits
# index_resids <- array(NA, dim = c(length(years),2))
# for(i in 1:(length(years))){
#   index_resids[i,1] <- report$Survey_Index_Obs[1,i]-report$Survey_Index_Pred[1,i]
#   index_resids[i,2] <- report$Survey_Index_Obs[2,i]-report$Survey_Index_Pred[2,i]
# }
# 
# sigma_res1 <- sum(index_resids[,1])/51
# sd_res1 <- sqrt(sigma_res1)
# sigma_res2 <- sum(index_resids[,2])/51
# sd_res2 <- sqrt(sigma_res2)
# 
# index_resids[,1] <- index_resids[,1]/sd_res1
# index_resids[,2] <- index_resids[,2]/sd_res2
# 
# windows(10,10)
# par(mfrow=c(1,2))
# plot(years, index_resids[,1], main = "Index season 1 std residuals")
# lines(years, rep(0,length(years)))
# plot(years, index_resids[,2], main = "Index season 2 std residuals")
# lines(years, rep(0,length(years)))
