
setwd(sce.dir)
# source("reptoRlist.r"); source("ADMB_read.r")
file <- read.rep('testing-2020-1.rep')
# report <- read.rep(file) ;par_est <- std_read("NSLSAP01.std")
setwd(AM.code.dir)
source("PlotFuncs.r")

pdf("proj sq error=TRUE.pdf",width=8,height=6)

PlotWL() # plot length-weight relationships
PlotML() # plot maturity
Plot_GM_tb(1,10) # visualize growth matrices (user defines the number of year-specific time-blocks for growth and how many year to grow)
# Plot_sel_f() # plot fleet selectivity
Plot_f() # plot fishing mortality
PlotTC() # plot total catch
Plot_f_c() # plot fishery-specific catch
PlotCC(1,1) # plot fishery-specific catch composition; specify fishery and season
PlotCC(2,1)
PlotCC(3,1)
PlotCCA(3) # plot aggregated fishery catch composition
PlotSelS(2000,1) # plot survey selectivity; specify year and survey
# PlotSelS(2000,2)
#PlotSelS(2000,3)
PlotSLC() # plot survey size composition
PlotSLCA() # plot aggregated survey size composition
# PlotSI(3) # plot survey index;n= # of indices
PlotM(2000)
PlotR()
PlotSSB()
PlotAbun()
PlotSpB(1)
PlotSexComp(1,0.5)
PlotFfit()
PlotLfifty()
PlotER()

#PlotStockStatus()
#PlotRoSSB(1985,2002,2006,4)
#PlotRo(1984,2012,2016,2)
#PlotRo_scale(1984,2008,2013,2)
#Getlike()
dev.off()

# projection
# run_projection(report = report, par_est = par_est, select_id = c(5,3,4), M = 0.35,
# 			    n_proj_year = 5, N_rep = 10000, F_trawl_proj = c(0.6,0.6,0.6,0.6,0.6), F_trap_proj = c(0.01,0.01,0.01,0.01,0.01), R_proj = c(500,500,500,500,500), growth_id = c(1,2), R_rand = T, R_sd = 0.5)














