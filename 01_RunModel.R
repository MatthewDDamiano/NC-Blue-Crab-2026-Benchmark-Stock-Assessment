# NC DMF blue crab stock assessment
# M. Damiano & M. Zink
# Last modified: MDD 4-2-26

rm(list=ls())
require(PBSadmb)
require(psychTools)
library(here)
setwd(here())
# Run script for model functions to get unchanging values like LW, Maturity
source("01_Functions_old.R")

# Run script to generate data inputs
# Set up data
# Model dimension
first_year = 1982
terminal_year = 2024
years = seq(first_year,terminal_year,1)

# read in data
catch.data <- as.matrix(read.csv("BC_catch_2seas.csv"))
survey.data <- as.matrix(read.csv("BC_survey_2seas.csv"))
env.data <- as.matrix(read.csv("NAO index annual_1982_2024_1yrlag.csv")) # loaded but only for plotting if function not used

# Some settings

# set the number of years 
years = rep(first_year:terminal_year)

# set annual time step
n_t = length(years)

year_start = years[1]

year_end = tail(years, n = 1)

# set number of fleets
n_fleets = 1

# set seasonal timestep; 1=no seasonality
n_s = c(1,2)[2] 

# for retro
nsteps = n_t*n_s

# set blocks for time-varying quantities, parameters
n_block = n_s

# set number of surveys
n_survey = 2

# proportion of recruitment in each season
rec_seas_prop = c(0.9,0.1) # set to 1 if annual
# proportion of growth in each season
grow_seas_prop = c(0.5,0.5) # set to 1 if annual

# run data file setup
source("01_SetData.R")

# Create input files for each model run
# Set input file directory
AM.code.dir = here("Assessment model and figure code")
result.dir = here()
scenario         = c('Base model github test 10 9-22-2026') # give a new name for each run so it creates a separate results folder
sce.dir          = paste(result.dir,'/',scenario,sep='')

input.dir        = paste(result.dir,'/',scenario,'/InputFiles',sep='')
input.season.dir = paste(input.dir,'/Season',sep='')
input.year.dir   = paste(input.dir,'/Year',sep='')

dir.create(sce.dir)
dir.create(input.dir)
dir.create(input.season.dir)
dir.create(input.year.dir)

file.copy(from=paste(AM.code.dir,'/NSLSAP01.exe',sep=''),to=sce.dir) 
file.copy(from=paste(AM.code.dir,'/PlotFuncs.R',sep=''),to=sce.dir)
file.copy(from=paste(AM.code.dir,'/reptoRlist.R',sep=''),to=sce.dir)
file.copy(from=paste(AM.code.dir,'/OutputPlots.R',sep=''),to=sce.dir)
file.copy(from=paste(AM.code.dir,'/01_RunPlots.R',sep=''),to=sce.dir)
file.copy(from=paste(AM.code.dir,'/01_RunProj.R',sep=''),to=sce.dir)
file.copy(from=paste(AM.code.dir,'/01_ReadADMB.R',sep=''),to=sce.dir)
file.copy(from=paste(AM.code.dir,'/01_BRP.R',sep=''),to=sce.dir)

# Set model priors/parameters to estimate
# I took this out of a separate R file to make quick alternate model runs a bit easier

# Prior for M - cannot be estimated for blue crab
M_prior = matrix(c(0.87, 0.1, 1.2, -1, 0.2, 0.0, 1), nrow = n_s, ncol = 7, byrow = TRUE) 

# Lorenzen natural mortality b parameter prior - not used for blue crab
M_b_prior = matrix(c(0.25, 0.1, 1.0, -1, 0.2, 0.0, 1), nrow = n_s, ncol = 7, byrow = TRUE)

# cpue_power_prior
cpue_power_prior = c(1, 0.5, 1.5, -1, 0.5, 0, 1)
# fleet selectivity prior
f_sel_pars_prior = matrix(c(20, 0, 1000, -1, 1, 0, 1), nrow = length(sel_pars_v), ncol = 7, byrow = TRUE)
# initial F prior
f_first_prior1 = matrix(c(F_init_v[1], 0, 5, 1, 1, 0.0, 1), nrow = n_fleets, ncol = 7, byrow = TRUE)
f_first_prior2 = matrix(c(F_init_v[2], 0, 5, 1, 1, 0.0, 1), nrow = n_fleets, ncol = 7, byrow = TRUE)
if(n_s > 1){
  f_first_prior = rbind(f_first_prior1, f_first_prior2)
} else {
  f_first_prior = f_first_prior1
  F_init_v = -0.1
}
# F deviations prior
f_devs_prior = matrix(c(0, -1000, 1000, 1, -1, 0.0, 1), nrow = n_s*n_fleets, ncol = 7, byrow = TRUE)
# Growth priors (season one)
Linf_prior1 = c(100, 50, 400, -1, 0.5, 0, 1)
K_prior1 = c(K, 0, 0.9, 1, 0.5, 0, 1)
sd_Linf_prior1 = c(0.1, 0, 0.5, -1, 0.5, 0, 1)
sd_K_prior1 = c(SEK, 0, 0.5, 1, 0.5, 0, 1)
rho_prior1 = c(-0.4, -0.5, 0.5, 1, 0.5, 0, 1)
# Growth priors (season two)
Linf_prior2 = c(100, 50, 400, -1, 0.5, 0, 1)
K_prior2 = c(K, 0, 0.9, 1, 0.5, 0, 1)
sd_Linf_prior2 = c(0.1, 0, 0.5, -1, 0.5, 0, 1)
sd_K_prior2 = c(SEK, 0, 0.5, 1, 0.5, 0, 1)
rho_prior2 = c(-0.4, -0.5, 0.5, -1, 0.5, 0, 1)

# Initial N prior
N_init_prior = c(N_init_tot, 0, N_init_tot*5, 1, 10, 0, 1)
# R proportion prior
R_l_pro_prior = c(0.1, 0, 1.0, 1, 0.5, 0, 1)
# SRR priors
alpha_prior = c(alpha, alpha*0.001, alpha*5, 1, 0.5, 0, 1)
beta_prior = c(beta, beta*0.01, beta*100, -1, 0.2, 0, 1)
# environmental effects prior
R_env_prior1 = matrix(c(0.25, 0.1, 1.0, -1, 0.2, 0.0, 1), nrow = 1, ncol = 7, byrow = TRUE)
R_env_prior2 = matrix(c(0.25, 0.1, 1.0, -1, 0.2, 0.0, 1), nrow = 1, ncol = 7, byrow = TRUE)
R_env_prior3 = matrix(c(0.25, 0.1, 1.0, -1, 0.2, 0.0, 1), nrow = 1, ncol = 7, byrow = TRUE)
R_env_prior = rbind(R_env_prior1, R_env_prior2, R_env_prior3)
# rec autocorrelation prior
R_autocor_prior = c(0, 0, 1, -1, 1, 0, 1)
# # rec devs prior information
R_devs_v_prior = c(0, -10, 10, 1, 1, 0, 1)
# R_devs_v_prior = c(0.05, -0.4, 0.62, 1, -1, 0, 1)
sd_R_devs_prior = c(sd_Rdevs, 0, 1, -1, -1, 0, 1)

# # sex information
Rsex = c(0.1, 0, 100, -1, 0.5, 0, 1)
# # length at 50% selection prior
L50_prior = c(123, 80, 150, -1, 0.5, 0, 1)

# survey priors
# survey catchability power prior
q_power_prior = c(exp(-6), 0, 1.5, 1, 0.5, 0, 1)
# survey selectivity prior (season one)
survey_sel_prior1_1 = matrix(c(20, 0, 1000, 1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)
survey_sel_prior1_2 = matrix(c(20, 0, 1000, 1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)
survey_sel_prior1_3 = matrix(c(20, 0, 1000, -1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)
survey_sel_prior1_4 = matrix(c(20, 0, 1000, -1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)
# survey selectivity prior (season two)
survey_sel_prior2_1 = matrix(c(20, 0, 1000, -1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)
survey_sel_prior2_2 = matrix(c(20, 0, 1000, 1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)
survey_sel_prior2_3 = matrix(c(20, 0, 1000, 1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)
survey_sel_prior2_4 = matrix(c(20, 0, 1000, 1, 1, 0.0, 1), nrow = 1, ncol = 7, byrow=TRUE)

if(n_s > 1){
  survey_sel_prior = rbind(survey_sel_prior1_1, survey_sel_prior1_2, survey_sel_prior1_3, survey_sel_prior1_4, survey_sel_prior2_1,
                           survey_sel_prior2_2, survey_sel_prior2_3, survey_sel_prior2_4)
} else {
  survey_sel_prior = rbind(survey_sel_prior1_1, survey_sel_prior1_2)
}

# Number of years projected - not used, but required for data file
num.proj = 1

# Run script to generate .dat files for model input
source("01_DataGen.R")

# Run assessment model
setwd(sce.dir)
file.create("NSLSAP01.DAT")
runAD("NSLSAP01",verbose=T) 

# Run plots
source("01_RunPlots.R")
