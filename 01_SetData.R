
# block_array = array(1, dim=3)

# set size bins
# when you calculate sel, you're using the midpoint of each size class;
# anything size-based, you want to use the midpoint, not just sel, mat,
# growth

# define size bins (mm)
size_bins = seq(5,235,10) # only one crab in 255 and 295 bins for FI data, removed for simplicity

# size_bins = seq(10,230,20)
# midpoints
size_bm = c()
for(i in 2:length(size_bins)){
  size_bm[i-1] = size_bins[i]-5
}
size_bm = size_bm

# set length structure variable based on number of size bins
n_l = length(size_bm)

# set last year to year end to faciliate retrospective analysis
year_end_retro = year_end

# Biological specifications
########### Weight-Length #################
# set weight-length relationship parameters
WL_pars_m = matrix(NA, nrow = n_t, ncol = 2) # important to note whether these and size bins are in mm or cm (currently mm)
avg_fem_par1 <- log(mean(c(0.0034865,0.0000551)))
avg_fem_par2 <- mean(c(2.1165,1.8660))

for (n in 1:n_t)
{
  WL_pars_m[n,]=c(log(0.0034865),2.1165) # Chesbay F
  # WL_pars_m[n,]=c(log(0.0000551), 1.8660) # FL F
  # WL_pars_m[n,]=c(avg_fem_par1, avg_fem_par2) # avg of F pars from ChesBay and FL
}
weight_m <- get_weight(WL_pars_m,size_bm,n_t)
# plot(size_bm, weight_m[1,]) # check

########### Maturity #################
# Define maturity/sex change variables
L50 = 123 # old value 120 # can be used for female maturity even if there are no sex changes
F50 = 110 #length at 50% female maturity (for sex changing spp only) - not used for BC 
inter = quantile(size_bm, 0.75)-quantile(size_bm,0.25) # interquartile range for use in the maturity function
# set maturity/sex change parameters
fe_prop_pars_m = matrix(NA, nrow = n_t, ncol = 2)
for(n in 1:n_t)
{
  fe_prop_pars_m[n,]=c(L50, inter)
}
fe_prop_pars_m
# plot(fe_prop_pars_m[,1]) # check
# for AM data gen
L50_mat = matrix(c(L50), ncol = 1, nrow = n_t, byrow = TRUE)
# check
maturity_m <- get_matu_prop(fe_prop_pars_m,F50,size_bm,n_t)
# plot(size_bm, maturity_m[1,]) # check
# sex ratio
sex_ratio <- array(NA,dim=c(n_t,n_l))# assumes 50:50 split between males and females
sex_ratio[] <- 0.5

########### Natural Mortality #################
M =  0.87 # 0.87
M_year_v = matrix(M, nrow = n_t, ncol = 1) 
# matrix for size effect on M
M_size_v = matrix(1, nrow = 1, ncol = n_l) # no size effect
M_value_as = M
M_m <- get_M(M_size_v,M_year_v)

########### Selectivity #################
# choose selectivity parameters
# logistic selectivity curve parameters (alpha, beta); slope and A50 (typically) age/size at 50% selection
sel_pars_v = c(130, 0.15) # 127

dome_pars = c(10,0.1,150,0.05) # not used currently

# set gear selectivity curve type
# 2 = Logistic
sel_switch = 2
switch_v = matrix(data=rep(2,n_fleets), nrow = n_fleets, ncol=1) 
# set sel parameters 
sel_pars_m = matrix(NA, nrow = n_t, ncol = length(sel_pars_v))
for(n in 1:n_t)
{
sel_pars_m[n,]= sel_pars_v 
}
# check
sels <- cal_sels(sel_pars_v,sel_switch,dome_pars,size_bm)
# plot(size_bm, sels) # check

########### Fishing Mortality #################
# create empty array for F by year, season and fleet - neededed for initialization
F_array = array(NA,c(n_t,n_s,n_fleets)) 
log_F_devs = array(-0.1,dim=dim(F_array)) # set to -1 so assessment model defaults to very small number
# Set initial F values: also -1 to default to small number within the AM
F_init_v = rep(-1,2)

# selectivity blocks and parameters for data generation files
f_sel_n_bloc = 1 

f_sel_switch_v = matrix(NA, ncol = 1, nrow = n_fleets)
for(i in 1:n_fleets)
{
  f_sel_switch_v[i] = switch_v[i]
}
f_sel_block_array = array(1, dim=c(n_t,n_fleets,n_s))

########### Growth #################
# Von Bertalanffy growth function parameters
# Set up minimum and maximum size and increments 
Linc = 10
Lmin = size_bm[1]
Lmax = size_bm[23]
# Increment by which size changes

# P195 est pars: Linf = 217, K = 0.47
# from Eggleston's work in the 2004 SA: Linf = 238, K = 0.74
Linf =  217
# CV for VBGF
vbcv = 0.1
# Standard error of Linf
# CV = sd/mean, sd = CV*mean
SELinf = vbcv*Linf
#K 
K = 0.47 # 
# Standard error of K
SEK = vbcv*K
# correlation between Linf and K, we don't know this, but Linf and K tend to be negatively correlated for finfish
rhoLinfK = -0.4 # mark for exploring - 0.7 is default
# for time-varying growth, create matrices for VBGF parameters - not used
Linf_v = matrix(Linf,nrow = n_t, ncol = 1)
SELinf_v = matrix(SELinf, nrow = n_t, ncol = 1)
K_v = matrix(K, nrow = n_t, ncol = 1)
SEK_v = matrix(SEK, nrow=n_t, ncol = 1)
rhoLinfK_v = matrix(rhoLinfK, nrow = n_t, ncol=1)
# set growth blocks to time-varying block
n_growblock = n_block
# matrix of growblocks for data generation
growblock_m <- matrix(NA, nrow = n_t, ncol = n_growblock)

########### Recruitment #################
# Define recruitment parameters # need for initialization
alpha = 30 # scale is ~ 300 million for BC

beta = 0.005 # not used in mean recruitment model, value from 2018 NC assessment
# SR relationship, 2 = 2 parameter BH
indicator = 1 # 1 is mean recruitment model, 2 is Beverton-Holt, 3 is Ricker, 4 is Cushing, 5 is BH with env effect, 6 is Ricker with env effects
# Spawning month: 
ssb_month = 7

########### Population Dynamics #################
# Create initial proportion at age (PIA) - only needed if PIA option set to 1
N_init_tot = 30 # scale is ~ 300 million for BC

# Initial size structure based on 1982 season one length composition
N_init_comp = c(0.037830121,	0.25945753,	0.140256959,	0.09350464,	0.074232691,	0.061027837,
                0.057458958,	0.053533191,	0.04675232,	0.045681656,	0.034261242,	0.026052819,
                0.030692363,	0.016773733,	0.013204854,	0.007851535,	0.001427552,	0,	0,
                0,	0,	0, 0)

# Uncomment these initial size structures for sensitivity models

# based on initial proportions in FI survey data (10mm bins)
# 1982 annual init lencomps
# N_init_comp = c(0.022263798, 0.182413471,	0.115809167,	0.098409729,	0.081384471,	0.076894294,	0.064546305,
# 0.063049579,	0.052198316,	0.050327409,	0.050888681,	0.03536015,	0.031992516,
# 0.027689429,	0.023573433,	0.017960711,	0.003180543,	0.001122544,	0.000935454,
# 0,	0,	0,	0)

# 1995 annual init lencomps
# N_init_comp = c(0.009893154,	0.120696478,	0.08547685,	0.077166601,
#                 0.074198655,	0.078155916,	0.06747131,	0.059556787,
#                 0.068460625,	0.057578156,	0.061733281,	0.058765334,
#                 0.049663633,	0.043529877,	0.039968342,	0.027107242,
#                 0.012267511,	0.006133755,	0.001780768,	0.000197863,
#                 0.000197863,	0,	0)

# 1995 s1 init lencomps
# N_init_comp = c(0.008287293,	0.159607121,	0.114487416,	0.097298956,	0.090239411,
#                 0.091160221,	0.083486802,	0.066605279,	0.073971762,	0.054941682,
#                 0.042971148,	0.035911602,	0.02854512,	0.024861878,	0.013505218,
#                 0.007673419,	0.004297115,	0.001227747,	0.00092081,	0,	0,	0,	0)


# 1982 s1 init lencomps with p100 lencomps included
# N_init_comp = c(0.294626942,	0.193277311,	0.10516934,	0.06824548,	0.054494525,	0.04456328,
#                 0.040998217,	0.038961039,	0.033613445,	0.032849503,	0.024446142,	0.018843901,
#                 0.021899669,	0.011968424,	0.009421951,	0.005602241,	0.001018589,	0,	0,	0,
#                 0,	0,	0)


# check that it sums to 1
# sum(N_init_comp)
# create vector for recruitment deviations - needed for initialization
R_devs_v = matrix(0, nrow = n_t, ncol = 1)

# set proportion of recruits going to each length class, depends on how we define recruitment

# base assumption: 2 size bins
R_l_pro = c(0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# Uncomment these for sensitivity models
# avg proportions under 130 mm size bin 
# R_l_pro = c(0.029109829,	0.114426497,	0.118427848,	0.134587891,	0.124909423,
#             0.105958414,	0.082603822,	0.071400257,	0.06153243,	0.053166683,
#             0.051807888,	0.052069019, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# R_l_pro = c(0.4,0.4,0.1,0.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
# R_l_pro = c(1.0,0,0,0,0,0,0,0,0,0,0)
# R_l_pro = c(0.2,0.2,0.1,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0,0,0,0,0,0,0,0,0,0,0)
# R_l_pro = c(0.2,0.2,0.2,0.1,0.1,0.1,0.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


n_R_l = length(which(R_l_pro != 0))

# sd for rec devs
sd_Rdevs = 0.7

# for environmental effects
# fit environmental effects to rec devs? 1 = yes, 0 = no
R_fit_env = 0

# number of environmental data sets - default setting is 3
n_env = 3

# environmental data matrix (year x n_env)
env_data = matrix(0, nrow = n_t, ncol = n_env)
# for(i in 1:n_t){
#   env_data[i] <- env.data[i,2] # uncomment if R_fit_env = 1
# }

########### Survey #################
# set survey catchability 
survey_q = exp(-6)

# Survey selectivity parameters (first season or annual)
dome_pars_s1 = c(20,0.5,300,1)
survey_sel_pars =  dome_pars_s1

# select survey selectivity relationship; 2 = logistic, 3 = double logistic, 4 = 4 par double
survey_sel_switch = 3
survey_sels <- cal_sels(survey_sel_pars,survey_sel_switch,dome_pars_s1,size_bm)
# plot(size_bm, survey_sels) # check

dome_pars_s2 = c(150,0.02,175,0.05)
survey_sels_s2 <- cal_sels(survey_sel_pars,survey_sel_switch,dome_pars_s2,size_bm)
# plot(size_bm, survey_sels_s2) # check

# matrix of survey sel pars
if(n_s > 1){
survey_sel_pars_m = matrix(NA, ncol = length(survey_sel_pars), nrow = n_s)
survey_sel_pars_m = rbind(dome_pars_s1, dome_pars_s2)
} else {
  # survey_sel_switch = 2
  survey_sel_pars_m = matrix(NA, ncol = length(survey_sel_pars), nrow = n_t) # use if annual, will overwrite above
  for (n in 1:n_t)
  {
    survey_sel_pars_m[n,] = survey_sel_pars # currently assumes same selectivity by season
  }
}


