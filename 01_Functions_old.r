# Last update on 03/31/2023 by M. Damiano
# 3/31/2023 - removed weight from SSB equations

############## weight-length ###############
get_weight <- function(WL_pars_m,size_bm,n_t) #size_bm is a 1 row matrix of size bins, WL_pars_M is a 1x2 matrix of weight parameters, and n_t is the number of years
{
  weight_m = matrix(NA,nrow=n_t,ncol=length(size_bm))
  for (t in 1:n_t)
  {
    weight_m[t,] = exp(WL_pars_m[t,1]+WL_pars_m[t,2]*log(size_bm)) #should create a 40x20 matrix of weight at length (k) for each year, time invariant
  }
  return(weight_m)
}

############## maturity-sex change ###################
#Age at 50% transition is 3.83, logistic function for BSB 2018
# Introduce switch for species that sex change
get_matu_prop <- function(fe_prop_pars_m,F50,size_bm,n_t) #fe_prop_pars_m is a 1x2 matrix of the 2 maturity parameters
{
  maturity_m = matrix(NA,nrow=n_t,ncol=length(size_bm))
  per_fem = matrix(NA,nrow=n_t,ncol=length(size_bm))
  for (t in 1:n_t)
  {
    maturity_m[t,] = (1/(1+exp(-2*log(3)*(size_bm-fe_prop_pars_m[t,1])/fe_prop_pars_m[t,2])))# should create logistic maturity 40x20 matrix for use in popdy equation
    # per_fem[t,] = 1-(1/(1+exp(-2*log(3)*(size_bm-F50)/fe_prop_pars_m[t,2])))
  }
  maturity_m <- maturity_m # *per_fem # turns off dome-shape for female-only mature
  #print(maturity_m)
  return(maturity_m)
}

# get_fec <- function(fec_pars_m,size_bm,weight_m,n_t){
#   fec_mat = matrix(NA, nrow = n_t, ncol = length(size_bm))
#   for(t in 1:n_t){
#     fec_mat[t,] = exp(fec_pars_m[t,1]+fec_pars_m[t,2]*log(weight_m[t,]))
#   }
#   fec_mat <- fec_mat/weight_m # added to workaround SSB calcs issue / standardize fec
#   return(fec_mat)
# }

############## set natural mortality ############
get_M <- function(M_size_v,M_year_v) #vectors of mortality by season, size and year (40x1), can get rid of seasonality here too, I think; M=0.38 constant, but decreasing with age
{
  
  n_l=length(M_size_v)
  n_t=length(M_year_v)
  
  M_m = array(NA, c(n_t,n_l,n_s))
  #What you could do instead, if you wanted to, is change this to the Lorenzen method but a M at length matrix instead of M at age
  for (n in 1:n_t)
  {
    M_m[n,,] = M_year_v[n]*M_size_v #Should give is a 40x20 matrix of M by length and year
  }
  # M_m
  return(M_m) #Modified to remove seasonality (10-18-19), returns a constant mortality schedule by year and length class
}

############# selectivities ###############
cal_sels <- function(sel_pars_v,sel_switch,dome_pars,size_bm) #pars_v contains the two parameters for a logistic selectivity curve. This is fine given that fleets for BSB follow this pattern
{
    if (sel_switch==2&length(sel_pars_v)!=2) 
    stop("check # of parameters")
  
    if (sel_switch==3&length(dome_pars)!=4)
    stop("check # of parameters for dome")
  
    sels <- array(NA, n_l)
    
    if (sel_switch==2) { # if logistic/flat-topped
      sels <- 1.0/(1.0+exp((sel_pars_v[1]- size_bm)* sel_pars_v[2])); #Gives a selectivity by length, changes pars_v[2] to negative to reflect logistic curve. 
      #Think of slope as how fast sel increases with size
      sels <- sels/max(sels)
    }
    
    else { # if dome-shaped, use double logistic model (Methot 1990)
      sels_temp1 <- 1.0/(1.0+exp((dome_pars[1] - size_bm)* dome_pars[2]))
      sels_temp2 <- 1.0 - 1.0/(1.0+exp((dome_pars[3] - size_bm)* dome_pars[4]))
      sels <- (sels_temp1*sels_temp2)/max(sels_temp1*sels_temp2)
  } 
  sels 
  return(sels)
}

###############
get_sels <- function(sel_pars_v,sel_switch,switch_v,n_fleets,size_bm,n_s,n_t,n_l)
  #create a selectivity array of n_t x n_l matrices by fleet
{
  # dim pars_m - c(n_block,6)
  #block_array <- array(NA, c(n_t,n_s,n_block))
  #block_array
  # swith_v swith for each block
  sels_4d = array(NA,c(n_t,n_l,n_s,n_fleets))
  #sels_4d
  for (f in 1:n_fleets){
    for (s in 1:n_s){
      for (t in 1:n_t){ # needs the block back so it can switch between sel
        
          sels_4d[t,,s,f] = cal_sels(sel_pars_v,switch_v[f],dome_pars,size_bm) #now it just calls the same selectivity parameters (logistic, 2 param)
      }
    }
  }
  sels_4d
  return(sels_4d)  # problem is in selectivities, I think 4-26-22
}

############## fishing mortality ##############
get_F <- function(F_array,sel_pars_v,sel_switch_v,n_fleets,size_bm,n_s,n_t,n_l)
{
  # dim(F_m) - c(n_t,n_s,n_fleets)
  sels <- get_sels(sel_pars_v,sel_switch,switch_v,n_fleets,size_bm,n_s,n_t,n_l)
  
  F_4d = array(NA,c(n_t,n_l,n_s,n_fleets))
  
  for (f in 1:n_fleets){
    for (s in 1:n_s){
      for (t in 1:n_t){
        F_4d[t,,s,f] = F_array[t,s,f]*sels[t,,s,f] 
      }
    }
  }
  
  F_3d = array(0,c(n_t,n_l,n_s))
  temp = array(0,c(n_t,n_l,n_s))
  # for (f in 1:n_fleets){
  #   temp = F_4d[,,,f]
  #   F_3d = temp + F_3d
  # }
  
  for (f in 1:n_fleets){ # error - inflating F at t10 to insane size
    if (n_s == 1){ #Whether 4D or 3D depending on if season is set to 1 for annual time steps
      temp[,,1] = F_4d[,,,f]
    }else{
      temp = F_4d[,,,f]
    }
    F_3d = temp + F_3d
  }
  
  f_obj = list(F_3d=F_3d,F_4d=F_4d)
  return(f_obj)
}

get_F_inits <- function(F_array) 
{
  f_first_year_index      = matrix(0,nrow=n_fleets,ncol=n_s)       # dim - c(n_fleets,season)
  log_f_first_year        = matrix(0,nrow=n_fleets,ncol=n_s)
  log_F_devs              = array(0,dim=dim(F_array))
  
  for (f in 1:n_fleets)
  {
    for (s in 1:n_s)
    { 
      for (t in 1:n_t)
      {
        if (F_array[t,s,f]>0)
        {
          f_first_year_index[f,s] = t
          break
        }
      }
    }
  }
  f_first_year_index
  for (f in 1:n_fleets)
  {
    for (s in 1:n_s)
    {
      year.index.temp        = f_first_year_index[f,s]
      
      if (year.index.temp > 0){
        log_f_first_year[f,s]  = log(F_array[year.index.temp,s,f])
        
        for (t in 1:n_t)
        {
          if (t==f_first_year_index[f,s])
          {
            log_F_devs[t,s,f] = 0
          }
          if (t>f_first_year_index[f,s])
          {
            # log_F_devs[t,s,f] = F_devs_array[t,s,f] # For conditioned values provided
            log_F_devs[t,s,f] = log(F_array[t,s,f])-log(F_array[t-1,s,f]) # if no conditioned values provided
          }
        }
        break
      }
    }
  }
  log_F_devs
  f_inits.obj = list(log_f_first_year=log_f_first_year,log_F_devs=log_F_devs)
  f_inits.obj  
  return(f_inits.obj)
  
}

############## growth transition matrix #############
#Linf = 502.1 mm
#K = 0.173
#t0 = -0.97

cal_GM <- function(Lmin,Lmax,Linc,Linf,SELinf,K,SEK,rhoLinfK)
{
  growmat <- NULL
  rhoLinfK <- abs(rhoLinfK)
  #7-9-2020: Here's the issue, this needs to measure distance between Linf and midpoint
  #7-18-2020: I think the problem is with the midpoints themselves. Even with alternative calcs
  #we never get the same result as the GTM in the AM built using the VBGF params.
  Ln <- seq(Lmin, Lmax+Linc, Linc)
  Ln
  COV <- rhoLinfK * SELinf * SEK
  #Delta L, change in length
  #The growth increment. 
  #7-9-2020: Here's the issue, this needs to measure distance between Linf and midpoint
  DL <- (Linf - Ln) * (1 - exp(-K))
  DL
  # DL <- DL[1:(length(which(DL >= 0)) + 1)] #7-17-2020: to prevent "shrinkage"
  # Ln <- Ln[1:(length(which(DL >= 0)) + 1)]
  DL <- ifelse(DL < 0, 0, DL)
  # DL[is.na(DL)] <- 0
  # Ln[is.na(Ln)] <- 0
  #Variance of delta L, 11-1-19 changed - 2*COV to + 2*COV per equation in 2017 paper
  # 6/26/2020: changed it back to be consistent with AM
  VL <- SELinf^2 * (1 - exp(-K))^2 + (Linf - Ln)^2 * SEK^2 * 
    exp(-K)*exp(-K) - 2 * COV * (1 - exp(-K)) * (Linf - Ln) * 
    exp(-K)
  VL
  sqrt(VL)
  Ln+DL
  #7-9-2020: math is consistent with AM for calculating variance
  growmat <- matrix(0, nrow = length(Ln) - 1, ncol = length(Ln) - 1)
  growmat
  for (L in 1:as.numeric(length(Ln) - 1)) {
    #print(Ln[L]+DL[L])
    # 7-18-2020: Needed to add -0.5*Linc to get appropriate mean; matches Excel GTM
    for (m in L:as.numeric(length(Ln) - 1)) {
      growmat[L, m] <- pnorm(Ln[m + 1]-(0.5*Linc), mean = Ln[L]+  # add -0.5*Linc to get midpoint size
                               DL[L], sqrt(VL[m])) - pnorm(Ln[m]-(0.5*Linc), mean = Ln[L]+ # add -0.5*Linc
                                                             DL[L], sqrt(VL[m]))
    }
  }
  growmat
  #11-1-19 need to replace row and column 20 with zeros b/c 0 probability of growing if they are as large as they can be
  growmat[is.na(growmat)] <- 0 
  growmat
  #11-1-19 every row should sum to 1 by rules of probability
  
  growmat <- growmat/rowSums(growmat, na.rm = T) 
  n.dim=(Lmax-Lmin)+Linc
  growmat_ful = matrix(0,nrow=n.dim,ncol=n.dim)
  
  if (nrow(growmat)<n.dim){
    growmat_ful[1:nrow(growmat),1:ncol(growmat)]=growmat
    for(nn in 1:(n.dim-nrow(growmat))){
      growmat_ful[nrow(growmat)+nn,ncol(growmat)+nn]=1
    }
  }else{
    growmat_ful=growmat[1:n.dim,1:n.dim]
  }
  
  return(growmat)
  
} 

get_GM <- function(n_l,n_growblock,Lmin,Lmax,Linc,Linf_v,SELinf_v,K_v,SEK_v,rhoLinfK_v)
{
  
  grow_3d = array(0,c(n_l,n_l,n_growblock))
  grow_3d
  for (t in 1:n_growblock){
    grow_3d[,,t] = cal_GM(Lmin,Lmax,Linc,Linf_v[t],SELinf_v[t],K_v[t],SEK_v[t],rhoLinfK_v[t])
  }
  return(grow_3d)
}

############## Stock Recruitment ################
get_R <- function(SSB,alpha,beta,indicator)
{
  if (indicator==1)R=alpha #no functional relationship
  if (indicator==2)R=alpha*SSB/(beta+SSB) #BH
  if (indicator==3)R=alpha*SSB*exp(-beta*SSB) #Ricker
  return(R)
}

get_ssb.fracs <- function(ssb_month,n_s) #This creates periods during the year when spawning biomass is contributed
{
  SSB_frac_v    = c()
  SSB_frac_v[1] = (ssb_month-1)/12
  temp=SSB_frac_v[1] *12;temp1=0
  for(i in 1:n_s)
  { 
    temp1=temp1+12/n_s
    if(temp<=temp1)
    {
      SSB_frac_v[2] = i
      SSB_frac_v[3] = (temp-(temp1-12/n_s))/(12/n_s)
      break
    }
  }
  return(SSB_frac_v)
}


############### pop dyn ##########################
pop_dyn <- function(N_init_tot,N_init_comp,F_3d,M_m,SSB_frac_v,maturity_m,weight_m,alpha,beta,indicator,R_devs_v,R_l_pro,n_t,n_s,n_growblock,growth_array,sex_ratio)
{ #10-25-19 changed inputs to alpha, beta and indicator instead of alpha_RS and switch_RS etc. to be consistent with get_R function and for simplicity, removed seasonal recruit component
  N_3d = array(0,c(n_t,n_l,n_s))
  SSB  = c()
  R    = c()
  
  for (t in 1:n_t){
    
    if(t>1)
      R[t] = get_R(SSB[t-1],alpha,beta,indicator)*exp(R_devs_v[t])
    
    for (s in 1:n_s){
      
      if(t==1&s==1){
        
        # F_temp_year = apply(F_3d[1,,1:n_s],1,sum)
        # M_temp_year = apply(M_m[1,,1:n_s],1,sum)
        
        F_temp_year = F_3d[1,,1] #10-25-19 error, dimensions are wrong. Currently, F_3d is n_t, n_l, n_fleets, would be n_t, n_l, n_s, n_fleets if 4D
        M_temp_year = M_m[1,,] #M_m is n_t, n_l dimensions
        #Checks
        M_m
        F_temp_year
        M_temp_year
        F_4d
        exp(-SSB_frac_v[3]*(F_temp_year+M_temp_year))
        
        N_3d[1,,1] = N_init_tot*N_init_comp #Starting N numbers by length
        N_3d
        #SSB[1] = sum(N_3d[1,,1]*exp(-SSB_frac_v[1]*(F_temp_year+M_temp_year))*PropFemale_m[1,]*weight_m[1,])
        SSB[1] = sum(N_3d[1,,1]*exp(-SSB_frac_v[3]*(F_temp_year+M_temp_year))*sex_ratio*maturity_m[1,]*weight_m[1,]) #Problem here is that SSB is a vector; I think it's trying to apply the sum of the first year of mortality and apply it to N to get SSB
        # -- to match up with assessment model (SSB1 will only be correct when spawning time is within season1)
        SSB
        R[1] = get_R(SSB[1],alpha,beta,indicator)*exp(R_devs_v[1])
        R
        N_3d[1,,1] = N_3d[1,,1] + R[1]*R_l_pro
        N_3d
        
      }
      
      if(t>1&s==1){
        N_3d[t,,s]=(N_3d[t-1,,n_s]*exp(-(F_3d[t-1,,n_s]+M_m[t-1,,1])))%*%growth_array[,,t-1]+R[t]*R_l_pro #%*% operator is for matrix mult
        N_3d
        n_growblock[t-1]
        if(SSB_frac_v[2]==1&n_s>1){
          SSB[t] = sum(N_3d[t,,1]*exp(-SSB_frac_v[3]*(F_3d[t,,]+M_m[t,]))*sex_ratio*maturity_m[t,]*weight_m[t,])
        }
      }
      
      if(n_s==1&t>1){
        SSB[t] = sum(N_3d[t,,n_s]*exp(-SSB_frac_v[1]*(F_3d[t,,n_s]+M_m[t,,1]))*sex_ratio*maturity_m[t,]*weight_m[t,])
      }
    }
  }
  #If there is seasonality in the model config. 
  if(s>1){
    N_3d[t,,s]=(N_3d[t,,s-1]*exp(-(F_3d[t,,s-1]+M_m[t,,s-1])))%*%grow_3d[,,growblock]+R[t]*R_l_pro
    
    if(s==SSB_frac_v[2]&n_s>1){
      SSB[t] = sum(N_3d[t,,s]*exp(-SSB_frac_v[3]*(F_3d[t,,s]+M_m[t,,s]))*sex_ratio*maturity_m[t,]*weight_m[t,])
    }
  }
  
  pop.dyn.obj = list(Abundance=N_3d,SSB=SSB,Recruits=R,F_3d=F_3d,M_2d=M_m)
  return(pop.dyn.obj)
}


############### calulate catch ##########################
cal_catch <- function(f_obj,M_m,N_3d,n_t,n_l,n_s,n_fleets)
{
  catch_4d  = array(NA,c(n_t,n_l,n_s,n_fleets))
  catch_3d  = array(NA,c(n_t,n_fleets,n_s))
  catch_tot = matrix(NA,nrow=n_t,ncol=n_fleets)
  
  F_3d = f_obj$F_3d
  F_4d = f_obj$F_4d
  
  
  for (f in 1:n_fleets)
  {
    for (t in 1:n_t)
    {
      for (s in 1:n_s)
      {
        catch_4d[t,,s,f] = N_3d[t,,s]*(1-exp(-(F_3d[t,,s]+M_m[t,,s])))*(F_4d[t,,s,f]/(F_3d[t,,s]+M_m[t,,s]))#Just the baranov catch equation, easy
        catch_3d[t,f,s] = sum(catch_4d[t,,s,f])
      }
    }
  }
  
  catch.obj = list(catch_4d=catch_4d,catch_3d=catch_3d)
  return(catch.obj)   
}

#10-28-2019 - also runs, but need to address NaNs before it's functional 
#11-11-2019 - Main issue here is that somewhere, something is driving the population too far low. The original problem was that the sex change pars were not on the log scale,
#so no SSB was being produced ergo no recruits. By changing that, we've gotten closer, but the survey data generation with error still doesn't work, and I'm pretty sure
#that's because popdy numbers are too close to zero in some places. So carefully review each of the functions contributing to popdy and see that all pars are on the right scale.

############### generate survey data (no error) ##########################
get_survey.fracs <- function(survey_months,n_survey,n_s) 
{
  survey_frac_m    = matrix(NA,nrow=n_survey,ncol=3) 
  for (n in 1:n_survey)
  {
    survey_frac_m[n,1] = survey_months[n]/12
    
    temp=(survey_months[n]/12)*12
    temp1=0
    
    for(i in 1:n_s)
    {
      temp1=temp1+12/n_s
      if (temp<=temp1)
      {
        survey_frac_m[n,2] = i
        survey_frac_m[n,3] = (temp-(temp1-12/n_s))/(12/n_s)
        break
      }
    }
    
  }
  return(survey_frac_m)
}

cal_survey_inds <- function(N_3d,F_3d,M_m,survey_frac_m,survey_q,survey_sel_switch,survey_sel_pars,dome_pars,survey_periods,size_bm,n_survey,n_s,n_t,n_l)
{  
  survey.temp  = array(0,c(n_t,n_l,n_survey))
  survey.sels  = matrix(NA,nrow=n_survey,ncol=n_l)
  survey.obj   = list()
  year_start = survey_periods[,1]
  year_end = survey_periods[,2]
  year_start
  year_end
  
  for (n in 1:n_survey)
  {
    survey.sels[n,] = cal_sels(survey_sel_pars,survey_sel_switch,dome_pars,size_bm)
    survey.sels #looks good
    n1 = survey_periods[n,1] - year_start + 1
    n2 = n_t - (year_end-year_start)+1
    n1
    n2
    
    for (t in 1:n_t)
    {
      if(n_s <= 1) 
      {
        survey.temp[t,,n] = N_3d[t,,1]*exp(-survey_frac_m[n,1]*(F_3d[t,,1]+M_m[t,,1]))*survey_q[n]*survey.sels[n,]
        survey.temp 
        survey_q
        N_3d[t,,1]
      }
      if(n_s > 1)
      {
        temp = survey.fracs[n,2]
        survey.temp[t,,n] = N_3d[t,,temp]*exp(-survey_frac_m[n,3]*(F_3d[t,,temp]+M_m[t,,temp]))*survey_q[n]*survey.sels[n,]
      }
      
    }
    survey.obj[[n]] = cbind(apply(survey.temp[,,n],1,sum),survey.temp[,,n]) 
  }
  return(survey.obj)
  
}

############### generate data (with random error) ##########################
add.error.lognormal <- function(data,cv)
{
  #dim.data   = dim(data)  # row - # of year; col - # of data set (fleet or survey)
  #data.error = matrix(NA,nrow=dim.data[1],ncol=dim.data[2])
  
  #for (n in 1:dim.data[2])
  #{
  sd               = sqrt(log((cv^2+1)))
  data.error       = rlnorm(length(data),log(data),sd)
  #}
  return(data.error)
}

add.error.multino <- function(data,ess)
{
  dim.data      = dim(data)  # row - # of year; col - # of data set (fleet or survey)
  data.error    = matrix(-1,nrow=dim.data[1],ncol=dim.data[2])
  sample.temp   = matrix(-1,nrow=dim.data[1],ncol=dim.data[2])
  
  original.com=sweep(data,1,rowSums(data),"/")
  
  for (i in 1:dim.data[1])
  {
    if (sum(data[i,])>0)
    {
      sample.temp[i,] = rmultinom(1, ess[i], prob=original.com[i,])
      data.error[i,]  = sample.temp[i,]/sum(sample.temp[i,])
    }  
  }
  return(data.error)
}

catch.data.gen <- function(year_start,year_end,catch.obj,catch.cv,catch.ess,add.error)
{
  
  n.season       = lapply(catch.obj, dim)$catch_4d[3]
  n.fleets       = lapply(catch.obj, dim)$catch_4d[4]
  n.lengthbin    = lapply(catch.obj, dim)$catch_4d[2]
  n.year         = lapply(catch.obj, dim)$catch_4d[1]
  
  n_s = n.season
  n_t = n.year
  n_fleets = n.fleets
  
  year_start = year_start
  year_end = year_end
  catch.data.m   = matrix(-1,nrow=n_t*n_s*n_fleets,ncol=9+length(size_bm)) #changed n.lengthbins to what I've been using to count them 'size_bm'
  catch.data.m
  length(size_bm)
  
  year      = rep(seq(year_start,year_end,1),n_s*n_fleets)
  season    = rep(rep(1:n_s,each=n_t),n_fleets)
  fleet     = rep(rep(1:n_fleets,each=n_t*n_s))
  catch.tot = c()
  cv        = c()
  cpue      = rep(-1,n_t*n_s*n_fleets)
  cpue.flag = rep(1,n_t*n_s*n_fleets)
  cpue.cv   = rep(-1,n_t*n_s*n_fleets)
  ess       = c()
  comp      = data.frame()
  
  for (f in 1:n_fleets)
  {
    for (s in 1:n_s)
    { 
      
      cv.temp          = catch.cv[,f,s]
      cv               = c(cv,cv.temp)
      cv
      catch.tot.temp   = catch.obj$catch_3d[,f,s]   
      catch.tot        = c(catch.tot,catch.tot.temp)
      catch.tot
      ess.temp         = catch.ess[,f,s]
      ess              = c(ess,ess.temp)
      ess
      comp.temp        = catch.obj$catch_4d[,,s,f]
      comp             = rbind(comp,comp.temp)
      comp
      
      if(add.error==TRUE)
      {
        catch.tot    = add.error.lognormal(catch.tot,cv)
        comp         = add.error.multino(comp,ess)
      }
      
    }
  }
  
  catch.data.m[,1]                   = year
  catch.data.m[,2]                   = season
  catch.data.m[,3]                   = fleet
  catch.data.m[,4]                   = catch.tot
  catch.data.m[,5]                   = cv
  catch.data.m[,6]                   = cpue
  catch.data.m[,7]                   = cpue.flag
  catch.data.m[,8]                   = cpue.cv
  catch.data.m[,9]                   = ess
  catch.data.m[,10:(9+n.lengthbin)]  = as.matrix(comp)
  catch.data.m
  return(catch.data.m)
  
}

#11-12-19 - the problem with this function appears to be with the nrows of survey.ess. I've made some adjustments below
#so that the survey ess nrows will reflect whatever we specify as the number of years. This will have to be fixed to accomodate
#multiple surveys, but for now, this does the job. 

# survey.data.gen <- function(survey_periods,survey_months,n_survey,survey.obj,survey.cv,survey.ess,add.error)
# {
#   nrow           = 0
#   nrow.temp      = 0
#   year           = c()
#   indexN         = c()
#   month          = c()
#   index          = c()
#   cv             = c()
#   ess            = c()
#   comp           = data.frame()
#   n_sizeb        = lapply(survey.obj, dim)[[1]][2] #1 gives dimensions, 2 calls the columns' dimensions
#   n_sizeb
#   for (i in 1:n_survey)
#   {
#     nrow.temp    = lapply(survey.obj, dim)[[i]][1]
#     nrow         = nrow + nrow.temp
#     nrow
#     year.temp    = seq(survey_periods[i,1],survey_periods[i,2],1)
#     year         = c(year,year.temp)
#     year
#     length(year)
#     indexN.temp  = rep(i,nrow.temp)
#     indexN       = c(indexN,indexN.temp)
#     indexN
#     length(indexN)
#     month.temp   = rep(survey_months[i],nrow.temp)
#     month        = c(month,month.temp)
#     month
#     length(month)
#     index.temp   = survey.obj[[i]][,1]
#     index        = c(index,index.temp)
#     index
#     length(index)
#     cv.temp      = survey.cv[[i]]
#     cv           = c(cv,cv.temp)
#     cv
#     ess.temp     = survey.ess[,i,] #Hot fix gets it to work. I think the main error here has to do with the array being n_t rows
#     ess          = c(ess,ess.temp)
#     ess
#     comp.temp    = survey.obj[[i]][,2:n_sizeb]
#     comp         = rbind(comp,comp.temp)
#     comp
#     colSums(comp)
#     dim(comp)
#     
#     if(add.error==TRUE)
#     {
#       index      = add.error.lognormal(index,cv)
#       comp       = add.error.multino(comp,ess)
#       
#     } 
#     
#   }
#   
#   survey.data.m  = matrix(-1,nrow=length(index),ncol=(6+n_sizeb-1))
#   survey.data.m
#   survey.data.m[,1]                 = year
#   survey.data.m[,2]                 = indexN
#   survey.data.m[,3]                 = month
#   survey.data.m[,4]                 = index
#   survey.data.m[,5]                 = cv
#   survey.data.m[,6]                 = ess
#   survey.data.m[,7:(6+n_sizeb-1)]   = as.matrix(comp)
#   survey.data.m
#   return(survey.data.m)  
#   
# }

survey.data.gen <- function(survey_periods,survey_months,n_survey,survey.obj,survey.cv,survey.ess,add.error)
{
  nrow           = 0
  nrow.temp      = 0
  year           = c()
  indexN         = c()
  month          = c()
  index          = c()
  cv             = c()
  ess            = c()
  comp           = data.frame()
  n_sizeb        = lapply(survey.obj, dim)[[1]][2]
  
  
  for (i in 1:n_survey)
  {
    nrow.temp    = lapply(survey.obj, dim)[[i]][1]
    nrow         = nrow + nrow.temp
    year.temp    = seq(survey_periods[i,1],survey_periods[i,2],1)
    year         = c(year,year.temp)
    indexN.temp  = rep(i,nrow.temp)
    indexN       = c(indexN,indexN.temp)
    month.temp   = rep(survey_months[i],nrow.temp)
    month        = c(month,month.temp)
    index.temp   = survey.obj[[i]][,1]
    index        = c(index,index.temp)
    cv.temp      = survey.cv[[i]]
    cv           = c(cv,cv.temp)
    ess.temp     = survey.ess[,i,]
    ess          = c(ess,ess.temp)
    comp.temp    = survey.obj[[i]][,2:n_sizeb]
    comp         = rbind(comp,comp.temp)
    
    if(add.error==TRUE)
    {
      index      = add.error.lognormal(index,cv)
      comp       = add.error.multino(comp,ess)
    } 
    
  }
  nrow
  year
  indexN
  month
  index
  cv
  ess
  comp
  
  
  survey.data.m  = matrix(-1,nrow=length(index),ncol=(6+n_sizeb-1))
  
  survey.data.m[,1]                 = year
  survey.data.m[,2]                 = indexN
  survey.data.m[,3]                 = month
  survey.data.m[,4]                 = index
  survey.data.m[,5]                 = cv
  survey.data.m[,6]                 = ess
  survey.data.m[,7:(6+n_sizeb-1)]   = as.matrix(comp)
  
  return(survey.data.m)  
  
}

############### retro pattern ##########################

run.retro <- function(n_peels)
  {
    fileNames=c("NSLSAP01.STD")
    file.name.full=paste("NSLSAP01-",year_start,"-",year_end,sep='')
    file.copy(from=paste('NSLSAP01.rep'),to=paste(paste(file.name.full),".rep",sep=''),overwrite=TRUE)
    file.copy(from=paste('NSLSAP01.par'),to=paste(paste(file.name.full),".par",sep=''),overwrite=TRUE)
    file.copy(from=paste('NSLSAP01.cor'),to=paste(paste(file.name.full),".cor",sep=''),overwrite=TRUE)

    SSB_retro = list(); SSB_rho_temp = c(); SSB_rho  = c()
    R_retro   = list(); R_rho_temp   = c(); R_rho    = c()
    F_retro   = list(); F_rho_temp   = c(); F_rho    = c()
    retros    = c()

    index.year = year_end-year_start+1-n_peels
    SSB_retro_ref <- read.admb(paste(file.name.full))$Spawning_stock_Biomass_input
    R_retro_ref   <- read.admb(paste(file.name.full))$recruitment_Pred
    F1_retro_ref   <- read.admb(paste(file.name.full))$Fishing_Mortality[1,]
    F2_retro_ref   <- read.admb(paste(file.name.full))$Fishing_Mortality[2,]

      for (i in 1:n_peels)
        {
          year_end_retro = year_end - i
          unlink("NSLSAP01.STD")

          source("/Users/jiecao/Desktop/Google Drive/work/Postdoc research/Natural_mortality/SSAP-SIM/write.control.r")
          runAD("NSLSAP01",verbose=F)

          if (file.exists(fileNames)==TRUE)
            {
              retro_name = paste0("SSAP01-",year_start,"-",year_end_retro,sep='')
              retro_name_rep = paste0(retro_name,".rep",sep="")
              retro_name_par = paste0(retro_name,".par",sep="")
              retro_name_cor = paste0(retro_name,".cor",sep="")
              file.copy(from=paste0('SSAP01.rep'),to=paste0(retro_name_rep),overwrite=TRUE)
              file.copy(from=paste0('SSAP01.par'),to=paste0(retro_name_par),overwrite=TRUE)
              file.copy(from=paste0('SSAP01.cor'),to=paste0(retro_name_cor),overwrite=TRUE)

              print(paste(retro_name,'converged'))
            }

          SSB_retro[[i]] <- read.admb(paste(retro_name))$Spawning_stock_Biomass_input
          R_retro[[i]]   <- read.admb(paste(retro_name))$recruitment_Pred
          F_retro[[i]]   <- read.admb(paste(retro_name))$Fishing_Mortality

          SSB_rho_temp[i] = (SSB_retro[[i]][index.year] - SSB_retro_ref[index.year])/SSB_retro_ref[index.year]
          R_rho_temp[i] = (R_retro[[i]][index.year] - R_retro_ref[index.year])/R_retro_ref[index.year]
          F_rho_temp[i] = (F_retro[[i]][index.year] - F_retro_ref[index.year])/F_retro_ref[index.year]

        }

    SSB_rho    = mean(SSB_rho_temp)
    R_rho      = mean(R_rho_temp)
    F_rho      = mean(F_rho_temp)

    retros     = c(SSB_rho,R_rho,F_rho)

    return(retros)
    #return(list(SSB_rho=SSB_rho,R_rho=R_rho,F_rho=F_rho))
  }


