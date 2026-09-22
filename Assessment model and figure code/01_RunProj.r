
# short-term projection model
#\code{run_projection} run the projection model with specified inputs 
#@param report, a .rep file contains outputs from the assessment model
#@param par_est, a .std file contains standard error estimates 
#@param select_id, a vector - number of selectivity time blocks, which one for trawl and which one for trap used in projection
#@param M, natrual mortality 
#@param n_proj_year, number of years (projection)
#@param N_rep, number of iterations
#@param F_trawl_proj, a vector of fishing mortality (trawl)
#@param F_trap_proj, a vector of fishing mortality (trap)
#@param growth_id, a vector of indicators (growth block)
#@param SFrac, fraction of year prior to spawning

run_projection <- function(report = report, par_est = par_est, select_id = select_id, M = M,
					n_proj_year = n_proj_year, N_rep = N_rep, F_trawl_proj = F_trawl_proj, F_trap_proj = F_trap_proj,
					R_proj = R_proj, growth_id = growth_id, SFrac = 1/3, R_rand = T, R_logsd = R_logsd){
		
		#selectivity
		n_sel_blocks = select_id[1]
		sel_trawl_blocl = select_id[2]
		sel_trap_block = select_id[3]
		
		GTM_block = growth_id
								
		recr_lbin <- report$Recrut_Sizebin_ratio
		recr_lbin_n <- length(recr_lbin[recr_lbin>0])
		lbin_n <- length(report$Size_bins)-1
		lbin <- report$Size_bins[1:lbin_n]+0.5

		PropFemale = report$Maturity_length[report$N_year,]
		weight = report$Weight_length[report$N_year,]
		
		GTM_1 = report$Growth[((GTM_block[1]-1)*lbin_n+1):((GTM_block[1]-1)*lbin_n+lbin_n),]
		GTM_2 = report$Growth[((GTM_block[2]-1)*lbin_n+1):((GTM_block[2]-1)*lbin_n+lbin_n),]

		sel_trawl_pars = c(par_est$est[(sel_trawl_blocl-1)*2+1][[1]],par_est$est[(sel_trawl_blocl-1)*2+2][[1]])
		sel_trap_pars = c(par_est$est[(sel_trap_block-1)*2+1][[1]],par_est$est[(sel_trap_block-1)*2+2][[1]])
		sel_trawl_pars_sd = c(par_est$std[(sel_trawl_blocl-1)*2+1][[1]],par_est$std[(sel_trawl_blocl-1)*2+2][[1]])
		sel_trap_pars_sd = c(par_est$std[(sel_trap_block-1)*2+1][[1]],par_est$std[(sel_trap_block-1)*2+2][[1]])

		N1_proj = matrix(NA,nrow=lbin_n,ncol=2)
		N1_proj[(recr_lbin_n+1):lbin_n,1] = par_est$est$MvNASPJD_1
		N1_proj[(recr_lbin_n+1):lbin_n,2] = par_est$std$MvNASPJD_1
		N1_proj_temp = c(); 

		SSB = matrix(NA,nrow=N_rep,ncol=n_proj_year)
		C_trawl_tot=matrix(NA,nrow=N_rep,ncol=n_proj_year)
		C_trap_tot=matrix(NA,nrow=N_rep,ncol=n_proj_year)
		N_proj_3d=array(NA,dim=c(N_rep,lbin_n,n_proj_year))

		for (y in 1:n_proj_year){
 			 for(i in 1:N_rep)
  					{
    						# define recruitment
    						if(R_rand){
    							R_proj_temp = rlnorm(1,log(R_proj[y])-R_logsd^2,R_logsd)
    						}else{
    							R_proj_temp = R_proj[y]
    							}
    						R_at_size_proj = R_proj_temp*recr_lbin
    
    						# randomly draw selectivety parameter
   							sel_trawl_pars_t=c(rnorm(1,sel_trawl_pars[1],sel_trawl_pars_sd[1]),
   													rnorm(1,sel_trawl_pars[2],sel_trawl_pars_sd[2]))
    						sel_trap_pars_t=c(rnorm(1,sel_trap_pars[1],sel_trap_pars_sd[1]),
    												rnorm(1,sel_trap_pars[2],sel_trap_pars_sd[2]))
    						sel_trawl = 1/(1+exp(-sel_trawl_pars_t[2]*(lbin-sel_trawl_pars_t[1])))
    						sel_trawl = sel_trawl/max(sel_trawl)
    						sel_trap  = 1/(1+exp(-sel_trap_pars_t[2]*(lbin-sel_trap_pars_t[1])));
    						sel_trap = sel_trap/max(sel_trap)
    						
    						# define F  
    						F_trawl = sel_trawl*F_trawl_proj[y]
    						F_trap = sel_trap*F_trap_proj[y]
    
    						# randomly draw N1_proj
    						if (y==1){
      								N1_proj_temp = R_at_size_proj
      										for (ii in (recr_lbin_n+1):lbin_n)
        											N1_proj_temp[ii] = rnorm(1,N1_proj[ii,1],N1_proj[ii,2])
      												N_proj_3d[i,,y] = N1_proj_temp
    								}else{
      									N_proj_3d[i,,y] = (((N_proj_3d[sample(1:N_rep,1,FALSE),,y-1] * 
      									exp(-sel_trawl*F_trawl_proj[y-1]-sel_trap*F_trap_proj[y-1]-M))%*%GTM_1) *
      									exp(-M))%*%GTM_2 + R_at_size_proj					
    									}
    
    		
    						C_trawl = N_proj_3d[i,,y]*(1-exp(-F_trawl-F_trap-M))*(F_trawl/(F_trawl+F_trap+M))
    						C_trap  = N_proj_3d[i,,y]*(1-exp(-F_trawl-F_trap-M))*(F_trap/(F_trawl+F_trap+M))
    
    						SSB[i,y] = sum(N_proj_3d[i,,y]*exp(-SFrac*(-F_trawl-F_trap-M))*PropFemale*weight)
    						C_trawl_tot[i,y]  = sum(C_trawl*weight)
    						C_trap_tot[i,y]   = sum(C_trap*weight)
  						}
					}

				#visulization 
				library("fanplot")
				par(mfrow=c(2,1))
				par(mar = c(4, 4.5, 2, 0.5))
				Year=seq(report$BeginYear,report$EndYear,1)
					if(report$SexAtSizeLamda>0)
						{SSB_estimates = report$Spawning_stock_Biomass
							}else
  								SSB_estimates=report$Spawning_stock_Biomass_input

				plot(Year,SSB_estimates,xlab="Year",ylab="Spawning biomass (mt)",type='l',
				xlim=c(report$BeginYear-1,report$EndYear+n_proj_year))
				fan(data = SSB, start = report$EndYear+1,rlab=F)

				plot(NULL,xlab="Year",ylab="Spawning biomass (mt)",
				type='l',xlim=c(report$EndYear,report$EndYear+n_proj_year),ylim=c(min(SSB),max(SSB)))
				fan(data = SSB, start = report$EndYear+1, ln=c(5, 50, 95))

				Table1 = array(NA, dim=c(n_proj_year,3), 
						dimnames=list(c(as.character(c(seq(1,n_proj_year,1)))),c("SSB","Catch_trawl","Catch_trap")) )

				Table1[,1] = apply(SSB,2,median)
				Table1[,2] = apply(C_trawl_tot,2,median)
				Table1[,3] = apply(C_trap_tot,2,median)

				print(Table1)
				
						
		}

