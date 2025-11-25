#ODE model function
model <- function(t, X, params){
  with(as.list(params),{
    DX <- array(0, length(as.vector(matind))+n_extra_states)
    
    names.nat.hist <- rownames(matind)
    names.type <- c('Sus','Res') #placeholder for a "resistant" subgroup - not used (pop set to 0)
    names.age <- c('0-14','15+')
    
    # births
    for(js in 1:length(names.type)){ 
      DX[matind['U',js, 1]] <- DX[matind['U',js, 1]] + birth.rate*prop.births[js]*sum(X[as.vector(matind)]) 
    }
    # mortality
    for(jn in 1:length(names.nat.hist)){
      for(js in 1:length(names.type)){
        for(ja in 1:length(names.age)){
          DX[matind[jn,js,ja]] <- DX[matind[jn,js,ja]] - mu[ja]*X[matind[jn,js,ja]]
          if(names.nat.hist[jn] %in% paste0('AA', 1:n.active)){
            DX[matind[jn,js,ja]] <- DX[matind[jn,js,ja]] - mu.AA*X[matind[jn,js,ja]]
          }
          else if	(names.nat.hist[jn] %in% paste0('SA', 1:n.active)){
            DX[matind[jn,js,ja]] <- DX[matind[jn,js,ja]] - mu.SA*X[matind[jn,js,ja]]
          }		
        }
      }
    }
    
    # aging
    for(jn in 1:length(names.nat.hist)){
      for(js in 1:length(names.type)){
        DX[matind[jn,js,'0-14']] <- DX[matind[jn,js,'0-14']] - 1/15*X[matind[jn,js,'0-14']] 
        DX[matind[jn,js,'15+']] <- DX[matind[jn,js,'15+']] + 1/15*X[matind[jn,js,'0-14']] 
      }
    }
    
    # transmissions
    foi <- c(0,0)
    tb <- rbind(which(names.nat.hist %in% paste0('AA', 1:n.active)), 
                which(names.nat.hist %in% paste0('SA', 1:n.active)))
    for(jn in c(1:2)){ #loop over infectious groups (pre-care-seeking AA and care-seeking SA)
      for(ja in 1:length(names.age)){ #loop over age groups
        foi[1] <- foi[1] + beta.stage[jn]*beta.age[ja]*sum(X[as.vector(matind[tb[jn,],1,ja])]) 		
        foi[2] <- foi[2] + beta.stage[jn]*beta.age[ja]*sum(X[as.vector(matind[tb[jn,],2,ja])]) 
      }
    }
    lambda <- c(0,0) #rate at which a given type is infected - based on rel. FOIs, decline in beta over time, beta by type, type mixing, beta itself
    lambda[1] <- rel.inf.type[1]*beta*exp(-t*decline.beta)*(foi[1] + sigma*foi[2])/(sum(X[as.vector(matind[,1,])]) + sigma*sum(X[as.vector(matind[,2,])]))
    lambda[2] <- rel.inf.type[2]*beta*exp(-t*decline.beta)*(foi[2] + sigma*foi[1])/(sum(X[as.vector(matind[,2,])]) + sigma*sum(X[as.vector(matind[,1,])]))
    
    for(js in 1:length(names.type)){
      for(ja in 1:length(names.age)){
        
        # new infections
        DX[matind['U',js,ja]] <- DX[matind['U',js,ja]] - rel.inf.age[ja]*lambda[js]*X[matind['U',js,ja]]
        DX[matind['L1',js,ja]] <- DX[matind['L1',js,ja]] + rel.inf.age[ja]*lambda[js]*X[matind['U',js,ja]]
        
        # reinfections (from recovered, cleared, and non-early latent state)
        DX[matind['R',js,ja]] <- DX[matind['R',js,ja]] - xi*rel.inf.age[ja]*lambda[js]*X[matind['R',js,ja]]
        DX[matind['L1',js,ja]] <- DX[matind['L1',js,ja]] + xi*rel.inf.age[ja]*lambda[js]*X[matind['R',js,ja]]
        
        DX[matind['C',js,ja]] <- DX[matind['C',js,ja]] - xi*rel.inf.age[ja]*lambda[js]*X[matind['C',js,ja]]
        DX[matind['L1',js,ja]] <- DX[matind['L1',js,ja]] + xi*rel.inf.age[ja]*lambda[js]*X[matind['C',js,ja]]
        
        DX[matind[paste0('L', 2:n.ltbi),js,ja]] <- DX[matind[paste0('L', 2:n.ltbi),js,ja]] - xi*rel.inf.age[ja]*lambda[js]*X[matind[paste0('L', 2:n.ltbi),js,ja]]
        DX[matind['L1',js,ja]] <- DX[[matind['L1',js,ja]]] + sum(xi*rel.inf.age[ja]*lambda[js]*X[matind[paste0('L', 2:n.ltbi),js,ja]])
        
        # rapid progression
        DX[matind['L1',js,ja]] <- DX[matind['L1',js,ja]] - p[ja]*X[matind['L1',js,ja]]
        DX[matind[paste0('AA', 1:n.active),js,ja]] <- DX[matind[paste0('AA', 1:n.active),js,ja]] + p[ja]*X[matind['L1',js,ja]]
        
        # stabilization
        for(i in 1:(n.ltbi-1)) {
          DX[matind[paste0('L', i), js, ja]] <-  DX[matind[paste0('L', i), js, ja]] - s[i]*X[matind[paste0('L', i), js, ja]]
          DX[matind[paste0('L', i+1), js, ja]] <- DX[matind[paste0('L', i+1), js, ja]] + s[i]*X[matind[paste0('L', i), js, ja]]
        }
        
        #clearance
        for(i in 1:n.ltbi) {
          DX[matind[paste0('L', i), js, ja]] <- DX[matind[paste0('L', i), js, ja]] - c[i]*X[matind[paste0('L', i), js, ja]]
          DX[matind['C',js,ja]] <- DX[matind['C',js,ja]] + c[i]*X[matind[paste0('L', i), js, ja]]
        }
        
        # reactivation
        DX[matind[paste0('L', 2:n.ltbi),js,ja]] <- DX[matind[paste0('L', 2:n.ltbi),js,ja]] - phi[,ja]*X[matind[paste0('L', 2:n.ltbi),js,ja]]
        DX[matind[paste0('AA', 1:n.active),js,ja]] <- DX[matind[paste0('AA', 1:n.active),js,ja]] + sum(phi[,ja]*X[matind[paste0('L', 2:n.ltbi),js,ja]])
        
        # progression from pre-care-seeking to care-seeking
        DX[matind[paste0('AA', 1:n.active),js,ja]] <- DX[matind[paste0('AA', 1:n.active),js,ja]] - r1*X[matind[paste0('AA', 1:n.active),js,ja]]
        DX[matind[paste0('SA', 1:n.active),js,ja]] <- DX[matind[paste0('SA', 1:n.active),js,ja]] + r1*X[matind[paste0('AA', 1:n.active),js,ja]]
        
        # regression from care-seeking back to pre-care-seeking
        DX[matind[paste0('SA', 1:n.active),js,ja]] <- DX[matind[paste0('SA', 1:n.active),js,ja]] - r2*X[matind[paste0('SA', 1:n.active),js,ja]]
        DX[matind[paste0('AA', 1:n.active),js,ja]] <- DX[matind[paste0('AA', 1:n.active),js,ja]] + r2*X[matind[paste0('SA', 1:n.active),js,ja]]
        
        # self-resolution
        DX[matind[paste0('AA', 1:n.active),js,ja]] <- DX[matind[paste0('AA', 1:n.active),js,ja]] - w*X[matind[paste0('AA', 1:n.active),js,ja]]
        DX[matind['R',js,ja]] <- DX[[matind['R',js,ja]]] + sum(w*X[matind[paste0('AA', 1:n.active),js,ja]])
        
        # treatment
        DX[matind[paste0('SA', 1:n.active),js,ja]] <- DX[matind[paste0('SA', 1:n.active),js,ja]] - 
          k*omega*(1+t*increase.omega)*X[matind[paste0('SA', 1:n.active),js,ja]]
        DX[matind['R',js,ja]] <- DX[[matind['R',js,ja]]] + 
          sum(k*omega*(1+t*increase.omega)*X[matind[paste0('SA', 1:n.active),js,ja]])
        
      }
    }
    
    #calculate additional tracked outputs
    for(ja in 1:length(names.age)){
      # new infections (among uninfected only), by type
      DX[length(as.vector(matind))+1] <- DX[length(as.vector(matind))+1] + lambda[1]*rel.inf.age[ja]*X[matind['U',1,ja]] 
      DX[length(as.vector(matind))+2] <- DX[length(as.vector(matind))+2] + lambda[2]*rel.inf.age[ja]*X[matind['U',2,ja]] 
      
      # new cases, by type
      DX[length(as.vector(matind))+3] <- DX[length(as.vector(matind))+3] + p[ja]*X[matind['L1',1,ja]] + sum(phi[,ja]*X[matind[paste0('L', 2:n.ltbi),1,ja]])
      DX[length(as.vector(matind))+4] <- DX[length(as.vector(matind))+4] + p[ja]*X[matind['L1',2,ja]] + sum(phi[,ja]*X[matind[paste0('L', 2:n.ltbi),2,ja]])
      
      # new TB-related deaths, by type
      DX[length(as.vector(matind))+5] <- DX[[length(as.vector(matind))+5]] + sum(mu.AA*X[matind[paste0('AA', 1:n.active),1,ja]]) +  sum(mu.SA*X[matind[paste0('SA', 1:n.active),1,ja]])
      DX[length(as.vector(matind))+6] <- DX[[length(as.vector(matind))+6]] + sum(mu.AA*X[matind[paste0('AA', 1:n.active),2,ja]]) +  sum(mu.SA*X[matind[paste0('SA', 1:n.active),2,ja]])
      
      # new treatment initiations, by type
      DX[length(as.vector(matind))+7] <- DX[[length(as.vector(matind))+7]] + 
        sum(omega*(1+t*increase.omega)*X[matind[paste0('SA', 1:n.active),1,ja]])
      DX[length(as.vector(matind))+8] <- DX[[length(as.vector(matind))+8]] + 
        sum(omega*(1+t*increase.omega)*X[matind[paste0('SA', 1:n.active),2,ja]])
      
      # new cases via recent transmission, by type
      DX[length(as.vector(matind))+9] <- DX[length(as.vector(matind))+9] + p[ja]*X[matind['L1',1,ja]] 
      DX[length(as.vector(matind))+10] <- DX[length(as.vector(matind))+10] + p[ja]*X[matind['L1',2,ja]] 
      
      # new cases via reactivation, by type
      DX[length(as.vector(matind))+11] <- DX[[length(as.vector(matind))+11]] + sum(phi[,ja]*X[matind[paste0('L', 2:n.ltbi),1,ja]]) 
      DX[length(as.vector(matind))+12] <- DX[[length(as.vector(matind))+12]] + sum(phi[,ja]*X[matind[paste0('L', 2:n.ltbi),2,ja]])
    }
    
    # new cases among < 15s, by type
    DX[length(as.vector(matind))+13] <- DX[[length(as.vector(matind))+13]] + p[1]*X[[matind['L1',1,1]]] + sum(phi[,1]*X[matind[paste0('L', 2:n.ltbi),1,1]])
    DX[length(as.vector(matind))+14] <- DX[[length(as.vector(matind))+14]] + p[1]*X[[matind['L1',2,1]]] + sum(phi[,1]*X[matind[paste0('L', 2:n.ltbi),2,1]])
    
    # new infections (among uninfected only) among < 15s, by type
    DX[length(as.vector(matind))+15] <- DX[length(as.vector(matind))+15] + lambda[1]*rel.inf.age[1]*X[matind['U',1,1]] 
    DX[length(as.vector(matind))+16] <- DX[length(as.vector(matind))+16] + lambda[2]*rel.inf.age[1]*X[matind['U',2,1]] 
    
    list(DX)
  })
}

#wrapper to main model function 
ode.model <- function (h, tmax, pars, ystart){
  sol <- as.data.frame(
    ode(
      ystart,
      times=seq(0,tmax,by=h),
      func=model,
      method="lsoda",
      parms=pars
    )
  )  
}

#function to simulate progression among a cohort of early LTBI adults over time, without reinfection
#this is used to compare against the 2 progression targets
sim_progression <- function(t, X, params) {
  DX <- array(0, 1 + n.ltbi) #early, late, and progressed states
  DX[1] <- -1*X[1]*(params$p[[2]] + 0.014 + params$s[1] + params$c[1])
  DX[2] <- X[1]*params$s[1] - X[2]*(params$phi[1, 2] + 0.014 + params$c[2]) #0.014= historical mortality of approx 1/70, rather than model-calibrated mortality
  DX[3] <- X[1]*params$p[[2]] + X[2:n.ltbi]%*%params$phi[,2] #don't subtract out deaths because these ppl progressed before dying
  list(DX)
}

#calculate modeled outputs used in calibration
calc_outputs <- function(y, pars, nt) {
  ns <- max(matind) #number of compartments/states
  ll <- nrow(y) #number of timesteps
  pop <- y[,-c(1,(ns+2):ncol(y))] #only the compartment sizes
  
  #incidence per 100,000 - cases over 1 yr/pop midpoint
  inc_all <-  100000*diff((y[,1+ns+3] + y[,1+ns+4]), nt)/
    ((rowSums(pop[1:(ll-nt),]) + rowSums(pop[(nt+1):ll,]))/2)
  inc_annual <- inc_all[seq(from=1, to=(ll-nt), by=nt)]
  inc <- inc_annual[[length(inc_annual)]]
  
  #annual average decline in incidence
  inc_decline <- -1*((inc_annual[[length(inc_annual)]]/tail(inc_annual, 6)[[1]])^(1/5) - 1)
  
  #mortality to incidence ratio
  mort_to_inc <- tail(diff((y[,1+ns+5] + y[,1+ns+6]), nt)/diff((y[,1+ns+3] + y[,1+ns+4]), nt), 1)
  
  #percent of incidence among < 15 age group
  prop_inc_15 <- tail(diff(y[,1+ns+13] + y[,1+ns+14], nt)/diff(y[,1+ns+3] + y[,1+ns+4], nt), 1)

  #prevalence of active TB 
  prev <- 100000*sum(pop[ll, as.vector(matind[c(paste0('AA', 1:n.active),paste0('SA', 1:n.active)),,2])])/
    sum(pop[ll,as.vector(matind[,,2])])
  
  #percent of prevalence that is pre-care-seeking 
  prop_prev_pre <- sum(pop[ll, as.vector(matind[paste0('AA', 1:n.active),,2])])/
    sum(pop[ll, as.vector(matind[c(paste0('AA', 1:n.active),paste0('SA', 1:n.active)),,2])])
  
  #prevalence of latent TB infection - includes latent, recovered, specified % of cleared
  prev_ltbi <- (sum(pop[ll, as.vector(matind[c('R', paste0('L', n.ltbi)),,2])]) +
                  pars$p_clear_react*sum(pop[ll, as.vector(matind['C',,2])]))/
    sum(pop[ll,as.vector(matind[,,2])])

  #prevalence-to-notification ratio 
  prev_to_notif <- sum(pop[ll, as.vector(matind[c(paste0('AA', 1:n.active),paste0('SA', 1:n.active)),,2])])/
    diff((y[,1+ns+7] + y[,1+ns+8]), nt)[[ll-nt]]
  
  #ARI: annual risk of infection (all ages, to match adjusted data)
  ari <- tail(diff(y[, 1+ns+1] + y[, 1+ns+2], nt), 1)/sum(pop[ll-nt, as.vector(matind['U',,])])
  
  #progression risk over time - function of parameters - use ODE solver (needed anyway when more than 2 LTBI compartments)
  t_seq <- 0:10
  y_start <- c(1, rep(0, n.ltbi))
  prog_all <- ode(y=y_start, times=t_seq, func=sim_progression, parms=pars)
  prog_early <- prog_all[[3, 2+n.ltbi]]
  prog_cum <- prog_all[[11, 2+n.ltbi]]

  #adjust for scenario where nobody has TB so we don't get divide by 0 errors
  if(prev==0 & inc==0) {
    inc_decline <- 0
    mort_to_inc <- 0
    mort_to_pop <- 0
    prop_inc_15 <- 0
    prop_prev_pre <- 0
    prev_to_notif <- 0
    ari <- 0
    cdr <- 0
  }
  output <- list("inc"=inc,
                 "inc_decline"=inc_decline,
                 "mort_to_inc"=mort_to_inc,
                 "prop_inc_15"=prop_inc_15,
                 "prev"=prev,
                 "prop_prev_pre"=prop_prev_pre,
                 "prev_ltbi"=prev_ltbi, 
                 "prev_to_notif"=prev_to_notif,
                 "ari"=ari,
                 "prog_early"=prog_early,
                 "prog_cum"=prog_cum)
  return(output)
}

#calculate likelihoods
calc_like <- function(output) {
  #incidence per 100,000
  inc <- dnorm(output$inc, mean=210, sd=(244-178)/(1.96*2), log=T)
  
  #annual average decline in incidence
  inc_decline <- dnorm(output$inc_decline, mean=0.029, sd=0.01/1.96, log=T)

  #mortality to incidence ratio
  mort_to_inc_dist <- beta_params(mean=0.17, sigma=(0.21-0.14)/(2*1.96))
  mort_to_inc <- dbeta(output$mort_to_inc, shape1=mort_to_inc_dist$alpha, shape2=mort_to_inc_dist$beta, log=T)

  #prevalence among 15+ age group per 100,000
  prev <- dnorm(output$prev, mean=316, sd=(342-290)/(1.96*2), log=T)
  
  #percent of prevalence that is sub-clinical
  prop_prev_pre <- dbeta(output$prop_prev_pre, shape1=67.3, shape2=16.8, log=T)
  
  #prevalence of latent TB infection (among 15+ age group)
  prev_ltbi_dist <- beta_params(mean=0.314, sigma=(0.335-0.272)) #"intermediate" version of the target
  prev_ltbi <- dbeta(output$prev_ltbi, shape1=prev_ltbi_dist$alpha, shape2=prev_ltbi_dist$beta, log=T)
  
  #prevalence to notification ratio (among 15+ age group)
  prev_to_notif <- dnorm(output$prev_to_notif, mean=2.84, sd=(3.10-2.61)/(2*1.96), log=T)
  
  #2 progression targets (2 yrs and 10 yrs)
  prog_early <- dbeta(output$prog_early, shape1=(214+373)/10, shape2=(7744+2672)/10 - (214+373)/10, log=T) #4-7%
  prog_cum <- dbeta(output$prog_cum, shape1=7.5*5, shape2=92.5*5, log=T) #5-10%
  prog_late <- dbeta((output$prog_cum - output$prog_early), shape1=(106+43)/20, shape2=(7371+2472)/20 - (106+43)/20, log=T) #0.6-2.7%
  
  #5 targets version
  log_like5 <- inc + prev_ltbi + prog_early + prog_late + inc_decline
  like5 <- exp(log_like5)
  
  #8 targets version
  log_like8 <- log_like5 + prev + prev_to_notif + mort_to_inc
  like8 <- exp(log_like8)
  
  #9 targets version (USED IN MAIN ANALYSIS)
  log_like9 <- log_like8 + prop_prev_pre
  like9 <- exp(log_like9)
  
  out_like <- list("inc_like"=inc,
                   "inc_decline_like"=inc_decline,
                   "mort_to_inc_like"=mort_to_inc,
                   "prev_like"=prev,
                   "prop_prev_pre_like"=prop_prev_pre,
                   "prev_to_notif_like"=prev_to_notif,
                   "prev_ltbi_like"=prev_ltbi,
                   "prog_early_like"=prog_early,
                   "prog_cum_like"=prog_cum,
                   "prog_late_like"=prog_late,
                   "log_like5"=log_like5,
                   "likelihood5"=like5,
                   "log_like8"=log_like8,
                   "likelihood8"=like8,
                   "log_like9"=log_like9,
                   "likelihood9"=like9
                   )
  return(out_like)
}

#function used in calibration - required for IMIS package - sample parameters from prior distributions
sample.prior <- function(nsmpl, par.range=par_range) {
  priors <- randomLHS(nsmpl, ncol(par.range))
  for(j in 1:ncol(par.range)){
    priors[,j] <- as.numeric(par.range['low',j]) + as.numeric(par.range['diff',j])*priors[,j]
  }
  priors <- as.matrix(priors)
  colnames(priors) <- colnames(par.range)
  return(priors)
}

#used in calibration - required for IMIS package - calculate prior 
prior <- function(pars, par.range=par_range) {
  like <- sapply(names(par.range), function(x)
    dunif(pars[,x], as.numeric(par.range['low', x]), as.numeric(par.range['hi',x])), simplify=F, USE.NAMES=T)
  like <- bind_cols(like)
  like <- rowProds(as.matrix(like))
  return(like)
}

#function to run the model, calculate outputs, calculate likelihood
model_output_like <- function(pars_use, pars.baseline) {
  ## set model parameters - all annual
  pars <- pars.baseline
  
  # update directly from sampled parameters
  pars[names(pars_use)] <- pars_use
  pars$decline.beta <- 0
  pars$increase.omega <- 0

  # additional updates to some params
  pars$beta.stage <- c(pars$beta.stage.1, pars$beta.stage.2) #rel transmissibility of pre-care-seeking vs. care-seeking TB
  pars$phi[,ages_u15] <- pars$phi.2*pars$rel.p.age.1
  pars$phi[,ages_o15] <- pars$phi.2
  pars$p[ages_u15] <- pars$p.2*pars$rel.p.age.1
  pars$p[ages_o15] <- pars$p.2
  pars$c <- rep(pars$c.2, n.ltbi)

  ## run transience
  nt <- 1 #1 yr timestep
  h <- 1/nt
  tmax <- 500
  ll <- tmax/h + 1
  y.trans <- ode.model(h=h, tmax=tmax, pars=pars, ystart=c(as.vector(ystart),rep(0, n_extra_states)))
  
  ## population size at the end of transience
  N.trans <- rowSums(y.trans[ll,-c(1,max(matind):ncol(y.trans))])
  
  ## run final 21 years (during which incidence declines are induced)
  nt <- 1 #still 1 yr timestep
  h <- 1/nt
  tmax <- 21
  
  pars$decline.beta <- ifelse("decline.beta" %in% names(pars_use), pars_use[['decline.beta']], pars$decline.beta)
  pars$increase.omega <- ifelse("increase.omega" %in% names(pars_use), pars_use[['increase.omega']], pars$increase.omega)

  ## adjust pop size - assuming 0.5% growth rate in the population over the last 21 years, and the population is 2 mil at the end
  ystart.cont <- c(as.numeric(y.trans[ll,-c(1,(2+max(matind)):ncol(y.trans))])* 1/(exp(21*0.0054)/2000000)/N.trans,
                   rep(0, n_extra_states))
  ystart.cont <- round(ystart.cont)
  y <- ode.model(h=h, tmax=tmax, pars=pars, ystart=ystart.cont)
  
  ## calculate outputs and likelihood
  output <- calc_outputs(y, pars, nt)
  like <- calc_like(output)

  out_all <- list(y=y,
                  ystart=ystart.cont,
                  pars=pars,
                  outputs=c(output, like))
  return(out_all)
  
}

#version of IMIS package function after initial sampling
#this is adapted from the IMIS R package from Raftery & Bao (from https://rdrr.io/cran/IMIS/man/IMIS.html)
#see also: Raftery AE, Bao L. Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Biometrics 2010
IMIS <- function(B=1000, B.re=3000, number_k=100, D=0, targets_set="9"){
  B0 = B*10
  X_all = X_k = sample.prior(B0)				# Draw initial samples from the prior distribution
  if (is.vector(X_all))	Sig2_global = var(X_all)	# the prior covariance
  if (is.matrix(X_all))	Sig2_global = cov(X_all)	# the prior covariance
  stat_all = matrix(NA, 6, number_k)				# 6 diagnostic statistics at each iteration
  center_all = prior_all = out_all = ystart_all = yend_all = like_all = NULL			# centers of Gaussian components, prior densities, and likelihoods
  sigma_all = list()						# covariance matrices of Gaussian components
  if (D>=1)	option.opt = 1					# use optimizer
  if (D==0){	option.opt = 0; D=1	}			# NOT use optimizer
  
  #allow likelihood to be a variant,based on targets_set
  like_name <- paste0("likelihood", targets_set)
  
  for (k in 1:number_k ){
    
    ptm.like = proc.time()
    prior_all = c(prior_all, prior(X_k))		# Calculate the prior densities
    #KEEP TRACK OF ADDITIONAL OUTPUT BEYOND THE LIKELIHOODS
    out_like_k <- gen_out_like(X_k) #includes outputs & likelihoods, ystart, yend
    out_all <- rbind(out_all, cbind(out_like_k$out, "sample"=1:nrow(out_like_k$out), "round"=k))
    ystart_all <- rbind(ystart_all, cbind(out_like_k$ystart, "sample"=1:nrow(out_like_k$ystart), "round"=k))
    yend_all <- rbind(yend_all, cbind(out_like_k$yend, "sample"=rep(1:nrow(out_like_k$ystart), each=2), "round"=k))
    like_all = c(like_all, out_like_k$out[,like_name]) # Calculate the likelihoods
    ptm.use = (proc.time() - ptm.like)[3]
    if (k==1)	print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
    
    if (k==1)	envelop_all = prior_all			# envelop stores the sampling densities
    if (k>1)	envelop_all = apply( rbind(prior_all*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2))
    Weights = prior_all*like_all / envelop_all	# importance weight is determined by the posterior density divided by the sampling density
    stat_all[1,k] = log(mean(Weights))			# the raw marginal likelihood
    Weights = Weights / sum(Weights)			
    stat_all[2,k] = sum(1-(1-Weights)^B.re)		# the expected number of unique points
    stat_all[3,k] = max(Weights)				# the maximum weight
    stat_all[4,k] = 1/sum(Weights^2)			# the effictive sample size
    stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))	# the entropy relative to uniform
    stat_all[6,k] = var(Weights/mean(Weights))	# the variance of scaled weights
    if (k==1)	print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
    print(c(k, round(stat_all[1:4,k], 3)))
    
    if (k==1 & option.opt==1){
      if (is.matrix(X_all))	Sig2_global = cov(X_all[which(like_all>min(like_all)),])
      X_k = which_exclude = NULL					# exclude the neighborhood of the local optima 
      label_weight = sort(Weights, decreasing = TRUE, index=TRUE)
      which_remain = which(Weights>label_weight$x[B0]) 	# the candidate inputs for the starting points
      size_remain = length(which_remain)
      for (i in 1:D){
        important = NULL
        if (length(which_remain)>0)
          important = which_remain[which(Weights[which_remain]==max(Weights[which_remain]))]
        if (length(important)>1)	important = sample(important,1)	
        if (is.vector(X_all))	X_imp = X_all[important]
        if (is.matrix(X_all))	X_imp = X_all[important,]
        # Remove the selected input from candidates
        which_exclude = union( which_exclude, important )
        which_remain = setdiff(which_remain, which_exclude)
        posterior = function(theta){	-log(prior(theta))-log(likelihood(theta)) } 
        
        if (is.vector(X_all)){
          if (length(important)==0)	X_imp = center_all[1]
          optimizer = optim(X_imp, posterior, method="BFGS", hessian=TRUE, 
                            control=list(parscale=sqrt(Sig2_global)/10,maxit=5000))
          print(paste("maximum posterior=", round(-optimizer$value,2), ", likelihood=", round(log(likelihood(optimizer$par)),2), 
                      ", prior=", round(log(prior(optimizer$par)),2), ", time used=", round(ptm.use/60,2), "minutes, convergence=", optimizer$convergence))
          center_all = c(center_all, optimizer$par)
          sigma_all[[i]] = solve(optimizer$hessian)
          X_k = c(X_k, rnorm(B, optimizer$par, sqrt(sigma_all[[i]])) )			# Draw new samples
          distance_remain = abs(X_all[which_remain]-optimizer$par)
        }
        if (is.matrix(X_all)){	
          # The rough optimizer uses the Nelder-Mead algorithm.
          if (length(important)==0)	X_imp = center_all[1,]
          ptm.opt = proc.time()
          optimizer = optim(X_imp, posterior, method="Nelder-Mead", 
                            control=list(maxit=1000, parscale=sqrt(diag(Sig2_global))) )
          theta.NM = optimizer$par
          
          # The more efficient optimizer uses the BFGS algorithm 
          optimizer = optim(theta.NM, posterior, method="BFGS", hessian=TRUE,
                            control=list(parscale=sqrt(diag(Sig2_global)), maxit=1000))
          ptm.use = (proc.time() - ptm.opt)[3]
          print(paste("maximum posterior=", round(-optimizer$value,2), ", likelihood=", round(log(likelihood(optimizer$par)),2), 
                      ", prior=", round(log(prior(optimizer$par)),2), ", time used=", round(ptm.use/60,2), "minutes, convergence=", optimizer$convergence))
          center_all = rbind(center_all, optimizer$par)						# the center of new samples
          if (min(eigen(optimizer$hessian)$values)>0)
            sigma_all[[i]] = solve(optimizer$hessian)						# the covariance of new samples
          if (min(eigen(optimizer$hessian)$values)<=0){						# If the hessian matrix is not positive definite, we define the covariance as following
            eigen.values = eigen(optimizer$hessian)$values
            eigen.values[which(eigen.values<0)] = 0
            hessian = eigen(optimizer$hessian)$vectors %*% diag(eigen.values) %*% t(eigen(optimizer$hessian)$vectors)
            sigma_all[[i]] = solve(hessian + diag(1/diag(Sig2_global)) )
          }
          X_k = rbind(X_k, rmvnorm(B, optimizer$par, sigma_all[[i]]) )			# Draw new samples
          distance_remain = mahalanobis(X_all[which_remain,], optimizer$par, diag(diag(Sig2_global)) )
        }
        # exclude the neighborhood of the local optima 
        label_dist = sort(distance_remain, decreasing = FALSE, index=TRUE)
        which_exclude = union( which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
        which_remain = setdiff(which_remain, which_exclude)
      }
      if (is.matrix(X_all))	X_all = rbind(X_all, X_k)
      if (is.vector(X_all))	X_all = c(X_all, X_k)
    }
    
    if (k>1 | option.opt==0){
      important = which(Weights == max(Weights))
      if (length(important)>1)	important = important[1]
      if (is.matrix(X_all))	X_imp = X_all[important,]				# X_imp is the maximum weight input
      if (is.vector(X_all))	X_imp = X_all[important]
      if (is.matrix(X_all))	center_all = rbind(center_all, X_imp)
      if (is.vector(X_all))	center_all = c(center_all, X_imp)
      if (is.matrix(X_all))	distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)) )
      if (is.vector(X_all))	distance_all = abs(X_all-X_imp)			# Calculate the distances to X_imp
      label_nr = sort(distance_all, decreasing = FALSE, index=TRUE)		# Sort the distances
      which_var = label_nr$ix[1:B]								# Pick B inputs for covariance calculation
      if (is.matrix(X_all))	Sig2 = cov.wt(X_all[which_var,], wt = Weights[which_var]+1/length(Weights), cor = FALSE, center = X_imp, method = "unbias")$cov
      if (is.vector(X_all)){
        Weights_var = Weights[which_var]+1/length(X_all)
        Weights_var = Weights_var/sum(Weights_var)
        Sig2 = (X_all[which_var]-X_imp)^2 %*% Weights_var
      }
      #fix Sig2 if it's not positive semi-definite (required for rtmvnorm, but not rmvnorm)
      Sig2_orig <- Sig2
      evalues <- eigen(Sig2)$values
      evectors <- eigen(Sig2)$vectors
      index_to_fix <- which(evalues<1e-15)
      Sig2 <- Sig2+evectors[,index_to_fix]%*%t(evectors[,index_to_fix])*
        (.Machine$double.eps-min(evalues[index_to_fix],0))
      sigma_all[[D+k-1]] = Sig2
      if (is.matrix(X_all))	X_k = rtmvnorm(B, X_imp, Sig2, lower=rep(0, length(X_imp))) #UPDATED SO NO NEGATIVE SAMPLES
      if (is.vector(X_all))	X_k = rnorm(B, X_imp, sqrt(Sig2))			# Draw new samples
      if (is.matrix(X_all))	X_all = rbind(X_all, X_k)
      if (is.vector(X_all))	X_all = c(X_all, X_k)
      #adding colnames to X_k
      colnames(X_k) <- colnames(X_all)
    }
    
    if (k==1){
      gaussian_all = matrix(NA, D, B0+D*B)
      for (i in 1:D){
        if (is.matrix(X_all))	gaussian_all[i,] = dmvnorm(X_all, center_all[i,], sigma_all[[i]])
        if (is.vector(X_all))	gaussian_all[i,] = dnorm(X_all, center_all[i], sqrt(sigma_all[[i]]))
      }
    }
    if (k>1){
      if (is.vector(X_all))	gaussian_new = matrix(0, D+k-1, length(X_all) )
      if (is.matrix(X_all))	gaussian_new = matrix(0, D+k-1, dim(X_all)[1] )
      if (is.matrix(X_all)){
        gaussian_new[1:(D+k-2), 1:(dim(X_all)[1]-B)] = gaussian_all
        gaussian_new[D+k-1, ] = dmvnorm(X_all, X_imp, sigma_all[[D+k-1]])
        for (j in 1:(D+k-2))	gaussian_new[j, (dim(X_all)[1]-B+1):dim(X_all)[1] ] = dmvnorm(X_k, center_all[j,], sigma_all[[j]])
      }
      if (is.vector(X_all)){
        gaussian_new[1:(D+k-2), 1:(length(X_all)-B)] = gaussian_all
        gaussian_new[D+k-1, ] = dnorm(X_all, X_imp, sqrt(sigma_all[[D+k-1]]))
        for (j in 1:(D+k-2))	gaussian_new[j, (length(X_all)-B+1):length(X_all) ] = dnorm(X_k, center_all[j], sqrt(sigma_all[[j]]))
      }
      gaussian_all = gaussian_new
    }
    if (stat_all[2,k] > (1-exp(-1))*B.re)	break
  } # end of k
  
  #sample and round get added to X_all
  X_all <- cbind(X_all, 
                 "sample"=c(1:(B*10), rep(1:B, (number_k))),
                 "round"=c(rep(1, B*10), rep(2:(number_k+1), each=B)))
  
  nonzero = which(Weights>0)
  which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
  if (is.matrix(X_all))	resample_X = X_all[which_X,]
  if (is.vector(X_all))	resample_X = X_all[which_X]
  
  #MORE OUTPUTS GET SAVED
  resample_out <- out_all[which_X, ]
  return(list(stat = t(stat_all), resample = resample_X, center = center_all,
              out = resample_out, params_all=X_all, out_all=out_all, weights=Weights,
              ystart_all=ystart_all, yend_all=yend_all, gaussian_all=gaussian_all))
} # end of IMIS

#wrapper function used in projections
model_output_proj <- function(pars_use, pars.baseline, ystart, tmax=21) {
  ## set model parameters - all annual
  pars <- pars.baseline
  
  # update directly from sampled parameters
  pars[names(pars_use)] <- pars_use
  
  # additional updates to some params
  pars$beta.stage <- c(pars$beta.stage.1, pars$beta.stage.2) #rel transmissibility of pre-care-seeking vs. care-seeking TB
  pars$phi[,ages_u15] <- pars$phi.2*pars$rel.p.age.1
  pars$phi[,ages_o15] <- pars$phi.2
  pars$p[ages_u15] <- pars$p.2*pars$rel.p.age.1
  pars$p[ages_o15] <- pars$p.2
  pars$c <- rep(pars$c.2, n.ltbi)
  
  ## run for 20 years, monthly timestep
  nt <- 12
  h <- 1/nt
  
  y <- ode.model(h=h, tmax=tmax, pars=pars, ystart=ystart)
  
  ## calculate outputs 
  output <- calc_outputs(y, pars, nt)
  
  out_all <- list(y=y,
                  pars=pars,
                  outputs=output)
  return(out_all)
}

#function for calculating additional outputs over time, used in projections
calc_outputs_time <- function(y, nt) {
  time <- 1:nrow(y)
  ns <- max(matind)
  ll <- nrow(y)
  
  #population
  pop <- rowSums(y[,-c(1,(ns+2):ncol(y))]) #only the compartment sizes
  pop_u15 <- rowSums(y[,1+as.vector(matind[,,ages_u15])])
  pop_o15 <- rowSums(y[,1+as.vector(matind[,,ages_o15])])
  pop_u <- rowSums(y[,1+as.vector(matind['U',,])])
  pop_u_u15 <- rowSums(y[,1+as.vector(matind['U',,ages_u15])])
  
  #cases and deaths
  cases <- c(0, diff((y[,1+ns+3] + y[,1+ns+4]), 1))
  deaths <- c(0, diff((y[,1+ns+5] + y[,1+ns+6]), 1))
  
  #incidence
  inc <-  c(rep(0, nt), 100000*diff((y[,1+ns+3] + y[,1+ns+4]), nt)/
              ((pop[1:(ll-nt)] + pop[(nt+1):ll])/2))
  inc_u15 <- c(rep(0, nt), 100000*diff((y[,1+ns+13] + y[,1+ns+14]), nt)/
                 ((pop_u15[1:(ll-nt)] + pop_u15[(nt+1):ll])/2))
  inc_o15 <- c(rep(0, nt), 100000*diff((y[,1+ns+3] + y[,1+ns+4] - y[,1+ns+13] - y[,1+ns+14]), nt)/
                 ((pop_o15[1:(ll-nt)] + pop_o15[(nt+1):ll])/2))
  inc_recent <- c(rep(0, nt), 100000*diff(y[,1+ns+9] + y[,1+ns+10], nt)/
                    ((pop[1:(ll-nt)] + pop[(nt+1):ll])/2))
  inc_react <- c(rep(0, nt), 100000*diff(y[,1+ns+11] + y[,1+ns+12], nt)/
                   ((pop[1:(ll-nt)] + pop[(nt+1):ll])/2))
  
  #ARI
  ari_all <- c(rep(NA, nt), diff((y[,1+ns+1] + y[,1+ns+2]), nt)/
                 ((pop_u[1:(ll-nt)] + pop_u[(nt+1):ll])/2))
  ari_u15 <- c(rep(NA, nt), diff((y[,1+ns+15] + y[,1+ns+16]), nt)/
                 ((pop_u_u15[1:(ll-nt)] + pop_u_u15[(nt+1):ll])/2))
  
  #mortality
  mort <- c(rep(NA, nt), 100000*diff((y[,1+ns+5] + y[,1+ns+6]), nt)/
              ((pop[1:(ll-nt)] + pop[(nt+1):ll])/2))
   
  
  #prevalence
  prev <- 100000*rowSums(y[, 1+as.vector(matind[c(paste0('AA', 1:n.active),paste0('SA', 1:n.active)),,])])/pop
  prev_u15 <- 100000*rowSums(y[, 1+as.vector(matind[c(paste0('AA', 1:n.active),paste0('SA', 1:n.active)),,ages_u15])])/pop_u15
  prev_o15 <- 100000*rowSums(y[, 1+as.vector(matind[c(paste0('AA', 1:n.active),paste0('SA', 1:n.active)),,ages_o15])])/pop_o15
  prev_sub <- 100000*rowSums(y[, 1+as.vector(matind[paste0('AA', 1:n.active),,])])/pop

  #LTBI prevalence
  prev_ltbi <- (rowSums(y[, 1+as.vector(matind[c(paste0('L', 1:n.ltbi), 'R'),,])]) +
    p_clear_react*rowSums(y[, 1+as.vector(matind['C',,])]))/pop
  prev_ltbi_u15 <- (rowSums(y[, 1+as.vector(matind[c(paste0('L', 1:n.ltbi), 'R'),,ages_u15])]) +
                  p_clear_react*rowSums(y[, 1+as.vector(matind['C',,ages_u15])]))/pop_u15
  prev_ltbi_o15 <- (rowSums(y[, 1+as.vector(matind[c(paste0('L', 1:n.ltbi), 'R'),,ages_o15])]) +
                  p_clear_react*rowSums(y[, 1+as.vector(matind['C',,ages_o15])]))/pop_o15
  prev_ltbi_early <- rowSums(y[, 1+as.vector(matind[c(paste0('L', 1)),,])])/pop

  #notifications
  notif <- c(rep(NA, nt), 100000*diff((y[,1+max(matind)+7] + y[,1+max(matind)+8]), nt)/
               ((pop[1:(ll-nt)] + pop[(nt+1):ll])/2))
  
  out <- data.frame(time, pop, pop_u15, pop_o15, cases, deaths,
                    inc, inc_u15, inc_o15, inc_recent, inc_react,
                    ari_all, ari_u15, mort, 
                    prev, prev_u15, prev_o15, prev_sub,
                    prev_ltbi, prev_ltbi_u15, prev_ltbi_o15, prev_ltbi_early, 
                    notif)
  return(out)
  
}

#additional output: avg duration of active disease among adults, by active1 vs. active2
#note that this changes over time if omega or k increase, which this function doesn't account for
calc_duration_active <- function(pars) {
  #asymptomatic active (AA) - now pre-care-seeking active: outflows = non-TB deaths, TB deaths, spontaneous resolution, symptom progression
  durs_AA <- 1/(pars$mu[[2]] + pars$mu.AA + pars$w + pars$r1)
  #symptomatic active (AA) - now care-seeking active: outflows = non-TB deaths, TB deaths, symptom regression (if modeling), treatment
  durs_SA <- 1/(pars$mu[[2]] + pars$mu.SA + pars$r2 + pars$omega*pars$k)
  #total active
  durs_active <- durs_AA + durs_SA*(pars$r1/(pars$mu[[2]] + pars$mu.AA + pars$w + pars$r1))
  
  return(list("AA"=durs_AA,
              "SA"=durs_SA,
              "Total"=durs_active))
}
