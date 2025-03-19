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
        
        # reinfections 
        DX[matind['R',js,ja]] <- DX[matind['R',js,ja]] - xi*rel.inf.age[ja]*lambda[js]*X[matind['R',js,ja]]
        DX[matind['L1',js,ja]] <- DX[matind['L1',js,ja]] + xi*rel.inf.age[ja]*lambda[js]*X[matind['R',js,ja]]
        
        DX[matind[paste0('L', 2:n.ltbi),js,ja]] <- DX[matind[paste0('L', 2:n.ltbi),js,ja]] - xi*rel.inf.age[ja]*lambda[js]*X[matind[paste0('L', 2:n.ltbi),js,ja]]
        DX[matind['L1',js,ja]] <- DX[[matind['L1',js,ja]]] + sum(xi*rel.inf.age[ja]*lambda[js]*X[matind[paste0('L', 2:n.ltbi),js,ja]])
        
        # rapid progression
        DX[matind['L1',js,ja]] <- DX[matind['L1',js,ja]] - p[ja]*X[matind['L1',js,ja]]
        DX[matind[paste0('AA', 1:n.active),js,ja]] <- DX[matind[paste0('AA', 1:n.active),js,ja]] + p[ja]*X[matind['L1',js,ja]]
        
        # stabilization (and clearance, if modeling it)
        for(i in 1:(n.ltbi-1)) {
          DX[matind[paste0('L', i), js, ja]] <-  DX[matind[paste0('L', i), js, ja]] - s[i]*X[matind[paste0('L', i), js, ja]]
          DX[matind[paste0('L', i+1), js, ja]] <- DX[matind[paste0('L', i+1), js, ja]] + s[i]*X[matind[paste0('L', i), js, ja]]
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
    lsoda(
      ystart,
      times=seq(0,tmax,by=h),
      func=model,
      parms=pars
    )
  )  
}

#function to simulate progression among a cohort of early LTBI adults over time, without reinfection
#this is used to compare against the 2 progression targets
sim_progression <- function(t, X, params) {
  DX <- array(0, 1 + n.ltbi)
  DX[1] <- -1*X[1]*(params$p[[2]] + 0.014 + params$s[1])
  if(n.ltbi>2) {
    for(i in 2:(n.ltbi-1)) {
      DX[i] <- params$s[i-1]*X[i-1] - 1*X[i]*(params$phi[i-1, 2] + 0.014 + params$s[i]) #0.014= historical mortality of approx 1/70, rather than model-calibrated mortality
    }
  }
  DX[n.ltbi] <- params$s[n.ltbi-1]*X[(n.ltbi-1)] - 1*X[n.ltbi]*(params$phi[(n.ltbi-1),2] + 0.014)
  DX[n.ltbi+1] <- X[1]*params$p[[2]] + X[2:n.ltbi]%*%params$phi[,2]
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
  
  #prevalence of latent TB infection 
  prev_ltbi <- sum(pop[ll, as.vector(matind[c('R', paste0('L', 1:2)),,2])])/sum(pop[ll,as.vector(matind[,,2])])

  #prevalence-to-notification ratio 
  prev_to_notif <- sum(pop[ll, as.vector(matind[c(paste0('AA', 1:n.active),paste0('SA', 1:n.active)),,2])])/
    diff((y[,1+ns+7] + y[,1+ns+8]), nt)[[ll-nt]]
  
  #ARI: annual risk of infection (all ages, to match adjusted data)
  ari <- tail(diff(y[, 1+ns+1] + y[, 1+ns+2], nt), 1)/sum(pop[ll-nt, as.vector(matind['U',,])])
  
  #progression risk over time per 1000 - function of parameters - use ODE solver (when more than 2 LTBI compartments)
  t_seq <- 0:10
  y_start <- c(1000, rep(0, n.ltbi))
  prog_all <- ode(y=y_start, times=t_seq, func=sim_progression, parms=pars)
  prog_early <- prog_all[[3, 2+n.ltbi]]/sum(y_start)
  prog_cum <- prog_all[[11, 2+n.ltbi]]/sum(y_start)

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
  
  #annual average decline in incidence - 2 different versions (narrow used in sensitivity analysis)
  inc_decline <- dnorm(output$inc_decline, mean=0.029, sd=0.01/1.96, log=T)
  inc_decline_narrow <- dnorm(output$inc_decline, mean=0.029, sd=0.0005/1.96, log=T)
  
  #mortality to incidence ratio
  mort_to_inc_dist <- beta_params(mean=0.17, sigma=(0.21-0.14)/(2*1.96))
  mort_to_inc <- dbeta(output$mort_to_inc, shape1=mort_to_inc_dist$alpha, shape2=mort_to_inc_dist$beta, log=T)
  
  #percent of incidence among < 15 age group
  prop_inc_15_dist <- beta_params(mean=0.12, sigma=(0.15-0.09)/(2*1.96))
  prop_inc_15 <- dbeta(output$prop_inc_15, shape1=prop_inc_15_dist$alpha, shape2=prop_inc_15_dist$beta, log=T)
  
  #prevalence among 15+ age group per 100,000
  prev <- dnorm(output$prev, mean=316, sd=(342-290)/(1.96*2), log=T)
  
  #percent of prevalence that is subclinical
  prop_prev_pre <- dbeta(output$prop_prev_pre, shape1=67.3, shape2=16.8, log=T)
  
  #prevalence of latent TB infection (among 15+ age group)
  prev_ltbi_dist <- beta_params(mean=0.314, sigma=(0.335-0.272)) #"intermediate" version of the target
  prev_ltbi <- dbeta(output$prev_ltbi, shape1=prev_ltbi_dist$alpha, shape2=prev_ltbi_dist$beta, log=T)
  
  #prevalence to notification ratio (among 15+ age group)
  prev_to_notif <- dnorm(output$prev_to_notif, mean=2.84, sd=(3.10-2.61)/(2*1.96), log=T)
  
  #ARI (annual risk of TB infection)
  ari <- dbeta(output$ari, shape1=3.9, shape2=191.1, log=T) #2% [0.05-4.3%]
  
  #2 progression targets (2 yrs and 10 yrs)
  prog_early <- dbeta(output$prog_early, shape1=(214+373)/10, shape2=(7744+2672)/10 - (214+373)/10, log=T) #4-7%
  prog_cum <- dbeta(output$prog_cum, shape1=7.5*5, shape2=92.5*5, log=T) #5-10%
  
  #sensitivity analysis version - only 4 targets
  log_like1 <- inc + prev_ltbi + prog_early + prog_cum
  like1 <- exp(log_like1)
  
  #sensitivity analysis version - leaving out incidence decline
  log_like2 <- log_like1 + mort_to_inc + prop_inc_15 + prev + prop_prev_pre + prev_to_notif + ari
  like2 <- exp(log_like2)
  
  #MAIN ANALSYIS VERSION
  log_like3 <- log_like2 + inc_decline
  like3 <- exp(log_like3)
  
  #sensitivity analysis version: narrower incidence decline
  log_like4 <- log_like2 + inc_decline_narrow
  like4 <- exp(log_like4)
  
  out_like <- list("inc_like"=inc,
                   "inc_decline_like"=inc_decline,
                   "inc_decline_narrow_like"=inc_decline_narrow,
                   "mort_to_inc_like"=mort_to_inc,
                   "prop_inc_15_like"=prop_inc_15,
                   "prev_like"=prev,
                   "prop_prev_pre_like"=prop_prev_pre,
                   "prev_to_notif_like"=prev_to_notif,
                   "ari_like"=ari,
                   "prev_ltbi_like"=prev_ltbi,
                   "prog_early_like"=prog_early,
                   "prog_cum_like"=prog_cum,
                   "log_like1"=log_like1,
                   "likelihood1"=like1,
                   "log_like2"=log_like2,
                   "likelihood2"=like2,
                   "log_like3"=log_like3,
                   "likelihood3"=like3,
                   "log_like4"=log_like4,
                   "likelihood4"=like4
                   )
  return(out_like)
}

#function used in calibration - required for IMIS package - sample parameters from prior distributions
sample.prior <- function(nsmpl, par.range) {
  priors <- randomLHS(nsmpl, ncol(par.range))
  for(j in 1:ncol(par.range)){
    priors[,j] <- as.numeric(par.range['low',j]) + as.numeric(par.range['diff',j])*priors[,j]
  }
  priors <- as.matrix(priors)
  colnames(priors) <- colnames(par.range)
  return(priors)
}

#used in calibration - required for IMIS package - calculate prior 
prior <- function(pars, par.range) {
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
  pars$phi[1,ages_u15] <- pars$phi.2*pars$rel.p.age.1
  pars$phi[1,ages_o15] <- pars$phi.2
  pars$p[ages_u15] <- pars$p.2*pars$rel.p.age.1
  pars$p[ages_o15] <- pars$p.2
  #only adjust first row of phi (progression during L2) - if clearance variation, 2nd row is fixed (at 0)
  
  if(scenario=="clearance") {
    pars$s[n.ltbi-1] <- pars_use["s2"]
  }
  
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

#just the random sampling and calculation of the likelihood - pt 1 of IMIS (can run more samples in parallel)
samples_like <- function(B, par.range, pars.baseline) {
  #sample params
  B0 = B * 10
  names_vary <- names(par.range)[par.range["diff",]!=0]
  
  X_all = X_k = sample.prior(B0, par.range[,names_vary])
  count <- 0
  
  #calculate prior
  prior_k <- prior(X_k, par.range[,names_vary])
  
  #run model, calculate outputs and likelihood
  out_all <- data.frame()
  pars_all <- data.frame()
  ystart_all <- data.frame()
  yend_all <- data.frame()
  for(js in 1:nrow(X_k)) {
    print(js)
    pars_use <- X_k[js,]
    out_js <- model_output_like(pars_use, pars.baseline)
    out_all <- rbind(out_all, unlist(out_js$outputs))
    pars_all <- rbind(pars_all, unlist(out_js$pars))
    ystart_all <- rbind(ystart_all, unlist(out_js$ystart))
    yend_all <- rbind(yend_all, tail(out_js$y, 2))
  }
  
  new_names <- str_replace(names(unlist(out_js$outputs)), "\\..*$", "")
  names(out_all) <- new_names
  names(pars_all) <- names(unlist(out_js$pars))
  names(ystart_all) <- paste0("X", 1:ncol(ystart_all))
  names(yend_all) <- c("t", names(ystart_all))
  
  out_all <- cbind(out_all, "prior"=prior_k)
  out_all <- cbind(out_all, "array"=part, "sample"=1:nrow(out_all))
  pars_all <- cbind(pars_all, "array"=part, "sample"=1:nrow(out_all))
  ystart_all <- cbind(ystart_all, "array"=part, "sample"=1:nrow(out_all))
  yend_all <- cbind(yend_all, "array"=part, "samples"=rep(1:nrow(out_all), each=2))
  
  #add priors as a column
  
  return(list("out"=out_all,
              "pars"=pars_all,
              "ystart"=ystart_all,
              "yend"=yend_all))
}

#version of IMIS package function after initial sampling (used in run_IMIS_pt2)
#this is adapted from the IMIS R package from Raftery & Bao (from https://rdrr.io/cran/IMIS/man/IMIS.html)
#see also: Raftery AE, Bao L. Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Biometrics 2010
IMIS_rounds <- function(B, B.re, number_k, par.range, pars.baseline, 
                        pars_all, outputs_all, ystart_all, yend_all, targets_set) {
  names_vary <- names(par.range)[par.range["diff",]!=0]
  
  X_all = X_k = pars_all[, names_vary]
  Sig2_global = cov(X_all)
  stat_all = matrix(NA, 6, number_k)
  center_all = prior_all = out_all = log_like_all = like_all = NULL
  sigma_all = list()
  for (k in 1:number_k) {
    ptm.like = proc.time()
    if(k==1) { #skip some steps
      prior_k <- outputs_all$prior
      prior_all = c(prior_all, prior_k)
      out_all <- outputs_all
      out_round <- outputs_all
      B0 <- nrow(out_all)
    } else {
      prior_k <- prior(X_k, par.range[,names_vary])
      prior_all = c(prior_all, prior_k)
      out_round <- data.frame()
      ystart_round <- data.frame()
      yend_round <- list()
      for(js in 1:B) {
        print(js)
        if(prior_k[js]==0) {
          out_js <- rep(-Inf, ncol(outputs_all)-4)
          ystart_js <- rep(0, ncol(ystart_all)-3)
          yend_js <- as.data.frame(rbind(rep(0, ncol(yend_all)- 3),
                                         rep(0, ncol(yend_all)- 3)))
        } else {
          pars_use <- X_k[js,]
          out_js_all <- model_output_like(pars_use, pars.baseline)
          out_js <- unlist(out_js_all$out)
          ystart_js <- out_js_all$ystart
          yend_js <- tail(out_js_all$y, 2)
        }
        out_round <- rbind(out_round, out_js)
        ystart_round <- rbind(ystart_round, ystart_js)
        yend_round[[js]] <- yend_js #so we don't get col name errors
      }
      #obtain correct column names for outputs (in case first entry has prior=0 and we don't run the model for it)
      out_names <- names(model_output_like(center_all[1,], pars.baseline)$outputs)
      names(out_round) <- out_names
      out_round <- cbind(out_round, "prior"=prior_k, "array"=part, "sample"=1:B, "round"=k)
      out_all <- rbind(out_all, out_round)
      
      ystart_round <- cbind(ystart_round, "array"=part, "sample"=1:B, "round"=k)
      names(ystart_round) <- names(ystart_all)
      ystart_all <- rbind(ystart_all, ystart_round)
      
      yend_round <- rbindlist(yend_round, use.names=FALSE)
      yend_round <- cbind(yend_round, "array"=part, "sample"=rep(1:B, each=2), "round"=k)
      names(yend_round) <- names(yend_all)
      yend_all <- rbind(yend_all, yend_round)
    }
    #specify which likelihood to use based on which set of targets we are using
    if(targets_set==1) {
      log_like_all = c(log_like_all, out_round$log_like1)
    } else if(targets_set==2) {
      log_like_all = c(log_like_all, out_round$log_like2)
    } else if(targets_set==3) {
      log_like_all = c(log_like_all, out_round$log_like3)
    } else if(targets_set==4) {
      log_like_all = c(log_like_all, out_round$log_like4)
    } 
    
    scale_fac = 650 - max(log_like_all) #to minimize rounding to 0 - doesn't affect weights otherwise
    like_all <- exp(log_like_all + scale_fac)
    ptm.use = (proc.time() - ptm.like)[3]
    if (k == 1) 
      envelop_all = prior_all
    if (k > 1) 
      envelop_all = apply(rbind(prior_all * B0/B, gaussian_all), 
                          2, sum)/(B0/B + k - 1)
    Weights = prior_all * like_all/envelop_all
    Weights[is.na(Weights)] <- 0
    stat_all[1, k] = log(mean(Weights))
    Weights = Weights/sum(Weights)
    stat_all[2, k] = sum(1 - (1 - Weights)^B.re)
    stat_all[3, k] = max(Weights)
    stat_all[4, k] = 1/sum(Weights^2)
    stat_all[5, k] = -sum(Weights * log(Weights), na.rm = TRUE)/log(length(Weights))
    stat_all[6, k] = var(Weights/mean(Weights))
    if (k == 1) 
      print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
    print(c(k, round(stat_all[1:4, k], 3)))
    important = which(Weights == max(Weights))
    if (length(important) > 1) 
      important = important[1]
    X_imp = unlist(X_all[important, ])
    center_all = rbind(center_all, X_imp)
    distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)))
    label_nr = sort(distance_all, decreasing = FALSE,  index = TRUE)
    which_var = label_nr$ix[1:100] #take "closest" B samples (based on mahalanobis distance)
    Sig2 = cov.wt(X_all[which_var, ], wt = Weights[which_var] + 
                    1/length(Weights), cor = FALSE, center = X_imp, 
                  method = "unbias")$cov
    sigma_all[[k]] = Sig2
    X_k = rmvnorm(B, X_imp, Sig2)
    X_all = rbind(X_all, X_k)
    if (k == 1) {
      gaussian_all = matrix(NA, 1, B0 + 1 * B)
      gaussian_all[1, ] = dmvnorm(X_all, center_all[1, ], sigma_all[[1]])
    }
    if (k > 1) {
      gaussian_new = matrix(0, k, dim(X_all)[1])
      gaussian_new[1:(k - 1), 1:(dim(X_all)[1] - B)] = gaussian_all
      gaussian_new[k, ] = dmvnorm(X_all, X_imp, sigma_all[[k]])
      for (j in 1:(k - 1)) 
        gaussian_new[j, (dim(X_all)[1] - B + 1):dim(X_all)[1]] = dmvnorm(X_k, center_all[j, ], sigma_all[[j]])
      gaussian_all = gaussian_new
      
    }
    out_temp <- list(stat = t(stat_all), center = center_all,
                     params_all=X_all, out_all=out_all, weights=Weights,
                     ystart_all=ystart_all, yend_all=yend_all)
    save(out_temp, file=paste(path_out, 'imis_temp',k, "_", part, '_', scenario, '.rda',sep=''))
  }
  nonzero = which(Weights > 0)
  which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
  
  #add info to output
  X_all <- head(X_all, -B) #last B samples don't get used in IMIS
  X_all <- cbind(X_all, array=c(pars_all$array, rep(part, B*(number_k-1))), 
                 sample=c(pars_all$sample, rep(1:B, (number_k-1))),
                 round=c(pars_all$round, rep(2:number_k, each=B)))
  
  resample_X = X_all[which_X, ]
  resample_out = out_all[which_X,]
  
  return(list(stat = t(stat_all), resample = resample_X, center = center_all,
              out = resample_out, params_all=X_all, out_all=out_all, weights=Weights,
              ystart_all=ystart_all, yend_all=yend_all))
}

#wrapper function used in projections
model_output_proj <- function(pars_use, pars.baseline, ystart, tmax=21) {
  ## set model parameters - all annual
  pars <- pars.baseline
  
  # update directly from sampled parameters
  pars[names(pars_use)] <- pars_use
  
  # additional updates to some params
  pars$beta.stage <- c(pars$beta.stage.1, pars$beta.stage.2) #rel transmissibility of pre-care-seeking vs. care-seeking TB
  pars$phi[1,ages_u15] <- pars$phi.2*pars$rel.p.age.1
  pars$phi[1,ages_o15] <- pars$phi.2
  pars$p[ages_u15] <- pars$p.2*pars$rel.p.age.1
  pars$p[ages_o15] <- pars$p.2
  #only adjust first row of phi (progression during L2) - if progressions variation, 2nd row is fixed (at 0)
  
   if(scenario=="clearance") {
    pars$s[n.ltbi-1] <- pars_use["s2"]
  }
  pars$s <- unlist(pars$s)
  
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
  prev_ltbi <- rowSums(y[, 1+as.vector(matind[c(paste0('L', 1:2)),,])])/pop
  prev_ltbi_u15 <- rowSums(y[, 1+as.vector(matind[c(paste0('L', 1:2)),,ages_u15])])/pop_u15
  prev_ltbi_o15 <- rowSums(y[, 1+as.vector(matind[c(paste0('L', 1:2)),,ages_o15])])/pop_o15
  prev_ltbi_early <- rowSums(y[, 1+as.vector(matind[c(paste0('L', 1)),,])])/pop

  #notifications
  notif <- c(rep(NA, nt), 100000*diff((y[,1+max(matind)+7] + y[,1+max(matind)+8]), nt)/
               ((pop[1:(ll-nt)] + pop[(nt+1):ll])/2))
  
  out <- data.frame(time, pop, pop_u15, pop_o15, 
                    inc, inc_u15, inc_o15, inc_recent, inc_react,
                    ari_all, ari_u15, mort, 
                    prev, prev_u15, prev_o15, prev_sub,
                    prev_ltbi, prev_ltbi_u15, prev_ltbi_o15, prev_ltbi_early, 
                    notif)
  return(out)
  
}
