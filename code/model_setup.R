n.active <- 1 #1 set of active states (each set includes subclinical and symptomatic)
n.ltbi <- ifelse(scenario=="clearance", 3, 2) #cleared & recovered is included in the model as a third LTBI state with no progression risk
n.ages <- 2 #2 age groups
n.types <- 2 #resistors and susceptibles - not using (resistor pops to 0)

## dummy array with indices
matind <- array(1:((2 + n.active*2 + n.ltbi)*n.types*n.ages),c((2 + n.active*2 + n.ltbi),n.types,n.ages))
dimnames(matind)[[1]] <- c('U', paste0('L', 1:n.ltbi), paste0('AA', 1:n.active) , paste0('SA', 1:n.active), 'R')
dimnames(matind)[[2]] <- c('Sus','Res')
dimnames(matind)[[3]] <- c('0-14','15+')
ages_u15 <- 1
ages_o15 <- 2

## starter params - most are replaced during simulation 

# birth and mortality rates
birth.rate <- 0.0165 # birth rate
mu.1 <- 0.00270 #background mortality
mu.2 <- 0.0169 #background mortality
mu <- c(mu.1, mu.2)   
mu.AA <- 0 #asymptomatic TB mortality
mu.SA <- 0.2 #symptomatic TB mortality

# aging rates
aging <- c(1/15, 0)

# population
N <- 2000000
prop.susceptible <- 1 ## no resistant pop
prop.resistant <- 1-prop.susceptible 
prop.births <- c(prop.susceptible, prop.resistant) #everyone born susceptible since no resistant pop

# transmission rates
beta <- 20

# rates of change in parameters
decline.beta <- 0 #transmission rate
increase.omega <- 0 #treatment rate

# age-specific transmission scalar (relative infectiousness)
beta.age.1 <- 0 #children aren't infectious
beta.age.2 <- 1
beta.age <- c(beta.age.1, beta.age.2) 

# stage-specific transmission scalar (relative infectiousness)
beta.stage.1 <- 1
beta.stage.2 <- 1
beta.stage <- c(beta.stage.1, beta.stage.2) #relative infectiousness of subclinical, symptomatic

# mixing (between resistant and susceptible types - not using)
sigma <- 1

# relative risk of infection by age
rel.inf.age <- 1 #doesn't vary
rel.inf.age <- rep(rel.inf.age, 2)
# relative risk of infection by type
rel.inf.type <- c(1, 0) #resistant can't get infected (not using)

# protection against reinfection
xi <- 0.4

# rapid progression by age
rel.p.age.1 <- 0.725
p_adults <- (0.0555+0.0476)/2
p <- c(rel.p.age.1*p_adults, p_adults)

# stabilization - starting value
s <- rep(0.5, n.ltbi-1)
if(scenario=="clearance") {
  s[n.ltbi-1] <- 1/(10 - 1/0.5) #10 yrs on average before clearance, subtract off first 2 years when in early LTBI
}

# reactivation rate - by age and LTBI type
phi_full <- colMeans(rbind(c(24.5, 9.5, 4.2, 1.1, 1.1, 2.2, 1.2, 0, 0)/1000,
                      c(2.5, 3.7, 2.9, 2.8, 1.8, 1.9, 1.3, 0.7, 0.5)/1000))
phi_full <- cbind(replicate(1, phi_full*rel.p.age.1), replicate(1, phi_full))
phi <- array(data=0, dim=c(n.ltbi-1, 2))
phi[,] <- phi_full[1:(n.ltbi-1),]
if(scenario=="clearance") {
  phi[n.ltbi-1, ] <- 0 #no progression from 3rd LTBI state ("cleared" state)
}

# progression to care-seeking
r1 <- 2.05

# regression to pre-care-seeking
r2 <- 0

# self-resolution (only among asymptomatic AA)
w <- 0.85

# treatment rates (only among symptomatic SA)
omega <- 0.8

# treatment success proportions
k <- 0.9

pars <- list(birth.rate=birth.rate,
             prop.births=prop.births,
             mu=mu,
             aging=aging,
             mu.AA=mu.AA,
             mu.SA=mu.SA,
             beta=beta,
             decline.beta=decline.beta,
             increase.omega=increase.omega,
             beta.stage.1=beta.stage.1,
             beta.stage.2=beta.stage.2,
             beta.stage=beta.stage,
             beta.age=beta.age,
             rel.inf.age=rel.inf.age,
             rel.inf.type=rel.inf.type,
             sigma=sigma,
             xi=xi,
             p=p,
             s=s,
             phi=phi,
             r1=r1,
             r2=r2,
             w=w,
             k=k,
             omega=omega
)
			
ystart <- array(1:((2 + n.active*2 + n.ltbi)*n.types*n.ages),c((2 + n.active*2 + n.ltbi),n.types,n.ages))


ystart.epi <- matrix(c(c(0.3, 0.1, rep(0.5/(n.ltbi-1), n.ltbi-1), rep(0.002/(n.active*2), n.active*2), 0.1-0.002),
                       c((rel.inf.type[1]) + 0.3*rel.inf.type[2], 
                         rel.inf.type[2]*c(0.1,  rep(0.5/(n.ltbi-1), n.ltbi-1), 
                                                  rep(0.002/(n.active*2), n.active*2), 0.1-0.002))), 
                     nrow=n.types, byrow=T)
ystart.type <- c(prop.susceptible,prop.resistant)
ystart.age <- c(0.37,0.63) 
			
for(je in 1:(2 + n.active*2 + n.ltbi)){
	for(js in 1:n.types){
		for(ja in 1:n.ages){
			ystart[je,js,ja] <- N*ystart.epi[js,je]*ystart.type[js]*ystart.age[ja]
		}
	}
}

#matrices = ages
#columns = types (susceptible vs resistant)
#rows = TB states

