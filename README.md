# tb-model-clearance
Code from "Modeling the impact of case finding for tuberculosis: The role of infection dynamics"

All of the code is organized to run a on a high-performance computing cluster (HPCC). 
The order to run each script is as follows:

#PART 1: MODEL CALIBRATION

1. Run run_IMIS_pt1.R. This script randomly samples from prior distributions and calculates modeled outputs and likelihoods.
This was run with 400 arrays on an HPCC with 2500 samples each (B=250)
We ran this for each model variation (scenario="base" and scenario="clearance")
In the main analysis, we set version="omega_increase", so that incidence declines were induced by increases in the treatment initiation rate (omega)
In sensitivity analysis, we set version="beta_decline" (declines in the transmission rate), or, for certain analyses, version="" (no declines in incidence)
The prior bounds of each model parameters are specified in CSV files in the "data" folder. 

2. Run run_IMIS_pt2.R. This script conducts targeted sampling in the areas of high likelihood idenfied by the random sampling in run_IMIS_pt1.R
This script was run with 20 arrays on an HPCC with XX samples each (B=XX/10). Each of the 20 arrays combines YY files from part 1.
We ran this for each model variation (scenario="base" and scenario="clearance")
In the main analysis, we set version="omega_increase" and targets_set=3, so that incidence declines were induced by increases in the treatment initiation rate (omega)
In sensitivity analysis, we set version="beta_decline" (declines in the transmission rate) or "", and used alternative targets_set (1, 2, or 4).

3. Run combine_IMIS.R. This script combines all of the outputs from run_IMIS_pt2.R and conducts likelihood-based weighted sampling
The outputs from this script are files of posterior parameter samples, along with their corresponding modeled outputs and starting and ending compartment sizes.

#PART 2: INTERVENTION PROJECTIONS (to be run after model calibration)

1. Run run_projections.R. 
This script projects outputs over 10 years with no interventions, mass case finding (ACF), and ACF with TPT.
Intervention-related parameters are specified in the script. 
This was run (for each model variation and analysis) on an HPCC with 20 arrays, each covering 100 or fewer unique posterior parameter samples.

2. Run combine_projections.R
This script combines all 20 sets of output from run_projections.R and calculates additional outputs, such as averted incidence and mortality. 

#HELPER SCRIPTS

1. model_functions.R 
This script contains the main model ODE function and other functions used in calibration and intervention projections

2. model_setup.R
This script sets up starting parameter values and a starting population/compartment sizes. 
It is used in both calibration and intervention projections.
