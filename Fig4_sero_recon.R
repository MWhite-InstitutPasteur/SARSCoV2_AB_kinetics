library(binom)

library(MASS)
library(ROCR)
library(randomForest)
library(pROC)


N_tree <- 10000

N_rep      <- 200
train_prop <- 2/3

pos_target_spec <- 0.99

###############################################
###############################################
##          ##                               ##
##   ####   ##  ####    ####  ######  ####   ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ##  ##  ##  ## ######   ##   ######  ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##   ####   ##  ####   ##  ##   ##   ##  ##  ##
##          ##                               ##
###############################################
###############################################





###################
###             ###
###  Spike IgG  ###
###             ###
###################

MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\Spike_IgG\\Spike_IgG_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)


load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\Spike_IgG.RData")

N_part <- dim(Spike_IgG_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(Spike_IgG_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

Spike_IgG_mod_val <- exp( log(Spike_IgG_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

Spike_IgG_mod_app <- exp( log(Spike_IgG_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

Spike_IgG_neg_data <- Spike_IgG_neg_data[[1]]

Spike_IgG_neg_data <- Spike_IgG_neg_data[ sample( length(Spike_IgG_neg_data), N_neg ) ]


###################
###             ###
###  Spike IgM  ###
###             ###
###################

MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\Spike_IgM\\Spike_IgM_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)



load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\Spike_IgM.RData")

N_part <- dim(Spike_IgM_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(Spike_IgM_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

Spike_IgM_mod_val <- exp( log(Spike_IgM_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

Spike_IgM_mod_app <- exp( log(Spike_IgM_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

Spike_IgM_neg_data <- Spike_IgM_neg_data[[1]]

Spike_IgM_neg_data <- Spike_IgM_neg_data[ sample( length(Spike_IgM_neg_data), N_neg ) ]



###################
###             ###
###  Spike IgA  ###
###             ###
###################

MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\Spike_IgA\\Spike_IgA_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)



load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\Spike_IgA.RData")

N_part <- dim(Spike_IgA_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(Spike_IgA_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

Spike_IgA_mod_val <- exp( log(Spike_IgA_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

Spike_IgA_mod_app <- exp( log(Spike_IgA_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

Spike_IgA_neg_data <- Spike_IgA_neg_data[[1]]

Spike_IgA_neg_data <- Spike_IgA_neg_data[ sample( length(Spike_IgA_neg_data), N_neg ) ]



###################
###             ###
###  RBD IgG    ###
###             ###
###################


MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\RBD_IgG\\RBD_IgG_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)




load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\RBD_IgG.RData")

N_part <- dim(RBD_IgG_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(RBD_IgG_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

RBD_IgG_mod_val <- exp( log(RBD_IgG_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

RBD_IgG_mod_app <- exp( log(RBD_IgG_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

RBD_IgG_neg_data <- RBD_IgG_neg_data[[1]]

RBD_IgG_neg_data <- RBD_IgG_neg_data[ sample( length(RBD_IgG_neg_data), N_neg ) ]


###################
###             ###
###  RBD IgM    ###
###             ###
###################


MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\RBD_IgM\\RBD_IgM_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)



load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\RBD_IgM.RData")

N_part <- dim(RBD_IgM_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(RBD_IgM_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

RBD_IgM_mod_val <- exp( log(RBD_IgM_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

RBD_IgM_mod_app <- exp( log(RBD_IgM_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

RBD_IgM_neg_data <- RBD_IgM_neg_data[[1]]

RBD_IgM_neg_data <- RBD_IgM_neg_data[ sample( length(RBD_IgM_neg_data), N_neg ) ]


###################
###             ###
###  RBD IgA    ###
###             ###
###################


MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\RBD_IgA\\RBD_IgA_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)



load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\RBD_IgA.RData")

N_part <- dim(RBD_IgA_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(RBD_IgA_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

RBD_IgA_mod_val <- exp( log(RBD_IgA_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

RBD_IgA_mod_app <- exp( log(RBD_IgA_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

RBD_IgA_neg_data <- RBD_IgA_neg_data[[1]]

RBD_IgA_neg_data <- RBD_IgA_neg_data[ sample( length(RBD_IgA_neg_data), N_neg ) ]


###################
###             ###
###  NP IgG     ###
###             ###
###################


MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\NP_IgG\\NP_IgG_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)



load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\NP_IgG.RData")

N_part <- dim(NP_IgG_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(NP_IgG_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

NP_IgG_mod_val <- exp( log(NP_IgG_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

NP_IgG_mod_app <- exp( log(NP_IgG_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

NP_IgG_neg_data <- NP_IgG_neg_data[[1]]

NP_IgG_neg_data <- NP_IgG_neg_data[ sample( length(NP_IgG_neg_data), N_neg ) ]




###################
###             ###
###  NP IgM     ###
###             ###
###################



MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\NP_IgM\\NP_IgM_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)



load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\NP_IgM.RData")

N_part <- dim(NP_IgM_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(NP_IgM_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

NP_IgM_mod_val <- exp( log(NP_IgM_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

NP_IgM_mod_app <- exp( log(NP_IgM_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

NP_IgM_neg_data <- NP_IgM_neg_data[[1]]

NP_IgM_neg_data <- NP_IgM_neg_data[ sample( length(NP_IgM_neg_data), N_neg ) ]





###################
###             ###
###  NP IgA     ###
###             ###
###################


MCMC <- read.table("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\NP_IgA\\NP_IgA_global_10m.txt")

N_cohort <- (ncol(MCMC) - 10 - 2)/5

sig_obs_vec <- MCMC[,4*N_cohort + 10 + 1]
sig_obs_vec <- sig_obs_vec[(0.5*length(sig_obs_vec)):length(sig_obs_vec)]
sig_obs <- median(sig_obs_vec)





load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\NP_IgA.RData")

N_part <- dim(NP_IgA_mod_Pasteur_draw)[1]  ## Number of positive participants
N_neg  <- length(NP_IgA_neg_data[[1]])     ## Number of negative controls

N_tt <- length(tt_plot)

NP_IgA_mod_val <- exp( log(NP_IgA_mod_Pasteur_draw[,,1]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

NP_IgA_mod_app <- exp( log(NP_IgA_mod_Pasteur_draw[,,2]) + 
                          matrix( rnorm(n=N_part*N_tt, mean=0, sd=sig_obs), nrow=N_part, ncol=N_tt ) )

NP_IgA_neg_data <- NP_IgA_neg_data[[1]]

NP_IgA_neg_data <- NP_IgA_neg_data[ sample( length(NP_IgA_neg_data), N_neg ) ]



#####################################
#####################################
##        ##                       ##  
##   ##   ##   ####  #### #     #  ##
##  ###   ##  ##      ##  ##   ##  ##
##   ##   ##   ####   ##  #######  ##
##   ##   ##      ##  ##  ## # ##  ## 
##  ####  ##   ####  #### ##   ##  ##
##        ##                       ##
#####################################
#####################################

########################################
##       ##  Simulate a data set for  ##
##  1.1  ##  algorithm VALidation     ##
##       ##                           ##
########################################
  
N_neg_val  <- 300  ## Simulated number of negative controls
N_pos1_val <- 300  ## Simulated number of infections (0-3 months)
N_pos2_val <- 300  ## Simulated number of infections (4-6 months)
N_pos3_val <- 300  ## Simulated number of infections (7-12 months)


N_tot_val <- N_neg_val + N_pos1_val + N_pos2_val + N_pos3_val

N_pos_val <- N_pos1_val + N_pos2_val + N_pos3_val

VAL <- matrix( NA, nrow = N_tot_val, ncol=10 )
colnames(VAL) <- c("group", "Spike_IgG", "RBD_IgG", "NP_IgG",
                            "Spike_IgM", "RBD_IgM", "NP_IgM", 
                            "Spike_IgA", "RBD_IgA", "NP_IgA" )

VAL <- as.data.frame(VAL)

VAL[1:N_neg_val,1]                                                                                            <- "neg"
VAL[(N_neg_val+1):(N_neg_val+N_pos1_val),1]                                                                   <- "pos1"
VAL[(N_neg_val+N_pos1_val+1):(N_neg_val+N_pos1_val+N_pos2_val),1]                                             <- "pos2"
VAL[(N_neg_val+N_pos1_val+N_pos2_val+1):(N_neg_val+N_pos1_val+N_pos2_val+N_pos3_val),1]                       <- "pos3"


###########################
## Simulate negatve controls

VAL[1:N_neg_val,2]  <- sample( x = Spike_IgG_neg_data, size=N_neg_val )
VAL[1:N_neg_val,3]  <- sample( x = RBD_IgG_neg_data,   size=N_neg_val )
VAL[1:N_neg_val,4]  <- sample( x = NP_IgG_neg_data,    size=N_neg_val )
VAL[1:N_neg_val,5]  <- sample( x = Spike_IgM_neg_data, size=N_neg_val )
VAL[1:N_neg_val,6]  <- sample( x = RBD_IgM_neg_data,   size=N_neg_val )
VAL[1:N_neg_val,7]  <- sample( x = NP_IgM_neg_data,    size=N_neg_val )
VAL[1:N_neg_val,8]  <- sample( x = Spike_IgA_neg_data, size=N_neg_val )
VAL[1:N_neg_val,9]  <- sample( x = RBD_IgA_neg_data,   size=N_neg_val )
VAL[1:N_neg_val,10] <- sample( x = RBD_IgA_neg_data,   size=N_neg_val )


###########################
## Simulate infections (0-3 months)

for( i in (N_neg_val+1):(N_neg_val+N_pos1_val) )
{ 
	tt_sample <- round(runif( 1, min=15, max=89 ))

	VAL[i,2]  <- Spike_IgG_mod_val[sample( nrow(Spike_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,3]  <- RBD_IgG_mod_val[sample( nrow(RBD_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,4]  <- NP_IgG_mod_val[sample( nrow(NP_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,5]  <- Spike_IgM_mod_val[sample( nrow(Spike_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,6]  <- RBD_IgM_mod_val[sample( nrow(RBD_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,7]  <- NP_IgM_mod_val[sample( nrow(NP_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,8]  <- Spike_IgA_mod_val[sample( nrow(Spike_IgA_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,9]  <- RBD_IgA_mod_val[sample( nrow(RBD_IgA_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,10] <- NP_IgA_mod_val[sample( nrow(NP_IgA_mod_val), 1),which(tt_plot==tt_sample)]
}

###########################
## Simulate infections (4-6 months)

for( i in (N_neg_val+N_pos1_val+1):(N_neg_val+N_pos1_val+N_pos2_val) )
{ 
	tt_sample <- round(runif( 1, min=91, max=179 ))

	VAL[i,2]  <- Spike_IgG_mod_val[sample( nrow(Spike_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,3]  <- RBD_IgG_mod_val[sample( nrow(RBD_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,4]  <- NP_IgG_mod_val[sample( nrow(NP_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,5]  <- Spike_IgM_mod_val[sample( nrow(Spike_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,6]  <- RBD_IgM_mod_val[sample( nrow(RBD_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,7]  <- NP_IgM_mod_val[sample( nrow(NP_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,8]  <- Spike_IgA_mod_val[sample( nrow(Spike_IgA_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,9]  <- RBD_IgA_mod_val[sample( nrow(RBD_IgA_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,10] <- NP_IgA_mod_val[sample( nrow(NP_IgA_mod_val), 1),which(tt_plot==tt_sample)]
}



###########################
## Simulate infections (7-12 months)


for( i in (N_neg_val+N_pos1_val+N_pos2_val+1):(N_neg_val+N_pos1_val+N_pos2_val+N_pos3_val) )
{ 
	tt_sample <- round(runif( 1, min=181, max=360 ))

	VAL[i,2]  <- Spike_IgG_mod_val[sample( nrow(Spike_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,3]  <- RBD_IgG_mod_val[sample( nrow(RBD_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,4]  <- NP_IgG_mod_val[sample( nrow(NP_IgG_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,5]  <- Spike_IgM_mod_val[sample( nrow(Spike_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,6]  <- RBD_IgM_mod_val[sample( nrow(RBD_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,7]  <- NP_IgM_mod_val[sample( nrow(NP_IgM_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,8]  <- Spike_IgA_mod_val[sample( nrow(Spike_IgA_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,9]  <- RBD_IgA_mod_val[sample( nrow(RBD_IgA_mod_val), 1),which(tt_plot==tt_sample)]
	VAL[i,10] <- NP_IgA_mod_val[sample( nrow(NP_IgA_mod_val), 1),which(tt_plot==tt_sample)]
}


###########################
## Transform data to log scale

VAL[,2:10] <- log(VAL[,2:10])


########################################
##       ##  Simulate a data set for  ##
##  1.2  ##  algorithm APPlication    ##
##       ##                           ##
########################################

APP_data_gen <- function( N_neg_app,     ## Number of negative controls
                          N_pos_waveA,   ## Number of positive samples from wave A
                          t_waveA,       ## Timing of wave A
                          sig_waveA,     ## Std dev in timing of wave A
                          N_pos_waveB,   ## Number of positive samples from wave B 
                          t_waveB,       ## Timing of wave B
                          sig_waveB )    ## Std dev in timing of wave B
{  
	###########################
	## Sample numbers

	N_tot_app <- N_neg_app + N_pos_waveA + N_pos_waveB
	
	N_pos_app <- N_pos_waveA + N_pos_waveB


	###########################
	## Simulate time since infection from wave A

	t_inf_waveA <- round(rnorm( n=N_pos_waveA, mean=t_waveA, sd=sig_waveA ))

	againdex <- which( t_inf_waveA <= 14 | t_inf_waveA > 360 )

	while( length(againdex) > 0 )
	{
		t_inf_waveA[againdex] <- round(rnorm( n=length(againdex), mean=t_waveA, sd=sig_waveA ))
	
		againdex <- which( t_inf_waveA <= 14 | t_inf_waveA > 360 )
	}


	###########################
	## Simulate time since infection from wave B

	t_inf_waveB <- round(rnorm( n=N_pos_waveB, mean=t_waveB, sd=sig_waveB ))
	
	againdex <- which( t_inf_waveB <= 14 | t_inf_waveB > 360 )

	while( length(againdex) > 0 )
	{
		t_inf_waveB[againdex] <- round(rnorm( n=length(againdex), mean=t_waveB, sd=sig_waveB ))
	
		againdex <- which( t_inf_waveB <= 14 | t_inf_waveB > 360 )
	}
	

	t_inf <- c( t_inf_waveA, t_inf_waveB )



	APP <- matrix( NA, nrow = N_tot_app, ncol=11 )
	colnames(APP) <- c("group", "Spike_IgG", "RBD_IgG", "NP_IgG",
	                            "Spike_IgM", "RBD_IgM", "NP_IgM", 
      	                      "Spike_IgA", "RBD_IgA", "NP_IgA", "t_inf" )

	APP <- as.data.frame(APP)

	APP[(N_neg_app+1):(N_neg_app+N_pos_app),11] <- t_inf


	APP[1:N_neg_app,1]                              <- "neg"
	
	APP[which( APP[,11] >14 & APP[,11] <= 90 )  ,1] <- "pos1"
	APP[which( APP[,11] >90 & APP[,11] <= 180 ) ,1] <- "pos2"
	APP[which( APP[,11] >180 & APP[,11] <= 360 ),1] <- "pos3"


	###########################
	## Simulate negatve controls

	APP[1:N_neg_app,2]  <- sample( x = Spike_IgG_neg_data, size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,3]  <- sample( x = RBD_IgG_neg_data,   size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,4]  <- sample( x = NP_IgG_neg_data,    size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,5]  <- sample( x = Spike_IgM_neg_data, size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,6]  <- sample( x = RBD_IgM_neg_data,   size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,7]  <- sample( x = NP_IgM_neg_data,    size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,8]  <- sample( x = Spike_IgA_neg_data, size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,9]  <- sample( x = RBD_IgA_neg_data,   size=N_neg_app, replace=TRUE)
	APP[1:N_neg_app,10] <- sample( x = RBD_IgA_neg_data,   size=N_neg_app, replace=TRUE)


	###########################
	## Simulate infections (0-3 months)

	for( i in (N_neg_app+1):(N_neg_app+N_pos_app) )
	{ 
		tt_sample <- APP[i,11]

		APP[i,2]  <- Spike_IgG_mod_app[sample( nrow(Spike_IgG_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,3]  <- RBD_IgG_mod_app[sample( nrow(RBD_IgG_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,4]  <- NP_IgG_mod_app[sample( nrow(NP_IgG_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,5]  <- Spike_IgM_mod_app[sample( nrow(Spike_IgM_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,6]  <- RBD_IgM_mod_app[sample( nrow(RBD_IgM_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,7]  <- NP_IgM_mod_app[sample( nrow(NP_IgM_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,8]  <- Spike_IgA_mod_app[sample( nrow(Spike_IgA_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,9]  <- RBD_IgA_mod_app[sample( nrow(RBD_IgA_mod_app), 1),which(tt_plot==tt_sample)]
		APP[i,10] <- NP_IgA_mod_app[sample( nrow(NP_IgA_mod_app), 1),which(tt_plot==tt_sample)]
	}


	###########################
	## Transform to log scale
	
	APP[,2:10] <- log(APP[,2:10])

	APP
}

APP_1 <- APP_data_gen( N_neg_app   = 600,     ## Number of negative controls
                       N_pos_waveA = 0,       ## Number of positive samples from wave A
                       t_waveA     = 250,     ## Timing of wave A
                       sig_waveA   = 30,      ## Std dev in timing of wave A
                       N_pos_waveB = 400,     ## Number of positive samples from wave B 
                       t_waveB     = 50,      ## Timing of wave B
                       sig_waveB   = 30 )     ## Std dev in timing of wave B


APP_2 <- APP_data_gen( N_neg_app   = 600,     ## Number of negative controls
                       N_pos_waveA = 400,     ## Number of positive samples from wave A
                       t_waveA     = 250,     ## Timing of wave A
                       sig_waveA   = 30,      ## Std dev in timing of wave A
                       N_pos_waveB = 0,       ## Number of positive samples from wave B 
                       t_waveB     = 50,      ## Timing of wave B
                       sig_waveB   = 30 )     ## Std dev in timing of wave B


APP_3 <- APP_data_gen( N_neg_app   = 600,     ## Number of negative controls
                       N_pos_waveA = 300,     ## Number of positive samples from wave A
                       t_waveA     = 250,     ## Timing of wave A
                       sig_waveA   = 30,      ## Std dev in timing of wave A
                       N_pos_waveB = 100,     ## Number of positive samples from wave B 
                       t_waveB     = 50,      ## Timing of wave B
                       sig_waveB   = 30 )     ## Std dev in timing of wave B



APP_4 <- APP_data_gen( N_neg_app   = 600,     ## Number of negative controls
                       N_pos_waveA = 200,     ## Number of positive samples from wave A
                       t_waveA     = 250,     ## Timing of wave A
                       sig_waveA   = 50,      ## Std dev in timing of wave A
                       N_pos_waveB = 200,     ## Number of positive samples from wave B 
                       t_waveB     = 70,     ## Timing of wave B
                       sig_waveB   = 70 )     ## Std dev in timing of wave B



########################################
##       ##  Plot an overview of the  ##
##  1.3  ##  VALidation and           ##
##       ##  APPlication data         ##
########################################

APP_data_plotter <- function( APP )
{

	par(mfrow=c(3,3))

	###########################
	## Spike IgG

	boxplot( VAL$Spike_IgG[which(VAL$group=="neg")],
	         APP$Spike_IgG[which(APP$group=="neg")],
	         VAL$Spike_IgG[which(VAL$group=="pos1")],
	         APP$Spike_IgG[which(APP$group=="pos1")],
	         VAL$Spike_IgG[which(VAL$group=="pos2")],
	         APP$Spike_IgG[which(APP$group=="pos2")],
	         VAL$Spike_IgG[which(VAL$group=="pos3")],
	         APP$Spike_IgG[which(APP$group=="pos3")],
	         main="Spike IgG", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )


	###########################
	## RBD IgG

	boxplot( VAL$RBD_IgG[which(VAL$group=="neg")],
	         APP$RBD_IgG[which(APP$group=="neg")],
	         VAL$RBD_IgG[which(VAL$group=="pos1")],
	         APP$RBD_IgG[which(APP$group=="pos1")],
	         VAL$RBD_IgG[which(VAL$group=="pos2")],
	         APP$RBD_IgG[which(APP$group=="pos2")],
	         VAL$RBD_IgG[which(VAL$group=="pos3")],
	         APP$RBD_IgG[which(APP$group=="pos3")],
	         main="RBD IgG", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )

	###########################
	## NP IgG

	boxplot( VAL$NP_IgG[which(VAL$group=="neg")],
	         APP$NP_IgG[which(APP$group=="neg")],
	         VAL$NP_IgG[which(VAL$group=="pos1")],
	         APP$NP_IgG[which(APP$group=="pos1")],
	         VAL$NP_IgG[which(VAL$group=="pos2")],
	         APP$NP_IgG[which(APP$group=="pos2")],
	         VAL$NP_IgG[which(VAL$group=="pos3")],
	         APP$NP_IgG[which(APP$group=="pos3")],
	         main="NP IgG", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )


	###########################
	## Spike IgM

	boxplot( VAL$Spike_IgM[which(VAL$group=="neg")],
	         APP$Spike_IgM[which(APP$group=="neg")],
	         VAL$Spike_IgM[which(VAL$group=="pos1")],
	         APP$Spike_IgM[which(APP$group=="pos1")],
	         VAL$Spike_IgM[which(VAL$group=="pos2")],
	         APP$Spike_IgM[which(APP$group=="pos2")],
	         VAL$Spike_IgM[which(VAL$group=="pos3")],
	         APP$Spike_IgM[which(APP$group=="pos3")],
	         main="Spike IgM", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )


	###########################
	## RBD IgM
	
	boxplot( VAL$RBD_IgM[which(VAL$group=="neg")],
	         APP$RBD_IgM[which(APP$group=="neg")],
	         VAL$RBD_IgM[which(VAL$group=="pos1")],
	         APP$RBD_IgM[which(APP$group=="pos1")],
	         VAL$RBD_IgM[which(VAL$group=="pos2")],
	         APP$RBD_IgM[which(APP$group=="pos2")],
	         VAL$RBD_IgM[which(VAL$group=="pos3")],
	         APP$RBD_IgM[which(APP$group=="pos3")],
	         main="RBD IgM", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
     	                 "4-6m", "4-6m", "7-12m", "7-12m")  )

	###########################
	## NP IgM

	boxplot( VAL$NP_IgM[which(VAL$group=="neg")],
	         APP$NP_IgM[which(APP$group=="neg")],
	         VAL$NP_IgM[which(VAL$group=="pos1")],
	         APP$NP_IgM[which(APP$group=="pos1")],
	         VAL$NP_IgM[which(VAL$group=="pos2")],
	         APP$NP_IgM[which(APP$group=="pos2")],
	         VAL$NP_IgM[which(VAL$group=="pos3")],
	         APP$NP_IgM[which(APP$group=="pos3")],
      	   main="NP IgM", ylab="log(antibody level)",
 	        col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )


	###########################
	## Spike IgA

	boxplot( VAL$Spike_IgA[which(VAL$group=="neg")],
	         APP$Spike_IgA[which(APP$group=="neg")],
	         VAL$Spike_IgA[which(VAL$group=="pos1")],
	         APP$Spike_IgA[which(APP$group=="pos1")],
	         VAL$Spike_IgA[which(VAL$group=="pos2")],
	         APP$Spike_IgA[which(APP$group=="pos2")],
	         VAL$Spike_IgA[which(VAL$group=="pos3")],
	         APP$Spike_IgA[which(APP$group=="pos3")],
	         main="Spike IgA", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )


	###########################
	## RBD IgA

	boxplot( VAL$RBD_IgA[which(VAL$group=="neg")],
		   APP$RBD_IgA[which(APP$group=="neg")],
		   VAL$RBD_IgA[which(VAL$group=="pos1")],
	         APP$RBD_IgA[which(APP$group=="pos1")],
	         VAL$RBD_IgA[which(VAL$group=="pos2")],
	         APP$RBD_IgA[which(APP$group=="pos2")],
	         VAL$RBD_IgA[which(VAL$group=="pos3")],
	         APP$RBD_IgA[which(APP$group=="pos3")],
	         main="RBD IgA", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )

	###########################
	## NP IgA

	boxplot( VAL$NP_IgA[which(VAL$group=="neg")],
	         APP$NP_IgA[which(APP$group=="neg")],
	         VAL$NP_IgA[which(VAL$group=="pos1")],
	         APP$NP_IgA[which(APP$group=="pos1")],
	         VAL$NP_IgA[which(VAL$group=="pos2")],
	         APP$NP_IgA[which(APP$group=="pos2")],
	         VAL$NP_IgA[which(VAL$group=="pos3")],
	         APP$NP_IgA[which(APP$group=="pos3")],
	         main="NP IgA", ylab="log(antibody level)",
	         col = rep( c("limegreen", "orangered"), 4),
	         names=c("neg", "neg", "0-3m", "0-3m", 
	                 "4-6m", "4-6m", "7-12m", "7-12m")  )
}


APP_data_plotter( APP_1 )

APP_data_plotter( APP_2 )

APP_data_plotter( APP_3 )

APP_data_plotter( APP_4 )


#####################################################
#####################################################
##          ##                                     ##  
##   ####   ##   ####  ##     ####   ####   ####   ##
##  ##  ##  ##  ##  ## ##    ##  ## ##     ##      ##
##     ##   ##  ##     ##    ######  ####   ####   ##
##    ##    ##  ##  ## ##    ##  ##     ##     ##  ##
##   #####  ##   ####  ##### ##  ##  ####   ####   ##
##          ##                                     ##
#####################################################
#####################################################

#####################################################
##                                                 ##
## 2.1 ALGORITHM 1v1 (IgG)                         ##
##                                                 ##
#####################################################

#########################################
## Random forests 4 way classifier: IgG only

RF_4way_IgG = randomForest( factor(group) ~ Spike_IgG + RBD_IgG + NP_IgG,  
                              data = VAL, importance = TRUE, ntree = N_tree)
 


#########################################
## RF 4 way prediction object: IgG only

ROC_4way_IgG <- multiclass.roc( VAL$group, RF_4way_IgG$votes )

RF_4way_IgG_pred <- ROC_4way_IgG$predictor

RF_4way_IgG_pred <- RF_4way_IgG_pred[,c( which(colnames(RF_4way_IgG_pred) == "neg"),
                                         which(colnames(RF_4way_IgG_pred) == "pos1"),                                         
						     which(colnames(RF_4way_IgG_pred) == "pos2"),
                                         which(colnames(RF_4way_IgG_pred) == "pos3") )]

RF_4way_IgG_pred <- as.data.frame(RF_4way_IgG_pred)

RF_4way_IgG_pred <- cbind( VAL$group, RF_4way_IgG_pred )
colnames(RF_4way_IgG_pred)[1] <- "group"



#####################################################
##                                                 ##
## 2.2 ALGORITHM 1v2 (IgG + IgM + IgA)             ##
##                                                 ##
#####################################################

#########################################
## Random forests 4 way classifier: IgG + IgM + IgA

RF_4way_IgGMA = randomForest( factor(group) ~ Spike_IgG + RBD_IgG + NP_IgG + 
                                              Spike_IgM + RBD_IgM + NP_IgM + 
                                              Spike_IgA + RBD_IgA + NP_IgA,  
                              data = VAL, importance = TRUE, ntree = N_tree)

 

#########################################
## RF 4 way prediction object: IgG + IgM + IgA

ROC_4way_IgGMA <- multiclass.roc( VAL$group, RF_4way_IgGMA$votes )

RF_4way_IgGMA_pred <- ROC_4way_IgGMA$predictor

RF_4way_IgGMA_pred <- RF_4way_IgGMA_pred[,c( which(colnames(RF_4way_IgGMA_pred) == "neg"),
                                             which(colnames(RF_4way_IgGMA_pred) == "pos1"),
                                             which(colnames(RF_4way_IgGMA_pred) == "pos2"),
                                             which(colnames(RF_4way_IgGMA_pred) == "pos3") )]

RF_4way_IgGMA_pred <- as.data.frame(RF_4way_IgGMA_pred)

RF_4way_IgGMA_pred <- cbind( VAL$group, RF_4way_IgGMA_pred )
colnames(RF_4way_IgGMA_pred)[1] <- "group"



#####################################################
##                                                 ##
## 2.3 ROC curve 1: negative vs positive           ##
##                                                 ##
#####################################################

#########################################
## ROC curve: binary infection (IgG only)

SS_cut_neg_pos <- sort( unique(c(0, RF_4way_IgGMA_pred[,2], 1)), decreasing=FALSE )
SS_cut_neg_pos <- 0.5*( SS_cut_neg_pos[-1] + SS_cut_neg_pos[-length(SS_cut_neg_pos)] )
SS_cut_neg_pos <- c(0, SS_cut_neg_pos, 1)

N_cut_neg_pos <- length(SS_cut_neg_pos)

ROC1_neg_pos <- matrix( NA, nrow=N_cut_neg_pos, ncol=3)
colnames(ROC1_neg_pos) <- c("cut", "sens", "spec")

ROC1_neg_pos[,1] <- SS_cut_neg_pos


index_pos <- which( RF_4way_IgGMA_pred$group %in% c("pos1", "pos2", "pos3") )
index_neg <- which( RF_4way_IgGMA_pred$group == "neg" )

for(i in 1:N_cut_neg_pos)
{
	## sensitivity (identifying positive)
	ROC1_neg_pos[i,2] <- length(which( RF_4way_IgGMA_pred[index_pos,2] <= SS_cut_neg_pos[i] ))/(N_pos_val) 	
 
	## specificity (identifying negative)
	ROC1_neg_pos[i,3] <- length(which( RF_4way_IgGMA_pred[index_neg,2] > SS_cut_neg_pos[i] ))/N_neg_val 	
}


#########################################
## Select cut off for classification of positive individuals
## based on IgG only algorithm

pos_target_sens <- ROC1_neg_pos[ max(which(ROC1_neg_pos[,3] >= pos_target_spec)),2 ]

pos_target_cut <- ROC1_neg_pos[ max(which(ROC1_neg_pos[,3] >= pos_target_spec)),1 ]


#########################################
## Make predictions

RF_4way_IgGMA_pred$pred <- rep( NA, nrow(RF_4way_IgGMA_pred) )

RF_4way_IgGMA_pred$pred[ which( RF_4way_IgGMA_pred[,2] > pos_target_cut ) ] <- "neg"

for(i in 1:nrow(RF_4way_IgGMA_pred))
{
	if( RF_4way_IgGMA_pred[i,2] > pos_target_cut )
	{
		RF_4way_IgGMA_pred$pred[i] <- "neg"
	}else{
		RF_4way_IgGMA_pred$pred[i] <- names(which.max(RF_4way_IgGMA_pred[i,3:5]))
	}
}





#########################################
## Confusion matrix

CMAT <- table( RF_4way_IgGMA_pred$group, RF_4way_IgGMA_pred$pred )
rownames(CMAT) <- c("cat_neg", "cat_pos1", "cat_pos2", "cat_pos3")
colnames(CMAT) <- c("class_neg", "class_pos1", "class_pos2", "class_pos3")

CMAT[1,] <- CMAT[1,]/N_neg_val
CMAT[2,] <- CMAT[2,]/N_pos1_val
CMAT[3,] <- CMAT[3,]/N_pos2_val
CMAT[4,] <- CMAT[4,]/N_pos3_val



#########################################
## Plot ROC curves


par(mfrow=c(3,1))

#########################################
## Panel 1: Previous infection

plot( x = 1 - ROC1_neg_pos[,3], y = ROC1_neg_pos[,2],
      type = 'S',
      xlim=c(0,1), ylim=c(0,1),
      main = "Classificying previous infection",
      xlab = "1 - specificity (negatives correctly classified)", 
      ylab = "sensitivity (previous infections identified)" ) 

points( x = 1 - pos_target_spec, y = pos_target_sens,
        pch=19, col="red")

points( x = c(0,1), y = c(0,1),
        type='l', lty="dashed")




#########################################
## Panel 2: Older than 3 months

plot( x = 1 - ROC_4way_IgGMA$rocs$`pos1/pos2`[[1]]$specificities, 
      y = ROC_4way_IgGMA$rocs$`pos1/pos2`[[1]]$sensitivities,
      type='S',
      xlim=c(0,1), ylim=c(0,1),
      main = "Classificying recent (<3 mths) vs older infections",
      xlab = "1 - specificity (old infection identified)", 
      ylab = "sensitivity (recent infection identified)" )

points( x = 1 - ROC_4way_IgGMA$rocs$`pos1/pos3`[[1]]$specificities, 
      y = ROC_4way_IgGMA$rocs$`pos1/pos3`[[1]]$sensitivities,
      type='S', col="forestgreen")

points( x = c(0,1), y = c(0,1),
        type='l', lty="dashed")

legend( x = "bottomright",
legend = c("<3 mths vs 3-6 mths", "<3 mths vs 7-12 mths"),
fill = c("black", "forestgreen"),
border = c("black", "forestgreen") )


#########################################
## Panel 3: Older than 6 months

plot( x = 1 - ROC_4way_IgGMA$rocs$`pos2/pos3`[[1]]$specificities, 
      y = ROC_4way_IgGMA$rocs$`pos2/pos3`[[1]]$sensitivities,
      type='S',
      xlim=c(0,1), ylim=c(0,1),
      main = "Classificying recent (3-6 mths) vs (7-12 mths) infections",
      xlab = "1 - specificity (old infection identified)", 
      ylab = "sensitivity (recent infection identified)" )

points( x = c(0,1), y = c(0,1),
        type='l', lty="dashed")








################################################
################################################
##          ##                                ##  
##   ####   ##  ##   ## ##   ##  ####  ##     ##
##  ##  ##  ##   ## ##  ##   ## ##  ## ##     ##
##     ##   ##    ###    ## ##  ###### ##     ##
##  ##  ##  ##   ## ##    ###   ##  ## ##     ##
##   ####   ##  ##   ##    #    ##  ## #####  ##
##          ##                                ##
################################################
################################################


xVAL_class <- function( APP )
{
	####################################
	## Cross-validated Random Forests algorithm

	ROC_test_neg_pos_xVAL_rep <- list()
	pos_target_sens_xVAL_rep  <- rep(NA, N_rep)


	CMAT_xVAL_rep <- array( NA, dim=c(N_rep,4,4) )

	RF_4way_IgGMA_pred_APP_rep <- list()


	for(n in 1:N_rep)
	{
		#################################################
		###       ## Create training and testing      ###
	  	###  3.1  ## data sets using validation data  ### 
		###       ##                                  ###
		#################################################
	
		index_train <- sample( N_tot_val, train_prop*N_tot_val )	
		index_test  <- setdiff( 1:N_tot_val, index_train )
	
		VAL_train <- VAL[index_train,]
		VAL_test  <- VAL[index_test,]
	
		N_neg_train  <- length(which( VAL_train$group == "neg" ))
		N_pos1_train <- length(which( VAL_train$group == "pos1" ))
		N_pos2_train <- length(which( VAL_train$group == "pos2" ))
		N_pos3_train <- length(which( VAL_train$group == "pos3" ))
	
		N_neg_test  <- length(which( VAL_test$group == "neg" ))
		N_pos1_test <- length(which( VAL_test$group == "pos1" ))
		N_pos2_test <- length(which( VAL_test$group == "pos2" ))
		N_pos3_test <- length(which( VAL_test$group == "pos3" ))
	
	
		####################################################
		###       ## Train RF algorithm IgG + IgM + IgA  ###
	  	###  3.2  ## on training data                    ### 
		###       ##                                     ###
		####################################################

		#################################################
		## 3.2.1. RF 5 way algorithm with IgG + IgA + IgM

		RF_4way_IgGMA_xVAL = randomForest( factor(group) ~ Spike_IgG + RBD_IgG + NP_IgG + 
		                                                   Spike_IgM + RBD_IgM + NP_IgM + 
		                                                   Spike_IgA + RBD_IgA + NP_IgA,  
		                             data = VAL_train, importance = TRUE, ntree = N_tree)

 
		#################################################
		## 3.2.2. RF 5 way prediction object

		ROC_4way_IgGMA_xVAL <- multiclass.roc( VAL_train$group, RF_4way_IgGMA_xVAL$votes )

		RF_4way_IgGMA_pred_xVAL <- ROC_4way_IgGMA_xVAL$predictor

		RF_4way_IgGMA_pred_xVAL <- RF_4way_IgGMA_pred_xVAL[,c( which(colnames(RF_4way_IgGMA_pred_xVAL) == "neg"),
	                                                             which(colnames(RF_4way_IgGMA_pred_xVAL) == "pos1"),
	                                                             which(colnames(RF_4way_IgGMA_pred_xVAL) == "pos2"),
	                                                             which(colnames(RF_4way_IgGMA_pred_xVAL) == "pos3") )]

		RF_4way_IgGMA_pred_xVAL <- as.data.frame(RF_4way_IgGMA_pred_xVAL)

		RF_4way_IgGMA_pred_xVAL <- cbind( VAL_train$group, RF_4way_IgGMA_pred_xVAL )
		colnames(RF_4way_IgGMA_pred_xVAL)[1] <- "group"


		#################################################
		###       ## Evaluate RF algorithm IgG only   ###
	  	###  3.3  ## on training data                 ### 
		###       ##                                  ###
		#################################################

		#################################################
		## 3.3.1. ROC curve on training data: binary 
		##        infection (IgG only)

		SS_cut_neg_pos_xVAL <- sort( unique(c(0, RF_4way_IgGMA_pred_xVAL[,2], 1)), decreasing=FALSE )
		SS_cut_neg_pos_xVAL <- 0.5*( SS_cut_neg_pos_xVAL[-1] + SS_cut_neg_pos_xVAL[-length(SS_cut_neg_pos_xVAL)] )
		SS_cut_neg_pos_xVAL <- c(0, SS_cut_neg_pos_xVAL, 1)

		N_cut_neg_pos_xVAL <- length(SS_cut_neg_pos_xVAL)

		ROC_neg_pos_xVAL <- matrix( NA, nrow=N_cut_neg_pos_xVAL, ncol=3)
		colnames(ROC_neg_pos_xVAL) <- c("cut", "sens", "spec")

		ROC_neg_pos_xVAL[,1] <- SS_cut_neg_pos_xVAL


		index_pos <- which( RF_4way_IgGMA_pred_xVAL$group %in% c("pos1", "pos2", "pos3") )
		index_neg <- which( RF_4way_IgGMA_pred_xVAL$group == "neg" )


		for(i in 1:N_cut_neg_pos_xVAL)
		{
			## sensitivity (identifying positive)
			ROC_neg_pos_xVAL[i,2] <- length(which( RF_4way_IgGMA_pred_xVAL[index_pos,2] <= SS_cut_neg_pos_xVAL[i] ))/length(index_pos) 	
 
			## specificity (identifying negative)
			ROC_neg_pos_xVAL[i,3] <- length(which( RF_4way_IgGMA_pred_xVAL[index_neg,2] > SS_cut_neg_pos_xVAL[i] ))/length(index_neg) 	
		}

		if( ROC_neg_pos_xVAL[1,2] > 0 )
		{
			ROC_neg_pos_xVAL <- rbind( c(0,0,1), ROC_neg_pos_xVAL )
		}

		#################################################
		## 3.3.2. Select cut off for classification of 
		##        positive individuals 

		pos_target_sens_xVAL <- ROC_neg_pos_xVAL[max(which(ROC_neg_pos_xVAL[,3] > pos_target_spec)),2]
	
		pos_target_cut_xVAL <- ROC_neg_pos_xVAL[max(which(ROC_neg_pos_xVAL[,3] > pos_target_spec)),1]


		#################################################
		###       ## Evaluate RF algorithm on         ###
	  	###  3.4  ## testing data                     ### 
		###       ##                                  ###
		#################################################

		#################################################
		## 3.4.1. Apply trained RF algorithm to  
		##        test data (IgG + IgM + IgM)

		RF_4way_IgGMA_test_pred_obj_xVAL = predict( RF_4way_IgGMA_xVAL, newdata=VAL_test, 
	                                                  predict.all=TRUE)

		RF_4way_IgGMA_test_pred_xVAL <- matrix(NA, nrow=nrow(VAL_test), ncol=5)
		colnames(RF_4way_IgGMA_test_pred_xVAL) <- c("group", "neg", "pos1", "pos2", "pos3")

		RF_4way_IgGMA_test_pred_xVAL <- as.data.frame(RF_4way_IgGMA_test_pred_xVAL)

		RF_4way_IgGMA_test_pred_xVAL$group <- VAL_test$group

		RF_4way_IgGMA_test_pred_xVAL$neg  <- rowSums( RF_4way_IgGMA_test_pred_obj_xVAL$individual=="neg")/N_tree
		RF_4way_IgGMA_test_pred_xVAL$pos1 <- rowSums( RF_4way_IgGMA_test_pred_obj_xVAL$individual=="pos1")/N_tree
		RF_4way_IgGMA_test_pred_xVAL$pos2 <- rowSums( RF_4way_IgGMA_test_pred_obj_xVAL$individual=="pos2")/N_tree
		RF_4way_IgGMA_test_pred_xVAL$pos3 <- rowSums( RF_4way_IgGMA_test_pred_obj_xVAL$individual=="pos3")/N_tree


		#########################################
		## 3.4.2. Make predictions

		RF_4way_IgGMA_test_pred_xVAL$pred <- rep( NA, nrow(RF_4way_IgGMA_test_pred_xVAL) )
	
		for(i in 1:nrow(RF_4way_IgGMA_test_pred_xVAL))
		{
			if( RF_4way_IgGMA_test_pred_xVAL[i,2] > pos_target_cut_xVAL )
			{
				RF_4way_IgGMA_test_pred_xVAL$pred[i] <- "neg"
			}else{
				RF_4way_IgGMA_test_pred_xVAL$pred[i] <- names(which.max(RF_4way_IgGMA_test_pred_xVAL[i,3:6]))
			}
		}


		#################################################
		## 3.4.3. ROC curve on testing data:  
		##        binary infection (IgG + IgM + IgA)

		SS_cut_neg_pos_xVAL <- sort( unique(c(0, RF_4way_IgGMA_test_pred_xVAL[,2], 1)), decreasing=FALSE )
		SS_cut_neg_pos_xVAL <- 0.5*( SS_cut_neg_pos_xVAL[-1] + SS_cut_neg_pos_xVAL[-length(SS_cut_neg_pos_xVAL)] )
		SS_cut_neg_pos_xVAL <- c(0, SS_cut_neg_pos_xVAL, 1)

		N_cut_neg_pos_xVAL <- length(SS_cut_neg_pos_xVAL)

		ROC_test_neg_pos_xVAL <- matrix( NA, nrow=N_cut_neg_pos_xVAL, ncol=3)
		colnames(ROC_test_neg_pos_xVAL) <- c("cut", "sens", "spec")

		ROC_test_neg_pos_xVAL[,1] <- SS_cut_neg_pos_xVAL


		index_pos <- which( RF_4way_IgGMA_test_pred_xVAL$group %in% c("pos1", "pos2", "pos3") )
		index_neg <- which( RF_4way_IgGMA_test_pred_xVAL$group == "neg" )

		for(i in 1:N_cut_neg_pos_xVAL)
		{
			## sensitivity
			ROC_test_neg_pos_xVAL[i,2] <- length(which( RF_4way_IgGMA_test_pred_xVAL[index_pos,2] <= SS_cut_neg_pos_xVAL[i] ))/length(index_pos) 	

			## spcificity
			ROC_test_neg_pos_xVAL[i,3] <- length(which( RF_4way_IgGMA_test_pred_xVAL[index_neg,2] > SS_cut_neg_pos_xVAL[i] ))/length(index_neg) 	
		}

		if( ROC_test_neg_pos_xVAL[1,2] > 0 )
		{
			ROC_test_neg_pos_xVAL <- rbind( c(0,0,1), ROC_test_neg_pos_xVAL )
		}


		#################################################
		## 3.4.4. Confusion matrix: trained algorithm 
		##        applied to testing data set
		##        Implemented the long way because 'table'
		##        doesn't work when variables are missing

		CMAT_xVAL <- matrix(NA, nrow=4, ncol=4)
		rownames(CMAT_xVAL) <- c("cat_neg", "cat_pos1", "cat_pos2", "cat_pos3")
		colnames(CMAT_xVAL) <- c("class_neg", "class_pos1", "class_pos2", "class_pos3")

		CMAT_xVAL[1,1] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "neg")] == "neg"))
		CMAT_xVAL[1,2] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "neg")] == "pos1"))
		CMAT_xVAL[1,3] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "neg")] == "pos2"))
		CMAT_xVAL[1,4] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "neg")] == "pos3"))

		CMAT_xVAL[2,1] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos1")] == "neg"))
		CMAT_xVAL[2,2] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos1")] == "pos1"))
		CMAT_xVAL[2,3] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos1")] == "pos2"))
		CMAT_xVAL[2,4] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos1")] == "pos3"))

		CMAT_xVAL[3,1] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos2")] == "neg"))
		CMAT_xVAL[3,2] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos2")] == "pos1"))
		CMAT_xVAL[3,3] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos2")] == "pos2"))
		CMAT_xVAL[3,4] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos2")] == "pos3"))

		CMAT_xVAL[4,1] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos3")] == "neg"))
		CMAT_xVAL[4,2] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos3")] == "pos1"))
		CMAT_xVAL[4,3] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos3")] == "pos2"))
		CMAT_xVAL[4,4] <- length(which(RF_4way_IgGMA_test_pred_xVAL$pred[which(RF_4way_IgGMA_test_pred_xVAL$group == "pos3")] == "pos3"))

		CMAT_xVAL[1,] <- CMAT_xVAL[1,]/N_neg_test
		CMAT_xVAL[2,] <- CMAT_xVAL[2,]/N_pos1_test
		CMAT_xVAL[3,] <- CMAT_xVAL[3,]/N_pos2_test
		CMAT_xVAL[4,] <- CMAT_xVAL[4,]/N_pos3_test


		#################################################
		###       ## Evaluate RF algorithm on         ###
	  	###  3.5  ## application data                 ### 
		###       ##                                  ###
		#################################################


		#################################################
		## 3.5.1. Apply trained RF algorithm to  
		##        test data (IgG + IgM + IgM)

		RF_4way_IgGMA_pred_obj_APP = predict( RF_4way_IgGMA_xVAL, newdata=APP, 
	                                            predict.all=TRUE)

		RF_4way_IgGMA_pred_APP <- matrix(NA, nrow=nrow(APP), ncol=5)
		colnames(RF_4way_IgGMA_pred_APP) <- c("group", "neg", "pos1", "pos2", "pos3")

		RF_4way_IgGMA_pred_APP <- as.data.frame(RF_4way_IgGMA_pred_APP)

		RF_4way_IgGMA_pred_APP$group <- APP$group

		RF_4way_IgGMA_pred_APP$neg  <- rowSums( RF_4way_IgGMA_pred_obj_APP$individual=="neg")/N_tree
		RF_4way_IgGMA_pred_APP$pos1 <- rowSums( RF_4way_IgGMA_pred_obj_APP$individual=="pos1")/N_tree
		RF_4way_IgGMA_pred_APP$pos2 <- rowSums( RF_4way_IgGMA_pred_obj_APP$individual=="pos2")/N_tree
		RF_4way_IgGMA_pred_APP$pos3 <- rowSums( RF_4way_IgGMA_pred_obj_APP$individual=="pos3")/N_tree


		#########################################
		## 3.5.2. Make predictions

		RF_4way_IgGMA_pred_APP$pred <- rep( NA, nrow(RF_4way_IgGMA_pred_APP) )
	
		for(i in 1:nrow(RF_4way_IgGMA_pred_APP))
		{
			if( RF_4way_IgGMA_pred_APP[i,2] > pos_target_cut_xVAL )
			{
				RF_4way_IgGMA_pred_APP$pred[i] <- "neg"
			}else{
				RF_4way_IgGMA_pred_APP$pred[i] <- names(which.max(RF_4way_IgGMA_pred_APP[i,3:5]))
			}
		}


		#################################################
		###       ## Store cross-validated output     ###
	  	###  3.6  ##                                  ### 
		###       ##                                  ###	
		#################################################

		ROC_test_neg_pos_xVAL_rep[[n]]  <- ROC_test_neg_pos_xVAL
		pos_target_sens_xVAL_rep[n]     <- pos_target_sens_xVAL

		CMAT_xVAL_rep[n,,]              <- CMAT_xVAL

		RF_4way_IgGMA_pred_APP_rep[[n]] <- RF_4way_IgGMA_pred_APP

	}

	output_list <- list()

	output_list[[1]] <- ROC_test_neg_pos_xVAL_rep
	output_list[[2]] <- pos_target_sens_xVAL_rep
	output_list[[3]] <- RF_4way_IgGMA_pred_APP_rep
	output_list[[4]] <- CMAT_xVAL_rep

	output_list
}


APP_1_pred <- xVAL_class( APP_1 )

APP_2_pred <- xVAL_class( APP_2 )

APP_3_pred <- xVAL_class( APP_3 )

APP_4_pred <- xVAL_class( APP_4 )


###############################################
###############################################
##          ##                               ## 
##  ##      ##  #####  #####   ####   ####   ##
##  ## ##   ##  ##  ## ##  ## ##  ## ##  ##  ##
##  ######  ##  #####  #####  ##  ## ##      ## 
##     ##   ##  ##     ## ##  ##  ## ##  ##  ## 
##     ##   ##  ##     ##  ##  ####   ####   ##
##          ##                               ##
###############################################
###############################################

####################################
####################################
##                                
##  Collate Xval results          

N_SS_cut <- 1000 

SS_seq <- seq(from=0, to=1, length=N_SS_cut)


APP_proc <- function( APP_pred )
{
	ROC_test_neg_pos_xVAL_rep  <- APP_pred[[1]]
	pos_target_sens_xVAL_rep   <- APP_pred[[2]]
	RF_4way_IgGMA_pred_APP_rep <- APP_pred[[3]]
	CMAT_xVAL_rep              <- APP_pred[[4]]


	#########################################
	## Binary infection: negative versus positive
	## Cross-validated variation in sensitivity
	## (i.e. holding specificity fixed)

	ROC_test_neg_pos_xVAL_sens_vary <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

	for(n in 1:N_rep)
	{
		for(i in 1:(N_SS_cut-1))
		{
			index_spec <- max(which(ROC_test_neg_pos_xVAL_rep[[n]][,3] > SS_seq[i]))

			ROC_test_neg_pos_xVAL_sens_vary[i,n] <- ROC_test_neg_pos_xVAL_rep[[n]][index_spec,2]
		}
	}
	ROC_test_neg_pos_xVAL_sens_vary[N_SS_cut,] <- 0


	ROC_test_neg_pos_xVAL_sens_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
	for(i in 1:N_SS_cut)
	{
		ROC_test_neg_pos_xVAL_sens_quant[i,] <- quantile( ROC_test_neg_pos_xVAL_sens_vary[i,], prob=c(0.5, 0.025, 0.975) )
	}


	#########################################
	## Binary infection: negative versus positive
	## Cross-validated variation in specificty
	## (i.e. holding sensitivity fixed)

	ROC_test_neg_pos_xVAL_spec_vary <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

	for(n in 1:N_rep)
	{
		for(i in 1:(N_SS_cut-1))
		{
			index_sens <- max(which(ROC_test_neg_pos_xVAL_rep[[n]][,2] <= SS_seq[i]))

			ROC_test_neg_pos_xVAL_spec_vary[i,n] <- ROC_test_neg_pos_xVAL_rep[[n]][index_sens,3]
		}
	}
	ROC_test_neg_pos_xVAL_spec_vary[N_SS_cut,] <- 0


	ROC_test_neg_pos_xVAL_spec_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
	for(i in 1:N_SS_cut)
	{
		ROC_test_neg_pos_xVAL_spec_quant[i,] <- quantile( ROC_test_neg_pos_xVAL_spec_vary[i,], prob=c(0.5, 0.025, 0.975) )
	}



	#########################################
	## Confusion matrix (with uncertainty)
	## based on cross-validation

	CMAT_xVAL_med <- matrix(NA, nrow=4, ncol=4)
	rownames(CMAT_xVAL_med) <- c("cat_neg", "cat_pos1", "cat_pos2", "cat_pos3")
	colnames(CMAT_xVAL_med) <- c("class_neg", "class_pos1", "class_pos2", "class_pos3")

	CMAT_xVAL_mean <- CMAT_xVAL_med
	CMAT_xVAL_lwr  <- CMAT_xVAL_med
	CMAT_xVAL_upr  <- CMAT_xVAL_med

	for(i in 1:4)
	{
		for(j in 1:4)
		{
			CMAT_xVAL_med[i,j] <- median( CMAT_xVAL_rep[,i,j] )

			CMAT_xVAL_mean[i,j] <- mean( CMAT_xVAL_rep[,i,j] )

			CMAT_xVAL_lwr[i,j] <- quantile( CMAT_xVAL_rep[,i,j], prob=0.025 )

			CMAT_xVAL_upr[i,j] <- quantile( CMAT_xVAL_rep[,i,j], prob=0.975 )
		}
	}
		


	#########################################
	## Confusion matrix (with uncertainty)
	## based on cross-validation

	APP_meas <- matrix(NA, nrow=N_rep, ncol=4)
	colnames(APP_meas) <- c("neg", "pos1", "pos2", "pos3")

	for(n in 1:N_rep)
	{
		N_tot_app <- nrow(RF_4way_IgGMA_pred_APP_rep[[n]])

		APP_meas[n,1] <- length(which(RF_4way_IgGMA_pred_APP_rep[[n]]$pred == "neg"))/N_tot_app
		APP_meas[n,2] <- length(which(RF_4way_IgGMA_pred_APP_rep[[n]]$pred == "pos1"))/N_tot_app
		APP_meas[n,3] <- length(which(RF_4way_IgGMA_pred_APP_rep[[n]]$pred == "pos2"))/N_tot_app
		APP_meas[n,4] <- length(which(RF_4way_IgGMA_pred_APP_rep[[n]]$pred == "pos3"))/N_tot_app
	}


	APP_meas_quant <- matrix(NA, nrow=3, ncol=4)
	colnames(APP_meas_quant) <- c("neg", "pos1", "pos2", "pos3")
	rownames(APP_meas_quant) <- c("med", "lwr", "upr")

	for(j in 1:4)
	{
		APP_meas_quant[,j] <- quantile( APP_meas[,j], prob=c(0.5, 0.025, 0.975) ) 
	}


	APP_adj <- matrix(NA, nrow=N_rep, ncol=4)
	colnames(APP_adj) <- c("neg", "pos1", "pos2", "pos3")

	for(n in 1:N_rep)
	{
		APP_adj[n,] <- c( solve( t(CMAT_xVAL_rep[n,,]) ) %*% APP_meas[n,] )

		index_neg <- 1 + which( APP_adj[n,2:4] < 0 )
		index_pos <- 1 + which( APP_adj[n,2:4] >= 0 )
	
		if( length(index_neg) > 0 )
		{
			APP_adj[n,index_pos] <- APP_adj[n,index_pos] + 
		                                  ( APP_adj[n,index_pos]/sum(APP_adj[n,index_pos]) )*( sum(APP_adj[n,index_neg]) ) 

			APP_adj[n,index_neg] <- 0
		}
	}


	APP_adj_quant <- matrix(NA, nrow=3, ncol=4)
	colnames(APP_adj_quant) <- c("neg", "pos1", "pos2", "pos3")
	rownames(APP_adj_quant) <- c("med", "lwr", "upr")

	for(j in 1:4)
	{
		APP_adj_quant[,j] <- quantile( APP_adj[,j], prob=c(0.5, 0.025, 0.975) ) 
	}

	output_list <- list()

	output_list[[1]] <- APP_meas_quant
	output_list[[2]] <- APP_adj_quant

	output_list
}

APP_1_proc <- APP_proc( APP_1_pred )

APP_2_proc <- APP_proc( APP_2_pred )

APP_3_proc <- APP_proc( APP_3_pred )

APP_4_proc <- APP_proc( APP_4_pred )


#############################################
#############################################
##         ##                              ## 
##  #####  ##  #####  ##     ####  ######  ##
##  ##     ##  ##  ## ##    ##  ##   ##    ##
##  ####   ##  #####  ##    ##  ##   ##    ## 
##     ##  ##  ##     ##    ##  ##   ##    ## 
##  ####   ##  ##     #####  ####    ##    ##
##         ##                              ##
#############################################
#############################################

t_bins <- seq(from=0, to=370, by=10)

N_bins <- length(t_bins)


#############################################
## Scenario 1

APP_1_meas_quant <- APP_1_proc[[1]]
APP_1_adj_quant  <- APP_1_proc[[2]]

N_neg_app_1 <- length(which( APP_1$group == "neg" ))
N_tot_app_1 <- nrow(APP_1)


inf_bins_1 <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	inf_bins_1[i] <- length(which( APP_1[,11] >= t_bins[i] & APP_1[,11] < t_bins[i+1] ))
}


#############################################
## Scenario 2

APP_2_meas_quant <- APP_2_proc[[1]]
APP_2_adj_quant  <- APP_2_proc[[2]]

N_neg_app_2 <- length(which( APP_2$group == "neg" ))
N_tot_app_2 <- nrow(APP_2)


inf_bins_2 <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	inf_bins_2[i] <- length(which( APP_2[,11] >= t_bins[i] & APP_2[,11] < t_bins[i+1] ))
}


#############################################
## Scenario 3

APP_3_meas_quant <- APP_3_proc[[1]]
APP_3_adj_quant  <- APP_3_proc[[2]]

N_neg_app_3 <- length(which( APP_3$group == "neg" ))
N_tot_app_3 <- nrow(APP_3)


inf_bins_3 <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	inf_bins_3[i] <- length(which( APP_3[,11] >= t_bins[i] & APP_3[,11] < t_bins[i+1] ))
}



#############################################
## Scenario 4

APP_4_meas_quant <- APP_4_proc[[1]]
APP_4_adj_quant  <- APP_4_proc[[2]]

N_neg_app_4 <- length(which( APP_4$group == "neg" ))
N_tot_app_4 <- nrow(APP_4)


inf_bins_4 <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	inf_bins_4[i] <- length(which( APP_4[,11] >= t_bins[i] & APP_4[,11] < t_bins[i+1] ))
}








neg_col <- "dodgerblue"
pos1_col <- "orangered"
pos2_col <- "forestgreen"
pos3_col <- "slateblue4"


tiff(file="Fig4_sero_recon.tif", width=40, height=20, units="cm", res=500)


lab.size   = 1.5
axis.size  = 1.25	
main.size  = 1.5
line.size  = 2
dash.line.size  = 1

lay.mat <- rbind( c( 1, 1, 2, 2 ), 
                  c( 3, 4, 5, 6 ),
                  c( 7, 7, 8, 8 ), 
                  c( 9,10,11,12 ), 
                  c(13,13,13,13 ) )
layout(lay.mat, heights=c(1.5,10,1.5,10,2), widths=c(2,8,2,8))
layout.show(13)


############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))

plot.new()
title( "(A) Recent wave scenario", 
        cex.main=3.0, line=-3)

plot.new()
title( "(B) Older wave scenario", 
        cex.main=3.0, line=-3)



####################
####################
###              ###
###  PANEL 1     ###
###  Scenario 1  ###
###  Negatives   ###
###              ###
####################
####################

par(mar = c(5,5,0.5,0.5))

line_seq_y <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

plot(x=1000, y=1000, 
xlim=c(0,100), ylim=c(0,1000),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

polygon( x = c(10,90,90,10),
         y = c(0,0,N_neg_app_1,N_neg_app_1),
         col=rgb(211/255,211/255,211/255,0.5), border=NA )




polygon( x = c(0,100,100,0), 
         y = c(APP_1_adj_quant[2,1], APP_1_adj_quant[2,1], APP_1_adj_quant[3,1], APP_1_adj_quant[3,1] )*N_tot_app_1,
         col=rgb(30/256,144/256,255/256,0.25), border=NA )

points( x=c(0,100), y=c(APP_1_adj_quant[1,1], APP_1_adj_quant[1,1])*N_tot_app_1,
        type='l', col=neg_col) 


axis(1, at = c(0, 100), 
        label=c("", ""), 
        cex.axis=axis.size )


axis(2, at = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), 
        label = c("0", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3.33, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="negatives")


####################
####################
###              ###
###  PANEL 2     ###
###  Scenario 1  ###
###  Positives   ###
###              ###
####################
####################


par(mar = c(5,5,0.5,2))
par(mgp=c(3,1.25,0))


line_seq_x <- -30*c(0, 3, 6, 9, 12)
line_seq_y <- c(0, 10, 20, 30, 40, 50, 60, 70)

plot(x=1000, y=1000, 
xlim=c(-360,30), ylim=c(0,70),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}



for(i in 1:N_bins)
{
	polygon( x=-t_bins[c(i,i+1,i+1,i)], y=c(0,0,inf_bins_1[i],inf_bins_1[i]),
               col=rgb(211/255,211/255,211/255,0.5), border=NA )
}

axis(1, at = -30*c(-1, 0, 3, 6, 9, 12), 
        label=c("", "0", "3", "6", "9", "12"), 
        cex.axis=1.5*axis.size )


axis(2, at = c(0, 10, 20, 30, 40, 50, 60, 70), 
        label = c("0", "10", "20", "30", "40", "50", "60", "70"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="positives: time before sampling (months)")

arrows( x0=0, y0=40, x1=0, y1=25,
        angle=30, lwd=2, length=0.1)

text( x=0, y=43, 
      labels="sampling", 
      cex=2)

###################
## 0-3 months

polygon( x = c(0,-90,-90,0), 
         y = c(APP_1_adj_quant[2,2], APP_1_adj_quant[2,2], APP_1_adj_quant[3,2], APP_1_adj_quant[3,2] )*N_tot_app_1/9,
         col=rgb(255/256,69/256,0/256,0.25), border=NA )

points( x = c(0,-90), 
        y = c(APP_1_adj_quant[1,2], APP_1_adj_quant[1,2] )*N_tot_app_1/9,
        type='l', col=pos1_col, lwd=2) 


###################
## 3-6 months

polygon( x = c(-90,-180,-180,-90), 
         y = c(APP_1_adj_quant[2,3], APP_1_adj_quant[2,3], APP_1_adj_quant[3,3], APP_1_adj_quant[3,3] )*N_tot_app_1/9,
         col=rgb(34/256,139/256,34/256,0.25), border=NA )

points( x = c(-90,-180), 
        y = c(APP_1_adj_quant[1,3], APP_1_adj_quant[1,3] )*N_tot_app_1/9,
        type='l', col=pos2_col, lwd=2) 

###################
## 7-12 months

polygon( x = c(-180,-360,-360,-180), 
         y = c(APP_1_adj_quant[2,4], APP_1_adj_quant[2,4], APP_1_adj_quant[3,4], APP_1_adj_quant[3,4] )*N_tot_app_1/18,
         col=rgb(71/256,60/256,139/256,0.25), border=NA )

points( x = c(-180,-360), 
        y = c(APP_1_adj_quant[1,4], APP_1_adj_quant[1,4] )*N_tot_app_1/18,
        type='l', col=pos3_col, lwd=2) 





####################
####################
###              ###
###  PANEL 3     ###
###  Scenario 1  ###
###  Negatives   ###
###              ###
####################
####################

par(mar = c(5,5,0.5,0.5))

line_seq_y <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

plot(x=1000, y=1000, 
xlim=c(0,100), ylim=c(0,1000),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

polygon( x = c(10,90,90,10),
         y = c(0,0,N_neg_app_2,N_neg_app_2),
         col=rgb(211/255,211/255,211/255,0.5), border=NA )




polygon( x = c(0,100,100,0), 
         y = c(APP_2_adj_quant[2,1], APP_2_adj_quant[2,1], APP_2_adj_quant[3,1], APP_2_adj_quant[3,1] )*N_tot_app_2,
         col=rgb(30/256,144/256,255/256,0.25), border=NA )

points( x=c(0,100), y=c(APP_2_adj_quant[1,1], APP_2_adj_quant[1,1])*N_tot_app_2,
        type='l', col=neg_col) 


axis(1, at = c(0, 100), 
        label=c("", ""), 
        cex.axis=axis.size )


axis(2, at = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), 
        label = c("0", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3.33, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="negatives")


####################
####################
###              ###
###  PANEL 4     ###
###  Scenario 2  ###
###  Positives   ###
###              ###
####################
####################


par(mar = c(5,5,0.5,2))
par(mgp=c(3,1.25,0))


line_seq_x <- -30*c(0, 3, 6, 9, 12)
line_seq_y <- c(0, 10, 20, 30, 40, 50, 60, 70)

plot(x=1000, y=1000, 
xlim=c(-360,30), ylim=c(0,70),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}



for(i in 1:N_bins)
{
	polygon( x=-t_bins[c(i,i+1,i+1,i)], y=c(0,0,inf_bins_2[i],inf_bins_2[i]),
               col=rgb(211/255,211/255,211/255,0.5), border=NA )
}

axis(1, at = -30*c(-1, 0, 3, 6, 9, 12), 
        label=c("", "0", "3", "6", "9", "12"), 
        cex.axis=1.5*axis.size )


axis(2, at = c(0, 10, 20, 30, 40, 50, 60, 70), 
        label = c("0", "10", "20", "30", "40", "50", "60", "70"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="positives: time before sampling (months)")

arrows( x0=0, y0=40, x1=0, y1=25,
        angle=30, lwd=2, length=0.1)

text( x=0, y=43, 
      labels="sampling", 
      cex=2)

###################
## 0-3 months

polygon( x = c(0,-90,-90,0), 
         y = c(APP_2_adj_quant[2,2], APP_2_adj_quant[2,2], APP_2_adj_quant[3,2], APP_2_adj_quant[3,2] )*N_tot_app_2/9,
         col=rgb(255/256,69/256,0/256,0.25), border=NA )

points( x = c(0,-90), 
        y = c(APP_2_adj_quant[1,2], APP_2_adj_quant[1,2] )*N_tot_app_2/9,
        type='l', col=pos1_col, lwd=2) 


###################
## 3-6 months

polygon( x = c(-90,-180,-180,-90), 
         y = c(APP_2_adj_quant[2,3], APP_2_adj_quant[2,3], APP_2_adj_quant[3,3], APP_2_adj_quant[3,3] )*N_tot_app_2/9,
         col=rgb(34/256,139/256,34/256,0.25), border=NA )

points( x = c(-90,-180), 
        y = c(APP_2_adj_quant[1,3], APP_2_adj_quant[1,3] )*N_tot_app_2/9,
        type='l', col=pos2_col, lwd=2) 

###################
## 7-12 months

polygon( x = c(-180,-360,-360,-180), 
         y = c(APP_2_adj_quant[2,4], APP_2_adj_quant[2,4], APP_2_adj_quant[3,4], APP_2_adj_quant[3,4] )*N_tot_app_2/18,
         col=rgb(71/256,60/256,139/256,0.25), border=NA )

points( x = c(-180,-360), 
        y = c(APP_2_adj_quant[1,4], APP_2_adj_quant[1,4] )*N_tot_app_2/18,
        type='l', col=pos3_col, lwd=2) 




############################
## Labels in middle       ## 
############################

par(mar = c(0,0,0,0))

plot.new()
title( "(C) Double wave scenario", 
        cex.main=3.0, line=-3)

plot.new()
title( "(D) Constant transmission scenario", 
        cex.main=3.0, line=-3)



####################
####################
###              ###
###  PANEL 5     ###
###  Scenario 3  ###
###  Negatives   ###
###              ###
####################
####################

par(mar = c(5,5,0.5,0.5))

line_seq_y <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

plot(x=1000, y=1000, 
xlim=c(0,100), ylim=c(0,1000),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

polygon( x = c(10,90,90,10),
         y = c(0,0,N_neg_app_3,N_neg_app_3),
         col=rgb(211/255,211/255,211/255,0.5), border=NA )




polygon( x = c(0,100,100,0), 
         y = c(APP_3_adj_quant[2,1], APP_3_adj_quant[2,1], APP_3_adj_quant[3,1], APP_3_adj_quant[3,1] )*N_tot_app_3,
         col=rgb(30/256,144/256,255/256,0.25), border=NA )

points( x=c(0,100), y=c(APP_3_adj_quant[1,1], APP_3_adj_quant[1,1])*N_tot_app_3,
        type='l', col=neg_col) 


axis(1, at = c(0, 100), 
        label=c("", ""), 
        cex.axis=axis.size )


axis(2, at = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), 
        label = c("0", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3.33, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="negatives")


####################
####################
###              ###
###  PANEL 6     ###
###  Scenario 3  ###
###  Positives   ###
###              ###
####################
####################


par(mar = c(5,5,0.5,2))
par(mgp=c(3,1.25,0))


line_seq_x <- -30*c(0, 3, 6, 9, 12)
line_seq_y <- c(0, 10, 20, 30, 40, 50)

plot(x=1000, y=1000, 
xlim=c(-360,30), ylim=c(0,50),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}



for(i in 1:N_bins)
{
	polygon( x=-t_bins[c(i,i+1,i+1,i)], y=c(0,0,inf_bins_3[i],inf_bins_3[i]),
               col=rgb(211/255,211/255,211/255,0.5), border=NA )
}

axis(1, at = -30*c(-1, 0, 3, 6, 9, 12), 
        label=c("", "0", "3", "6", "9", "12"), 
        cex.axis=1.5*axis.size )


axis(2, at = c(0, 10, 20, 30, 40, 50), 
        label = c("0", "10", "20", "30", "40", "50"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="positives: time before sampling (months)")

arrows( x0=0, y0=40, x1=0, y1=25,
        angle=30, lwd=2, length=0.1)

text( x=0, y=43, 
      labels="sampling", 
      cex=2)

###################
## 0-3 months

polygon( x = c(0,-90,-90,0), 
         y = c(APP_3_adj_quant[2,2], APP_3_adj_quant[2,2], APP_3_adj_quant[3,2], APP_3_adj_quant[3,2] )*N_tot_app_3/9,
        col=rgb(255/256,69/256,0/256,0.25), border=NA )

points( x = c(0,-90), 
        y = c(APP_3_adj_quant[1,2], APP_3_adj_quant[1,2] )*N_tot_app_3/9,
        type='l', col=pos1_col, lwd=2) 


###################
## 3-6 months

polygon( x = c(-90,-180,-180,-90), 
         y = c(APP_3_adj_quant[2,3], APP_3_adj_quant[2,3], APP_3_adj_quant[3,3], APP_3_adj_quant[3,3] )*N_tot_app_3/9,
         col=rgb(34/256,139/256,34/256,0.25), border=NA )

points( x = c(-90,-180), 
        y = c(APP_3_adj_quant[1,3], APP_3_adj_quant[1,3] )*N_tot_app_3/9,
        type='l', col=pos2_col, lwd=2) 

###################
## 7-12 months

polygon( x = c(-180,-360,-360,-180), 
         y = c(APP_3_adj_quant[2,4], APP_3_adj_quant[2,4], APP_3_adj_quant[3,4], APP_3_adj_quant[3,4] )*N_tot_app_3/18,
         col=rgb(71/256,60/256,139/256,0.25), border=NA )

points( x = c(-180,-360), 
        y = c(APP_3_adj_quant[1,4], APP_3_adj_quant[1,4] )*N_tot_app_3/18,
        type='l', col=pos3_col, lwd=2) 





####################
####################
###              ###
###  PANEL 7     ###
###  Scenario 4  ###
###  Negatives   ###
###              ###
####################
####################

par(mar = c(5,5,0.5,0.5))

line_seq_y <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

plot(x=1000, y=1000, 
xlim=c(0,100), ylim=c(0,1000),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

polygon( x = c(10,90,90,10),
         y = c(0,0,N_neg_app_4,N_neg_app_4),
         col=rgb(211/255,211/255,211/255,0.5), border=NA )




polygon( x = c(0,100,100,0), 
         y = c(APP_4_adj_quant[2,1], APP_4_adj_quant[2,1], APP_4_adj_quant[3,1], APP_4_adj_quant[3,1] )*N_tot_app_4,
         col=rgb(30/256,144/256,255/256,0.25), border=NA )

points( x=c(0,100), y=c(APP_4_adj_quant[1,1], APP_4_adj_quant[1,1])*N_tot_app_4,
        type='l', col=neg_col) 


axis(1, at = c(0, 100), 
        label=c("", ""), 
        cex.axis=axis.size )


axis(2, at = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), 
        label = c("0", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3.33, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="negatives")


####################
####################
###              ###
###  PANEL 8     ###
###  Scenario 4  ###
###  Positives   ###
###              ###
####################
####################


par(mar = c(5,5,0.5,2))
par(mgp=c(3,1.25,0))


line_seq_x <- -30*c(0, 3, 6, 9, 12)
line_seq_y <- c(0, 10, 20, 30, 40, 50)

plot(x=1000, y=1000, 
xlim=c(-360,30), ylim=c(0,50),
yaxt='n', xaxt='n', bty='n', xaxs="i", yaxs="i",
ylab="", xlab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}



for(i in 1:N_bins)
{
	polygon( x=-t_bins[c(i,i+1,i+1,i)], y=c(0,0,inf_bins_4[i],inf_bins_4[i]),
               col=rgb(211/255,211/255,211/255,0.5), border=NA )
}

axis(1, at = -30*c(-1, 0, 3, 6, 9, 12), 
        label=c("", "0", "3", "6", "9", "12"), 
        cex.axis=1.5*axis.size )


axis(2, at = c(0, 10, 20, 30, 40, 50), 
        label = c("0", "10", "20", "30", "40", "50"),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="positives: time before sampling (months)")

arrows( x0=0, y0=40, x1=0, y1=25,
        angle=30, lwd=2, length=0.1)

text( x=0, y=43, 
      labels="sampling", 
      cex=2)

###################
## 0-3 months

polygon( x = c(0,-90,-90,0), 
         y = c(APP_4_adj_quant[2,2], APP_4_adj_quant[2,2], APP_4_adj_quant[3,2], APP_4_adj_quant[3,2] )*N_tot_app_4/9,
         col=rgb(255/256,69/256,0/256,0.25), border=NA )

points( x = c(0,-90), 
        y = c(APP_4_adj_quant[1,2], APP_4_adj_quant[1,2] )*N_tot_app_4/9,
        type='l', col=pos1_col, lwd=2) 


###################
## 3-6 months

polygon( x = c(-90,-180,-180,-90), 
         y = c(APP_4_adj_quant[2,3], APP_4_adj_quant[2,3], APP_4_adj_quant[3,3], APP_4_adj_quant[3,3] )*N_tot_app_4/9,
         col=rgb(34/256,139/256,34/256,0.25), border=NA )

points( x = c(-90,-180), 
        y = c(APP_4_adj_quant[1,3], APP_4_adj_quant[1,3] )*N_tot_app_4/9,
        type='l', col=pos2_col, lwd=2) 

###################
## 7-12 months

polygon( x = c(-180,-360,-360,-180), 
         y = c(APP_4_adj_quant[2,4], APP_4_adj_quant[2,4], APP_4_adj_quant[3,4], APP_4_adj_quant[3,4] )*N_tot_app_4/18,
         col=rgb(71/256,60/256,139/256,0.25), border=NA )

points( x = c(-180,-360), 
        y = c(APP_4_adj_quant[1,4], APP_4_adj_quant[1,4] )*N_tot_app_4/18,
        type='l', col=pos3_col, lwd=2) 



###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar=c(0,0,0,0))

plot.new()

legend(x='center', 
       legend = c("data (simulated)", "predicted: negative", "predicted: 0-3 months", "predicted: 3-6 months", "predicted: 6-12 months"), 
       fill = c(rgb(211/255,211/255,211/255,0.5), neg_col, pos1_col, pos2_col, pos3_col), 
       border = c(rgb(211/255,211/255,211/255,0.5), neg_col, pos1_col, pos2_col, pos3_col), 
       ncol=5, cex=2, bty="n" )



dev.off()




c( length(which( is.na(APP_1[,11])==TRUE )),
   length(which( APP_1[,11] > 0 & APP_1[,11] <= 90 )),
   length(which( APP_1[,11] > 90 & APP_1[,11] <= 180 )),
   length(which( APP_1[,11] > 180  ))  )/nrow(APP_1)

APP_1_meas_quant

APP_1_adj_quant




c( length(which( is.na(APP_2[,11])==TRUE )),
   length(which( APP_2[,11] > 0 & APP_2[,11] <= 90 )),
   length(which( APP_2[,11] > 90 & APP_2[,11] <= 180 )),
   length(which( APP_2[,11] > 180  ))  )/nrow(APP_2)

APP_2_meas_quant

APP_2_adj_quant



c( length(which( is.na(APP_3[,11])==TRUE )),
   length(which( APP_3[,11] > 0 & APP_3[,11] <= 90 )),
   length(which( APP_3[,11] > 90 & APP_3[,11] <= 180 )),
   length(which( APP_3[,11] > 180  ))  )/nrow(APP_3)

APP_3_meas_quant

APP_3_adj_quant





c( length(which( is.na(APP_4[,11])==TRUE )),
   length(which( APP_4[,11] > 0 & APP_4[,11] <= 90 )),
   length(which( APP_4[,11] > 90 & APP_4[,11] <= 180 )),
   length(which( APP_4[,11] > 180  ))  )/nrow(APP_1)

APP_4_meas_quant

APP_4_adj_quant

save.image("Fig4_sero_recon.RData")





