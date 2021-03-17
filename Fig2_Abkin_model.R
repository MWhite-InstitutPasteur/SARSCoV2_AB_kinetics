


cohort_names <- c("Pelleau et al.", "Iyer et al.", "WangTo et al.", "Tang et al.",
                  "Isho et al.", "Seow et al.", "Roltgen et al.", "Dan et al.")

cohort_cols  <- c("royalblue", "gold", "firebrick", "orangered", 
                  "forestgreen", "yellowgreen", "darkmagenta", "pink1")

cohort_rgb <- c( rgb(65/255,  105/255, 225/255, 1/3),
                 rgb(255/255, 215/255, 0/255,   1/3),
                 rgb(178/255, 34/255,  34/255,  1/3),
                 rgb(255/255, 69/255,  0/255,   1/3),
                 rgb(34/255,  139/255, 34/255,  1/3),
                 rgb(154/255, 205/255, 50/255, 1/3),
                 rgb(139/255, 0/255,   139/255, 1/3),
                 rgb(255/255, 181/255, 197/255, 1/3) )
 

###################################
###################################
##                               ##
##  ####    ####  ######  ####   ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ## ######   ##   ######  ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ####   ##  ##   ##   ##  ##  ##
##                               ##
###################################
###################################

tt_rdc <- 30*c(3, 6, 12, 18, 24)

###################
###             ###
###  Spike IgG  ###
###             ###
###################

Spike_IgG_cols <- cohort_cols[c(1,2,4,5,6,7,8)]
Spike_IgG_rgb  <- cohort_rgb[c(1,2,4,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\Spike_IgG.RData")

N_cohort_Spike_IgG <- length(Spike_IgG_N_part)

Spike_IgG_rdc <- list()

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_rdc[[c]] <- matrix(NA, nrow=dim(Spike_IgG_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		Spike_IgG_rdc[[c]][,j] <- (Spike_IgG_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - Spike_IgG_mod[[c]][,1,2]) / 
                                      (Spike_IgG_mod[[c]][,which(tt_plot==14),2] - Spike_IgG_mod[[c]][,1,2]) 

	}

}


Spike_IgG_mod_median <- matrix(NA, nrow=N_cohort_Spike_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_mod_median[c,] <- apply( X = (Spike_IgG_mod[[c]][,,2] - Spike_IgG_mod[[c]][,1,2]) / 
                                             (Spike_IgG_mod[[c]][,which(tt_plot==14),2] - Spike_IgG_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}


###################
###             ###
###  Spike IgM  ###
###             ###
###################

Spike_IgM_cols <- cohort_cols[c(1,2,5,6,7)]
Spike_IgM_rgb  <- cohort_rgb[c(1,2,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\Spike_IgM.RData")

N_cohort_Spike_IgM <- length(Spike_IgM_N_part)

Spike_IgM_rdc <- list()

for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_rdc[[c]] <- matrix(NA, nrow=dim(Spike_IgM_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		Spike_IgM_rdc[[c]][,j] <- (Spike_IgM_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - Spike_IgM_mod[[c]][,1,2]) / 
                                      (Spike_IgM_mod[[c]][,which(tt_plot==14),2] - Spike_IgM_mod[[c]][,1,2]) 

	}

}


Spike_IgM_mod_median <- matrix(NA, nrow=N_cohort_Spike_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_mod_median[c,] <- apply( X = (Spike_IgM_mod[[c]][,,2] - Spike_IgM_mod[[c]][,1,2]) / 
                                             (Spike_IgM_mod[[c]][,which(tt_plot==14),2] - Spike_IgM_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}



###################
###             ###
###  Spike IgA  ###
###             ###
###################

Spike_IgA_cols <- cohort_cols[c(1,2,5,6,7,8)]
Spike_IgA_rgb  <- cohort_rgb[c(1,2,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\Spike_IgA.RData")

N_cohort_Spike_IgA <- length(Spike_IgA_N_part)

Spike_IgA_rdc <- list()

for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_rdc[[c]] <- matrix(NA, nrow=dim(Spike_IgA_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		Spike_IgA_rdc[[c]][,j] <- (Spike_IgA_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - Spike_IgA_mod[[c]][,1,2]) / 
                                      (Spike_IgA_mod[[c]][,which(tt_plot==14),2] - Spike_IgA_mod[[c]][,1,2]) 

	}

}



Spike_IgA_mod_median <- matrix(NA, nrow=N_cohort_Spike_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_mod_median[c,] <- apply( X = (Spike_IgA_mod[[c]][,,2] - Spike_IgA_mod[[c]][,1,2]) / 
                                             (Spike_IgA_mod[[c]][,which(tt_plot==14),2] - Spike_IgA_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}



###################
###             ###
###  RBD IgG    ###
###             ###
###################

RBD_IgG_cols <- cohort_cols[c(1,3,4,5,6,7,8)]
RBD_IgG_rgb  <- cohort_rgb[c(1,3,4,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\RBD_IgG.RData")

N_cohort_RBD_IgG <- length(RBD_IgG_N_part)

RBD_IgG_rdc <- list()

for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_rdc[[c]] <- matrix(NA, nrow=dim(RBD_IgG_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		RBD_IgG_rdc[[c]][,j] <- (RBD_IgG_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - RBD_IgG_mod[[c]][,1,2]) / 
                                      (RBD_IgG_mod[[c]][,which(tt_plot==14),2] - RBD_IgG_mod[[c]][,1,2]) 

	}

}



RBD_IgG_mod_median <- matrix(NA, nrow=N_cohort_RBD_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_mod_median[c,] <- apply( X = (RBD_IgG_mod[[c]][,,2] - RBD_IgG_mod[[c]][,1,2]) / 
                                           (RBD_IgG_mod[[c]][,which(tt_plot==14),2] - RBD_IgG_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}



###################
###             ###
###  RBD IgM    ###
###             ###
###################

RBD_IgM_cols <- cohort_cols[c(1,3,5,6,7)]
RBD_IgM_rgb  <- cohort_rgb[c(1,3,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\RBD_IgM.RData")

N_cohort_RBD_IgM <- length(RBD_IgM_N_part)

RBD_IgM_rdc <- list()

for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_rdc[[c]] <- matrix(NA, nrow=dim(RBD_IgM_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		RBD_IgM_rdc[[c]][,j] <- (RBD_IgM_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - RBD_IgM_mod[[c]][,1,2]) / 
                                      (RBD_IgM_mod[[c]][,which(tt_plot==14),2] - RBD_IgM_mod[[c]][,1,2]) 

	}

}


RBD_IgM_mod_median <- matrix(NA, nrow=N_cohort_RBD_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_mod_median[c,] <- apply( X = (RBD_IgM_mod[[c]][,,2] - RBD_IgM_mod[[c]][,1,2]) / 
                                           (RBD_IgM_mod[[c]][,which(tt_plot==14),2] - RBD_IgM_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}


###################
###             ###
###  RBD IgA    ###
###             ###
###################

RBD_IgA_cols <- cohort_cols[c(1,5,6,7,8)]
RBD_IgA_rgb  <- cohort_rgb[c(1,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\RBD_IgA.RData")

N_cohort_RBD_IgA <- length(RBD_IgA_N_part)

RBD_IgA_rdc <- list()

for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_rdc[[c]] <- matrix(NA, nrow=dim(RBD_IgA_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		RBD_IgA_rdc[[c]][,j] <- (RBD_IgA_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - RBD_IgA_mod[[c]][,1,2]) / 
                                      (RBD_IgA_mod[[c]][,which(tt_plot==14),2] - RBD_IgA_mod[[c]][,1,2]) 

	}

}


RBD_IgA_mod_median <- matrix(NA, nrow=N_cohort_RBD_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_mod_median[c,] <- apply( X = (RBD_IgA_mod[[c]][,,2] - RBD_IgA_mod[[c]][,1,2]) / 
                                           (RBD_IgA_mod[[c]][,which(tt_plot==14),2] - RBD_IgA_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}


###################
###             ###
###  NP IgG     ###
###             ###
###################

NP_IgG_cols <- cohort_cols[c(1,3,5,6,7,8)]
NP_IgG_rgb  <- cohort_rgb[c(1,3,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\NP_IgG.RData")

N_cohort_NP_IgG <- length(NP_IgG_N_part)

NP_IgG_rdc <- list()

for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_rdc[[c]] <- matrix(NA, nrow=dim(NP_IgG_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		NP_IgG_rdc[[c]][,j] <- (NP_IgG_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - NP_IgG_mod[[c]][,1,2]) / 
                                      (NP_IgG_mod[[c]][,which(tt_plot==14),2] - NP_IgG_mod[[c]][,1,2]) 

	}

}


NP_IgG_mod_median <- matrix(NA, nrow=N_cohort_NP_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_mod_median[c,] <- apply( X = (NP_IgG_mod[[c]][,,2] - NP_IgG_mod[[c]][,1,2]) / 
                                           (NP_IgG_mod[[c]][,which(tt_plot==14),2] - NP_IgG_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}




###################
###             ###
###  NP IgM     ###
###             ###
###################

NP_IgM_cols <- cohort_cols[c(1,3,5,6,7)]
NP_IgM_rgb  <- cohort_rgb[c(1,3,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\NP_IgM.RData")

N_cohort_NP_IgM <- length(NP_IgM_N_part)


NP_IgM_rdc <- list()

for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_rdc[[c]] <- matrix(NA, nrow=dim(NP_IgM_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		NP_IgM_rdc[[c]][,j] <- (NP_IgM_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - NP_IgM_mod[[c]][,1,2]) / 
                                      (NP_IgM_mod[[c]][,which(tt_plot==14),2] - NP_IgM_mod[[c]][,1,2]) 

	}

}


NP_IgM_mod_median <- matrix(NA, nrow=N_cohort_NP_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_mod_median[c,] <- apply( X = (NP_IgM_mod[[c]][,,2] - NP_IgM_mod[[c]][,1,2]) / 
                                           (NP_IgM_mod[[c]][,which(tt_plot==14),2] - NP_IgM_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}





###################
###             ###
###  NP IgA     ###
###             ###
###################

NP_IgA_cols <- cohort_cols[c(1,5,6,7)]
NP_IgA_rgb  <- cohort_rgb[c(1,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\NP_IgA.RData")

N_cohort_NP_IgA <- length(NP_IgA_N_part)


NP_IgA_rdc <- list()

for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_rdc[[c]] <- matrix(NA, nrow=dim(NP_IgA_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		NP_IgA_rdc[[c]][,j] <- (NP_IgA_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - NP_IgA_mod[[c]][,1,2]) / 
                                      (NP_IgA_mod[[c]][,which(tt_plot==14),2] - NP_IgA_mod[[c]][,1,2]) 

	}

}



NP_IgA_mod_median <- matrix(NA, nrow=N_cohort_NP_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_mod_median[c,] <- apply( X = (NP_IgA_mod[[c]][,,2] - NP_IgA_mod[[c]][,1,2]) / 
                                           (NP_IgA_mod[[c]][,which(tt_plot==14),2] - NP_IgA_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}



###################
###             ###
###  S2 IgG     ###
###             ###
###################

S2_IgG_cols <- cohort_cols[c(1,3,5,6,7,8)]
S2_IgG_rgb  <- cohort_rgb[c(1,3,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\S2_IgG.RData")

N_cohort_S2_IgG <- length(S2_IgG_N_part)

S2_IgG_rdc <- list()

for(c in 1:N_cohort_S2_IgG)
{
	S2_IgG_rdc[[c]] <- matrix(NA, nrow=dim(S2_IgG_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		S2_IgG_rdc[[c]][,j] <- (S2_IgG_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - S2_IgG_mod[[c]][,1,2]) / 
                                      (S2_IgG_mod[[c]][,which(tt_plot==14),2] - S2_IgG_mod[[c]][,1,2]) 

	}

}


S2_IgG_mod_median <- matrix(NA, nrow=N_cohort_S2_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_S2_IgG)
{
	S2_IgG_mod_median[c,] <- apply( X = (S2_IgG_mod[[c]][,,2] - S2_IgG_mod[[c]][,1,2]) / 
                                           (S2_IgG_mod[[c]][,which(tt_plot==14),2] - S2_IgG_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}




###################
###             ###
###  S2IgM      ###
###             ###
###################

S2_IgM_cols <- cohort_cols[c(1,3,5,6,7)]
S2_IgM_rgb  <- cohort_rgb[c(1,3,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\S2_IgM.RData")

N_cohort_S2_IgM <- length(S2_IgM_N_part)


S2_IgM_rdc <- list()

for(c in 1:N_cohort_S2_IgM)
{
	S2_IgM_rdc[[c]] <- matrix(NA, nrow=dim(S2_IgM_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		S2_IgM_rdc[[c]][,j] <- (S2_IgM_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - S2_IgM_mod[[c]][,1,2]) / 
                                      (S2_IgM_mod[[c]][,which(tt_plot==14),2] - S2_IgM_mod[[c]][,1,2]) 

	}

}


S2_IgM_mod_median <- matrix(NA, nrow=N_cohort_S2_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_S2_IgM)
{
	S2_IgM_mod_median[c,] <- apply( X = (S2_IgM_mod[[c]][,,2] - S2_IgM_mod[[c]][,1,2]) / 
                                           (S2_IgM_mod[[c]][,which(tt_plot==14),2] - S2_IgM_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}





###################
###             ###
###  S2 IgA     ###
###             ###
###################

S2_IgA_cols <- cohort_cols[c(1,5,6,7)]
S2_IgA_rgb  <- cohort_rgb[c(1,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\S2_IgA.RData")

N_cohort_S2_IgA <- length(S2_IgA_N_part)


S2_IgA_rdc <- list()

for(c in 1:N_cohort_S2_IgA)
{
	S2_IgA_rdc[[c]] <- matrix(NA, nrow=dim(S2_IgA_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		S2_IgA_rdc[[c]][,j] <- (S2_IgA_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - S2_IgA_mod[[c]][,1,2]) / 
                                      (S2_IgA_mod[[c]][,which(tt_plot==14),2] - S2_IgA_mod[[c]][,1,2]) 

	}

}



S2_IgA_mod_median <- matrix(NA, nrow=N_cohort_S2_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_S2_IgA)
{
	S2_IgA_mod_median[c,] <- apply( X = (S2_IgA_mod[[c]][,,2] - S2_IgA_mod[[c]][,1,2]) / 
                                           (S2_IgA_mod[[c]][,which(tt_plot==14),2] - S2_IgA_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}



###################
###             ###
###  ME IgG     ###
###             ###
###################

ME_IgG_cols <- cohort_cols[c(1,3,5,6,7,8)]
ME_IgG_rgb  <- cohort_rgb[c(1,3,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\ME_IgG.RData")

N_cohort_ME_IgG <- length(ME_IgG_N_part)

ME_IgG_rdc <- list()

for(c in 1:N_cohort_ME_IgG)
{
	ME_IgG_rdc[[c]] <- matrix(NA, nrow=dim(ME_IgG_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		ME_IgG_rdc[[c]][,j] <- (ME_IgG_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - ME_IgG_mod[[c]][,1,2]) / 
                                      (ME_IgG_mod[[c]][,which(tt_plot==14),2] - ME_IgG_mod[[c]][,1,2]) 

	}

}


ME_IgG_mod_median <- matrix(NA, nrow=N_cohort_ME_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_ME_IgG)
{
	ME_IgG_mod_median[c,] <- apply( X = (ME_IgG_mod[[c]][,,2] - ME_IgG_mod[[c]][,1,2]) / 
                                           (ME_IgG_mod[[c]][,which(tt_plot==14),2] - ME_IgG_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}




###################
###             ###
###  ME IgM     ###
###             ###
###################

ME_IgM_cols <- cohort_cols[c(1,3,5,6,7)]
ME_IgM_rgb  <- cohort_rgb[c(1,3,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\ME_IgM.RData")

N_cohort_ME_IgM <- length(ME_IgM_N_part)


ME_IgM_rdc <- list()

for(c in 1:N_cohort_ME_IgM)
{
	ME_IgM_rdc[[c]] <- matrix(NA, nrow=dim(ME_IgM_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		ME_IgM_rdc[[c]][,j] <- (ME_IgM_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - ME_IgM_mod[[c]][,1,2]) / 
                                      (ME_IgM_mod[[c]][,which(tt_plot==14),2] - ME_IgM_mod[[c]][,1,2]) 

	}

}


ME_IgM_mod_median <- matrix(NA, nrow=N_cohort_ME_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_ME_IgM)
{
	ME_IgM_mod_median[c,] <- apply( X = (ME_IgM_mod[[c]][,,2] - ME_IgM_mod[[c]][,1,2]) / 
                                           (ME_IgM_mod[[c]][,which(tt_plot==14),2] - ME_IgM_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}





###################
###             ###
###  ME IgA     ###
###             ###
###################

ME_IgA_cols <- cohort_cols[c(1,5,6,7)]
ME_IgA_rgb  <- cohort_rgb[c(1,5,6,7)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\ME_IgA.RData")

N_cohort_ME_IgA <- length(ME_IgA_N_part)


ME_IgA_rdc <- list()

for(c in 1:N_cohort_ME_IgA)
{
	ME_IgA_rdc[[c]] <- matrix(NA, nrow=dim(ME_IgA_mod[[c]])[1], ncol=length(tt_rdc))	

	for(j in 1:length(tt_rdc))
	{
		ME_IgA_rdc[[c]][,j] <- (ME_IgA_mod[[c]][,which(tt_plot==tt_rdc[j]),2] - ME_IgA_mod[[c]][,1,2]) / 
                                      (ME_IgA_mod[[c]][,which(tt_plot==14),2] - ME_IgA_mod[[c]][,1,2]) 

	}

}



ME_IgA_mod_median <- matrix(NA, nrow=N_cohort_ME_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_ME_IgA)
{
	ME_IgA_mod_median[c,] <- apply( X = (ME_IgA_mod[[c]][,,2] - ME_IgA_mod[[c]][,1,2]) / 
                                           (ME_IgA_mod[[c]][,which(tt_plot==14),2] - ME_IgA_mod[[c]][,1,2]), 
                                         FUN=median, MARGIN=2, na.rm=TRUE)
}





#########################################
#########################################
##                                     ##  
##  #####  ##### ####   ##  ##  ####   ##
##  ##  ## ##    ## ##  ##  ## ##  ##  ##
##  #####  ####  ##  ## ##  ## ##      ##
##  ## ##  ##    ## ##  ##  ## ##  ##  ## 
##  ##  ## ##### ####    ####   ####   ##
##                                     ##
#########################################
#########################################


Spike_IgG_j <- c()

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_j <- c(Spike_IgG_j, Spike_IgG_rdc[[c]][,2] )
}

Spike_IgG_j <- Spike_IgG_j[which(is.na(Spike_IgG_j)==FALSE)] 
Spike_IgG_j <- Spike_IgG_j[which(Spike_IgG_j<2)] 

quantile( Spike_IgG_j, prob=c(0.5, 0.025, 0.975) )




RBD_IgG_j <- c()

for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_j <- c(RBD_IgG_j, RBD_IgG_rdc[[c]][,2] )
}

RBD_IgG_j <- RBD_IgG_j[which(is.na(RBD_IgG_j)==FALSE)] 
RBD_IgG_j <- RBD_IgG_j[which(RBD_IgG_j<2)] 

quantile( RBD_IgG_j, prob=c(0.5, 0.025, 0.975) )



NP_IgG_j <- c()

for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_j <- c(NP_IgG_j, NP_IgG_rdc[[c]][,2] )
}

NP_IgG_j <- NP_IgG_j[which(is.na(NP_IgG_j)==FALSE)] 
NP_IgG_j <- NP_IgG_j[which(NP_IgG_j<2)] 

quantile( NP_IgG_j, prob=c(0.5, 0.025, 0.975) )





Spike_IgM_j <- c()

for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_j <- c(Spike_IgM_j, Spike_IgM_rdc[[c]][,2] )
}

Spike_IgM_j <- Spike_IgM_j[which(is.na(Spike_IgM_j)==FALSE)] 
Spike_IgM_j <- Spike_IgM_j[which(Spike_IgM_j<2)] 

quantile( Spike_IgM_j, prob=c(0.5, 0.025, 0.975) )




RBD_IgM_j <- c()

for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_j <- c(RBD_IgM_j, RBD_IgM_rdc[[c]][,2] )
}

RBD_IgM_j <- RBD_IgM_j[which(is.na(RBD_IgM_j)==FALSE)] 
RBD_IgM_j <- RBD_IgM_j[which(RBD_IgM_j<2)] 

quantile( RBD_IgM_j, prob=c(0.5, 0.025, 0.975) )



NP_IgM_j <- c()

for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_j <- c(NP_IgM_j, NP_IgM_rdc[[c]][,2] )
}

NP_IgM_j <- NP_IgM_j[which(is.na(NP_IgM_j)==FALSE)] 
NP_IgM_j <- NP_IgM_j[which(NP_IgM_j<2)] 

quantile( NP_IgM_j, prob=c(0.5, 0.025, 0.975) )




Spike_IgA_j <- c()

for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_j <- c(Spike_IgA_j, Spike_IgA_rdc[[c]][,2] )
}

Spike_IgA_j <- Spike_IgA_j[which(is.na(Spike_IgA_j)==FALSE)] 
Spike_IgA_j <- Spike_IgA_j[which(Spike_IgA_j<2)] 

quantile( Spike_IgA_j, prob=c(0.5, 0.025, 0.975) )





RBD_IgA_j <- c()

for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_j <- c(RBD_IgA_j, RBD_IgA_rdc[[c]][,2] )
}

RBD_IgA_j <- RBD_IgA_j[which(is.na(RBD_IgA_j)==FALSE)] 
RBD_IgA_j <- RBD_IgA_j[which(RBD_IgA_j<2)] 

quantile( RBD_IgA_j, prob=c(0.5, 0.025, 0.975) )




NP_IgA_j <- c()

for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_j <- c(NP_IgA_j, NP_IgA_rdc[[c]][,2] )
}

NP_IgA_j <- NP_IgA_j[which(is.na(NP_IgA_j)==FALSE)] 
NP_IgA_j <- NP_IgA_j[which(NP_IgA_j<2)] 

quantile( NP_IgA_j, prob=c(0.5, 0.025, 0.975) )








S2_IgG_j <- c()

for(c in 1:N_cohort_S2_IgG)
{
	S2_IgG_j <- c(S2_IgG_j, S2_IgG_rdc[[c]][,5] )
}

S2_IgG_j <- S2_IgG_j[which(is.na(S2_IgG_j)==FALSE)] 
S2_IgG_j <- S2_IgG_j[which(S2_IgG_j<2)] 

quantile( S2_IgG_j, prob=c(0.5, 0.025, 0.975) )



S2_IgM_j <- c()

for(c in 1:N_cohort_S2_IgM)
{
	S2_IgM_j <- c(S2_IgM_j, S2_IgM_rdc[[c]][,5] )
}

S2_IgM_j <- S2_IgM_j[which(is.na(S2_IgM_j)==FALSE)] 
S2_IgM_j <- S2_IgM_j[which(S2_IgM_j<2)] 

quantile( S2_IgM_j, prob=c(0.5, 0.025, 0.975) )



S2_IgA_j <- c()

for(c in 1:N_cohort_S2_IgA)
{
	S2_IgA_j <- c(S2_IgA_j, S2_IgA_rdc[[c]][,5] )
}

S2_IgA_j <- S2_IgA_j[which(is.na(S2_IgA_j)==FALSE)] 
S2_IgA_j <- S2_IgA_j[which(S2_IgA_j<2)] 

quantile( S2_IgA_j, prob=c(0.5, 0.025, 0.975) )






ME_IgG_j <- c()

for(c in 1:N_cohort_ME_IgG)
{
	ME_IgG_j <- c(ME_IgG_j, ME_IgG_rdc[[c]][,5] )
}

ME_IgG_j <- ME_IgG_j[which(is.na(ME_IgG_j)==FALSE)] 
ME_IgG_j <- ME_IgG_j[which(ME_IgG_j<2)] 

quantile( ME_IgG_j, prob=c(0.5, 0.025, 0.975) )



ME_IgM_j <- c()

for(c in 1:N_cohort_ME_IgM)
{
	ME_IgM_j <- c(ME_IgM_j, ME_IgM_rdc[[c]][,5] )
}

ME_IgM_j <- ME_IgM_j[which(is.na(ME_IgM_j)==FALSE)] 
ME_IgM_j <- ME_IgM_j[which(ME_IgM_j<2)] 

quantile( ME_IgM_j, prob=c(0.5, 0.025, 0.975) )



ME_IgA_j <- c()

for(c in 1:N_cohort_ME_IgA)
{
	ME_IgA_j <- c(ME_IgA_j, ME_IgA_rdc[[c]][,5] )
}

ME_IgA_j <- ME_IgA_j[which(is.na(ME_IgA_j)==FALSE)] 
ME_IgA_j <- ME_IgA_j[which(ME_IgA_j<2)] 

quantile( ME_IgA_j, prob=c(0.5, 0.025, 0.975) )








###################################
###################################
##                               ##  
##  ######  ####  #   ##  ####   ##
##    ##   ##  ## ##  ## ##      ##
##    ##   ###### ###### ## ###  ##
##    ##   ##  ## ## ### ##  ##  ## 
##    ##   ##  ## ##  ##  ####   ##
##                               ##
###################################
###################################

batch_file <- read.table( "C:/U/CoronaVirus/PELLEAU_et_al/5_kinetics/2_model/BATCH_files/RBD_IgG.bat" )

N_cohort <- batch_file[[3]]

if( length(batch_file) != 6*N_cohort + 4 )
{
	print( "Check batch file entries" )
}

prior_file <- paste( batch_file[[2]] )
pop_file   <- paste( batch_file[[6*N_cohort + 4]] )


N_tt         <- rep(NA, N_cohort)
N_part       <- rep(NA, N_cohort)
N_neg        <- rep(NA, N_cohort)
ind_pos_file <- rep(NA, N_cohort)
ind_neg_file <- rep(NA, N_cohort)
ind_out_file <- rep(NA, N_cohort)

for( c in 1:N_cohort )
{
	N_tt[c]         <-        batch_file[[4 + 6*(c-1)]] 
	N_part[c]       <-        batch_file[[5 + 6*(c-1)]] 
	N_neg[c]        <-        batch_file[[6 + 6*(c-1)]] 
	ind_pos_file[c] <- paste( batch_file[[7 + 6*(c-1)]] )
	ind_neg_file[c] <- paste( batch_file[[8 + 6*(c-1)]] )
	ind_out_file[c] <- paste( batch_file[[9 + 6*(c-1)]] )
}


#####################################
## Cohort names

glob_prior_read <- read.table( prior_file )

cohort_names <- glob_prior_read[seq(from=1, by=4, length=N_cohort),2]
cohort_names <- as.vector(cohort_names)

for(c in 1:N_cohort)
{
	cohort_names[c] <- substr( cohort_names[c], start= 15, stop=nchar(cohort_names[c]) )
}


tt_plot_long <- seq(from=-30, to=2500, by=10)

N_tt_plot_long <- length(tt_plot_long)

AB_mod_quant <- array(NA, dim=c(N_part[3], length(tt_plot_long), 3) )


ind_out_file[3] <- "C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\2_model\\OUTPUT\\RBD_IgG\\Tang_RBD_IgG_local_10m.txt"

MCMC_ind <- read.table( ind_out_file[3] )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]


for(n in 1:N_part[3] )
{
	###################################
	## Model prediction for participant i
	## Posterior projections

	AB_mod_n <- matrix(NA, nrow=nrow(MCMC_ind), ncol=N_tt_plot_long)
	
	for(i in 1:nrow(AB_mod_n))
	{
		Ab_0    = MCMC_ind[i,(n-1)*8+1]
		beta    = MCMC_ind[i,(n-1)*8+2]
		tau     = MCMC_ind[i,(n-1)*8+3]		
		t_short = MCMC_ind[i,(n-1)*8+4]
		t_long  = MCMC_ind[i,(n-1)*8+5]
		t_IgG   = MCMC_ind[i,(n-1)*8+6]
		rho     = MCMC_ind[i,(n-1)*8+7]

		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		##AB_tt = Ab_0*exp(-r_cl*tt_plot_long)
		AB_tt =  rep( Ab_0, length(tt_plot_long) )

		tt_temp = (tt_plot_long - tau)[which(tt_plot_long >= tau )]
		AB_tt[which(tt_plot_long >= tau)] = AB_tt[which(tt_plot_long >= tau)] +
          	                                      beta*(         (rho/(r_a - r_cs))*( exp(-r_cs*tt_temp) - exp(-r_a*tt_temp) ) +
		                                       ((1.0 - rho)/(r_a - r_cl))*( exp(-r_cl*tt_temp) - exp(-r_a*tt_temp) ) )
		AB_mod_n[i,] = AB_tt
	}

	for(j in 1:N_tt_plot_long)
	{
		AB_mod_quant[n,j,] <- quantile( AB_mod_n[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}

AB_Tang_lwr <- AB_mod_quant[1,,1]
AB_Tang_med <- AB_mod_quant[1,,2]
AB_Tang_upr <- AB_mod_quant[1,,3]


##################################
##################################
##                              ##  
##  #####  ##     ####  ######  ##
##  ##  ## ##    ##  ##   ##    ##
##  #####  ##    ##  ##   ##    ##
##  ##     ##    ##  ##   ##    ## 
##  ##     #####  ####    ##    ##
##                              ##
##################################
##################################

lab.size1   = 1.5
axis.size1  = 1.5	
main.size1  = 1.5

line.size1  = 2
dash.line.size1 = 0.75


pos.point.size = 0.3
neg.point.size = 0.3
lab.size2   = 1.5
axis.size2  = 1.25	
main.size2  = 1.5

line.size2  = 2
dash.line.size2 = 0.75


ww <- 6



n_Pelleau <- 387
n_Iyer <- 316
n_WangTo <- 1
n_Tang <- 1
n_Isho <- 7
n_Seow <- 57
n_Roltgen <- 77
n_Dan <- 8






tiff(file="Fig2_Abkin_model.tif", width=40, height=26, units="cm", res=500)


lay.mat <- rbind( c( 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,  7, 7, 7,  8, 8, 8, 9, 9, 9 ), 
                  c( 1, 10,10,10, 10,10,10, 10,10,11, 11,11,11, 11,11,11, 11,12,12, 12,12,12, 12,12,12 ),
                  c(13, 13,13,13, 13,13,13, 13,13,14, 14,14,14, 14,14,14, 14,15,15, 15,15,15, 15,15,15 ), 
                  c(16, 16,16,16, 16,16,16, 16,16,17, 17,17,17, 17,17,17, 17,18,18, 18,18,18, 18,18,18 ), 
                  c(19, 19,19,19, 19,19,19, 19,19,20, 20,20,20, 20,20,20, 20,21,21, 21,21,21, 21,21,21 ), 
                  c(22, 22,22,22, 22,22,22, 22,22,22, 22,22,22, 22,22,22, 22,22,22, 22,22,22, 22,22,22 ) )
layout(lay.mat, heights = c(6,3,10,10,11,1.5), 
widths = c(1, 3,3,3, 3,3,3, 3,3,3, 3,3,3, 3,3,3, 3,3,3, 3,3,3, 3,3,3))
layout.show(22)


############################
## Empty corner panel     ## 
############################

par(mar=c(0,0,0,0))

plot.new()


###################
###################
###             ###
###  PANEL 1    ###
###  Pelleau    ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- Spike_IgG_data[[1]][n_Pelleau,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- Spike_IgG_data[[1]][n_Pelleau,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- Spike_IgG_mod[[1]][n_Pelleau,,2]
AB_mod_plot_lwr <- Spike_IgG_mod[[1]][n_Pelleau,,1]
AB_mod_plot_upr <- Spike_IgG_mod[[1]][n_Pelleau,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[1], 
xlim=c(0,270), ylim=c(9e-6,0.03), log="y",
yaxt='n', xaxt='n', bty='n',
ylab="Spike IgG", xlab="",
main="Pelleau et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)


line_seq_x <- 30*c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), 
        label=c("0", "", "", "3", "", "", "6", "", "", "9"), 
        cex.axis=0.7*axis.size1 )

line_seq_y <- c(9e-6,0.03)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )

mtext(side = 1, line = 2.5, 
cex=0.8*lab.size1, 
text="time (months)")


points(x=tt_plot, y=AB_mod_plot_med, 
       type='l', col=cohort_cols[1])

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( AB_mod_plot_lwr, rev(AB_mod_plot_upr) ),
        col=cohort_rgb[1], border=NA)


###################
###################
###             ###
###  PANEL 2    ###
###  Iyer       ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- Spike_IgG_data[[2]][n_Iyer,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- Spike_IgG_data[[2]][n_Iyer,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- Spike_IgG_mod[[2]][n_Iyer,,2]
AB_mod_plot_lwr <- Spike_IgG_mod[[2]][n_Iyer,,1]
AB_mod_plot_upr <- Spike_IgG_mod[[2]][n_Iyer,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[2], 
xlim=c(0,100), ylim=c(0.01,100), log="y",
yaxt='n', xaxt='n', bty='n',
ylab="Spike IgG", xlab="",
main="Iyer et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)



line_seq_x <- 30*c(0, 1, 2, 3)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 1, 2, 3), 
        label=c(0, 1, 2, 3), 
        cex.axis=0.8*axis.size1 )

line_seq_y <- c(0.01,100)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )

mtext(side = 1, line = 2.5, 
cex=0.8*lab.size1, 
text="time (months)")

points(x=tt_plot, y=AB_mod_plot_med, 
       type='l', col=cohort_cols[2])

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( AB_mod_plot_lwr, rev(AB_mod_plot_upr) ),
        col=cohort_rgb[2], border=NA)




###################
###################
###             ###
###  PANEL 3    ###
###  WangTo     ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- RBD_IgG_data[[2]][n_WangTo,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- RBD_IgG_data[[2]][n_WangTo,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- RBD_IgG_mod[[2]][n_WangTo,,2]
AB_mod_plot_lwr <- RBD_IgG_mod[[2]][n_WangTo,,1]
AB_mod_plot_upr <- RBD_IgG_mod[[2]][n_WangTo,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[3], 
xlim=c(0,30), ylim=c(0.001,0.1), log="y",
yaxt='n', xaxt='n', bty='n',
ylab="RBD IgG", xlab="",
main="WangTo et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)




line_seq_x <- 30*c(0, 1/3, 2/3, 1)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 1/3, 2/3, 1), 
        label=c("0", "", "", "1"), 
        cex.axis=0.8*axis.size1 )

line_seq_y <- c(0.001,0.1)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )


mtext(side = 1, line = 2.5, 
cex=0.8*lab.size1, 
text="time (months)")

points(x=tt_plot, y=AB_mod_plot_med, 
       type='l', col=cohort_cols[3])

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( AB_mod_plot_lwr, rev(AB_mod_plot_upr) ),
        col=cohort_rgb[3], border=NA)



###################
###################
###             ###
###  PANEL 4    ###
###  Tang       ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- Spike_IgG_data[[3]][n_Tang,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- Spike_IgG_data[[3]][n_Tang,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- Spike_IgG_mod[[3]][n_Tang,,2]
AB_mod_plot_lwr <- Spike_IgG_mod[[3]][n_Tang,,1]
AB_mod_plot_upr <- Spike_IgG_mod[[3]][n_Tang,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[4], 
xlim=c(0,2200), ylim=c(1,1000), log="y",
yaxt='n', xaxt='n', bty='n',
ylab="", xlab="",
main="Tang et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)



line_seq_x <- 30*c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72), 
        label=c("0", "", "", "", "24", "", "", "", "48", "", "", "", "72"), 
        cex.axis=0.8*axis.size1 )

line_seq_y <- c(1,1000)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )


mtext(side = 2, line = 1.8, 
cex=0.9, 
text="SARS-CoV-1")

mtext(side = 2, line = 0.7, 
cex=0.9, 
text="Spike IgG")


mtext(side = 1, line = 3, 
cex=0.8*lab.size1, 
text="time (months)")

points(x=tt_plot_long, y=AB_Tang_med, 
       type='l', col=cohort_cols[4])

polygon(x=c(tt_plot_long, rev(tt_plot_long)), 
	  y=c( AB_Tang_lwr, rev(AB_Tang_upr) ),
        col=cohort_rgb[4], border=NA)



###################
###################
###             ###
###  PANEL 5    ###
###  Isho       ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- Spike_IgG_data[[4]][n_Isho,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- Spike_IgG_data[[4]][n_Isho,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- Spike_IgG_mod[[4]][n_Isho,,2]
AB_mod_plot_lwr <- Spike_IgG_mod[[4]][n_Isho,,1]
AB_mod_plot_upr <- Spike_IgG_mod[[4]][n_Isho,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[5], 
xlim=c(0,150), ylim=c(0.01,10), log="y",,
yaxt='n', xaxt='n', bty='n',
ylab="Spike IgG", xlab="",
main="Isho et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)



line_seq_x <- 30*c(0, 1, 2, 3, 4, 5)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 1, 2, 3, 4, 5), 
        label=c(0, 1, 2, 3, 4, 5), 
        cex.axis=0.8*axis.size1 )

line_seq_y <- c(0.01,10)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )

mtext(side = 1, line = 2.5, 
cex=0.8*lab.size1, 
text="time (months)")


points(x=tt_plot, y=AB_mod_plot_med, 
       type='l', col=cohort_cols[5])

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( AB_mod_plot_lwr, rev(AB_mod_plot_upr) ),
        col=cohort_rgb[5], border=NA)


###################
###################
###             ###
###  PANEL 6    ###
###  Seow       ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- Spike_IgG_data[[5]][n_Seow,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- Spike_IgG_data[[5]][n_Seow,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- Spike_IgG_mod[[5]][n_Seow,,2]
AB_mod_plot_lwr <- Spike_IgG_mod[[5]][n_Seow,,1]
AB_mod_plot_upr <- Spike_IgG_mod[[5]][n_Seow,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[6], 
xlim=c(0,90), ylim=c(0.01,10), log="y",
yaxt='n', xaxt='n', bty='n',
ylab="Spike IgG", xlab="",
main="Seow et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)



line_seq_x <- 30*c(0, 1, 2, 3)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 1, 2, 3), 
        label=c(0, 1, 2, 3), 
        cex.axis=0.8*axis.size1 )

line_seq_y <- c(0.01,10)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )

mtext(side = 1, line = 2.5, 
cex=0.8*lab.size1, 
text="time (months)")

points(x=tt_plot, y=AB_mod_plot_med, 
       type='l', col=cohort_cols[6])

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( AB_mod_plot_lwr, rev(AB_mod_plot_upr) ),
        col=cohort_rgb[6], border=NA)



###################
###################
###             ###
###  PANEL 7    ###
###  Roltgen    ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- Spike_IgG_data[[6]][n_Roltgen,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- Spike_IgG_data[[6]][n_Roltgen,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- Spike_IgG_mod[[6]][n_Roltgen,,2]
AB_mod_plot_lwr <- Spike_IgG_mod[[6]][n_Roltgen,,1]
AB_mod_plot_upr <- Spike_IgG_mod[[6]][n_Roltgen,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[7], 
xlim=c(0,60), ylim=c(0.01,10), log="y",
yaxt='n', xaxt='n', bty='n',
ylab="Spike IgG", xlab="",
main="Roltgen et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)



line_seq_x <- 30*c(0, 0.5, 1, 1.5, 2)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 0.5, 1, 1.5, 2), 
        label=c("0", "", "1", "", "2"), 
        cex.axis=0.8*axis.size1 )

line_seq_y <- c(0.01,10)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )

mtext(side = 1, line = 2.5, 
cex=0.8*lab.size1, 
text="time (months)")

points(x=tt_plot, y=AB_mod_plot_med, 
       type='l', col=cohort_cols[7])

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( AB_mod_plot_lwr, rev(AB_mod_plot_upr) ),
        col=cohort_rgb[7], border=NA)



###################
###################
###             ###
###  PANEL 8    ###
###  Dan        ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,3,2,1))
par(mgp=c(1.5, 1, 0))

AB_data_plot <- Spike_IgG_data[[7]][n_Dan,,1]
AB_data_plot <- AB_data_plot[which(is.na(AB_data_plot)==FALSE)]

tt_data_plot <- Spike_IgG_data[[7]][n_Dan,,2]
tt_data_plot <- tt_data_plot[which(is.na(tt_data_plot)==FALSE)]

AB_mod_plot_med <- Spike_IgG_mod[[7]][n_Dan,,2]
AB_mod_plot_lwr <- Spike_IgG_mod[[7]][n_Dan,,1]
AB_mod_plot_upr <- Spike_IgG_mod[[7]][n_Dan,,3]


plot( x=tt_data_plot, y=AB_data_plot, 
pch=19, cex=1, col=cohort_cols[8], 
xlim=c(0,240), ylim=c(100,10000), log="y",
yaxt='n', xaxt='n', bty='n',
ylab="Spike IgG", xlab="",
main="Dan et al.",
cex.lab=lab.size1, cex.axis=axis.size1, cex.main=main.size1)



line_seq_x <- 30*c(0, 1, 2, 3, 4, 5, 6, 7, 8)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size1, col="grey", lty="dashed")
}



axis(1, at = 30*c(0, 1, 2, 3, 4, 5, 6, 7, 8), 
        label=c("0", "", "2", "", "4", "", "6", "", "8"), 
        cex.axis=0.8*axis.size1 )

line_seq_y <- c(100,10000)

axis(2, at = line_seq_y, 
        label=c("", ""), 
        las=2, cex.axis=axis.size1 )

mtext(side = 1, line = 2.5, 
cex=0.8*lab.size1, 
text="time (months)")

points(x=tt_plot, y=AB_mod_plot_med, 
       type='l', col=cohort_cols[8])

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( AB_mod_plot_lwr, rev(AB_mod_plot_upr) ),
        col=cohort_rgb[8], border=NA)






############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))

plot.new()
title( "Spike", 
        cex.main=3.0, line=-3.5)

plot.new()
title( "RBD", 
        cex.main=3.0, line=-3.5)

plot.new()
title( "Nucleocapsid", 
        cex.main=3.0, line=-3.5)


###################
###################
###             ###
###  PANEL 9    ###
###  Spike IgG  ###
###             ###
###################
###################

line_seq_x <- c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30
line_seq_y <- c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5)


par(mar=c(3,8,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )


axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )


mtext(side = 2, line = 6, 
cex=lab.size2, 
text="%IgG remaining")

mtext(side = 2, line = 4, 
cex=lab.size2, 
text="(relative to day 14)")


###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_Spike_IgG)
	{
		Spike_IgG_j <- Spike_IgG_rdc[[c]][,j]	
		Spike_IgG_j <- Spike_IgG_j[which(is.na(Spike_IgG_j)==FALSE)] 
		Spike_IgG_j <- Spike_IgG_j[which(Spike_IgG_j<2)] 


		DEN_Spike_IgG_j = density( Spike_IgG_j, cut=0, na.rm=TRUE )

		DEN_Spike_IgG_j$y <- DEN_Spike_IgG_j$y/max(DEN_Spike_IgG_j$y)

		polygon( x = tt_rdc[j] - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_Spike_IgG_j$y, rev(DEN_Spike_IgG_j$y) ),
			   y=c( DEN_Spike_IgG_j$x, rev(DEN_Spike_IgG_j$x) ), 
	         	   col=Spike_IgG_rgb[c], border=NA)

		for(n in 1:length(Spike_IgG_j))
		{
			index <- which.min(abs(Spike_IgG_j[n] - DEN_Spike_IgG_j$x))

			points( x = tt_rdc[j] - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_Spike_IgG_j$y[index], max=DEN_Spike_IgG_j$y[index]),
      	              y = Spike_IgG_j[n],
      	              pch=19, cex=pos.point.size, col=Spike_IgG_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_Spike_IgG)
{
	points( x = tt_plot,
              y = Spike_IgG_mod_median[c,],
              col=Spike_IgG_cols[c], type='l')
}







###################
###################
###             ###
###  PANEL 10   ###
###  RBD IgG    ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )


axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )



###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_RBD_IgG)
	{
		RBD_IgG_j <- RBD_IgG_rdc[[c]][,j]	
		RBD_IgG_j <- RBD_IgG_j[which(is.na(RBD_IgG_j)==FALSE)] 
		RBD_IgG_j <- RBD_IgG_j[which(RBD_IgG_j<2)] 


		DEN_RBD_IgG_j = density( RBD_IgG_j, cut=0, na.rm=TRUE )

		DEN_RBD_IgG_j$y <- DEN_RBD_IgG_j$y/max(DEN_RBD_IgG_j$y)

		polygon( x = tt_rdc[j] - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_RBD_IgG_j$y, rev(DEN_RBD_IgG_j$y) ),
			   y=c( DEN_RBD_IgG_j$x, rev(DEN_RBD_IgG_j$x) ), 
	         	   col=RBD_IgG_rgb[c], border=NA)

		for(n in 1:length(RBD_IgG_j))
		{
			index <- which.min(abs(RBD_IgG_j[n] - DEN_RBD_IgG_j$x))

			points( x = tt_rdc[j] - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_RBD_IgG_j$y[index], max=DEN_RBD_IgG_j$y[index]),
      	              y = RBD_IgG_j[n],
      	              pch=19, cex=pos.point.size, col=RBD_IgG_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_RBD_IgG)
{
	points( x = tt_plot,
              y = RBD_IgG_mod_median[c,],
              col=RBD_IgG_cols[c], type='l')
}






###################
###################
###             ###
###  PANEL 11   ###
###  NP IgG     ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )


axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )


###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_NP_IgG)
	{
		NP_IgG_j <- NP_IgG_rdc[[c]][,j]	
		NP_IgG_j <- NP_IgG_j[which(is.na(NP_IgG_j)==FALSE)] 
		NP_IgG_j <- NP_IgG_j[which(NP_IgG_j<2)] 


		DEN_NP_IgG_j = density( NP_IgG_j, cut=0, na.rm=TRUE )

		DEN_NP_IgG_j$y <- DEN_NP_IgG_j$y/max(DEN_NP_IgG_j$y)

		polygon( x = tt_rdc[j] - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_NP_IgG_j$y, rev(DEN_NP_IgG_j$y) ),
			   y=c( DEN_NP_IgG_j$x, rev(DEN_NP_IgG_j$x) ), 
	         	   col=NP_IgG_rgb[c], border=NA)

		for(n in 1:length(NP_IgG_j))
		{
			index <- which.min(abs(NP_IgG_j[n] - DEN_NP_IgG_j$x))

			points( x = tt_rdc[j] - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_NP_IgG_j$y[index], max=DEN_NP_IgG_j$y[index]),
      	              y = NP_IgG_j[n],
      	              pch=19, cex=pos.point.size, col=NP_IgG_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_NP_IgG)
{
	points( x = tt_plot,
              y = NP_IgG_mod_median[c,],
              col=NP_IgG_cols[c], type='l')
}





###################
###################
###             ###
###  PANEL 12   ###
###  Spike IgM  ###
###             ###
###################
###################

par(mar=c(3,8,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )


axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )


mtext(side = 2, line = 6, 
cex=lab.size2, 
text="%IgM remaining")

mtext(side = 2, line = 4, 
cex=lab.size2, 
text="(relative to day 14)")


###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_Spike_IgM)
	{
		Spike_IgM_j <- Spike_IgM_rdc[[c]][,j]	
		Spike_IgM_j <- Spike_IgM_j[which(is.na(Spike_IgM_j)==FALSE)] 
		Spike_IgM_j <- Spike_IgM_j[which(Spike_IgM_j<2)] 


		DEN_Spike_IgM_j = density( Spike_IgM_j, cut=0, na.rm=TRUE )

		DEN_Spike_IgM_j$y <- DEN_Spike_IgM_j$y/max(DEN_Spike_IgM_j$y)

		polygon( x = tt_rdc[j] - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_Spike_IgM_j$y, rev(DEN_Spike_IgM_j$y) ),
			   y=c( DEN_Spike_IgM_j$x, rev(DEN_Spike_IgM_j$x) ), 
	         	   col=Spike_IgM_rgb[c], border=NA)

		for(n in 1:length(Spike_IgM_j))
		{
			index <- which.min(abs(Spike_IgM_j[n] - DEN_Spike_IgM_j$x))

			points( x = tt_rdc[j] - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_Spike_IgM_j$y[index], max=DEN_Spike_IgM_j$y[index]),
      	              y = Spike_IgM_j[n],
      	              pch=19, cex=pos.point.size, col=Spike_IgM_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_Spike_IgM)
{
	points( x = tt_plot,
              y = Spike_IgM_mod_median[c,],
              col=Spike_IgM_cols[c], type='l')
}





###################
###################
###             ###
###  PANEL 13   ###
###  RBD IgM    ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )

axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )



###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_RBD_IgM)
	{
		RBD_IgM_j <- RBD_IgM_rdc[[c]][,j]	
		RBD_IgM_j <- RBD_IgM_j[which(is.na(RBD_IgM_j)==FALSE)] 
		RBD_IgM_j <- RBD_IgM_j[which(RBD_IgM_j<2)] 


		DEN_RBD_IgM_j = density( RBD_IgM_j, cut=0, na.rm=TRUE )

		DEN_RBD_IgM_j$y <- DEN_RBD_IgM_j$y/max(DEN_RBD_IgM_j$y)

		polygon( x = tt_rdc[j] - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_RBD_IgM_j$y, rev(DEN_RBD_IgM_j$y) ),
			   y=c( DEN_RBD_IgM_j$x, rev(DEN_RBD_IgM_j$x) ), 
	         	   col=RBD_IgM_rgb[c], border=NA)

		for(n in 1:length(RBD_IgM_j))
		{
			index <- which.min(abs(RBD_IgM_j[n] - DEN_RBD_IgM_j$x))

			points( x = tt_rdc[j] - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_RBD_IgM_j$y[index], max=DEN_RBD_IgM_j$y[index]),
      	              y = RBD_IgM_j[n],
      	              pch=19, cex=pos.point.size, col=RBD_IgM_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_RBD_IgM)
{
	points( x = tt_plot,
              y = RBD_IgM_mod_median[c,],
              col=RBD_IgM_cols[c], type='l')
}



###################
###################
###             ###
###  PANEL 14   ###
###  NP IgM     ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )

axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )


###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_NP_IgM)
	{
		NP_IgM_j <- NP_IgM_rdc[[c]][,j]	
		NP_IgM_j <- NP_IgM_j[which(is.na(NP_IgM_j)==FALSE)] 
		NP_IgM_j <- NP_IgM_j[which(NP_IgM_j<2)] 


		DEN_NP_IgM_j = density( NP_IgM_j, cut=0, na.rm=TRUE )

		DEN_NP_IgM_j$y <- DEN_NP_IgM_j$y/max(DEN_NP_IgM_j$y)

		polygon( x = tt_rdc[j] - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_NP_IgM_j$y, rev(DEN_NP_IgM_j$y) ),
			   y=c( DEN_NP_IgM_j$x, rev(DEN_NP_IgM_j$x) ), 
	         	   col=NP_IgM_rgb[c], border=NA)

		for(n in 1:length(NP_IgM_j))
		{
			index <- which.min(abs(NP_IgM_j[n] - DEN_NP_IgM_j$x))

			points( x = tt_rdc[j] - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_NP_IgM_j$y[index], max=DEN_NP_IgM_j$y[index]),
      	              y = NP_IgM_j[n],
      	              pch=19, cex=pos.point.size, col=NP_IgM_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_NP_IgM)
{
	points( x = tt_plot,
              y = NP_IgM_mod_median[c,],
              col=NP_IgM_cols[c], type='l')
}



###################
###################
###             ###
###  PANEL 15   ###
###  Spike IgA  ###
###             ###
###################
###################

par(mar=c(6,8,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )


axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )


mtext(side = 2, line = 6, 
cex=lab.size2, 
text="%IgA remaining")

mtext(side = 2, line = 4, 
cex=lab.size2, 
text="(relative to day 14)")

mtext(side = 1, line = 3, 
cex=lab.size2, 
text="time post symptoms (months)")


###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_Spike_IgA)
	{
		Spike_IgA_j <- Spike_IgA_rdc[[c]][,j]	
		Spike_IgA_j <- Spike_IgA_j[which(is.na(Spike_IgA_j)==FALSE)] 
		Spike_IgA_j <- Spike_IgA_j[which(Spike_IgA_j<2)] 


		DEN_Spike_IgA_j = density( Spike_IgA_j, cut=0, na.rm=TRUE )

		DEN_Spike_IgA_j$y <- DEN_Spike_IgA_j$y/max(DEN_Spike_IgA_j$y)

		polygon( x = tt_rdc[j] - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_Spike_IgA_j$y, rev(DEN_Spike_IgA_j$y) ),
			   y=c( DEN_Spike_IgA_j$x, rev(DEN_Spike_IgA_j$x) ), 
	         	   col=Spike_IgA_rgb[c], border=NA)

		for(n in 1:length(Spike_IgA_j))
		{
			index <- which.min(abs(Spike_IgA_j[n] - DEN_Spike_IgA_j$x))

			points( x = tt_rdc[j] - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_Spike_IgA_j$y[index], max=DEN_Spike_IgA_j$y[index]),
      	              y = Spike_IgA_j[n],
      	              pch=19, cex=pos.point.size, col=Spike_IgA_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_Spike_IgA)
{
	points( x = tt_plot,
              y = Spike_IgA_mod_median[c,],
              col=Spike_IgA_cols[c], type='l')
}




###################
###################
###             ###
###  PANEL 16   ###
###  RBD IgA    ###
###             ###
###################
###################

par(mar=c(6,3,1,1))
par(mgp=c(2.5, 1, 0))



###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )


axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )

mtext(side = 1, line = 3, 
cex=lab.size2, 
text="time post symptoms (months)")


###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_RBD_IgA)
	{
		RBD_IgA_j <- RBD_IgA_rdc[[c]][,j]	
		RBD_IgA_j <- RBD_IgA_j[which(is.na(RBD_IgA_j)==FALSE)] 
		RBD_IgA_j <- RBD_IgA_j[which(RBD_IgA_j<2)] 


		DEN_RBD_IgA_j = density( RBD_IgA_j, cut=0, na.rm=TRUE )

		DEN_RBD_IgA_j$y <- DEN_RBD_IgA_j$y/max(DEN_RBD_IgA_j$y)

		polygon( x = tt_rdc[j] - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_RBD_IgA_j$y, rev(DEN_RBD_IgA_j$y) ),
			   y=c( DEN_RBD_IgA_j$x, rev(DEN_RBD_IgA_j$x) ), 
	         	   col=RBD_IgA_rgb[c], border=NA)

		for(n in 1:length(RBD_IgA_j))
		{
			index <- which.min(abs(RBD_IgA_j[n] - DEN_RBD_IgA_j$x))

			points( x = tt_rdc[j] - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_RBD_IgA_j$y[index], max=DEN_RBD_IgA_j$y[index]),
      	              y = RBD_IgA_j[n],
      	              pch=19, cex=pos.point.size, col=RBD_IgA_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_RBD_IgA)
{
	points( x = tt_plot,
              y = RBD_IgA_mod_median[c,],
              col=RBD_IgA_cols[c], type='l')
}




###################
###################
###             ###
###  PANEL 17   ###
###  NP IgA     ###
###             ###
###################
###################

par(mar=c(6,3,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(0,750), ylim=c(0,1.5), 
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size2, cex.axis=axis.size2, cex.main=main.size2)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,1000), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1e10), 
             type='l', lwd=dash.line.size2, col="grey", lty="dashed")
}




axis(1, at = c(0, 3, 6, 9, 12, 15, 18, 21, 24)*30, 
        label=c(0, 3, 6, 9, 12, 15, 18, 21, 24), 
        cex.axis=axis.size2 )


axis(2, at = c(0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), 
        label = c("0%", "25%", "50%", "75%", "100%", "125%", "150%"),
        las=2, cex.axis=axis.size2 )


mtext(side = 1, line = 3, 
cex=lab.size2, 
text="time post symptoms (months)")


###################
## Positive points
## Month 1

##for(j in 1:length(tt_rdc))
for(j in c(2,3,4,5))
{
	for(c in 1:N_cohort_NP_IgA)
	{
		NP_IgA_j <- NP_IgA_rdc[[c]][,j]	
		NP_IgA_j <- NP_IgA_j[which(is.na(NP_IgA_j)==FALSE)] 
		NP_IgA_j <- NP_IgA_j[which(NP_IgA_j<2)] 


		DEN_NP_IgA_j = density( NP_IgA_j, cut=0, na.rm=TRUE )

		DEN_NP_IgA_j$y <- DEN_NP_IgA_j$y/max(DEN_NP_IgA_j$y)

		polygon( x = tt_rdc[j] - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
	                   ww*c( -DEN_NP_IgA_j$y, rev(DEN_NP_IgA_j$y) ),
			   y=c( DEN_NP_IgA_j$x, rev(DEN_NP_IgA_j$x) ), 
	         	   col=NP_IgA_rgb[c], border=NA)

		for(n in 1:length(NP_IgA_j))
		{
			index <- which.min(abs(NP_IgA_j[n] - DEN_NP_IgA_j$x))

			points( x = tt_rdc[j] - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
      	                  ww*runif(1, min=-DEN_NP_IgA_j$y[index], max=DEN_NP_IgA_j$y[index]),
      	              y = NP_IgA_j[n],
      	              pch=19, cex=pos.point.size, col=NP_IgA_cols[c])
		}
	}
}


###################
## Model prediction

for(c in 1:N_cohort_NP_IgA)
{
	points( x = tt_plot,
              y = NP_IgA_mod_median[c,],
              col=NP_IgA_cols[c], type='l')
}




###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar=c(0,0,0,0))

plot.new()

legend(x='center', 
       legend = c("Pelleau et al.", "Iyer et al.", "WangTo et al.", "Tang et al.",
                  "Isho et al.", "Seow et al.", "Roltgen et al.", "Dan et al."), 
       col = cohort_cols, 
       pch=rep(19,8),
       ncol=8, cex=2.2, bty="n" )


dev.off()






