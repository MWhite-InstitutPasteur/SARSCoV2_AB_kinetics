#####################################
#####################################
##                                 ## 
##  ##### #### ##    #####  ####   ## 
##  ##     ##  ##    ##    ##      ## 
##  ####   ##  ##    ####   ####   ## 
##  ##     ##  ##    ##        ##  ## 
##  ##    #### ##### #####  ####   ## 
##                                 ## 
#####################################
#####################################

rm(list=ls())

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




##########################################################################
##########################################################################
##                                                                      ## 
##  #####   ####  #####  ##  ## ##     ####  ###### ####  ####  #   ##  ##
##  ##  ## ##  ## ##  ## ##  ## ##    ##  ##   ##    ##  ##  ## ##  ##  ##
##  #####  ##  ## #####  ##  ## ##    ######   ##    ##  ##  ## ### ##  ##
##  ##     ##  ## ##     ##  ## ##    ##  ##   ##    ##  ##  ## ## ###  ##
##  ##      ####  ##      ####  ##### ##  ##   ##   ####  ####  ##  ##  ##
##                                                                      ## 
##########################################################################
##########################################################################

##########################################################################
## Read in and process global/population chains

MCMC_pop <- read.table( pop_file )

for( c in 1:N_cohort )
{
	colnames(MCMC_pop)[c]                   <- paste( "AB_bg_", cohort_names[c], sep="" )
	colnames(MCMC_pop)[N_cohort + c]        <- paste( "beta_", cohort_names[c], sep="" )
	colnames(MCMC_pop)[2*N_cohort + 5 + c]  <- paste( "sig_AB_bg_", cohort_names[c], sep="" )
	colnames(MCMC_pop)[3*N_cohort + 5 + c]  <- paste( "sig_beta_", cohort_names[c], sep="" )
	colnames(MCMC_pop)[4*N_cohort + 10 + c] <- paste( "sig_obs_", cohort_names[c], sep="" )
}

colnames(MCMC_pop)[(2*N_cohort+1):(2*N_cohort+5)]       <- c("tau", "t_short", "t_long", "t_IgG", "rho")

colnames(MCMC_pop)[(4*N_cohort+5+1):(4*N_cohort+10)]    <- c( "sig_tau",  "sig_t_short", "sig_t_long", "sig_t_IgG", "sig_rho")

colnames(MCMC_pop)[(5*N_cohort+10+1):(5*N_cohort+10+2)] <- c( "likelihood", "prior")


N_par = ncol(MCMC_pop) - 2

N_glob_par = 2*N_cohort + 5


##########################################################################
## Convert from precision to standard deviation

for(k in (N_glob_par + 1):(2*N_glob_par) )
{
	MCMC_pop[,k] = 1/sqrt(MCMC_pop[,k])
}

##########################################################################
## Convert from log to linear scale

for(k in 1:(N_glob_par - 1))
{
	MCMC_pop[,k]            = exp( MCMC_pop[,k] + 0.5*MCMC_pop[,N_glob_par+k]^2 )
	MCMC_pop[,N_glob_par+k] = sqrt( exp(MCMC_pop[,N_glob_par+k]^2) -1 )*MCMC_pop[,k]
}


##########################################################################
## Convert from logit to linear scale

for(k in 1:nrow(MCMC_pop))
{
	f_m <- function(x)
	{ 
		0.3989423*exp( -0.5*( ( log(x/(1-x))-MCMC_pop[k,N_glob_par] )/MCMC_pop[k,2*N_glob_par] )^2 )/( MCMC_pop[k,2*N_glob_par]*(1-x) )
	}

	f_m2 <- function(x)
	{ 
		x*0.3989423*exp( -0.5*( (log(x/(1-x))-MCMC_pop[k,N_glob_par])/MCMC_pop[k,2*N_glob_par] )^2 )/( MCMC_pop[k,2*N_glob_par]*(1-x) )
	}
	
	tryCatch(
	{
		moment_1 <- integrate(f_m, lower=0, upper=1)$value
		moment_2 <- integrate(f_m2, lower=0, upper=1)$value
	}, error=function(e){ NULL }
	)

	MCMC_pop[k,N_glob_par] <- moment_1
	MCMC_pop[k,2*N_glob_par] <- sqrt( moment_2 - moment_1^2 )
}


##########################################################################
## Drop the first 10% for burn-in

MCMC_pop_burn <- MCMC_pop[floor(0.1*nrow(MCMC_pop)):(nrow(MCMC_pop)-1),]



##########################################################################
## Read in file for global/population priors

glob_prior_read <- read.table( prior_file )

glob_prior_mat <- matrix( glob_prior_read[1:(4*N_glob_par),1], nrow=N_glob_par, ncol=4, byrow=TRUE)
colnames(glob_prior_mat) <- c("mean", "mean_CV", "SD", "SD_CV")

rownames(glob_prior_mat) <- 1:nrow(glob_prior_mat)

for( c in 1:N_cohort )
{
	rownames(glob_prior_mat)[c]             <- paste( "AB_bg_", cohort_names[c], sep="" )
	rownames(glob_prior_mat)[N_cohort + c]  <- paste( "beta_", cohort_names[c], sep="" )
}

rownames(glob_prior_mat)[(2*N_cohort+1):(2*N_cohort+5)] <- c( "tau", "t_short", 
                                                              "t_long", "t_IgG", "rho")

glob_prior_mat <- cbind( glob_prior_mat, matrix(NA, nrow=N_glob_par, ncol=8) )
colnames(glob_prior_mat)[5:12] <- c("logN_mu", "logN_mu_CV", "logN_sig", "logN_sig_CV", 
                                    "prior_mu", "prior_tau", "prior_k", "prior_theta" )


obs_prior_mat <- matrix( glob_prior_read[(4*N_glob_par+1):(4*N_glob_par+2*N_cohort),1],
                         nrow=N_cohort, ncol=2, byrow=TRUE )

colnames(obs_prior_mat) <- c("SIG", "SIG_CV")
rownames(obs_prior_mat) <-cohort_names
 

##########################################################################
## Calculate parameters on log scale

for( p in 1:(N_glob_par-1) )
{
	glob_prior_mat[p,5] <- log( glob_prior_mat[p,1] / sqrt(1 + (glob_prior_mat[p,3]/glob_prior_mat[p,1])^2 ) )
	glob_prior_mat[p,6] <- glob_prior_mat[p,2]
	glob_prior_mat[p,7] <- sqrt(log( 1.0 + (glob_prior_mat[p,3]/glob_prior_mat[p,1])^2 ))
	glob_prior_mat[p,8] <- glob_prior_mat[p,4]
}

glob_prior_mat[N_glob_par,5] <- 2.625251
glob_prior_mat[N_glob_par,6] <- 2.5*glob_prior_mat[N_glob_par,2]
glob_prior_mat[N_glob_par,7] <- 1.096994
glob_prior_mat[N_glob_par,8] <- 2.5*glob_prior_mat[N_glob_par,4]


##########################################################################
## Calculate parameters of prior distributions
## Normal conjugate prior on mean
## Gamma conjugate prior on precision

for( p in 1:N_glob_par )
{
	glob_prior_mat[p,9]  = glob_prior_mat[p,5]
	glob_prior_mat[p,10] = 1/( (glob_prior_mat[p,5]*glob_prior_mat[p,6])^2 )

	glob_prior_mat[p,11] = 1/( (2*glob_prior_mat[p,8])^2 )
	glob_prior_mat[p,12] = ( 2*glob_prior_mat[p,8]/glob_prior_mat[p,7] )^2
}


##########################################################################
##########################################################################
##                                                                      ##
##  MCMC chains for population-level means                              ##
##                                                                      ##
##########################################################################
##########################################################################

par(ask=TRUE)

if( N_cohort >= 6 )
{
	par(mfrow=c(3,N_cohort))
}
if( N_cohort < 6 )
{
	par(mfrow=c(4,N_cohort))
}

for(i in 1:N_glob_par)
{
	plot(x=1:nrow(MCMC_pop), y=MCMC_pop[,i], 
           pch=19, cex=0.01, col="grey", 
	     ylim=quantile(MCMC_pop_burn[,i], prob=c(0.001,0.999), na.rm=TRUE),
	     main=paste( colnames(MCMC_pop)[i] ), ylab="", xlab="")
}



##########################################################################
##########################################################################
##                                                                      ##
##  MCMC chains for population-level standard deviations                ##
##                                                                      ##
##########################################################################
##########################################################################

if( N_cohort >= 6 )
{
	par(mfrow=c(3,N_cohort))
}
if( N_cohort < 6 )
{
	par(mfrow=c(4,N_cohort))
}

for(i in (N_glob_par+1):(2*N_glob_par))
{
	plot(x=1:nrow(MCMC_pop), y=MCMC_pop[,i], 
           pch=19, cex=0.01, col="grey", 
	     ylim=quantile(MCMC_pop_burn[,i], prob=c(0.001,0.999), na.rm=TRUE),
	     main=paste( colnames(MCMC_pop)[i] ), ylab="", xlab="")
}



##########################################################################
##########################################################################
##                                                                      ##
##  MCMC chains for observational variance and likelihood               ##
##                                                                      ##
##########################################################################
##########################################################################

par(mfrow=c(floor(N_cohort/2),3))

for(i in (2*N_glob_par+1):(2*N_glob_par+N_cohort))
{
	plot(x=1:nrow(MCMC_pop), y=MCMC_pop[,i], 
           pch=19, cex=0.01, col="grey", 
	     ylim=quantile(MCMC_pop_burn[,i], prob=c(0.001,0.999), na.rm=TRUE),
	     main=paste( colnames(MCMC_pop)[i] ), ylab="", xlab="")
}


plot(x=1:nrow(MCMC_pop), y=MCMC_pop[,N_par+1], 
pch=19, cex=0.01, col="grey", 
ylim=quantile(MCMC_pop[,N_par+1], prob=c(0.0001,1), na.rm=TRUE),
main=paste( colnames(MCMC_pop)[N_par+1] ), ylab="", xlab="")



##########################################################################
##########################################################################
##                                                                      ##
##  Posterior and priors for population-level means                     ##
##                                                                      ##
##########################################################################
##########################################################################

N_prior_sam <- 10000
N_prior_rho_sam <- 1000


if( N_cohort >= 6 )
{
	par(mfrow=c(3,N_cohort))
}
if( N_cohort < 6 )
{
	par(mfrow=c(4,N_cohort))
}

for(i in 1:N_glob_par)
{
	DEN_post = density( MCMC_pop_burn[,i], na.rm=TRUE )

	DEN_post$x <- DEN_post$x

	QUANT_post = quantile( MCMC_pop_burn[,i], prob=c(0.025, 0.5, 0.975), na.rm=TRUE )

	plot_upper <- 1.5*quantile( MCMC_pop_burn[,i], prob=0.99, na.rm=TRUE )

	prior_upper <- 100*quantile( MCMC_pop_burn[,i], prob=0.99, na.rm=TRUE )


	###############################
	## Posterior

	plot(x=DEN_post$x, y=DEN_post$y, type='l', #log="x",
	xlim=c(0, plot_upper),
	main=paste( colnames(MCMC_pop_burn)[i] ), ylab=NULL, xlab=NULL)

	low_index  = which(DEN_post$x<QUANT_post[1])
	mid_index  = intersect( which(DEN_post$x>=QUANT_post[1]), which(DEN_post$x<=QUANT_post[3]) )
	high_index = which(DEN_post$x>QUANT_post[3])

	polygon( x=c( DEN_post$x[low_index], rev(DEN_post$x[low_index]) ),
		   y=c( rep(0,length(low_index)), rev(DEN_post$y[low_index]) ), 
               col="pink")

	polygon( x=c( DEN_post$x[mid_index], rev(DEN_post$x[mid_index]) ),
		   y=c( rep(0,length(mid_index)), rev(DEN_post$y[mid_index]) ), 
               col="grey")

	polygon( x=c( DEN_post$x[high_index], rev(DEN_post$x[high_index]) ),
		   y=c( rep(0,length(high_index)), rev(DEN_post$y[high_index]) ), 
               col="pink")

	points(x=rep(QUANT_post[2],2), y=c(0,max(DEN_post$y)), type='l', lty="dashed", lwd=2)


	######################################################
	## Prior: Normal conjugate applied on the log scale

	if( i < N_glob_par)
	{
		prior_mu    <- glob_prior_mat[i,9]
		prior_tau   <- glob_prior_mat[i,10]
		prior_k     <- glob_prior_mat[i,11]
		prior_theta <- glob_prior_mat[i,12]
	
		prior_sig <- 1/sqrt(prior_tau)

		meanlog_sam <- rnorm( N_prior_sam, mean=prior_mu, sd=prior_sig)
		taulog_sam  <- rgamma(N_prior_sam, shape = prior_k, scale=prior_theta)

		siglog_sam <- 1/sqrt(taulog_sam)


		mean_sam <- exp( meanlog_sam + 0.5*siglog_sam^2 )

		mean_sam <- mean_sam[which(mean_sam > 0 & mean_sam < prior_upper)]


 		DEN_prior = density( mean_sam, na.rm=TRUE )

		points(x=DEN_prior$x, y=DEN_prior$y, type='l', col="red")
	}

	######################################################
	## Prior: Normal conjugate applied on the logit scale

	if( i == N_glob_par)
	{
		prior_mu    <- glob_prior_mat[i,9]
		prior_tau   <- glob_prior_mat[i,10]
		prior_k     <- glob_prior_mat[i,11]
		prior_theta <- glob_prior_mat[i,12]
	
		prior_sig <- 1/sqrt(prior_tau)

		meanlogit_sam <- rnorm( N_prior_rho_sam, mean=prior_mu, sd=prior_sig)
		taulogit_sam  <- rgamma(N_prior_rho_sam, shape = prior_k, scale=prior_theta)

		siglogit_sam <- 1/sqrt(taulogit_sam)

		rho_sam <- rep(NA, N_prior_rho_sam)
	
		for(k in 1:N_prior_rho_sam)
		{
			f_m <- function(x)
			{ 
				0.3989423*exp( -0.5*( ( log(x/(1-x))-meanlogit_sam[k] )/siglogit_sam[k] )^2 )/( siglogit_sam[k]*(1-x) )
			}
	
			tryCatch(
			{
				moment_1 <- integrate(f_m, lower=0, upper=1)$value
			}, error=function(e){ NULL }
			)

			rho_sam[k] <- moment_1
		}


 		DEN_prior = density( rho_sam, na.rm=TRUE )

		points(x=DEN_prior$x, y=DEN_prior$y, type='l', col="red")
	}
}


##########################################################################
##########################################################################
##                                                                      ##
##  Posterior and priors for population-level standard deviations       ##
##                                                                      ##
##########################################################################
##########################################################################

if( N_cohort >= 6 )
{
	par(mfrow=c(3,N_cohort))
}
if( N_cohort < 6 )
{
	par(mfrow=c(4,N_cohort))
}

for(i in (N_glob_par+1):(2*N_glob_par))
{
	DEN_1 = density( MCMC_pop_burn[,i], na.rm=TRUE )

	DEN_1$x <- DEN_1$x

	QUANT_1 = quantile( MCMC_pop_burn[,i], prob=c(0.025, 0.5, 0.975), na.rm=TRUE )

	plot_upper <- 1.5*quantile( MCMC_pop_burn[,i], prob=0.99, na.rm=TRUE )

	prior_upper <- 100*quantile( MCMC_pop_burn[,i], prob=0.99, na.rm=TRUE )


	###############################
	## Posterior

	plot(x=DEN_1$x, y=DEN_1$y, type='l', #log="x",
	xlim=c(0, plot_upper),
	main=paste( colnames(MCMC_pop_burn)[i] ), ylab=NULL, xlab=NULL)

	low_index  = which(DEN_1$x<QUANT_1[1])
	mid_index  = intersect( which(DEN_1$x>=QUANT_1[1]), which(DEN_1$x<=QUANT_1[3]) )
	high_index = which(DEN_1$x>QUANT_1[3])

	polygon( x=c( DEN_1$x[low_index], rev(DEN_1$x[low_index]) ),
		   y=c( rep(0,length(low_index)), rev(DEN_1$y[low_index]) ), 
               col="pink")

	polygon( x=c( DEN_1$x[mid_index], rev(DEN_1$x[mid_index]) ),
		   y=c( rep(0,length(mid_index)), rev(DEN_1$y[mid_index]) ), 
               col="grey")

	polygon( x=c( DEN_1$x[high_index], rev(DEN_1$x[high_index]) ),
		   y=c( rep(0,length(high_index)), rev(DEN_1$y[high_index]) ), 
               col="pink")

	points(x=rep(QUANT_1[2],2), y=c(0,max(DEN_1$y)), type='l', lty="dashed", lwd=2)


	######################################################
	## Prior: Gamma conjugate applied on precision

	if( (i > N_glob_par) & (i < 2*N_glob_par) )
	{
		prior_mu    <- glob_prior_mat[i-N_glob_par,9]
		prior_tau   <- glob_prior_mat[i-N_glob_par,10]
		prior_k     <- glob_prior_mat[i-N_glob_par,11]
		prior_theta <- glob_prior_mat[i-N_glob_par,12]
	
		prior_sig <- 1/sqrt(prior_tau)

		meanlog_sam <- rnorm( N_prior_sam, mean=prior_mu, sd=prior_sig)
		taulog_sam  <- rgamma(N_prior_sam, shape = prior_k, scale=prior_theta)

		siglog_sam <- 1/sqrt(taulog_sam)

		sd_sam <- exp( meanlog_sam + 0.5*siglog_sam^2 )*sqrt( exp(siglog_sam^2) - 1 )

		sd_sam <- sd_sam[which(sd_sam > 0 & sd_sam < prior_upper)]


 		DEN_prior = density( sd_sam, na.rm=TRUE )

		points(x=DEN_prior$x, y=DEN_prior$y, type='l', col="red")
	}

	######################################################
	## Prior: Normal conjugate applied on the logit scale

	if( i == 2*N_glob_par)
	{
		prior_mu    <- glob_prior_mat[i-N_glob_par,9]
		prior_tau   <- glob_prior_mat[i-N_glob_par,10]
		prior_k     <- glob_prior_mat[i-N_glob_par,11]
		prior_theta <- glob_prior_mat[i-N_glob_par,12]
	
		prior_sig <- 1/sqrt(prior_tau)

		meanlogit_sam <- rnorm( N_prior_rho_sam, mean=prior_mu, sd=prior_sig)
		taulogit_sam  <- rgamma(N_prior_rho_sam, shape = prior_k, scale=prior_theta)

		siglogit_sam <- 1/sqrt(taulogit_sam)

		rho_sd_sam <- rep(NA, N_prior_rho_sam)
	
		for(k in 1:N_prior_rho_sam)
		{
			f_m <- function(x)
			{ 
				0.3989423*exp( -0.5*( ( log(x/(1-x))-meanlogit_sam[k] )/siglogit_sam[k] )^2 )/( siglogit_sam[k]*(1-x) )
			}

			f_m2 <- function(x)
			{ 
				x*0.3989423*exp( -0.5*( (log(x/(1-x))-meanlogit_sam[k])/siglogit_sam[k] )^2 )/( siglogit_sam[k]*(1-x) )
			}
	
			tryCatch(
			{
				moment_1 <- integrate(f_m, lower=0, upper=1)$value
				moment_2 <- integrate(f_m2, lower=0, upper=1)$value
			}, error=function(e){ NULL }
			)

			rho_sd_sam[k] <- sqrt( moment_2 - moment_1^2 )
		}


 		DEN_prior = density( rho_sd_sam, na.rm=TRUE )

		points(x=DEN_prior$x, y=DEN_prior$y, type='l', col="red")
	}
}




##########################################################################
##########################################################################
##                                                                      ##
##  Posterior and priors for observational variance                     ##
##                                                                      ##
##########################################################################
##########################################################################

N_prior_sam <- 10000



par(mfrow=c(2,2))

for(i in (2*N_glob_par+1):(2*N_glob_par+4))
{
	DEN_post = density( MCMC_pop_burn[,i], na.rm=TRUE )

	DEN_post$x <- DEN_post$x

	QUANT_post = quantile( MCMC_pop_burn[,i], prob=c(0.025, 0.5, 0.975), na.rm=TRUE )


	###############################
	## Posterior

	plot(x=DEN_post$x, y=DEN_post$y, type='l', #log="x",
	main=paste( colnames(MCMC_pop_burn)[i] ), ylab=NULL, xlab=NULL)

	low_index  = which(DEN_post$x<QUANT_post[1])
	mid_index  = intersect( which(DEN_post$x>=QUANT_post[1]), which(DEN_post$x<=QUANT_post[3]) )
	high_index = which(DEN_post$x>QUANT_post[3])

	polygon( x=c( DEN_post$x[low_index], rev(DEN_post$x[low_index]) ),
		   y=c( rep(0,length(low_index)), rev(DEN_post$y[low_index]) ), 
               col="pink")

	polygon( x=c( DEN_post$x[mid_index], rev(DEN_post$x[mid_index]) ),
		   y=c( rep(0,length(mid_index)), rev(DEN_post$y[mid_index]) ), 
               col="grey")

	polygon( x=c( DEN_post$x[high_index], rev(DEN_post$x[high_index]) ),
		   y=c( rep(0,length(high_index)), rev(DEN_post$y[high_index]) ), 
               col="pink")

	points(x=rep(QUANT_post[2],2), y=c(0,max(DEN_post$y)), type='l', lty="dashed", lwd=2)
}


##########################################################################
##########################################################################
##                                                                      ##
##  Output means and standard deviations                                ##
##                                                                      ##
##########################################################################
##########################################################################

for(i in 1:N_par)
{
	print( colnames(MCMC_pop_burn)[i] )

	print(paste( quantile( MCMC_pop_burn[,i], prob=c(0.5, 0.025, 0.975), na.rm=TRUE ) ))
}


#######################################################################
#######################################################################
##                                                                   ##
##  #### #   ## ####   #### ##   ## #### ####   ##  ##  ####  ##     ##
##   ##  ##  ## ## ##   ##  ##   ##  ##  ## ##  ##  ## ##  ## ##     ##
##   ##  ### ## ##  ##  ##   ## ##   ##  ##  ## ##  ## ###### ##     ##
##   ##  ## ### ## ##   ##    ###    ##  ## ##  ##  ## ##  ## ##     ## 
##  #### ##  ## ####   ####    #    #### ####    ####  ##  ## #####  ##
##                                                                   ## 
#######################################################################
#######################################################################

for( c in 1:N_cohort )
{
	AB_data = read.table( ind_pos_file[c] )

	AB_max = max(AB_data[,(N_tt[c]+1):(2*N_tt[c])])
	tt_max = max(AB_data[,1:N_tt[c]])


	AB_min <- AB_data[,(N_tt[c]+1):(2*N_tt[c])]
	AB_min[which(AB_min < -50, arr.ind=TRUE)] <- NA
	AB_min = min( AB_min, na.rm=TRUE)

	tt_min <- AB_data[,1:N_tt[c]]
	tt_min[which(tt_min < -50, arr.ind=TRUE)] <- NA
	tt_min = min( tt_min, na.rm=TRUE)

	tt_min <- tt_min - 10

	tt_upper <- 1.2*tt_max

	tt_plot <- seq(from=tt_min, to=tt_upper, by=1)
	log2 <- log(2)


	plot_order <- rowSums(AB_data[,1:(0.5*ncol(AB_data))] > -90)
	plot_order <- order(plot_order, decreasing=TRUE)

	MCMC_ind <- read.table( ind_out_file[c] )

	MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]


	par(ask=TRUE)

	lay.mat <- rbind( c( 1, 1, 2, 2, 3, 4, 5, 6 ),
	                  c( 1, 1, 2, 2, 7, 8, 9,10 ),
		            c(11,11,12,12,13,14,15,16 ),
	                  c(11,11,12,12,17,18,19,20 ) )
	layout(lay.mat)
	layout.show(20)


	##for(n in 1:N_part[c] )
	for(n in plot_order )
	{
		###################################
		## Data for participant n

		N_sam = length(which(AB_data[n,1:N_tt[c]] > -50))

		tt = as.numeric(AB_data[n,1:N_sam])
		AB = as.numeric(AB_data[n,(1+N_tt[c]):(N_tt[c]+N_sam)])

	
		###################################
		## Model prediction for participant i
		## Posterior projections

		AB_mod <- matrix(NA, nrow=nrow(MCMC_ind), ncol=length(tt_plot))
	
		for(i in 1:nrow(AB_mod))
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

			##AB_tt = Ab_0*exp(-r_cl*tt_plot)
			AB_tt =  rep( Ab_0, length(tt_plot) )

			tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
			AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*(         (rho/(r_a - r_cs))*( exp(-r_cs*tt_temp) - exp(-r_a*tt_temp) ) +
			                                       ((1.0 - rho)/(r_a - r_cl))*( exp(-r_cl*tt_temp) - exp(-r_a*tt_temp) ) )

			AB_mod[i,] = AB_tt
		}

		AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
		for(j in 1:length(tt_plot))
		{
			AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
		}


		########################
		## PANEL 1

		plot_title <- paste(cohort_names[c], " case #", n, " / ", N_part[c], sep="")

		plot( x=tt, y=AB, pch=19, cex=1, col="red", 
	            xlim=c(tt_min,tt_upper), ylim=c(0,2*max(AB)),
			xlab="time (days)", ylab="antibody titre",
	            main=plot_title )

		points(x=tt_plot, y=AB_quant[2,], type='l')

		polygon(x=c(tt_plot, rev(tt_plot)), 
			  y=c( AB_quant[1,], rev(AB_quant[3,]) ),
		        col=rgb(190/256,190/256,190/256,0.4), border=NA)

		##points(x=c(-100,10000), y=c(0.02,0.02), type='l', lty="dashed")


		########################
		## PANEL 2

		plot( x=tt, y=AB, pch=19, cex=1, col="red", 
      	      xlim=c(tt_min,tt_upper), ylim=c(0.5*AB_min,2*AB_max), log="y",
			xlab="time (days)", ylab="antibody titre",
      	      main=plot_title )


		points(x=tt_plot, y=AB_quant[2,], type='l')


		polygon(x=c(tt_plot, rev(tt_plot)), 
			y=c( AB_quant[1,], rev(AB_quant[3,]) ),
			col=rgb(190/256,190/256,190/256,0.4), border=NA)

		##points(x=c(-100,1000), y=c(0.02,0.02), type='l', lty="dashed")
		points(x=c(-100,10000), y=c(AB_min,AB_min), type='l', lty="dashed")



		#################################
		## PANEL 3: Ab_0 chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+1], 
      	     log="y", pch=19, cex=0.01, col="grey", 
      	     xlab="MCMC iteration", ylab="", main="Ab_0" )


		#################################
		## PANEL 4: beta chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+2], 
      	     log="y", pch=19, cex=0.01, col="grey", 
      	     xlab="MCMC iteration", ylab="", main="beta" )


		#################################
		## PANEL 5: tau chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+3], 
      	     pch=19, cex=0.01, col="grey", 
      	     xlab="MCMC iteration", ylab="", main="tau" )


		#################################
		## PANEL 6: t_short chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+4], 
      	     pch=19, cex=0.01, col="grey", 
      	     xlab="MCMC iteration", ylab="", main="t_short" )


		#################################
		## PANEL 7: t_long chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+5], 
      	     pch=19, cex=0.01, col="grey", log="y",
      	     xlab="MCMC iteration", ylab="", main="t_long" )


		#################################
		## PANEL 8: t_IgG chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+6], 
      	     pch=19, cex=0.01, col="grey", 
      	     xlab="MCMC iteration", ylab="", main="t_IgG" )


		#################################
		## PANEL 9: rho chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+7], 
      	     pch=19, cex=0.01, col="grey", ylim=c(0,1),
      	     xlab="MCMC iteration", ylab="", main="rho" )


		#################################
		## PANEL 10: lhood chain

		plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*8+8], 
      	     pch=19, cex=0.01, col="grey", 
      	     xlab="MCMC iteration", ylab="", main="likelihood" )

	}
}











