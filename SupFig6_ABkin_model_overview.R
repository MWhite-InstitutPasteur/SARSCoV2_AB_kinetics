

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

###################
###             ###
###  Spike IgG  ###
###             ###
###################

Spike_IgG_cols <- cohort_cols[c(1,2,4,5,6,7,8)]
Spike_IgG_rgb  <- cohort_rgb[c(1,2,4,5,6,7,8)]

load("C:\\U\\CoronaVirus\\PELLEAU_et_al\\5_kinetics\\3_processed_output\\Spike_IgG.RData")

N_cohort_Spike_IgG <- length(Spike_IgG_N_part)

Spike_IgG_norm <- rep(NA, N_cohort_Spike_IgG)

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_norm[c] <- exp(mean(log( Spike_IgG_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))

	##Spike_IgG_norm[c] <- median( Spike_IgG_mod[[c]][,which(tt_plot==14),2], na.rm=TRUE)
}


Spike_IgG_mod_median <- matrix(NA, nrow=N_cohort_Spike_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_mod_median[c,] <- apply(X=Spike_IgG_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/Spike_IgG_norm[c]
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

Spike_IgM_norm <- rep(NA, N_cohort_Spike_IgM)

for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_norm[c] <- exp(mean(log( Spike_IgM_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}



Spike_IgM_mod_median <- matrix(NA, nrow=N_cohort_Spike_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_mod_median[c,] <- apply(X=Spike_IgM_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/Spike_IgM_norm[c]
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

Spike_IgA_norm <- rep(NA, N_cohort_Spike_IgA)

for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_norm[c] <- exp(mean(log( Spike_IgA_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}



Spike_IgA_mod_median <- matrix(NA, nrow=N_cohort_Spike_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_mod_median[c,] <- apply(X=Spike_IgA_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/Spike_IgA_norm[c]
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

RBD_IgG_norm <- rep(NA, N_cohort_RBD_IgG)

for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_norm[c] <- exp(mean(log( RBD_IgG_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}



RBD_IgG_mod_median <- matrix(NA, nrow=N_cohort_RBD_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_mod_median[c,] <- apply(X=RBD_IgG_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/RBD_IgG_norm[c]
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

RBD_IgM_norm <- rep(NA, N_cohort_RBD_IgM)

for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_norm[c] <- exp(mean(log( RBD_IgM_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}



RBD_IgM_mod_median <- matrix(NA, nrow=N_cohort_RBD_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_mod_median[c,] <- apply(X=RBD_IgM_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/RBD_IgM_norm[c]
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

RBD_IgA_norm <- rep(NA, N_cohort_RBD_IgA)

for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_norm[c] <- exp(mean(log( RBD_IgA_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}


RBD_IgA_mod_median <- matrix(NA, nrow=N_cohort_RBD_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_mod_median[c,] <- apply(X=RBD_IgA_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/RBD_IgA_norm[c]
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

NP_IgG_norm <- rep(NA, N_cohort_NP_IgG)

for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_norm[c] <- exp(mean(log( NP_IgG_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}



NP_IgG_mod_median <- matrix(NA, nrow=N_cohort_NP_IgG, ncol=length(tt_plot))

for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_mod_median[c,] <- apply(X=NP_IgG_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/NP_IgG_norm[c]
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

NP_IgM_norm <- rep(NA, N_cohort_NP_IgM)

for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_norm[c] <- exp(mean(log( NP_IgM_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}


NP_IgM_mod_median <- matrix(NA, nrow=N_cohort_NP_IgM, ncol=length(tt_plot))

for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_mod_median[c,] <- apply(X=NP_IgM_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/NP_IgM_norm[c]
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

NP_IgA_norm <- rep(NA, N_cohort_NP_IgA)

for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_norm[c] <- exp(mean(log( NP_IgA_mod[[c]][,which(tt_plot==14),2] ), na.rm=TRUE))
}

NP_IgA_mod_median <- matrix(NA, nrow=N_cohort_NP_IgA, ncol=length(tt_plot))

for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_mod_median[c,] <- apply(X=NP_IgA_mod[[c]][,,2], FUN=median, MARGIN=2, na.rm=TRUE)/NP_IgA_norm[c]
}







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

line_seq_x <- c(-1, 0, 2, 3, 4, 5, 6, 7, 8)*(365/4)
line_seq_y <- c(0.001, 0.01, 0.1, 1, 10, 100)

par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))

pos.point.size = 0.3
neg.point.size = 0.3
lab.size   = 1.5
axis.size  = 1.25	
main.size  = 1.5
line.size  = 2
dash.line.size = 0.5

ww <- 6




tiff(file="SupFig6_ABkin_model_overview.tif", width=40, height=32, units="cm", res=500)


lay.mat <- rbind( c( 1, 2, 3 ), 
                  c( 4, 5, 6 ),
                  c( 7, 8, 9 ), 
                  c(10,11,12 ), 
                  c(13,13,13 ) )
layout(lay.mat, heights=c(1.5,10,10,11,1), widths=c(11,10,10))
layout.show(13)


############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))

plot.new()
title( "Spike", 
        cex.main=4.0, line=-3)

plot.new()
title( "RBD", 
        cex.main=4.0, line=-3)

plot.new()
title( "Nucleoprotein", 
        cex.main=4.0, line=-3)


###################
###################
###             ###
###  PANEL 1    ###
###  Spike IgG  ###
###             ###
###################
###################

par(mar=c(3,8,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 6, 
cex=lab.size, 
text="IgG antibody level")

mtext(side = 2, line = 4, 
cex=lab.size, 
text="(relative to day 14)")

###################
## Negative controls


for(c in 1:N_cohort_Spike_IgG)
{
	if( Spike_IgG_N_neg[c] > 0 )
	{
		DEN_Spike_IgG_neg = density( log(Spike_IgG_neg_data[[c]]/Spike_IgG_norm[c]), cut=0, na.rm=TRUE )

		DEN_Spike_IgG_neg$x <- exp(DEN_Spike_IgG_neg$x)
		DEN_Spike_IgG_neg$y <- DEN_Spike_IgG_neg$y/max(DEN_Spike_IgG_neg$y)

		polygon( x = -90 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) +
                         ww*c( -DEN_Spike_IgG_neg$y, rev(DEN_Spike_IgG_neg$y) ),
			   y=c( DEN_Spike_IgG_neg$x, rev(DEN_Spike_IgG_neg$x) ), 
      	   	   col=Spike_IgG_rgb[c], border=NA)

		for(n in 1:Spike_IgG_N_neg[c])
		{
			index <- which.min(abs(Spike_IgG_neg_data[[c]][n]/Spike_IgG_norm[c] - DEN_Spike_IgG_neg$x))

			points( x = -90 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_Spike_IgG_neg$y[index], max=DEN_Spike_IgG_neg$y[index]),
                          y = Spike_IgG_neg_data[[c]][n]/Spike_IgG_norm[c],
                          pch=19, cex=neg.point.size, col=Spike_IgG_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_1m <- Spike_IgG_mod[[c]][,which(tt_plot==14),2]/Spike_IgG_norm[c]	

	DEN_Spike_IgG_1m = density( log(Spike_IgG_1m), cut=0, na.rm=TRUE )

	DEN_Spike_IgG_1m$x <- exp(DEN_Spike_IgG_1m$x)
	DEN_Spike_IgG_1m$y <- DEN_Spike_IgG_1m$y/max(DEN_Spike_IgG_1m$y)

	polygon( x = 30 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_Spike_IgG_1m$y, rev(DEN_Spike_IgG_1m$y) ),
		   y=c( DEN_Spike_IgG_1m$x, rev(DEN_Spike_IgG_1m$x) ), 
         	   col=Spike_IgG_rgb[c], border=NA)

	for(n in 1:length(Spike_IgG_1m))
	{
		index <- which.min(abs(Spike_IgG_1m[n] - DEN_Spike_IgG_1m$x))

		points( x = 30 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgG_1m$y[index], max=DEN_Spike_IgG_1m$y[index]),
                    y = Spike_IgG_1m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgG_cols[c])
	}

}



###################
## Positive points
## Month 6

off_set <- -2*ww

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_6m <- Spike_IgG_mod[[c]][,which(tt_plot==180),2]/Spike_IgG_norm[c]	

	DEN_Spike_IgG_6m = density( log(Spike_IgG_6m), cut=0, na.rm=TRUE )

	DEN_Spike_IgG_6m$x <- exp(DEN_Spike_IgG_6m$x)
	DEN_Spike_IgG_6m$y <- DEN_Spike_IgG_6m$y/max(DEN_Spike_IgG_6m$y)

	polygon( x = 180 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) +
                   ww*c( -DEN_Spike_IgG_6m$y, rev(DEN_Spike_IgG_6m$y) ),
		   y=c( DEN_Spike_IgG_6m$x, rev(DEN_Spike_IgG_6m$x) ), 
         	   col=Spike_IgG_rgb[c], border=NA)

	for(n in 1:length(Spike_IgG_6m))
	{
		index <- which.min(abs(Spike_IgG_6m[n] - DEN_Spike_IgG_6m$x))

		points( x = 180 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgG_6m$y[index], max=DEN_Spike_IgG_6m$y[index]),
                    y = Spike_IgG_6m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgG_cols[c])
	}

	off_set <- off_set - 2*ww
}




###################
## Positive points
## Month 12

off_set <- -2*ww

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_12m <- Spike_IgG_mod[[c]][,which(tt_plot==365),2]/Spike_IgG_norm[c]	

	DEN_Spike_IgG_12m = density( log(Spike_IgG_12m), cut=0, na.rm=TRUE )

	DEN_Spike_IgG_12m$x <- exp(DEN_Spike_IgG_12m$x)
	DEN_Spike_IgG_12m$y <- DEN_Spike_IgG_12m$y/max(DEN_Spike_IgG_12m$y)

	polygon( x = 365 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_Spike_IgG_12m$y, rev(DEN_Spike_IgG_12m$y) ),
		   y=c( DEN_Spike_IgG_12m$x, rev(DEN_Spike_IgG_12m$x) ), 
         	   col=Spike_IgG_rgb[c], border=NA)

	for(n in 1:length(Spike_IgG_12m))
	{
		index <- which.min(abs(Spike_IgG_12m[n] - DEN_Spike_IgG_12m$x))

		points( x = 365 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_Spike_IgG_12m$y[index], max=DEN_Spike_IgG_12m$y[index]),
                    y = Spike_IgG_12m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgG_cols[c])
	}

	off_set <- off_set - 2*ww
}



###################
## Positive points
## Month 24

off_set <- -2*ww

for(c in 1:N_cohort_Spike_IgG)
{
	Spike_IgG_24m <- Spike_IgG_mod[[c]][,which(tt_plot==2*365),2]/Spike_IgG_norm[c]	

	DEN_Spike_IgG_24m = density( log(Spike_IgG_24m), cut=0, na.rm=TRUE )

	DEN_Spike_IgG_24m$x <- exp(DEN_Spike_IgG_24m$x)
	DEN_Spike_IgG_24m$y <- DEN_Spike_IgG_24m$y/max(DEN_Spike_IgG_24m$y)

	polygon( x = 2*365 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) +
                   ww*c( -DEN_Spike_IgG_24m$y, rev(DEN_Spike_IgG_24m$y) ),
		   y=c( DEN_Spike_IgG_24m$x, rev(DEN_Spike_IgG_24m$x) ), 
         	   col=Spike_IgG_rgb[c], border=NA)

	for(n in 1:length(Spike_IgG_24m))
	{
		index <- which.min(abs(Spike_IgG_24m[n] - DEN_Spike_IgG_24m$x))

		points( x = 2*365 - N_cohort_Spike_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgG_24m$y[index], max=DEN_Spike_IgG_24m$y[index]),
                    y = Spike_IgG_24m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgG_cols[c])
	}

	off_set <- off_set - 2*ww
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
###  PANEL 2    ###
###  RBD IgG    ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


###################
## Negative controls


for(c in 1:N_cohort_RBD_IgG)
{
	if( RBD_IgG_N_neg[c] > 0 )
	{
		DEN_RBD_IgG_neg = density( log(RBD_IgG_neg_data[[c]]/RBD_IgG_norm[c]), cut=0, na.rm=TRUE )

		DEN_RBD_IgG_neg$x <- exp(DEN_RBD_IgG_neg$x)
		DEN_RBD_IgG_neg$y <- DEN_RBD_IgG_neg$y/max(DEN_RBD_IgG_neg$y)

		polygon( x = -90 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) +
                         ww*c( -DEN_RBD_IgG_neg$y, rev(DEN_RBD_IgG_neg$y) ),
			   y=c( DEN_RBD_IgG_neg$x, rev(DEN_RBD_IgG_neg$x) ), 
      	   	   col=RBD_IgG_rgb[c], border=NA)

		for(n in 1:RBD_IgG_N_neg[c])
		{
			index <- which.min(abs(RBD_IgG_neg_data[[c]][n]/RBD_IgG_norm[c] - DEN_RBD_IgG_neg$x))

			points( x = -90 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_RBD_IgG_neg$y[index], max=DEN_RBD_IgG_neg$y[index]),
                          y = RBD_IgG_neg_data[[c]][n]/RBD_IgG_norm[c],
                          pch=19, cex=neg.point.size, col=RBD_IgG_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_1m <- RBD_IgG_mod[[c]][,which(tt_plot==14),2]/RBD_IgG_norm[c]	

	DEN_RBD_IgG_1m = density( log(RBD_IgG_1m), cut=0, na.rm=TRUE )

	DEN_RBD_IgG_1m$x <- exp(DEN_RBD_IgG_1m$x)
	DEN_RBD_IgG_1m$y <- DEN_RBD_IgG_1m$y/max(DEN_RBD_IgG_1m$y)

	polygon( x = 30 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_RBD_IgG_1m$y, rev(DEN_RBD_IgG_1m$y) ),
		   y=c( DEN_RBD_IgG_1m$x, rev(DEN_RBD_IgG_1m$x) ), 
         	   col=RBD_IgG_rgb[c], border=NA)

	for(n in 1:length(RBD_IgG_1m))
	{
		index <- which.min(abs(RBD_IgG_1m[n] - DEN_RBD_IgG_1m$x))

		points( x = 30 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgG_1m$y[index], max=DEN_RBD_IgG_1m$y[index]),
                    y = RBD_IgG_1m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgG_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_6m <- RBD_IgG_mod[[c]][,which(tt_plot==180),2]/RBD_IgG_norm[c]	

	DEN_RBD_IgG_6m = density( log(RBD_IgG_6m), cut=0, na.rm=TRUE )

	DEN_RBD_IgG_6m$x <- exp(DEN_RBD_IgG_6m$x)
	DEN_RBD_IgG_6m$y <- DEN_RBD_IgG_6m$y/max(DEN_RBD_IgG_6m$y)

	polygon( x = 180 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) +
                   ww*c( -DEN_RBD_IgG_6m$y, rev(DEN_RBD_IgG_6m$y) ),
		   y=c( DEN_RBD_IgG_6m$x, rev(DEN_RBD_IgG_6m$x) ), 
         	   col=RBD_IgG_rgb[c], border=NA)

	for(n in 1:length(RBD_IgG_6m))
	{
		index <- which.min(abs(RBD_IgG_6m[n] - DEN_RBD_IgG_6m$x))

		points( x = 180 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgG_6m$y[index], max=DEN_RBD_IgG_6m$y[index]),
                    y = RBD_IgG_6m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgG_cols[c])
	}

}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_12m <- RBD_IgG_mod[[c]][,which(tt_plot==365),2]/RBD_IgG_norm[c]	

	DEN_RBD_IgG_12m = density( log(RBD_IgG_12m), cut=0, na.rm=TRUE )

	DEN_RBD_IgG_12m$x <- exp(DEN_RBD_IgG_12m$x)
	DEN_RBD_IgG_12m$y <- DEN_RBD_IgG_12m$y/max(DEN_RBD_IgG_12m$y)

	polygon( x = 365 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_RBD_IgG_12m$y, rev(DEN_RBD_IgG_12m$y) ),
		   y=c( DEN_RBD_IgG_12m$x, rev(DEN_RBD_IgG_12m$x) ), 
         	   col=Spike_IgG_rgb[c], border=NA)

	for(n in 1:length(RBD_IgG_12m))
	{
		index <- which.min(abs(RBD_IgG_12m[n] - DEN_RBD_IgG_12m$x))

		points( x = 365 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_RBD_IgG_12m$y[index], max=DEN_RBD_IgG_12m$y[index]),
                    y = RBD_IgG_12m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgG_cols[c])
	}

}



###################
## Positive points
## Month 24

off_set <- -2*ww

for(c in 1:N_cohort_RBD_IgG)
{
	RBD_IgG_24m <- RBD_IgG_mod[[c]][,which(tt_plot==2*365),2]/RBD_IgG_norm[c]	

	DEN_RBD_IgG_24m = density( log(RBD_IgG_24m), cut=0, na.rm=TRUE )

	DEN_RBD_IgG_24m$x <- exp(DEN_RBD_IgG_24m$x)
	DEN_RBD_IgG_24m$y <- DEN_RBD_IgG_24m$y/max(DEN_RBD_IgG_24m$y)

	polygon( x = 2*365 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) +
                   ww*c( -DEN_RBD_IgG_24m$y, rev(DEN_RBD_IgG_24m$y) ),
		   y=c( DEN_RBD_IgG_24m$x, rev(DEN_RBD_IgG_24m$x) ), 
         	   col=RBD_IgG_rgb[c], border=NA)

	for(n in 1:length(RBD_IgG_24m))
	{
		index <- which.min(abs(RBD_IgG_24m[n] - DEN_RBD_IgG_24m$x))

		points( x = 2*365 - N_cohort_RBD_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgG_24m$y[index], max=DEN_RBD_IgG_24m$y[index]),
                    y = RBD_IgG_24m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgG_cols[c])
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
###  PANEL 3    ###
###  NP IgG     ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


###################
## Negative controls


for(c in 1:N_cohort_NP_IgG)
{
	if( NP_IgG_N_neg[c] > 0 )
	{
		DEN_NP_IgG_neg = density( log(NP_IgG_neg_data[[c]]/NP_IgG_norm[c]), cut=0, na.rm=TRUE )

		DEN_NP_IgG_neg$x <- exp(DEN_NP_IgG_neg$x)
		DEN_NP_IgG_neg$y <- DEN_NP_IgG_neg$y/max(DEN_NP_IgG_neg$y)

		polygon( x = -90 - N_cohort_NP_IgG*ww + 2*ww*(c-1) +
                         ww*c( -DEN_NP_IgG_neg$y, rev(DEN_NP_IgG_neg$y) ),
			   y=c( DEN_NP_IgG_neg$x, rev(DEN_NP_IgG_neg$x) ), 
      	   	   col=NP_IgG_rgb[c], border=NA)

		for(n in 1:NP_IgG_N_neg[c])
		{
			index <- which.min(abs(NP_IgG_neg_data[[c]][n]/NP_IgG_norm[c] - DEN_NP_IgG_neg$x))

			points( x = -90 - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_NP_IgG_neg$y[index], max=DEN_NP_IgG_neg$y[index]),
                          y = NP_IgG_neg_data[[c]][n]/NP_IgG_norm[c],
                          pch=19, cex=neg.point.size, col=NP_IgG_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_1m <- NP_IgG_mod[[c]][,which(tt_plot==14),2]/NP_IgG_norm[c]	

	DEN_NP_IgG_1m = density( log(NP_IgG_1m), cut=0, na.rm=TRUE )

	DEN_NP_IgG_1m$x <- exp(DEN_NP_IgG_1m$x)
	DEN_NP_IgG_1m$y <- DEN_NP_IgG_1m$y/max(DEN_NP_IgG_1m$y)

	polygon( x = 30 - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_NP_IgG_1m$y, rev(DEN_NP_IgG_1m$y) ),
		   y=c( DEN_NP_IgG_1m$x, rev(DEN_NP_IgG_1m$x) ), 
         	   col=NP_IgG_rgb[c], border=NA)

	for(n in 1:length(NP_IgG_1m))
	{
		index <- which.min(abs(NP_IgG_1m[n] - DEN_NP_IgG_1m$x))

		points( x = 30 - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgG_1m$y[index], max=DEN_NP_IgG_1m$y[index]),
                    y = NP_IgG_1m[n],
                    pch=19, cex=pos.point.size, col=NP_IgG_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_6m <- NP_IgG_mod[[c]][,which(tt_plot==180),2]/NP_IgG_norm[c]	

	DEN_NP_IgG_6m = density( log(NP_IgG_6m), cut=0, na.rm=TRUE )

	DEN_NP_IgG_6m$x <- exp(DEN_NP_IgG_6m$x)
	DEN_NP_IgG_6m$y <- DEN_NP_IgG_6m$y/max(DEN_NP_IgG_6m$y)

	polygon( x = 180 - N_cohort_NP_IgG*ww + 2*ww*(c-1) +
                   ww*c( -DEN_NP_IgG_6m$y, rev(DEN_NP_IgG_6m$y) ),
		   y=c( DEN_NP_IgG_6m$x, rev(DEN_NP_IgG_6m$x) ), 
         	   col=NP_IgG_rgb[c], border=NA)

	for(n in 1:length(NP_IgG_6m))
	{
		index <- which.min(abs(NP_IgG_6m[n] - DEN_NP_IgG_6m$x))

		points( x = 180 - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgG_6m$y[index], max=DEN_NP_IgG_6m$y[index]),
                    y = NP_IgG_6m[n],
                    pch=19, cex=pos.point.size, col=NP_IgG_cols[c])
	}

}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_12m <- NP_IgG_mod[[c]][,which(tt_plot==365),2]/NP_IgG_norm[c]	

	DEN_NP_IgG_12m = density( log(NP_IgG_12m), cut=0, na.rm=TRUE )

	DEN_NP_IgG_12m$x <- exp(DEN_NP_IgG_12m$x)
	DEN_NP_IgG_12m$y <- DEN_NP_IgG_12m$y/max(DEN_NP_IgG_12m$y)

	polygon( x = 365 - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_NP_IgG_12m$y, rev(DEN_NP_IgG_12m$y) ),
		   y=c( DEN_NP_IgG_12m$x, rev(DEN_NP_IgG_12m$x) ), 
         	   col=NP_IgG_rgb[c], border=NA)

	for(n in 1:length(NP_IgG_12m))
	{
		index <- which.min(abs(NP_IgG_12m[n] - DEN_NP_IgG_12m$x))

		points( x = 365 - N_cohort_NP_IgG*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_NP_IgG_12m$y[index], max=DEN_NP_IgG_12m$y[index]),
                    y = NP_IgG_12m[n],
                    pch=19, cex=pos.point.size, col=NP_IgG_cols[c])
	}

}



###################
## Positive points
## Month 24

off_set <- -2*ww

for(c in 1:N_cohort_NP_IgG)
{
	NP_IgG_24m <- NP_IgG_mod[[c]][,which(tt_plot==2*365),2]/NP_IgG_norm[c]	

	DEN_NP_IgG_24m = density( log(NP_IgG_24m), cut=0, na.rm=TRUE )

	DEN_NP_IgG_24m$x <- exp(DEN_NP_IgG_24m$x)
	DEN_NP_IgG_24m$y <- DEN_NP_IgG_24m$y/max(DEN_NP_IgG_24m$y)

	polygon( x = 2*365 - N_cohort_NP_IgG*ww + 2*ww*(c-1) +
                   ww*c( -DEN_NP_IgG_24m$y, rev(DEN_NP_IgG_24m$y) ),
		   y=c( DEN_NP_IgG_24m$x, rev(DEN_NP_IgG_24m$x) ), 
         	   col=NP_IgG_rgb[c], border=NA)

	for(n in 1:length(NP_IgG_24m))
	{
		index <- which.min(abs(NP_IgG_24m[n] - DEN_NP_IgG_24m$x))

		points( x = 2*365 - N_cohort_NP_IgG*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgG_24m$y[index], max=DEN_NP_IgG_24m$y[index]),
                    y = NP_IgG_24m[n],
                    pch=19, cex=pos.point.size, col=NP_IgG_cols[c])
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
###  PANEL 4    ###
###  Spike IgM  ###
###             ###
###################
###################

par(mar=c(3,8,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )

mtext(side = 2, line = 6, 
cex=lab.size, 
text="IgM antibody level")

mtext(side = 2, line = 4, 
cex=lab.size, 
text="(relative to day 14)")


###################
## Negative controls

for(c in 1:N_cohort_Spike_IgM)
{
	if( Spike_IgM_N_neg[c] > 0 )
	{
		DEN_Spike_IgM_neg = density( log(Spike_IgM_neg_data[[c]]/Spike_IgM_norm[c]), cut=0, na.rm=TRUE )

		DEN_Spike_IgM_neg$x <- exp(DEN_Spike_IgM_neg$x)
		DEN_Spike_IgM_neg$y <- DEN_Spike_IgM_neg$y/max(DEN_Spike_IgM_neg$y)

		polygon( x = -90 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) +
                         ww*c( -DEN_Spike_IgM_neg$y, rev(DEN_Spike_IgM_neg$y) ),
			   y=c( DEN_Spike_IgM_neg$x, rev(DEN_Spike_IgM_neg$x) ), 
      	   	   col=Spike_IgM_rgb[c], border=NA)

		for(n in 1:Spike_IgM_N_neg[c])
		{
			index <- which.min(abs(Spike_IgM_neg_data[[c]][n]/Spike_IgM_norm[c] - DEN_Spike_IgM_neg$x))

			points( x = -90 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_Spike_IgM_neg$y[index], max=DEN_Spike_IgM_neg$y[index]),
                          y = Spike_IgM_neg_data[[c]][n]/Spike_IgM_norm[c],
                          pch=19, cex=neg.point.size, col=Spike_IgM_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_1m <- Spike_IgM_mod[[c]][,which(tt_plot==14),2]/Spike_IgM_norm[c]	

	DEN_Spike_IgM_1m = density( log(Spike_IgM_1m), cut=0, na.rm=TRUE )

	DEN_Spike_IgM_1m$x <- exp(DEN_Spike_IgM_1m$x)
	DEN_Spike_IgM_1m$y <- DEN_Spike_IgM_1m$y/max(DEN_Spike_IgM_1m$y)

	polygon( x = 30 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_Spike_IgM_1m$y, rev(DEN_Spike_IgM_1m$y) ),
		   y=c( DEN_Spike_IgM_1m$x, rev(DEN_Spike_IgM_1m$x) ), 
         	   col=Spike_IgM_rgb[c], border=NA)

	for(n in 1:length(Spike_IgM_1m))
	{
		index <- which.min(abs(Spike_IgM_1m[n] - DEN_Spike_IgM_1m$x))

		points( x = 30 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgM_1m$y[index], max=DEN_Spike_IgM_1m$y[index]),
                    y = Spike_IgM_1m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgM_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_6m <- Spike_IgM_mod[[c]][,which(tt_plot==180),2]/Spike_IgM_norm[c]	

	DEN_Spike_IgM_6m = density( log(Spike_IgM_6m), cut=0, na.rm=TRUE )

	DEN_Spike_IgM_6m$x <- exp(DEN_Spike_IgM_6m$x)
	DEN_Spike_IgM_6m$y <- DEN_Spike_IgM_6m$y/max(DEN_Spike_IgM_6m$y)

	polygon( x = 180 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) +
                   ww*c( -DEN_Spike_IgM_6m$y, rev(DEN_Spike_IgM_6m$y) ),
		   y=c( DEN_Spike_IgM_6m$x, rev(DEN_Spike_IgM_6m$x) ), 
         	   col=Spike_IgM_rgb[c], border=NA)

	for(n in 1:length(Spike_IgM_6m))
	{
		index <- which.min(abs(Spike_IgM_6m[n] - DEN_Spike_IgM_6m$x))

		points( x = 180 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgM_6m$y[index], max=DEN_Spike_IgM_6m$y[index]),
                    y = Spike_IgM_6m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgM_cols[c])
	}

	off_set <- off_set - 2*ww
}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_12m <- Spike_IgM_mod[[c]][,which(tt_plot==365),2]/Spike_IgM_norm[c]	

	DEN_Spike_IgM_12m = density( log(Spike_IgM_12m), cut=0, na.rm=TRUE )

	DEN_Spike_IgM_12m$x <- exp(DEN_Spike_IgM_12m$x)
	DEN_Spike_IgM_12m$y <- DEN_Spike_IgM_12m$y/max(DEN_Spike_IgM_12m$y)

	polygon( x = 365 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_Spike_IgM_12m$y, rev(DEN_Spike_IgM_12m$y) ),
		   y=c( DEN_Spike_IgM_12m$x, rev(DEN_Spike_IgM_12m$x) ), 
         	   col=Spike_IgM_rgb[c], border=NA)

	for(n in 1:length(Spike_IgM_12m))
	{
		index <- which.min(abs(Spike_IgM_12m[n] - DEN_Spike_IgM_12m$x))

		points( x = 365 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_Spike_IgM_12m$y[index], max=DEN_Spike_IgM_12m$y[index]),
                    y = Spike_IgM_12m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgM_cols[c])
	}

}



###################
## Positive points
## Month 24



for(c in 1:N_cohort_Spike_IgM)
{
	Spike_IgM_24m <- Spike_IgM_mod[[c]][,which(tt_plot==2*365),2]/Spike_IgM_norm[c]	

	DEN_Spike_IgM_24m = density( log(Spike_IgM_24m), cut=0, na.rm=TRUE )

	DEN_Spike_IgM_24m$x <- exp(DEN_Spike_IgM_24m$x)
	DEN_Spike_IgM_24m$y <- DEN_Spike_IgM_24m$y/max(DEN_Spike_IgM_24m$y)

	polygon( x = 2*365 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) +
                   ww*c( -DEN_Spike_IgM_24m$y, rev(DEN_Spike_IgM_24m$y) ),
		   y=c( DEN_Spike_IgM_24m$x, rev(DEN_Spike_IgM_24m$x) ), 
         	   col=Spike_IgM_rgb[c], border=NA)

	for(n in 1:length(Spike_IgM_24m))
	{
		index <- which.min(abs(Spike_IgM_24m[n] - DEN_Spike_IgM_24m$x))

		points( x = 2*365 - N_cohort_Spike_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgM_24m$y[index], max=DEN_Spike_IgM_24m$y[index]),
                    y = Spike_IgM_24m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgM_cols[c])
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
###  PANEL 5    ###
###  RBD IgM    ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )



###################
## Negative controls

for(c in 1:N_cohort_RBD_IgM)
{
	if( RBD_IgM_N_neg[c] > 0 )
	{
		DEN_RBD_IgM_neg = density( log(RBD_IgM_neg_data[[c]]/RBD_IgM_norm[c]), cut=0, na.rm=TRUE )

		DEN_RBD_IgM_neg$x <- exp(DEN_RBD_IgM_neg$x)
		DEN_RBD_IgM_neg$y <- DEN_RBD_IgM_neg$y/max(DEN_RBD_IgM_neg$y)

		polygon( x = -90 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) +
                         ww*c( -DEN_RBD_IgM_neg$y, rev(DEN_RBD_IgM_neg$y) ),
			   y=c( DEN_RBD_IgM_neg$x, rev(DEN_RBD_IgM_neg$x) ), 
      	   	   col=RBD_IgM_rgb[c], border=NA)

		for(n in 1:RBD_IgM_N_neg[c])
		{
			index <- which.min(abs(RBD_IgM_neg_data[[c]][n]/RBD_IgM_norm[c] - DEN_RBD_IgM_neg$x))

			points( x = -90 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_RBD_IgM_neg$y[index], max=DEN_RBD_IgM_neg$y[index]),
                          y = RBD_IgM_neg_data[[c]][n]/RBD_IgM_norm[c],
                          pch=19, cex=neg.point.size, col=RBD_IgM_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_1m <- RBD_IgM_mod[[c]][,which(tt_plot==14),2]/RBD_IgM_norm[c]	

	DEN_RBD_IgM_1m = density( log(RBD_IgM_1m), cut=0, na.rm=TRUE )

	DEN_RBD_IgM_1m$x <- exp(DEN_RBD_IgM_1m$x)
	DEN_RBD_IgM_1m$y <- DEN_RBD_IgM_1m$y/max(DEN_RBD_IgM_1m$y)

	polygon( x = 30 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_RBD_IgM_1m$y, rev(DEN_RBD_IgM_1m$y) ),
		   y=c( DEN_RBD_IgM_1m$x, rev(DEN_RBD_IgM_1m$x) ), 
         	   col=RBD_IgM_rgb[c], border=NA)

	for(n in 1:length(RBD_IgM_1m))
	{
		index <- which.min(abs(RBD_IgM_1m[n] - DEN_RBD_IgM_1m$x))

		points( x = 30 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgM_1m$y[index], max=DEN_RBD_IgM_1m$y[index]),
                    y = RBD_IgM_1m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgM_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_6m <- RBD_IgM_mod[[c]][,which(tt_plot==180),2]/RBD_IgM_norm[c]	

	DEN_RBD_IgM_6m = density( log(RBD_IgM_6m), cut=0, na.rm=TRUE )

	DEN_RBD_IgM_6m$x <- exp(DEN_RBD_IgM_6m$x)
	DEN_RBD_IgM_6m$y <- DEN_RBD_IgM_6m$y/max(DEN_RBD_IgM_6m$y)

	polygon( x = 180 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) +
                   ww*c( -DEN_RBD_IgM_6m$y, rev(DEN_RBD_IgM_6m$y) ),
		   y=c( DEN_RBD_IgM_6m$x, rev(DEN_RBD_IgM_6m$x) ), 
         	   col=RBD_IgM_rgb[c], border=NA)

	for(n in 1:length(RBD_IgM_6m))
	{
		index <- which.min(abs(RBD_IgM_6m[n] - DEN_RBD_IgM_6m$x))

		points( x = 180 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgM_6m$y[index], max=DEN_RBD_IgM_6m$y[index]),
                    y = RBD_IgM_6m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgM_cols[c])
	}

	off_set <- off_set - 2*ww
}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_12m <- RBD_IgM_mod[[c]][,which(tt_plot==365),2]/RBD_IgM_norm[c]	

	DEN_RBD_IgM_12m = density( log(RBD_IgM_12m), cut=0, na.rm=TRUE )

	DEN_RBD_IgM_12m$x <- exp(DEN_RBD_IgM_12m$x)
	DEN_RBD_IgM_12m$y <- DEN_RBD_IgM_12m$y/max(DEN_RBD_IgM_12m$y)

	polygon( x = 365 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_RBD_IgM_12m$y, rev(DEN_RBD_IgM_12m$y) ),
		   y=c( DEN_RBD_IgM_12m$x, rev(DEN_RBD_IgM_12m$x) ), 
         	   col=RBD_IgM_rgb[c], border=NA)

	for(n in 1:length(RBD_IgM_12m))
	{
		index <- which.min(abs(RBD_IgM_12m[n] - DEN_RBD_IgM_12m$x))

		points( x = 365 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_RBD_IgM_12m$y[index], max=DEN_RBD_IgM_12m$y[index]),
                    y = RBD_IgM_12m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgM_cols[c])
	}

}



###################
## Positive points
## Month 24



for(c in 1:N_cohort_RBD_IgM)
{
	RBD_IgM_24m <- RBD_IgM_mod[[c]][,which(tt_plot==2*365),2]/RBD_IgM_norm[c]	

	DEN_RBD_IgM_24m = density( log(RBD_IgM_24m), cut=0, na.rm=TRUE )

	DEN_RBD_IgM_24m$x <- exp(DEN_RBD_IgM_24m$x)
	DEN_RBD_IgM_24m$y <- DEN_RBD_IgM_24m$y/max(DEN_RBD_IgM_24m$y)

	polygon( x = 2*365 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) +
                   ww*c( -DEN_RBD_IgM_24m$y, rev(DEN_RBD_IgM_24m$y) ),
		   y=c( DEN_RBD_IgM_24m$x, rev(DEN_RBD_IgM_24m$x) ), 
         	   col=RBD_IgM_rgb[c], border=NA)

	for(n in 1:length(RBD_IgM_24m))
	{
		index <- which.min(abs(RBD_IgM_24m[n] - DEN_RBD_IgM_24m$x))

		points( x = 2*365 - N_cohort_RBD_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgM_24m$y[index], max=DEN_RBD_IgM_24m$y[index]),
                    y = RBD_IgM_24m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgM_cols[c])
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
###  PANEL 6    ###
###  NP IgM     ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


###################
## Negative controls

for(c in 1:N_cohort_NP_IgM)
{
	if( NP_IgM_N_neg[c] > 0 )
	{
		DEN_NP_IgM_neg = density( log(NP_IgM_neg_data[[c]]/NP_IgM_norm[c]), cut=0, na.rm=TRUE )

		DEN_NP_IgM_neg$x <- exp(DEN_NP_IgM_neg$x)
		DEN_NP_IgM_neg$y <- DEN_NP_IgM_neg$y/max(DEN_NP_IgM_neg$y)

		polygon( x = -90 - N_cohort_NP_IgM*ww + 2*ww*(c-1) +
                         ww*c( -DEN_NP_IgM_neg$y, rev(DEN_NP_IgM_neg$y) ),
			   y=c( DEN_NP_IgM_neg$x, rev(DEN_NP_IgM_neg$x) ), 
      	   	   col=NP_IgM_rgb[c], border=NA)

		for(n in 1:NP_IgM_N_neg[c])
		{
			index <- which.min(abs(NP_IgM_neg_data[[c]][n]/NP_IgM_norm[c] - DEN_NP_IgM_neg$x))

			points( x = -90 - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_NP_IgM_neg$y[index], max=DEN_NP_IgM_neg$y[index]),
                          y = NP_IgM_neg_data[[c]][n]/NP_IgM_norm[c],
                          pch=19, cex=neg.point.size, col=NP_IgM_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_1m <- NP_IgM_mod[[c]][,which(tt_plot==14),2]/NP_IgM_norm[c]	

	DEN_NP_IgM_1m = density( log(NP_IgM_1m), cut=0, na.rm=TRUE )

	DEN_NP_IgM_1m$x <- exp(DEN_NP_IgM_1m$x)
	DEN_NP_IgM_1m$y <- DEN_NP_IgM_1m$y/max(DEN_NP_IgM_1m$y)

	polygon( x = 30 - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_NP_IgM_1m$y, rev(DEN_NP_IgM_1m$y) ),
		   y=c( DEN_NP_IgM_1m$x, rev(DEN_NP_IgM_1m$x) ), 
         	   col=NP_IgM_rgb[c], border=NA)

	for(n in 1:length(NP_IgM_1m))
	{
		index <- which.min(abs(NP_IgM_1m[n] - DEN_NP_IgM_1m$x))

		points( x = 30 - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgM_1m$y[index], max=DEN_NP_IgM_1m$y[index]),
                    y = NP_IgM_1m[n],
                    pch=19, cex=pos.point.size, col=NP_IgM_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_6m <- NP_IgM_mod[[c]][,which(tt_plot==180),2]/NP_IgM_norm[c]	

	DEN_NP_IgM_6m = density( log(NP_IgM_6m), cut=0, na.rm=TRUE )

	DEN_NP_IgM_6m$x <- exp(DEN_NP_IgM_6m$x)
	DEN_NP_IgM_6m$y <- DEN_NP_IgM_6m$y/max(DEN_NP_IgM_6m$y)

	polygon( x = 180 - N_cohort_NP_IgM*ww + 2*ww*(c-1) +
                   ww*c( -DEN_NP_IgM_6m$y, rev(DEN_NP_IgM_6m$y) ),
		   y=c( DEN_NP_IgM_6m$x, rev(DEN_NP_IgM_6m$x) ), 
         	   col=NP_IgM_rgb[c], border=NA)

	for(n in 1:length(NP_IgM_6m))
	{
		index <- which.min(abs(NP_IgM_6m[n] - DEN_NP_IgM_6m$x))

		points( x = 180 - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgM_6m$y[index], max=DEN_NP_IgM_6m$y[index]),
                    y = NP_IgM_6m[n],
                    pch=19, cex=pos.point.size, col=NP_IgM_cols[c])
	}

	off_set <- off_set - 2*ww
}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_12m <- NP_IgM_mod[[c]][,which(tt_plot==365),2]/NP_IgM_norm[c]	

	DEN_NP_IgM_12m = density( log(NP_IgM_12m), cut=0, na.rm=TRUE )

	DEN_NP_IgM_12m$x <- exp(DEN_NP_IgM_12m$x)
	DEN_NP_IgM_12m$y <- DEN_NP_IgM_12m$y/max(DEN_NP_IgM_12m$y)

	polygon( x = 365 - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_NP_IgM_12m$y, rev(DEN_NP_IgM_12m$y) ),
		   y=c( DEN_NP_IgM_12m$x, rev(DEN_NP_IgM_12m$x) ), 
         	   col=NP_IgM_rgb[c], border=NA)

	for(n in 1:length(NP_IgM_12m))
	{
		index <- which.min(abs(NP_IgM_12m[n] - DEN_NP_IgM_12m$x))

		points( x = 365 - N_cohort_NP_IgM*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_NP_IgM_12m$y[index], max=DEN_NP_IgM_12m$y[index]),
                    y = NP_IgM_12m[n],
                    pch=19, cex=pos.point.size, col=NP_IgM_cols[c])
	}

}



###################
## Positive points
## Month 24



for(c in 1:N_cohort_NP_IgM)
{
	NP_IgM_24m <- NP_IgM_mod[[c]][,which(tt_plot==2*365),2]/NP_IgM_norm[c]	

	DEN_NP_IgM_24m = density( log(NP_IgM_24m), cut=0, na.rm=TRUE )

	DEN_NP_IgM_24m$x <- exp(DEN_NP_IgM_24m$x)
	DEN_NP_IgM_24m$y <- DEN_NP_IgM_24m$y/max(DEN_NP_IgM_24m$y)

	polygon( x = 2*365 - N_cohort_NP_IgM*ww + 2*ww*(c-1) +
                   ww*c( -DEN_NP_IgM_24m$y, rev(DEN_NP_IgM_24m$y) ),
		   y=c( DEN_NP_IgM_24m$x, rev(DEN_NP_IgM_24m$x) ), 
         	   col=NP_IgM_rgb[c], border=NA)

	for(n in 1:length(NP_IgM_24m))
	{
		index <- which.min(abs(NP_IgM_24m[n] - DEN_NP_IgM_24m$x))

		points( x = 2*365 - N_cohort_NP_IgM*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgM_24m$y[index], max=DEN_NP_IgM_24m$y[index]),
                    y = NP_IgM_24m[n],
                    pch=19, cex=pos.point.size, col=NP_IgM_cols[c])
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
###  PANEL 7    ###
###  Spike IgA  ###
###             ###
###################
###################

par(mar=c(8,8,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )

mtext(side = 2, line = 6, 
cex=lab.size, 
text="IgA antibody level")

mtext(side = 2, line = 4, 
cex=lab.size, 
text="(relative to day 14)")

mtext(side = 1, line = 4, 
cex=lab.size, 
text="months after symptom onset")


###################
## Negative controls

for(c in 1:N_cohort_Spike_IgA)
{
	if( Spike_IgA_N_neg[c] > 0 )
	{
		DEN_Spike_IgA_neg = density( log(Spike_IgA_neg_data[[c]]/Spike_IgA_norm[c]), cut=0, na.rm=TRUE )

		DEN_Spike_IgA_neg$x <- exp(DEN_Spike_IgA_neg$x)
		DEN_Spike_IgA_neg$y <- DEN_Spike_IgA_neg$y/max(DEN_Spike_IgA_neg$y)

		polygon( x = -90 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) +
                         ww*c( -DEN_Spike_IgA_neg$y, rev(DEN_Spike_IgA_neg$y) ),
			   y=c( DEN_Spike_IgA_neg$x, rev(DEN_Spike_IgA_neg$x) ), 
      	   	   col=Spike_IgA_rgb[c], border=NA)

		for(n in 1:Spike_IgA_N_neg[c])
		{
			index <- which.min(abs(Spike_IgA_neg_data[[c]][n]/Spike_IgA_norm[c] - DEN_Spike_IgA_neg$x))

			points( x = -90 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_Spike_IgA_neg$y[index], max=DEN_Spike_IgA_neg$y[index]),
                          y = Spike_IgA_neg_data[[c]][n]/Spike_IgA_norm[c],
                          pch=19, cex=neg.point.size, col=Spike_IgA_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_1m <- Spike_IgA_mod[[c]][,which(tt_plot==14),2]/Spike_IgA_norm[c]	

	DEN_Spike_IgA_1m = density( log(Spike_IgA_1m), cut=0, na.rm=TRUE )

	DEN_Spike_IgA_1m$x <- exp(DEN_Spike_IgA_1m$x)
	DEN_Spike_IgA_1m$y <- DEN_Spike_IgA_1m$y/max(DEN_Spike_IgA_1m$y)

	polygon( x = 30 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_Spike_IgA_1m$y, rev(DEN_Spike_IgA_1m$y) ),
		   y=c( DEN_Spike_IgA_1m$x, rev(DEN_Spike_IgA_1m$x) ), 
         	   col=Spike_IgA_rgb[c], border=NA)

	for(n in 1:length(Spike_IgA_1m))
	{
		index <- which.min(abs(Spike_IgA_1m[n] - DEN_Spike_IgA_1m$x))

		points( x = 30 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgA_1m$y[index], max=DEN_Spike_IgA_1m$y[index]),
                    y = Spike_IgA_1m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgA_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_6m <- Spike_IgA_mod[[c]][,which(tt_plot==180),2]/Spike_IgA_norm[c]	

	DEN_Spike_IgA_6m = density( log(Spike_IgA_6m), cut=0, na.rm=TRUE )

	DEN_Spike_IgA_6m$x <- exp(DEN_Spike_IgA_6m$x)
	DEN_Spike_IgA_6m$y <- DEN_Spike_IgA_6m$y/max(DEN_Spike_IgA_6m$y)

	polygon( x = 180 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) +
                   ww*c( -DEN_Spike_IgA_6m$y, rev(DEN_Spike_IgA_6m$y) ),
		   y=c( DEN_Spike_IgA_6m$x, rev(DEN_Spike_IgA_6m$x) ), 
         	   col=Spike_IgA_rgb[c], border=NA)

	for(n in 1:length(Spike_IgA_6m))
	{
		index <- which.min(abs(Spike_IgA_6m[n] - DEN_Spike_IgA_6m$x))

		points( x = 180 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgA_6m$y[index], max=DEN_Spike_IgA_6m$y[index]),
                    y = Spike_IgA_6m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgA_cols[c])
	}

	off_set <- off_set - 2*ww
}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_12m <- Spike_IgA_mod[[c]][,which(tt_plot==365),2]/Spike_IgA_norm[c]	

	DEN_Spike_IgA_12m = density( log(Spike_IgA_12m), cut=0, na.rm=TRUE )

	DEN_Spike_IgA_12m$x <- exp(DEN_Spike_IgA_12m$x)
	DEN_Spike_IgA_12m$y <- DEN_Spike_IgA_12m$y/max(DEN_Spike_IgA_12m$y)

	polygon( x = 365 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_Spike_IgA_12m$y, rev(DEN_Spike_IgA_12m$y) ),
		   y=c( DEN_Spike_IgA_12m$x, rev(DEN_Spike_IgA_12m$x) ), 
         	   col=Spike_IgA_rgb[c], border=NA)

	for(n in 1:length(Spike_IgA_12m))
	{
		index <- which.min(abs(Spike_IgA_12m[n] - DEN_Spike_IgA_12m$x))

		points( x = 365 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_Spike_IgA_12m$y[index], max=DEN_Spike_IgA_12m$y[index]),
                    y = Spike_IgA_12m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgA_cols[c])
	}

}



###################
## Positive points
## Month 24



for(c in 1:N_cohort_Spike_IgA)
{
	Spike_IgA_24m <- Spike_IgA_mod[[c]][,which(tt_plot==2*365),2]/Spike_IgA_norm[c]	

	DEN_Spike_IgA_24m = density( log(Spike_IgA_24m), cut=0, na.rm=TRUE )

	DEN_Spike_IgA_24m$x <- exp(DEN_Spike_IgA_24m$x)
	DEN_Spike_IgA_24m$y <- DEN_Spike_IgA_24m$y/max(DEN_Spike_IgA_24m$y)

	polygon( x = 2*365 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) +
                   ww*c( -DEN_Spike_IgA_24m$y, rev(DEN_Spike_IgA_24m$y) ),
		   y=c( DEN_Spike_IgA_24m$x, rev(DEN_Spike_IgA_24m$x) ), 
         	   col=Spike_IgA_rgb[c], border=NA)

	for(n in 1:length(Spike_IgA_24m))
	{
		index <- which.min(abs(Spike_IgA_24m[n] - DEN_Spike_IgA_24m$x))

		points( x = 2*365 - N_cohort_Spike_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_Spike_IgA_24m$y[index], max=DEN_Spike_IgA_24m$y[index]),
                    y = Spike_IgA_24m[n],
                    pch=19, cex=pos.point.size, col=Spike_IgA_cols[c])
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
###  PANEL 8    ###
###  RBD IgA    ###
###             ###
###################
###################

par(mar=c(8,3,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )

mtext(side = 1, line = 4, 
cex=lab.size, 
text="months after symptom onset")


###################
## Negative controls

for(c in 1:N_cohort_RBD_IgA)
{
	if( RBD_IgA_N_neg[c] > 0 )
	{
		DEN_RBD_IgA_neg = density( log(RBD_IgA_neg_data[[c]]/RBD_IgA_norm[c]), cut=0, na.rm=TRUE )

		DEN_RBD_IgA_neg$x <- exp(DEN_RBD_IgA_neg$x)
		DEN_RBD_IgA_neg$y <- DEN_RBD_IgA_neg$y/max(DEN_RBD_IgA_neg$y)

		polygon( x = -90 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) +
                         ww*c( -DEN_RBD_IgA_neg$y, rev(DEN_RBD_IgA_neg$y) ),
			   y=c( DEN_RBD_IgA_neg$x, rev(DEN_RBD_IgA_neg$x) ), 
      	   	   col=RBD_IgA_rgb[c], border=NA)

		for(n in 1:RBD_IgA_N_neg[c])
		{
			index <- which.min(abs(RBD_IgA_neg_data[[c]][n]/RBD_IgA_norm[c] - DEN_RBD_IgA_neg$x))

			points( x = -90 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_RBD_IgA_neg$y[index], max=DEN_RBD_IgA_neg$y[index]),
                          y = RBD_IgA_neg_data[[c]][n]/RBD_IgA_norm[c],
                          pch=19, cex=neg.point.size, col=RBD_IgA_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_1m <- RBD_IgA_mod[[c]][,which(tt_plot==14),2]/RBD_IgA_norm[c]	

	DEN_RBD_IgA_1m = density( log(RBD_IgA_1m), cut=0, na.rm=TRUE )

	DEN_RBD_IgA_1m$x <- exp(DEN_RBD_IgA_1m$x)
	DEN_RBD_IgA_1m$y <- DEN_RBD_IgA_1m$y/max(DEN_RBD_IgA_1m$y)

	polygon( x = 30 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_RBD_IgA_1m$y, rev(DEN_RBD_IgA_1m$y) ),
		   y=c( DEN_RBD_IgA_1m$x, rev(DEN_RBD_IgA_1m$x) ), 
         	   col=RBD_IgA_rgb[c], border=NA)

	for(n in 1:length(RBD_IgA_1m))
	{
		index <- which.min(abs(RBD_IgA_1m[n] - DEN_RBD_IgA_1m$x))

		points( x = 30 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgA_1m$y[index], max=DEN_RBD_IgA_1m$y[index]),
                    y = RBD_IgA_1m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgA_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_6m <- RBD_IgA_mod[[c]][,which(tt_plot==180),2]/RBD_IgA_norm[c]	

	DEN_RBD_IgA_6m = density( log(RBD_IgA_6m), cut=0, na.rm=TRUE )

	DEN_RBD_IgA_6m$x <- exp(DEN_RBD_IgA_6m$x)
	DEN_RBD_IgA_6m$y <- DEN_RBD_IgA_6m$y/max(DEN_RBD_IgA_6m$y)

	polygon( x = 180 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) +
                   ww*c( -DEN_RBD_IgA_6m$y, rev(DEN_RBD_IgA_6m$y) ),
		   y=c( DEN_RBD_IgA_6m$x, rev(DEN_RBD_IgA_6m$x) ), 
         	   col=RBD_IgA_rgb[c], border=NA)

	for(n in 1:length(RBD_IgA_6m))
	{
		index <- which.min(abs(RBD_IgA_6m[n] - DEN_RBD_IgA_6m$x))

		points( x = 180 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgA_6m$y[index], max=DEN_RBD_IgA_6m$y[index]),
                    y = RBD_IgA_6m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgA_cols[c])
	}

}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_12m <- RBD_IgA_mod[[c]][,which(tt_plot==365),2]/RBD_IgA_norm[c]	

	DEN_RBD_IgA_12m = density( log(RBD_IgA_12m), cut=0, na.rm=TRUE )

	DEN_RBD_IgA_12m$x <- exp(DEN_RBD_IgA_12m$x)
	DEN_RBD_IgA_12m$y <- DEN_RBD_IgA_12m$y/max(DEN_RBD_IgA_12m$y)

	polygon( x = 365 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_RBD_IgA_12m$y, rev(DEN_RBD_IgA_12m$y) ),
		   y=c( DEN_RBD_IgA_12m$x, rev(DEN_RBD_IgA_12m$x) ), 
         	   col=RBD_IgA_rgb[c], border=NA)

	for(n in 1:length(RBD_IgA_12m))
	{
		index <- which.min(abs(RBD_IgA_12m[n] - DEN_RBD_IgA_12m$x))

		points( x = 365 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_RBD_IgA_12m$y[index], max=DEN_RBD_IgA_12m$y[index]),
                    y = RBD_IgA_12m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgA_cols[c])
	}

}



###################
## Positive points
## Month 24



for(c in 1:N_cohort_RBD_IgA)
{
	RBD_IgA_24m <- RBD_IgA_mod[[c]][,which(tt_plot==2*365),2]/RBD_IgA_norm[c]	

	DEN_RBD_IgA_24m = density( log(RBD_IgA_24m), cut=0, na.rm=TRUE )

	DEN_RBD_IgA_24m$x <- exp(DEN_RBD_IgA_24m$x)
	DEN_RBD_IgA_24m$y <- DEN_RBD_IgA_24m$y/max(DEN_RBD_IgA_24m$y)

	polygon( x = 2*365 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) +
                   ww*c( -DEN_RBD_IgA_24m$y, rev(DEN_RBD_IgA_24m$y) ),
		   y=c( DEN_RBD_IgA_24m$x, rev(DEN_RBD_IgA_24m$x) ), 
         	   col=RBD_IgA_rgb[c], border=NA)

	for(n in 1:length(RBD_IgA_24m))
	{
		index <- which.min(abs(RBD_IgA_24m[n] - DEN_RBD_IgA_24m$x))

		points( x = 2*365 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_RBD_IgA_24m$y[index], max=DEN_RBD_IgA_24m$y[index]),
                    y = RBD_IgA_24m[n],
                    pch=19, cex=pos.point.size, col=RBD_IgA_cols[c])
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
###  PANEL 9    ###
###  NP IgA     ###
###             ###
###################
###################

par(mar=c(8,3,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-120,750), ylim=c(0.001,100), log="y",
	yaxt='n', xaxt='n', bty='n',
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




axis(1, at = c(-90, 30, 180, 365, 2*365), 
        label=c("neg", "1", "6", "12", "24"), 
        cex.axis=0.8*axis.size )


axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )

mtext(side = 1, line = 4, 
cex=lab.size, 
text="months after symptom onset")


###################
## Negative controls

for(c in 1:N_cohort_NP_IgA)
{
	if( NP_IgA_N_neg[c] > 0 )
	{
		DEN_NP_IgA_neg = density( log(NP_IgA_neg_data[[c]]/NP_IgA_norm[c]), cut=0, na.rm=TRUE )

		DEN_NP_IgA_neg$x <- exp(DEN_NP_IgA_neg$x)
		DEN_NP_IgA_neg$y <- DEN_NP_IgA_neg$y/max(DEN_NP_IgA_neg$y)

		polygon( x = -90 - N_cohort_RBD_IgA*ww + 2*ww*(c-1) +
                         ww*c( -DEN_NP_IgA_neg$y, rev(DEN_NP_IgA_neg$y) ),
			   y=c( DEN_NP_IgA_neg$x, rev(DEN_NP_IgA_neg$x) ), 
      	   	   col=NP_IgA_rgb[c], border=NA)

		for(n in 1:NP_IgA_N_neg[c])
		{
			index <- which.min(abs(NP_IgA_neg_data[[c]][n]/NP_IgA_norm[c] - DEN_NP_IgA_neg$x))

			points( x = -90 - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
                              ww*runif(1, min=-DEN_NP_IgA_neg$y[index], max=DEN_NP_IgA_neg$y[index]),
                          y = NP_IgA_neg_data[[c]][n]/NP_IgA_norm[c],
                          pch=19, cex=neg.point.size, col=NP_IgA_cols[c])
		}
	}

}



###################
## Positive points
## Month 1


for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_1m <- NP_IgA_mod[[c]][,which(tt_plot==14),2]/NP_IgA_norm[c]	

	DEN_NP_IgA_1m = density( log(NP_IgA_1m), cut=0, na.rm=TRUE )

	DEN_NP_IgA_1m$x <- exp(DEN_NP_IgA_1m$x)
	DEN_NP_IgA_1m$y <- DEN_NP_IgA_1m$y/max(DEN_NP_IgA_1m$y)

	polygon( x = 30 - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_NP_IgA_1m$y, rev(DEN_NP_IgA_1m$y) ),
		   y=c( DEN_NP_IgA_1m$x, rev(DEN_NP_IgA_1m$x) ), 
         	   col=NP_IgA_rgb[c], border=NA)

	for(n in 1:length(NP_IgA_1m))
	{
		index <- which.min(abs(NP_IgA_1m[n] - DEN_NP_IgA_1m$x))

		points( x = 30 - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgA_1m$y[index], max=DEN_NP_IgA_1m$y[index]),
                    y = NP_IgA_1m[n],
                    pch=19, cex=pos.point.size, col=NP_IgA_cols[c])
	}

}



###################
## Positive points
## Month 6

for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_6m <- NP_IgA_mod[[c]][,which(tt_plot==180),2]/NP_IgA_norm[c]	

	DEN_NP_IgA_6m = density( log(NP_IgA_6m), cut=0, na.rm=TRUE )

	DEN_NP_IgA_6m$x <- exp(DEN_NP_IgA_6m$x)
	DEN_NP_IgA_6m$y <- DEN_NP_IgA_6m$y/max(DEN_NP_IgA_6m$y)

	polygon( x = 180 - N_cohort_NP_IgA*ww + 2*ww*(c-1) +
                   ww*c( -DEN_NP_IgA_6m$y, rev(DEN_NP_IgA_6m$y) ),
		   y=c( DEN_NP_IgA_6m$x, rev(DEN_NP_IgA_6m$x) ), 
         	   col=NP_IgA_rgb[c], border=NA)

	for(n in 1:length(NP_IgA_6m))
	{
		index <- which.min(abs(NP_IgA_6m[n] - DEN_NP_IgA_6m$x))

		points( x = 180 - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgA_6m$y[index], max=DEN_NP_IgA_6m$y[index]),
                    y = NP_IgA_6m[n],
                    pch=19, cex=pos.point.size, col=NP_IgA_cols[c])
	}

}




###################
## Positive points
## Month 12


for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_12m <- NP_IgA_mod[[c]][,which(tt_plot==365),2]/NP_IgA_norm[c]	

	DEN_NP_IgA_12m = density( log(NP_IgA_12m), cut=0, na.rm=TRUE )

	DEN_NP_IgA_12m$x <- exp(DEN_NP_IgA_12m$x)
	DEN_NP_IgA_12m$y <- DEN_NP_IgA_12m$y/max(DEN_NP_IgA_12m$y)

	polygon( x = 365 - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
                   ww*c( -DEN_NP_IgA_12m$y, rev(DEN_NP_IgA_12m$y) ),
		   y=c( DEN_NP_IgA_12m$x, rev(DEN_NP_IgA_12m$x) ), 
         	   col=NP_IgA_rgb[c], border=NA)

	for(n in 1:length(NP_IgA_12m))
	{
		index <- which.min(abs(NP_IgA_12m[n] - DEN_NP_IgA_12m$x))

		points( x = 365 - N_cohort_NP_IgA*ww + 2*ww*(c-1) +  
                        ww*runif(1, min=-DEN_NP_IgA_12m$y[index], max=DEN_NP_IgA_12m$y[index]),
                    y = NP_IgA_12m[n],
                    pch=19, cex=pos.point.size, col=NP_IgA_cols[c])
	}

}



###################
## Positive points
## Month 24



for(c in 1:N_cohort_NP_IgA)
{
	NP_IgA_24m <- NP_IgA_mod[[c]][,which(tt_plot==2*365),2]/NP_IgA_norm[c]	

	DEN_NP_IgA_24m = density( log(NP_IgA_24m), cut=0, na.rm=TRUE )

	DEN_NP_IgA_24m$x <- exp(DEN_NP_IgA_24m$x)
	DEN_NP_IgA_24m$y <- DEN_NP_IgA_24m$y/max(DEN_NP_IgA_24m$y)

	polygon( x = 2*365 - N_cohort_NP_IgA*ww + 2*ww*(c-1) +
                   ww*c( -DEN_NP_IgA_24m$y, rev(DEN_NP_IgA_24m$y) ),
		   y=c( DEN_NP_IgA_24m$x, rev(DEN_NP_IgA_24m$x) ), 
         	   col=NP_IgA_rgb[c], border=NA)

	for(n in 1:length(NP_IgA_24m))
	{
		index <- which.min(abs(NP_IgA_24m[n] - DEN_NP_IgA_24m$x))

		points( x = 2*365 - N_cohort_NP_IgA*ww + 2*ww*(c-1) + 
                        ww*runif(1, min=-DEN_NP_IgA_24m$y[index], max=DEN_NP_IgA_24m$y[index]),
                    y = NP_IgA_24m[n],
                    pch=19, cex=pos.point.size, col=NP_IgA_cols[c])
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
       legend = cohort_names, 
       col = cohort_cols, 
       pch=rep(19,length(cohort_names)),
       ncol=length(cohort_names), cex=2.0, bty="n" )


dev.off()









