

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

line_seq_x <- c(-30, 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390)
line_seq_y <- c(0.001, 0.01, 0.1, 1, 10, 100)

par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))

pos.point.size = 0.5
neg.point.size = 0.2
lab.size   = 1.5
axis.size  = 1.25	
main.size  = 1.5
line.size  = 2
dash.line.size = 0.5

ww <- 3




tiff(file="SupFig5_Abkin_data_overview.tif", width=40, height=32, units="cm", res=500)


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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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



axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
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
## Positive points

for(c in 1:N_cohort_Spike_IgG)
{
	for(n in 1:Spike_IgG_N_part[c])
	{
		points( x = Spike_IgG_data[[c]][n,,2],
                    y = Spike_IgG_data[[c]][n,,1]/Spike_IgG_norm[c],
			  pch=19, cex=pos.point.size, col=Spike_IgG_cols[c])   
	}
}


for(n in 1:Spike_IgG_N_part[3])
{
	index <- which( Spike_IgG_data[[3]][n,,2] > 330 )

	points( x = 360 + (Spike_IgG_data[[3]][n,index,2] - 360)*((420-360)/(2400-360)),
              y = Spike_IgG_data[[3]][n,index,1]/Spike_IgG_norm[3],
		  pch=19, cex=pos.point.size, col=Spike_IgG_cols[3])   
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_Spike_IgG)
{
	if( Spike_IgG_N_neg[c] > 0 )
	{
		DEN_Spike_IgG_neg = density( Spike_IgG_neg_data[[c]]/Spike_IgG_norm[c], cut=0, na.rm=TRUE )

		DEN_Spike_IgG_neg$y <- DEN_Spike_IgG_neg$y/max(DEN_Spike_IgG_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_Spike_IgG_neg$y, rev(DEN_Spike_IgG_neg$y) ),
			   y=c( DEN_Spike_IgG_neg$x, rev(DEN_Spike_IgG_neg$x) ), 
      	   	   col=Spike_IgG_rgb[c], border=NA)

		for(n in 1:Spike_IgG_N_neg[c])
		{
			index <- which.min(abs(Spike_IgG_neg_data[[c]][n]/Spike_IgG_norm[c] - DEN_Spike_IgG_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_Spike_IgG_neg$y[index], max=DEN_Spike_IgG_neg$y[index]),
                          y = Spike_IgG_neg_data[[c]][n]/Spike_IgG_norm[c],
                          pch=19, cex=neg.point.size, col=Spike_IgG_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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


axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


###################
## Positive points

for(c in 1:N_cohort_RBD_IgG)
{
	for(n in 1:RBD_IgG_N_part[c])
	{
		points( x = RBD_IgG_data[[c]][n,,2],
                    y = RBD_IgG_data[[c]][n,,1]/RBD_IgG_norm[c],
			  pch=19, cex=pos.point.size, col=RBD_IgG_cols[c])   
	}
}


for(n in 1:RBD_IgG_N_part[3])
{
	index <- which( RBD_IgG_data[[3]][n,,2] > 330 )

	points( x = 360 + (RBD_IgG_data[[3]][n,index,2] - 360)*((420-360)/(2400-360)),
              y = RBD_IgG_data[[3]][n,index,1]/RBD_IgG_norm[3],
		  pch=19, cex=pos.point.size, col=RBD_IgG_cols[3])   
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_RBD_IgG)
{
	if( RBD_IgG_N_neg[c] > 0 )
	{
		DEN_RBD_IgG_neg = density( RBD_IgG_neg_data[[c]]/RBD_IgG_norm[c], cut=0, na.rm=TRUE )

		DEN_RBD_IgG_neg$y <- DEN_RBD_IgG_neg$y/max(DEN_RBD_IgG_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_RBD_IgG_neg$y, rev(DEN_RBD_IgG_neg$y) ),
			   y=c( DEN_RBD_IgG_neg$x, rev(DEN_RBD_IgG_neg$x) ), 
      	   	   col=RBD_IgG_rgb[c], border=NA)

		for(n in 1:RBD_IgG_N_neg[c])
		{
			index <- which.min(abs(RBD_IgG_neg_data[[c]][n]/RBD_IgG_norm[c] - DEN_RBD_IgG_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_RBD_IgG_neg$y[index], max=DEN_RBD_IgG_neg$y[index]),
                          y = RBD_IgG_neg_data[[c]][n]/RBD_IgG_norm[c],
                          pch=19, cex=neg.point.size, col=RBD_IgG_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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

axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


###################
## Positive points

for(c in 1:N_cohort_NP_IgG)
{
	for(n in 1:NP_IgG_N_part[c])
	{
		points( x = NP_IgG_data[[c]][n,,2],
                    y = NP_IgG_data[[c]][n,,1]/NP_IgG_norm[c],
			  pch=19, cex=pos.point.size, col=NP_IgG_cols[c])   
	}
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_NP_IgG)
{
	if( NP_IgG_N_neg[c] > 0 )
	{
		DEN_NP_IgG_neg = density( NP_IgG_neg_data[[c]]/NP_IgG_norm[c], cut=0, na.rm=TRUE )

		DEN_NP_IgG_neg$y <- DEN_NP_IgG_neg$y/max(DEN_NP_IgG_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_NP_IgG_neg$y, rev(DEN_NP_IgG_neg$y) ),
			   y=c( DEN_NP_IgG_neg$x, rev(DEN_NP_IgG_neg$x) ), 
      	   	   col=NP_IgG_rgb[c], border=NA)

		for(n in 1:NP_IgG_N_neg[c])
		{
			index <- which.min(abs(NP_IgG_neg_data[[c]][n]/NP_IgG_norm[c] - DEN_NP_IgG_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_NP_IgG_neg$y[index], max=DEN_NP_IgG_neg$y[index]),
                          y = NP_IgG_neg_data[[c]][n]/NP_IgG_norm[c],
                          pch=19, cex=neg.point.size, col=NP_IgG_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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


axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
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
## Positive points

for(c in 1:N_cohort_Spike_IgM)
{
	for(n in 1:Spike_IgM_N_part[c])
	{
		points( x = Spike_IgM_data[[c]][n,,2],
                    y = Spike_IgM_data[[c]][n,,1]/Spike_IgM_norm[c],
			  pch=19, cex=pos.point.size, col=Spike_IgM_cols[c])   
	}
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_Spike_IgM)
{
	if( Spike_IgM_N_neg[c] > 0 )
	{
		DEN_Spike_IgM_neg = density( Spike_IgM_neg_data[[c]]/Spike_IgM_norm[c], cut=0, na.rm=TRUE )

		DEN_Spike_IgM_neg$y <- DEN_Spike_IgM_neg$y/max(DEN_Spike_IgM_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_Spike_IgM_neg$y, rev(DEN_Spike_IgM_neg$y) ),
			   y=c( DEN_Spike_IgM_neg$x, rev(DEN_Spike_IgM_neg$x) ), 
      	   	   col=Spike_IgM_rgb[c], border=NA)

		for(n in 1:Spike_IgM_N_neg[c])
		{
			index <- which.min(abs(Spike_IgM_neg_data[[c]][n]/Spike_IgM_norm[c] - DEN_Spike_IgM_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_Spike_IgM_neg$y[index], max=DEN_Spike_IgM_neg$y[index]),
                          y = Spike_IgM_neg_data[[c]][n]/Spike_IgM_norm[c],
                          pch=19, cex=neg.point.size, col=Spike_IgM_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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


axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


###################
## Positive points

for(c in 1:N_cohort_RBD_IgM)
{
	for(n in 1:RBD_IgM_N_part[c])
	{
		points( x = RBD_IgM_data[[c]][n,,2],
                    y = RBD_IgM_data[[c]][n,,1]/RBD_IgM_norm[c],
			  pch=19, cex=pos.point.size, col=RBD_IgM_cols[c])   
	}
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_RBD_IgM)
{
	if( RBD_IgM_N_neg[c] > 0 )
	{
		DEN_RBD_IgM_neg = density( RBD_IgM_neg_data[[c]]/RBD_IgM_norm[c], cut=0, na.rm=TRUE )

		DEN_RBD_IgM_neg$y <- DEN_RBD_IgM_neg$y/max(DEN_RBD_IgM_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_RBD_IgM_neg$y, rev(DEN_RBD_IgM_neg$y) ),
			   y=c( DEN_RBD_IgM_neg$x, rev(DEN_RBD_IgM_neg$x) ), 
      	   	   col=RBD_IgM_rgb[c], border=NA)

		for(n in 1:RBD_IgM_N_neg[c])
		{
			index <- which.min(abs(RBD_IgM_neg_data[[c]][n]/RBD_IgM_norm[c] - DEN_RBD_IgM_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_RBD_IgM_neg$y[index], max=DEN_RBD_IgM_neg$y[index]),
                          y = RBD_IgM_neg_data[[c]][n]/RBD_IgM_norm[c],
                          pch=19, cex=neg.point.size, col=RBD_IgM_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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


axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )


###################
## Positive points

for(c in 1:N_cohort_NP_IgM)
{
	for(n in 1:NP_IgM_N_part[c])
	{
		points( x = NP_IgM_data[[c]][n,,2],
                    y = NP_IgM_data[[c]][n,,1]/NP_IgM_norm[c],
			  pch=19, cex=pos.point.size, col=NP_IgM_cols[c])   
	}
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_NP_IgM)
{
	if( NP_IgM_N_neg[c] > 0 )
	{
		DEN_NP_IgM_neg = density( NP_IgM_neg_data[[c]]/NP_IgM_norm[c], cut=0, na.rm=TRUE )

		DEN_NP_IgM_neg$y <- DEN_NP_IgM_neg$y/max(DEN_NP_IgM_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_NP_IgM_neg$y, rev(DEN_NP_IgM_neg$y) ),
			   y=c( DEN_NP_IgM_neg$x, rev(DEN_NP_IgM_neg$x) ), 
      	   	   col=NP_IgM_rgb[c], border=NA)

		for(n in 1:NP_IgM_N_neg[c])
		{
			index <- which.min(abs(NP_IgM_neg_data[[c]][n]/NP_IgM_norm[c] - DEN_NP_IgM_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_NP_IgM_neg$y[index], max=DEN_NP_IgM_neg$y[index]),
                          y = NP_IgM_neg_data[[c]][n]/NP_IgM_norm[c],
                          pch=19, cex=neg.point.size, col=NP_IgM_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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


axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
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
## Positive points

for(c in 1:N_cohort_Spike_IgA)
{
	for(n in 1:Spike_IgA_N_part[c])
	{
		points( x = Spike_IgA_data[[c]][n,,2],
                    y = Spike_IgA_data[[c]][n,,1]/Spike_IgA_norm[c],
			  pch=19, cex=pos.point.size, col=Spike_IgA_cols[c])   
	}
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_Spike_IgA)
{
	if( Spike_IgA_N_neg[c] > 0 )
	{
		DEN_Spike_IgA_neg = density( Spike_IgA_neg_data[[c]]/Spike_IgA_norm[c], cut=0, na.rm=TRUE )

		DEN_Spike_IgA_neg$y <- DEN_Spike_IgA_neg$y/max(DEN_Spike_IgA_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_Spike_IgA_neg$y, rev(DEN_Spike_IgA_neg$y) ),
			   y=c( DEN_Spike_IgA_neg$x, rev(DEN_Spike_IgA_neg$x) ), 
      	   	   col=Spike_IgA_rgb[c], border=NA)

		for(n in 1:Spike_IgA_N_neg[c])
		{
			index <- which.min(abs(Spike_IgA_neg_data[[c]][n]/Spike_IgA_norm[c] - DEN_Spike_IgA_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_Spike_IgA_neg$y[index], max=DEN_Spike_IgA_neg$y[index]),
                          y = Spike_IgA_neg_data[[c]][n]/Spike_IgA_norm[c],
                          pch=19, cex=neg.point.size, col=Spike_IgA_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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


axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )

mtext(side = 1, line = 4, 
cex=lab.size, 
text="months after symptom onset")


###################
## Positive points

for(c in 1:N_cohort_RBD_IgA)
{
	for(n in 1:RBD_IgA_N_part[c])
	{
		points( x = RBD_IgA_data[[c]][n,,2],
                    y = RBD_IgA_data[[c]][n,,1]/RBD_IgA_norm[c],
			  pch=19, cex=pos.point.size, col=RBD_IgA_cols[c])   
	}
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_RBD_IgA)
{
	if( RBD_IgA_N_neg[c] > 0 )
	{
		DEN_RBD_IgA_neg = density( RBD_IgA_neg_data[[c]]/RBD_IgA_norm[c], cut=0, na.rm=TRUE )

		DEN_RBD_IgA_neg$y <- DEN_RBD_IgA_neg$y/max(DEN_RBD_IgA_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_RBD_IgA_neg$y, rev(DEN_RBD_IgA_neg$y) ),
			   y=c( DEN_RBD_IgA_neg$x, rev(DEN_RBD_IgA_neg$x) ), 
      	   	   col=RBD_IgA_rgb[c], border=NA)

		for(n in 1:RBD_IgA_N_neg[c])
		{
			index <- which.min(abs(RBD_IgA_neg_data[[c]][n]/RBD_IgA_norm[c] - DEN_RBD_IgA_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_RBD_IgA_neg$y[index], max=DEN_RBD_IgA_neg$y[index]),
                          y = RBD_IgA_neg_data[[c]][n]/RBD_IgA_norm[c],
                          pch=19, cex=neg.point.size, col=RBD_IgA_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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
      xlim=c(-60,420), ylim=c(0.001,100), log="y",
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


axis(1, at=c(-50, -30, -10), 
        label=c("", "neg", ""), 
        cex.axis=0.8*axis.size )

axis(1, at=c( 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330), 
        label=c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
        cex.axis=0.8*axis.size )

axis(1, at =    360 + (c(12, 46, 80) - 12)*((420-360)/(80-12)) , 
        label = c( "12", "46", "80"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100), 
        label = c(0.001, 0.01, 0.1, 1, 10, 100),
        las=2, cex.axis=axis.size )

mtext(side = 1, line = 4, 
cex=lab.size, 
text="months after symptom onset")


###################
## Positive points

for(c in 1:N_cohort_NP_IgA)
{
	for(n in 1:NP_IgA_N_part[c])
	{
		points( x = NP_IgA_data[[c]][n,,2],
                    y = NP_IgA_data[[c]][n,,1]/NP_IgA_norm[c],
			  pch=19, cex=pos.point.size, col=NP_IgA_cols[c])   
	}
}


###################
## Negative controls

off_set <- -2*ww

for(c in 1:N_cohort_NP_IgA)
{
	if( NP_IgA_N_neg[c] > 0 )
	{
		DEN_NP_IgA_neg = density( NP_IgA_neg_data[[c]]/NP_IgA_norm[c], cut=0, na.rm=TRUE )

		DEN_NP_IgA_neg$y <- DEN_NP_IgA_neg$y/max(DEN_NP_IgA_neg$y)

		polygon( x = -30 + off_set + ww*c( -DEN_NP_IgA_neg$y, rev(DEN_NP_IgA_neg$y) ),
			   y=c( DEN_NP_IgA_neg$x, rev(DEN_NP_IgA_neg$x) ), 
      	   	   col=NP_IgA_rgb[c], border=NA)

		for(n in 1:NP_IgA_N_neg[c])
		{
			index <- which.min(abs(NP_IgA_neg_data[[c]][n]/NP_IgA_norm[c] - DEN_NP_IgA_neg$x))

			points( x = -30 + off_set + ww*runif(1, min=-DEN_NP_IgA_neg$y[index], max=DEN_NP_IgA_neg$y[index]),
                          y = NP_IgA_neg_data[[c]][n]/NP_IgA_norm[c],
                          pch=19, cex=neg.point.size, col=NP_IgA_cols[c])
		}
	}

	off_set <- off_set - 2*ww
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









