library(binom)

AB_data <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Manuscript/Database/IDEA_db.csv")


AB_data <- AB_data[ which(AB_data$covid_status %in% c("flow-positive", "pcr-positive", "pre-epidemic-negative")),]
AB_data <- AB_data[,-(38:39)]


IMM_db <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Manuscript/Database/IMM_flw_db.csv")

IMM_db_part <- IMM_db$part_ID[ which( IMM_db$site == "IMM_r2" ) ]

IMM_db <- IMM_db[which( IMM_db[,2] %in% IMM_db_part ),]

IMM_db <- IMM_db[,c(1,2,3,4,5,6,
                    20,19,21,22, 23:49)]


AB_data <- rbind( AB_data, IMM_db )


AB_data$status <- rep( "pos", nrow(AB_data) )
AB_data$status[ which(AB_data$covid_status=="pre-epidemic-negative") ] <- "neg"



AB_data$colour <- rep("dodgerblue", nrow(AB_data))

AB_data$colour[ which(AB_data$site %in% c("bichat", "cochin", "crepy1", "crepy2", "strasbourg")) ] <- "gold" 

AB_data$colour[ which(AB_data$site %in% c("IMM_r1", "IMM_r2")) ] <- "slateblue4" 

AB_data$colour[ which(AB_data$site == "Trinity") ] <- "forestgreen"

AB_data$colour[ which(AB_data$site == "Strasbourg2") ] <- "tan1"




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

line_seq_x <- ((-1):12)*(365/12)

line_seq_y <- c(3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 0.001, 0.003, 0.01, 0.03)

par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))

pos.point.size = 0.5
neg.point.size = 0.2
lab.size   = 1.5
axis.size  = 1.25	
main.size  = 1.5
line.size  = 2
dash.line.size = 0.5
dash.line.size2 = 1

ww <- 10


min_AB <- min(AB_data$spike_IgG)
max_AB <- 0.02


AB_data$S2_IgG[which(AB_data$S2_IgG > 0.02)] <- 0.02
AB_data$S2_IgM[which(AB_data$S2_IgM > 0.02)] <- 0.02
AB_data$S2_IgA[which(AB_data$S2_IgA > 0.02)] <- 0.02

AB_data$ME_IgG[which(AB_data$ME_IgG > 0.02)] <- 0.02
AB_data$ME_IgM[which(AB_data$ME_IgM > 0.02)] <- 0.02
AB_data$ME_IgA[which(AB_data$ME_IgA > 0.02)] <- 0.02

AB_data$OC43_IgG[which(AB_data$OC43_IgG > 0.02)] <- 0.02
AB_data$OC43_IgM[which(AB_data$OC43_IgM > 0.02)] <- 0.02
AB_data$OC43_IgA[which(AB_data$OC43_IgA > 0.02)] <- 0.02


AB_data$X229E_IgG[which(AB_data$X229E_IgG > 0.02)] <- 0.02
AB_data$X229E_IgM[which(AB_data$X229E_IgM > 0.02)] <- 0.02
AB_data$X229E_IgA[which(AB_data$X229E_IgA > 0.02)] <- 0.02

AB_data$HKU1_IgG[which(AB_data$HKU1_IgG > 0.02)] <- 0.02
AB_data$HKU1_IgM[which(AB_data$HKU1_IgM > 0.02)] <- 0.02
AB_data$HKU1_IgA[which(AB_data$HKU1_IgA > 0.02)] <- 0.02

AB_data$NL63_IgG[which(AB_data$NL63_IgG > 0.02)] <- 0.02
AB_data$NL63_IgM[which(AB_data$NL63_IgM > 0.02)] <- 0.02
AB_data$NL63_IgA[which(AB_data$NL63_IgA > 0.02)] <- 0.02




tiff(file="SupFig2_Ab_data_S2_ME_OC43.tif", width=40, height=22, units="cm", res=500)


lay.mat <- rbind( c( 1, 2, 3 ), 
                  c( 4, 5, 6 ),
                  c( 7, 8, 9 ), 
                  c(10,11,12 ), 
                  c(13,13,13 ) )
layout(lay.mat, heights=c(1.5,10,10,11,1.5), widths=c(11,10,10))
layout.show(13)


############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))

plot.new()
title( "Spike subunit 2", 
        cex.main=3.5, line=-2.5)

plot.new()
title( "Membrane-Envelope", 
        cex.main=3.5, line=-2.5)

plot.new()
title( "OC43 Spike", 
        cex.main=3.5, line=-2.5)


###################
###################
###             ###
###  PANEL 1    ###
###  S2 IgG     ###
###             ###
###################
###################

par(mar=c(3,8,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )


mtext(side = 2, line = 4, 
cex=0.7*lab.size, 
text="IgG antibody level (RAU)")



###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$S2_IgG,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

Spike_IgG_N_neg <- AB_data$S2_IgG[ which( AB_data$status == "neg" ) ]

DEN_Spike_IgG_neg = density( log(Spike_IgG_N_neg), cut=0, na.rm=TRUE )

DEN_Spike_IgG_neg$y <- DEN_Spike_IgG_neg$y/max(DEN_Spike_IgG_neg$y)

polygon( x = -20 +  ww*c( -DEN_Spike_IgG_neg$y, rev(DEN_Spike_IgG_neg$y) ),
	   y=c( exp(DEN_Spike_IgG_neg$x), rev(exp(DEN_Spike_IgG_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(Spike_IgG_N_neg))
{
	index <- which.min(abs(log(Spike_IgG_N_neg[n]) - DEN_Spike_IgG_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_Spike_IgG_neg$y[index], max=DEN_Spike_IgG_neg$y[index]),
                  y = Spike_IgG_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
}







###################
###################
###             ###
###  PANEL 2    ###
###  ME IgG     ###
###             ###
###################
###################

par(mar=c(3,3,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )



###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$ME_IgG,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

RBD_IgG_N_neg <- AB_data$ME_IgG[ which( AB_data$status == "neg" ) ]

DEN_RBD_IgG_neg = density( log(RBD_IgG_N_neg), cut=0, na.rm=TRUE )

DEN_RBD_IgG_neg$y <- DEN_RBD_IgG_neg$y/max(DEN_RBD_IgG_neg$y)

polygon( x = -20 +  ww*c( -DEN_RBD_IgG_neg$y, rev(DEN_RBD_IgG_neg$y) ),
	   y=c( exp(DEN_RBD_IgG_neg$x), rev(exp(DEN_RBD_IgG_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(RBD_IgG_N_neg))
{
	index <- which.min(abs(log(RBD_IgG_N_neg[n]) - DEN_RBD_IgG_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_RBD_IgG_neg$y[index], max=DEN_RBD_IgG_neg$y[index]),
                  y = RBD_IgG_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
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
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )



###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$OC43_IgG,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

NP_IgG_N_neg <- AB_data$OC43_IgG[ which( AB_data$status == "neg" ) ]

DEN_NP_IgG_neg = density( log(NP_IgG_N_neg), cut=0, na.rm=TRUE )

DEN_NP_IgG_neg$y <- DEN_NP_IgG_neg$y/max(DEN_NP_IgG_neg$y)

polygon( x = -20 +  ww*c( -DEN_NP_IgG_neg$y, rev(DEN_NP_IgG_neg$y) ),
	   y=c( exp(DEN_NP_IgG_neg$x), rev(exp(DEN_NP_IgG_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(NP_IgG_N_neg))
{
	index <- which.min(abs(log(NP_IgG_N_neg[n]) - DEN_NP_IgG_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_NP_IgG_neg$y[index], max=DEN_NP_IgG_neg$y[index]),
                  y = NP_IgG_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
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
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )


mtext(side = 2, line = 4, 
cex=0.7*lab.size, 
text="IgM antibody level (RAU)")



###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$S2_IgM,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

Spike_IgM_N_neg <- AB_data$S2_IgM[ which( AB_data$status == "neg" ) ]
Spike_IgM_N_neg <- Spike_IgM_N_neg[which(is.na(Spike_IgM_N_neg)==FALSE)]


DEN_Spike_IgM_neg = density( log(Spike_IgM_N_neg), cut=0, na.rm=TRUE )

DEN_Spike_IgM_neg$y <- DEN_Spike_IgM_neg$y/max(DEN_Spike_IgM_neg$y)

polygon( x = -20 +  ww*c( -DEN_Spike_IgM_neg$y, rev(DEN_Spike_IgM_neg$y) ),
	   y=c( exp(DEN_Spike_IgM_neg$x), rev(exp(DEN_Spike_IgM_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(Spike_IgM_N_neg))
{
	index <- which.min(abs(log(Spike_IgM_N_neg[n]) - DEN_Spike_IgM_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_Spike_IgM_neg$y[index], max=DEN_Spike_IgM_neg$y[index]),
                  y = Spike_IgM_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
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
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )





###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$ME_IgM,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

RBD_IgM_N_neg <- AB_data$ME_IgM[ which( AB_data$status == "neg" ) ]
RBD_IgM_N_neg <- RBD_IgM_N_neg[which(is.na(RBD_IgM_N_neg)==FALSE)]

DEN_RBD_IgM_neg = density( log(RBD_IgM_N_neg), cut=0, na.rm=TRUE )

DEN_RBD_IgM_neg$y <- DEN_RBD_IgM_neg$y/max(DEN_RBD_IgM_neg$y)

polygon( x = -20 +  ww*c( -DEN_RBD_IgM_neg$y, rev(DEN_RBD_IgM_neg$y) ),
	   y=c( exp(DEN_RBD_IgM_neg$x), rev(exp(DEN_RBD_IgM_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(RBD_IgM_N_neg))
{
	index <- which.min(abs(log(RBD_IgM_N_neg[n]) - DEN_RBD_IgM_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_RBD_IgM_neg$y[index], max=DEN_RBD_IgM_neg$y[index]),
                  y = RBD_IgM_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
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
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )




###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$OC43_IgM,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

NP_IgM_N_neg <- AB_data$OC43_IgM[ which( AB_data$status == "neg" ) ]
NP_IgM_N_neg <- NP_IgM_N_neg[which(is.na(NP_IgM_N_neg)==FALSE)]

DEN_NP_IgM_neg = density( log(NP_IgM_N_neg), cut=0, na.rm=TRUE )

DEN_NP_IgM_neg$y <- DEN_NP_IgM_neg$y/max(DEN_NP_IgM_neg$y)

polygon( x = -20 +  ww*c( -DEN_NP_IgM_neg$y, rev(DEN_NP_IgM_neg$y) ),
	   y=c( exp(DEN_NP_IgM_neg$x), rev(exp(DEN_NP_IgM_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(NP_IgM_N_neg))
{
	index <- which.min(abs(log(NP_IgM_N_neg[n]) - DEN_NP_IgM_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_NP_IgM_neg$y[index], max=DEN_NP_IgM_neg$y[index]),
                  y = NP_IgM_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
}



###################
###################
###             ###
###  PANEL 7    ###
###  Spike IgA  ###
###             ###
###################
###################

par(mar=c(6,8,1,1))
par(mgp=c(2.5, 1, 0))

###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )


mtext(side = 2, line = 4, 
cex=0.7*lab.size, 
text="IgA antibody level (RAU)")

mtext(side = 1, line = 3, 
cex=lab.size, 
text="time post symptoms (months)")

###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$S2_IgA,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

Spike_IgA_N_neg <- AB_data$S2_IgA[ which( AB_data$status == "neg" ) ]
Spike_IgA_N_neg <- Spike_IgA_N_neg[which(is.na(Spike_IgA_N_neg)==FALSE)]

DEN_Spike_IgA_neg = density( log(Spike_IgA_N_neg), cut=0, na.rm=TRUE )

DEN_Spike_IgA_neg$y <- DEN_Spike_IgA_neg$y/max(DEN_Spike_IgA_neg$y)

polygon( x = -20 +  ww*c( -DEN_Spike_IgA_neg$y, rev(DEN_Spike_IgA_neg$y) ),
	   y=c( exp(DEN_Spike_IgA_neg$x), rev(exp(DEN_Spike_IgA_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(Spike_IgA_N_neg))
{
	index <- which.min(abs(log(Spike_IgA_N_neg[n]) - DEN_Spike_IgA_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_Spike_IgA_neg$y[index], max=DEN_Spike_IgA_neg$y[index]),
                  y = Spike_IgA_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
}



###################
###################
###             ###
###  PANEL 8    ###
###  RBD IgA    ###
###             ###
###################
###################

par(mar=c(6,3,1,1))
par(mgp=c(2.5, 1, 0))


###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )



mtext(side = 1, line = 3, 
cex=lab.size, 
text="time post symptoms (months)")

###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$ME_IgA,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

RBD_IgA_N_neg <- AB_data$ME_IgA[ which( AB_data$status == "neg" ) ]
RBD_IgA_N_neg <- RBD_IgA_N_neg[which(is.na(RBD_IgA_N_neg)==FALSE)]

DEN_RBD_IgA_neg = density( log(RBD_IgA_N_neg), cut=0, na.rm=TRUE )

DEN_RBD_IgA_neg$y <- DEN_RBD_IgA_neg$y/max(DEN_RBD_IgA_neg$y)

polygon( x = -20 +  ww*c( -DEN_RBD_IgA_neg$y, rev(DEN_RBD_IgA_neg$y) ),
	   y=c( exp(DEN_RBD_IgA_neg$x), rev(exp(DEN_RBD_IgA_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(RBD_IgA_N_neg))
{
	index <- which.min(abs(log(RBD_IgA_N_neg[n]) - DEN_RBD_IgA_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_RBD_IgA_neg$y[index], max=DEN_RBD_IgA_neg$y[index]),
                  y = RBD_IgA_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
}






###################
###################
###             ###
###  PANEL 9    ###
###  NP IgA     ###
###             ###
###################
###################

par(mar=c(6,3,1,1))
par(mgp=c(2.5, 1, 0))



###################
## Frame and axes

plot( x=1000, y=1000,
      xlim=c(-31,365), ylim=c(9e-6,0.03), log="y",
	yaxt='n', xaxt='n', bty='n',
	ylab="", xlab="",
	main="",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,0.03), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}

points( x=c(-1000,365), y=rep(min_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")

points( x=c(-1000,365), y=rep(max_AB,2), 
        type='l', lwd=dash.line.size2, col="black", lty="dashed")



axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*(365/12), 
        label=c( "", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
        cex.axis=axis.size )


axis(2, at =    c(1e-5, 1e-4, 0.001, 0.01), 
        label = c(1e-5, 1e-4, 0.001, 0.01),
        las=2, cex.axis=axis.size, tck=-0.03 )


axTicks <-  c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) )
axis(2, at    =c( 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:3) ), 
        label = rep("", length(axTicks)),
        las=2, cex.axis=axis.size, tck=-0.015 )


mtext(side = 1, line = 3, 
cex=lab.size, 
text="time post symptoms (months)")

###################
## Positive points

points( x = AB_data$days_pso, y = AB_data$OC43_IgA,
	  pch=19, cex=pos.point.size, 
        col=AB_data$colour )



###################
## Negative controls

NP_IgA_N_neg <- AB_data$OC43_IgA[ which( AB_data$status == "neg" ) ]
NP_IgA_N_neg <- NP_IgA_N_neg[which(is.na(NP_IgA_N_neg)==FALSE)]

DEN_NP_IgA_neg = density( log(NP_IgA_N_neg), cut=0, na.rm=TRUE )

DEN_NP_IgA_neg$y <- DEN_NP_IgA_neg$y/max(DEN_NP_IgA_neg$y)	

polygon( x = -20 +  ww*c( -DEN_NP_IgA_neg$y, rev(DEN_NP_IgA_neg$y) ),
	   y=c( exp(DEN_NP_IgA_neg$x), rev(exp(DEN_NP_IgA_neg$x)) ), 
         col=rgb(30/255, 144/255, 255/255, 0.5), border=NA)

for(n in 1:length(NP_IgA_N_neg))
{
	index <- which.min(abs(log(NP_IgA_N_neg[n]) - DEN_NP_IgA_neg$x))

	points( x = -20 + ww*runif(1, min=-DEN_NP_IgA_neg$y[index], max=DEN_NP_IgA_neg$y[index]),
                  y = NP_IgA_N_neg[n],
                  pch=19, cex=neg.point.size, col="dodgerblue")
}




###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar=c(0,0,0,0))

plot.new()

legend(x='center', 
       legend = c("negative panel", "positive panel", "Dublin", "Strasbourg", "Paris"), 
       col = c("dodgerblue", "gold", "forestgreen", "tan1", "slateblue4"), 
       pch = rep(19,5),
       ncol=5, cex=2.5, bty="n" )


dev.off()







