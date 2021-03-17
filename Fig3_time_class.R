library(binom)

library(MASS)
library(ROCR)
library(randomForest)
library(pROC)

set.seed(1212)

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


##################################
## Read in IDEA database

AB_data <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Manuscript/Database/IDEA_db.csv")



##################################
## Exclude unknown and post-2020 negative

AB_data <- AB_data[ which(AB_data$covid_status %in% c("flow-positive", "pcr-positive", "pre-epidemic-negative")),]
AB_data <- AB_data[,-(38:39)]



drop_NA <- which( AB_data$site %in% c("bichat", "cochin", "crepy1", "crepy2",
                                  "strasbourg", "Strasbourg2", "Trinity", "IMM_r1", "IMM_r2") &
                  is.na(AB_data$days_pso) == TRUE )

AB_data <- AB_data[-drop_NA,]


drop_early <- which( AB_data$site %in% c("bichat", "cochin", "crepy1", "crepy2",
                                  "strasbourg", "Strasbourg2", "Trinity", "IMM_r1", "IMM_r2") &
                     AB_data$days_pso <= 14 )

AB_data <- AB_data[-drop_early,]


##################################
## Read in IMM database

IMM_db <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Manuscript/Database/IMM_flw_db.csv")


IMM_db_part <- IMM_db$part_ID[ which( IMM_db$site == "IMM_r2" ) ]

IMM_db <- IMM_db[which( IMM_db[,2] %in% IMM_db_part ),]

IMM_db <- IMM_db[,c(1,2,3,4,5,6,
                    20,19,21,22, 23:49)]

IMM_db$covid_status <- "IMM_pos"


IMM_db$days_pso[which( is.na(IMM_db$days_pso)==TRUE )] <- median( IMM_db$days_pso[which(IMM_db$site=="IMM_r1")], na.rm=TRUE )



##################################
## Merge database

AB_data <- rbind( AB_data, IMM_db )

drop_AB_missing <- unique(which(is.na(AB_data[,c(11,12,13,20,21,22,29,30,31)])==TRUE, arr.ind=TRUE)[,1])


AB_data <- AB_data[-drop_AB_missing,]


##################################
## Add in binary infection category

AB_data$bin_status <- rep( "pos", nrow(AB_data) )
AB_data$bin_status[ which(AB_data$covid_status == "pre-epidemic-negative") ] <- "neg"


##################################
## Add in time infection category

AB_data$status <- rep( "neg", nrow(AB_data) )
AB_data$status[ which(AB_data$days_pso <= 90) ] <- "pos1"
AB_data$status[ which(AB_data$days_pso > 90 & AB_data$days_pso <= 180 ) ] <- "pos2"
AB_data$status[ which(AB_data$days_pso > 180 ) ] <- "pos3"



neg_col <- "dodgerblue"
pos1_col <- "orangered"
pos2_col <- "forestgreen"
pos3_col <- "slateblue4"


AB_data$colour <- rep(neg_col, nrow(AB_data))

AB_data$colour[ which(AB_data$status == "pos1") ] <- pos1_col

AB_data$colour[ which(AB_data$status == "pos2") ] <- pos2_col

AB_data$colour[ which(AB_data$status == "pos3") ] <- pos3_col

AB_data$status <- as.factor(AB_data$status)


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

RF_4way_IgGMA = randomForest( status ~ spike_IgG + RBD_IgG + NP_IgG +
                                       spike_IgM + RBD_IgM + NP_IgM +
                                       spike_IgA + RBD_IgA + NP_IgA,  
                              data = AB_data, importance = TRUE, ntree = 10000)
 

#RF_4way_IgGMA = randomForest( status ~ spike_IgG + RBD_IgG + NP_IgG + ME_IgG,  
#                              data = AB_data, importance = TRUE, ntree = 10000)
 


#########################################
## RF 4 way prediction object: IgG only

ROC_4way_IgGMA <- multiclass.roc( AB_data$status, RF_4way_IgGMA$votes )

RF_4way_IgGMA_pred <- ROC_4way_IgGMA$predictor

RF_4way_IgGMA_pred <- RF_4way_IgGMA_pred[,c( which(colnames(RF_4way_IgGMA_pred) == "neg"),
                                             which(colnames(RF_4way_IgGMA_pred) == "pos1"),                                         
						         which(colnames(RF_4way_IgGMA_pred) == "pos2"),
                                             which(colnames(RF_4way_IgGMA_pred) == "pos3") )]

RF_4way_IgGMA_pred <- as.data.frame(RF_4way_IgGMA_pred)

RF_4way_IgGMA_pred <- cbind( AB_data$status, RF_4way_IgGMA_pred )
colnames(RF_4way_IgGMA_pred)[1] <- "status"



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


index_pos <- which( RF_4way_IgGMA_pred$status %in% c("pos1", "pos2", "pos3") )
index_neg <- which( RF_4way_IgGMA_pred$status == "neg" )

N_pos1_val <- length(which( AB_data$status %in% c("pos1") ))
N_pos2_val <- length(which( AB_data$status %in% c("pos2") ))
N_pos3_val <- length(which( AB_data$status %in% c("pos3") ))
N_pos_val <- length(which( AB_data$status %in% c("pos1", "pos2", "pos3") ))
N_neg_val <- length(which( AB_data$status %in% c("neg") ))



for(i in 1:N_cut_neg_pos)
{
	## sensitivity (identifying positive)
	ROC1_neg_pos[i,2] <- length(which( RF_4way_IgGMA_pred[index_pos,2] <= SS_cut_neg_pos[i] ))/(N_pos_val) 	
 
	## specificity (identifying negative)
	ROC1_neg_pos[i,3] <- length(which( RF_4way_IgGMA_pred[index_neg,2] > SS_cut_neg_pos[i] ))/N_neg_val 	
}


ROC1_neg_pos <- rbind( c(0,0,1), ROC1_neg_pos )

#########################################
## Select cut off for classification of positive individuals
## based on IgG only algorithm

pos_target_spec <- 0.99

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

CMAT <- table( RF_4way_IgGMA_pred$status, RF_4way_IgGMA_pred$pred )
rownames(CMAT) <- c("cat_neg", "cat_pos1", "cat_pos2", "cat_pos3")
colnames(CMAT) <- c("class_neg", "class_pos1", "class_pos2", "class_pos3")

CMAT[1,] <- CMAT[1,]/N_neg_val
CMAT[2,] <- CMAT[2,]/N_pos1_val
CMAT[3,] <- CMAT[3,]/N_pos2_val
CMAT[4,] <- CMAT[4,]/N_pos3_val





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

RF_4way_IgGMA_pred$days_pso <- AB_data$days_pso


t_bins <- seq(from=0, to=370, by=10)

N_bins <- length(t_bins)




neg_bins <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	neg_bins[i] <- length(which( RF_4way_IgGMA_pred$pred == "neg" &
                                   RF_4way_IgGMA_pred$days_pso >= t_bins[i] & 
                                   RF_4way_IgGMA_pred$days_pso < t_bins[i+1] ))
}


pos1_bins <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	pos1_bins[i] <- length(which( RF_4way_IgGMA_pred$pred == "pos1" &
                                    RF_4way_IgGMA_pred$days_pso >= t_bins[i] & 
                                    RF_4way_IgGMA_pred$days_pso < t_bins[i+1] ))
}



pos2_bins <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	pos2_bins[i] <- length(which( RF_4way_IgGMA_pred$pred == "pos2" &
                                    RF_4way_IgGMA_pred$days_pso >= t_bins[i] & 
                                    RF_4way_IgGMA_pred$days_pso < t_bins[i+1] ))
}



pos3_bins <- rep(NA, N_bins)

for(i in 1:N_bins)
{
	pos3_bins[i] <- length(which( RF_4way_IgGMA_pred$pred == "pos3" &
                                    RF_4way_IgGMA_pred$days_pso >= t_bins[i] & 
                                    RF_4way_IgGMA_pred$days_pso < t_bins[i+1] ))
}




neg_pred_neg <- length(which(RF_4way_IgGMA_pred$pred[which( RF_4way_IgGMA_pred$status == "neg" )] == "neg"))

neg_pred_pos1 <- length(which(RF_4way_IgGMA_pred$pred[which( RF_4way_IgGMA_pred$status == "neg" )] == "pos1"))

neg_pred_pos2 <- length(which(RF_4way_IgGMA_pred$pred[which( RF_4way_IgGMA_pred$status == "neg" )] == "pos2"))

neg_pred_pos3 <- length(which(RF_4way_IgGMA_pred$pred[which( RF_4way_IgGMA_pred$status == "neg" )] == "pos3"))









lab.size   = 1.5
axis.size  = 1.25	
main.size  = 1.5
main.size.top  = 1.75
line.size  = 2
dash.line.size <- 1



tiff(file="Fig3_time_class.tif", width=30, height=18, units="cm", res=500)

lay.mat <- rbind( c( 1, 1, 2, 2, 3, 3 ),
                  c( 4, 4, 4, 4, 4, 4 ),
                  c( 5, 6, 6, 6, 6, 6 ),
                  c( 7, 7, 7, 7, 7, 7 ) )
layout(lay.mat, heights=c(20,2,30,2), widths=c(1,1,1))
layout.show(7)



########################
########################
##                    ##              
##  PANEL 1           ##            
##  Infection ROC     ##
##                    ##
########################
########################


par(mar=c(4,6,2,1))
par(mgp=c(2.5, 0.75,0))


line_seq_x <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
line_seq_y <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)



plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(A) Negative or positive?",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size.top)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}


mtext(side = 1, line = 3, 
cex=1.0, 
text="1 - specificity")

mtext(side = 2, line = 4, 
cex=1.0, 
text="sensitivity")




axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=1.5) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=1.5 ) 

points( x = 1 - ROC1_neg_pos[,3], y = ROC1_neg_pos[,2], 
		  type='S', lwd=2, col=neg_col )


#points( x = 1 - pos_target_spec, y = pos_target_sens,
#        pch=19, col="red")




########################
########################
##                    ##              
##  PANEL 2           ##            
##  3 month ROC       ##
##                    ##
########################
########################


par(mar=c(4,6,2,1))
par(mgp=c(2.5, 1.0,0))


line_seq_x <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
line_seq_y <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)



plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(B) Infected <3 months ago?",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size.top)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}


mtext(side = 1, line = 2.5, 
cex=1.0, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
cex=1.0, 
text="sensitivity")




axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=1.5) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=1.5 ) 

points( x = 1 - ROC_4way_IgGMA$rocs$`pos1/pos2`[[1]]$specificities, 
        y = ROC_4way_IgGMA$rocs$`pos1/pos2`[[1]]$sensitivities, 
		  type='S', lwd=2, col=pos2_col )


points( x = 1 - ROC_4way_IgGMA$rocs$`pos1/pos3`[[1]]$specificities, 
      y = ROC_4way_IgGMA$rocs$`pos1/pos3`[[1]]$sensitivities,
      type='S', lwd=2, col=pos3_col)

points( x = c(0,1), y = c(0,1),
        type='l', lty="dashed")

legend( x = "bottomright",
legend = c("<3 mths vs 3-6 mths", "<3 mths vs 6-12 mths"),
fill = c(pos2_col, pos3_col),
border = c(pos2_col, pos3_col), 
bty='n', cex=1.25 )




########################
########################
##                    ##              
##  PANEL 3           ##            
##  6 month ROC       ##
##                    ##
########################
########################


par(mar=c(4,6,2,1))
par(mgp=c(2.5, 1.0,0))


line_seq_x <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
line_seq_y <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)



plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(C) Infected <6 months ago?",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size.top)

for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}
for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(1e-10,1), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}


mtext(side = 1, line = 2.5, 
cex=1.0, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
cex=1.0, 
text="sensitivity")




axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=1.5) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=1.5 ) 

points( x = 1 - ROC_4way_IgGMA$rocs$`pos2/pos3`[[1]]$specificities, 
        y = ROC_4way_IgGMA$rocs$`pos2/pos3`[[1]]$sensitivities, 
		  type='S', lwd=2, col=pos3_col )


legend( x = "bottomright",
legend = c("3-6 mths vs 6-12 mths"),
fill = c( pos3_col),
border = c( pos3_col), 
bty='n', cex=1.25 )




############################
## Labels in middle       ## 
############################

par(mar = c(0,0,0,0))

plot.new()
title( "(D) Estimation of time since infection", 
        cex.main=2.5, line=-1.75)

####################
####################
###              ###
###  PANEL 4     ###
###  Scenario 1  ###
###  Negatives   ###
###              ###
####################
####################

par(mar = c(5,5,0.5,0.5))

line_seq_y <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)

plot(x=1000, y=1000, 
xlim=c(0,100), ylim=c(0,400),
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
         y = c(0,0,neg_pred_neg,neg_pred_neg),
         col=neg_col, border=NA )


polygon( x = c(10,90,90,10),
         y = neg_pred_neg + c(0,0,neg_pred_pos1,neg_pred_pos1),
         col=pos1_col, border=NA )

polygon( x = c(10,90,90,10),
         y = neg_pred_neg + neg_pred_pos1 + c(0,0,neg_pred_pos2,neg_pred_pos2),
         col=pos2_col, border=NA )

polygon( x = c(10,90,90,10),
         y = neg_pred_neg + neg_pred_pos1 + neg_pred_pos2 + c(0,0,neg_pred_pos3,neg_pred_pos3),
         col=pos3_col, border=NA )





axis(1, at = c(0, 100), 
        label=c("", ""), 
        cex.axis=axis.size )


axis(2, at = c(0, 50, 100, 150, 200, 250, 300, 350, 400), 
        label = c("0", "50", "100", "150", "200", "250", "300", "350", "400"),
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
###  PANEL 5     ###
###  Scenario 1  ###
###  Positives   ###
###              ###
####################
####################


par(mar = c(5,5,0.5,2))
par(mgp=c(3,1.25,0))


line_seq_x <- -30*c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
line_seq_y <- c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200)

plot(x=1000, y=1000, 
xlim=c(-360,30), ylim=c(0,200),
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
	polygon( x = -t_bins[c(i,i+1,i+1,i)], 
               y = c(0,0,neg_bins[i],neg_bins[i]),
               col=neg_col, border=NA )
}

for(i in 1:N_bins)
{
	polygon( x = -t_bins[c(i,i+1,i+1,i)], 
               y = neg_bins[i] + c(0,0,pos1_bins[i],pos1_bins[i]),
               col=pos1_col, border=NA )
}

for(i in 1:N_bins)
{
	polygon( x = -t_bins[c(i,i+1,i+1,i)], 
               y = neg_bins[i] + pos1_bins[i] + c(0,0,pos2_bins[i],pos2_bins[i]),
               col=pos2_col, border=NA )
}


for(i in 1:N_bins)
{
	polygon( x = -t_bins[c(i,i+1,i+1,i)], 
               y = neg_bins[i] + pos1_bins[i] + pos2_bins[i] + c(0,0,pos3_bins[i],pos3_bins[i]),
               col=pos3_col, border=NA )
}

axis(1, at = -30*c(-1, 0, 3, 6, 9, 12), 
        label=c("", "0", "3", "6", "9", "12"), 
        cex.axis=1.5*axis.size )


axis(2, at = c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200), 
        label = c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200),
        las=2, cex.axis=axis.size )


mtext(side = 2, line = 3, 
cex=lab.size, 
text="count")

mtext(side = 1, line = 2.75, 
cex=lab.size, 
text="positives: time before sampling (months)")

arrows( x0=0, y0=135, x1=0, y1=115,
        angle=30, lwd=2, length=0.1)

text( x=0, y=150, 
      labels="sampling", 
      cex=2)


#points(x=c(-360,-360), y=c(0,1000), 
#       type='l', lwd=2, col="black", lty="dashed")
#
#points(x=c(-180,-180), y=c(0,1000), 
#       type='l', lwd=2, col="black", lty="dashed")
#
#points(x=c(-90,-90), y=c(0,1000), 
#       type='l', lwd=2, col="black", lty="dashed")
#
#points(x=c(-0,-0), y=c(0,1000), 
#       type='l', lwd=2, col="black", lty="dashed")


###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar=c(0,0,0,0))

plot.new()

legend(x='center', 
       legend = c("negative controls", "0-3 months", "3-6 months", "6-12 months"), 
       fill = c(neg_col, pos1_col, pos2_col, pos3_col), 
       border = c(neg_col, pos1_col, pos2_col, pos3_col), 
       ncol=4, cex=2, bty="n" )



dev.off()




