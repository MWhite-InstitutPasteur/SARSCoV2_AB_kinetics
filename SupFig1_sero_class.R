library(binom)

library(MASS)
library(ROCR)
library(randomForest)
library(pROC)

set.seed(1234)


N_tree <- 10000

N_rep <- 50

train_prop <- 0.666

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

###################################
###################################
##                                
##  VALIDATION DATA               


VAL_data <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Manuscript/Database/IDEA_db.csv")


##################################
## Exclude unknown and post-2020 negative

VAL_data <- VAL_data[ which(VAL_data$covid_status %in% c("flow-positive", "pcr-positive", "pre-epidemic-negative")),]
VAL_data <- VAL_data[,-(38:39)]


VAL_data$status <- rep( "pos", nrow(VAL_data) )
VAL_data$status[ which(VAL_data$covid_status=="pre-epidemic-negative") ] <- "neg"



##################################
## Exclude individuals with missing days since symptom
## and with days <= 14

VAL_data <- VAL_data[-which( VAL_data$status == "pos" & is.na(VAL_data$days_pso)==TRUE ),]

VAL_data <- VAL_data[-which( VAL_data$status == "pos" & VAL_data$days_pso<=14 ),]

VAL_data <- VAL_data[,c(11:13,20:22,29:31,38)]

VAL_data <- VAL_data[-unique(which( is.na(VAL_data)==TRUE, arr.ind=TRUE )[,1]),]

VAL_data$status <- as.factor(VAL_data$status)


###################################
###################################
##                                
##  IMM DATA                      

IMM_data <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Manuscript/Database/IMM_db.csv")

IMM_data <- IMM_data[which(IMM_data$site == "IMM_r1"),]

IMM_data <- IMM_data[,c(23:25,32:34,41:43)]

IMM_data <- IMM_data[-unique(which( is.na(IMM_data)==TRUE, arr.ind=TRUE )[,1]),]


IMM_xVAL_data <- IMM_data





#########################################
#########################################
##                                     ##
##   ####  ##     ####   ####   ####   ##
##  ##  ## ##    ##  ## ##     ##      ##
##  ##     ##    ######  ####   ####   ##
##  ##  ## ##    ##  ##     ##     ##  ##
##   ####  ##### ##  ##  ####   ####   ##
##                                     ##
#########################################
#########################################

###################################
###################################
##                                
##  RANDOM FORESTS                
##  Based on all validation data  

####################################
## Train RF algorithm on validation data

RF_VAL = randomForest( status ~ ., data=VAL_data, 
                                      importance=TRUE, ntree=N_tree)


####################################
## Variable importance plot

varImpPlot( RF_VAL )


####################################
## Evaluate classification performance

RF_VAL_roc <- roc(VAL_data$status, RF_VAL$votes[,2])


RF_VAL_roc$specificities[min(which(RF_VAL_roc$specificities > 0.99))]

RF_VAL_roc$sensitivities[min(which(RF_VAL_roc$specificities > 0.99))]

RF_VAL_roc$threshold[min(which(RF_VAL_roc$specificities > 0.99))]

RF_VAL_roc$auc


####################################
## Apply RF validation algorithm to IMM data

RF_IMM_predict = predict( RF_VAL, newdata=IMM_data, predict.all=TRUE)

RF_IMM_score = rowSums(RF_IMM_predict$individual=="pos")/N_tree

IMM_data$RF_IMM_score <- RF_IMM_score


####################################
## Calculate sero-prevalence

RF_VAL_cut <- RF_VAL_roc$threshold

N_cut <- length(RF_VAL_cut)

RF_VAL_cut[1] <- 0
RF_VAL_cut[N_cut] <- 1

RF_IMM_meas_prev <- rep(NA, N_cut)

for(i in 1:N_cut)
{
	RF_IMM_meas_prev[i] <- length(which(IMM_data$RF_IMM_score >= RF_VAL_cut[i]))/nrow(IMM_data)
}

RF_IMM_adj_prev <- rep(NA, N_cut)

for(i in 1:N_cut)
{
	RF_IMM_adj_prev[i] <- ( RF_IMM_meas_prev[i] + RF_VAL_roc$specificities[i] - 1 )/
                            ( RF_VAL_roc$sensitivities[i] + RF_VAL_roc$specificities[i] - 1 )  

}

RF_IMM_adj_prev[which(RF_IMM_adj_prev < 0)] <- 0





plot( x=1-RF_VAL_roc$specificities, y=RF_VAL_roc$sensitivities)




plot( x=1-RF_VAL_roc$specificities, y=RF_IMM_meas_prev)



plot( x=1-RF_VAL_roc$specificities, y=RF_IMM_adj_prev)

####################################
####################################
##                                ##
##  ##   ## ##   ##  ####  ##     ##
##   ## ##  ##   ## ##  ## ##     ##
##    ###    ## ##  ###### ##     ##
##   ## ##    ###   ##  ## ##     ##
##  ##   ##    #    ##  ## #####  ##
##                                ##
####################################
####################################

####################################
####################################
##                                
##  RANDOM FORESTS                
##  Cross-validated               


RF_xVAL_sens <- list()
RF_xVAL_spec <- list()

RF_xIMM_meas_prev <- list()
RF_xIMM_adj_prev <- list()


for(n in 1:N_rep)
{
	####################################
	## Split the validation data into training
	## and testing sets

	index_train <- sample( length(VAL_data$status), train_prop*length(VAL_data$status) )	
	index_test  <- setdiff( 1:length(VAL_data$status), index_train )

	#status_train <- VAL_data$status[index_train]
	#status_test  <- VAL_data$status[index_test]

	VAL_data_train <- VAL_data[index_train,]
	VAL_data_test  <- VAL_data[index_test,]


	####################################
	## Train RF algorithm on validation data

	RF_xVAL = randomForest( status ~ ., data=VAL_data_train, 
      	                  importance=FALSE, ntree=N_tree)

	RF_xVAL_pred_obj = predict( RF_xVAL, newdata=VAL_data_test, predict.all=TRUE)

	RF_xVAL_votes = rowSums( RF_xVAL_pred_obj$individual=="pos")/N_tree

	RF_xVAL_roc <- roc(VAL_data_test$status, RF_xVAL_votes)


	####################################
	## Apply cross-validated RF algorithm to IMM data

	RF_IMM_xpredict = predict( RF_xVAL, newdata=IMM_xVAL_data, predict.all=TRUE)

	RF_IMM_xscore = rowSums(RF_IMM_xpredict$individual=="pos")/N_tree


	####################################
	## Calculate sero-prevalence

	RF_xVAL_sens[[n]] <- RF_xVAL_roc$sensitivities
	RF_xVAL_spec[[n]] <- RF_xVAL_roc$specificities

	RF_xVAL_cut <- RF_xVAL_roc$threshold

	N_cut <- length(RF_xVAL_cut)

	RF_xVAL_cut[1] <- 0
	RF_xVAL_cut[N_cut] <- 1

	RF_xIMM_meas_prev[[n]] <- rep(NA, N_cut)

	for(i in 1:N_cut)
	{
		RF_xIMM_meas_prev[[n]][i] <- length(which(RF_IMM_xscore >= RF_xVAL_cut[i]))/nrow(IMM_data)
	}

	RF_xIMM_adj_prev[[n]] <- rep(NA, N_cut)

	for(i in 1:N_cut)
	{
		RF_xIMM_adj_prev[[n]][i] <- ( RF_xIMM_meas_prev[[n]][i] +  RF_xVAL_roc$specificities[i] - 1 )/
                                        ( RF_xVAL_roc$sensitivities[i] + RF_xVAL_roc$specificities[i] - 1 )  

	}

	RF_xIMM_adj_prev[[n]][which(RF_xIMM_adj_prev[[n]] < 0)] <- 0
}


####################################
####################################
##                                
##  Collate Xval results          

N_SS_cut <- 10000 

SS_seq <- seq(from=0, to=1, length=N_SS_cut)


#########################################
## Cross-validated variation in sensitivity

RF_xVAL_sens_vary <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:(N_SS_cut-1))
	{
		RF_xVAL_sens_vary[i,n] <- RF_xVAL_sens[[n]][ min(which(RF_xVAL_spec[[n]] > SS_seq[i])) ]
	}
}
RF_xVAL_sens_vary[N_SS_cut,] <- 0


RF_xVAL_sens_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
for(i in 1:N_SS_cut)
{
	RF_xVAL_sens_quant[i,] <- quantile( RF_xVAL_sens_vary[i,], prob=c(0.5, 0.025, 0.975) )
}



#########################################
## Cross-validated variation in measured prevalence

RF_xIMM_meas_prev_vary <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:(N_SS_cut-1))
	{
		RF_xIMM_meas_prev_vary[i,n] <- RF_xIMM_meas_prev[[n]][ min(which(RF_xVAL_spec[[n]] > SS_seq[i])) ]
	}
}
RF_xIMM_meas_prev_vary[N_SS_cut,] <- 0


RF_xIMM_meas_prev_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
for(i in 1:N_SS_cut)
{
	RF_xIMM_meas_prev_quant[i,] <- quantile( RF_xIMM_meas_prev_vary[i,], prob=c(0.5, 0.025, 0.975) )
}




#########################################
## Cross-validated variation in adjusted prevalence

RF_xIMM_adj_prev_vary <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:(N_SS_cut-1))
	{
		RF_xIMM_adj_prev_vary[i,n] <- RF_xIMM_adj_prev[[n]][ min(which(RF_xVAL_spec[[n]] > SS_seq[i])) ]
	}
}
RF_xIMM_adj_prev_vary[N_SS_cut,] <- 0


RF_xIMM_adj_prev_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
for(i in 1:N_SS_cut)
{
	RF_xIMM_adj_prev_quant[i,] <- quantile( RF_xIMM_adj_prev_vary[i,], prob=c(0.5, 0.025, 0.975) )
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


tiff(file="SupFig1_sero_class.tif", width=12, height=24, units="cm", res=500)

par(mfrow=c(3,1))


lab.size   = 1.5
axis.size  = 1
axis.size2 = 1.35
main.size  = 1.75

par(mar=c(5,5,4,1))
par(mgp=c(1.5, 0.9, 0))


#####################################
#####################################
##                                 ## 
##  PANEL 1                        ##
##  Cross-validated ROC curve      ##
##  Variation in specificity       ##
##                                 ##
#####################################
##################################### 

line_seq_x <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1)
line_seq_y <- c(0.9, 0.92, 0.94, 0.96, 0.98, 1)


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,0.1), ylim=c(0.9,1),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Cross-validated ROC curve",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3, 
cex=lab.size, 
text="sensitivity")

for(i in 1:length(line_seq_y))
{
	points(x=c(-1,0.1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

axis(1, at = c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels = c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis = axis.size2) 

axis(2, at = c(0.9, 0.92, 0.94, 0.96, 0.98, 1), 
        labels=c("90%", "92%", "94%", "96%", "98%", "100%"), 
        las=2, cex.axis=axis.size ) 


points( x=1 - SS_seq, y=RF_xVAL_sens_quant[,1], 
    	  type='S', lwd=2, col="aquamarine4" )

polygon(x=c(1 - SS_seq, rev(1 - SS_seq)), 
	  y=c( RF_xVAL_sens_quant[,2], rev(RF_xVAL_sens_quant[,3]) ),
        col=rgb(69/256,139/256,116/256,0.2), border=NA)




#####################################
#####################################
##                                 ## 
##  PANEL 2                        ##
##  IMM measured prevalence        ##
##  Variation in specificity       ##
##                                 ##
#####################################
##################################### 

line_seq_x <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1)
line_seq_y <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,0.1), ylim=c(0.0,0.2),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Measured sero-prevalence (Paris, Apr 2020)",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3, 
cex=lab.size, 
text="sero-prevalence")

for(i in 1:length(line_seq_y))
{
	points(x=c(-1,0.1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

axis(1, at = c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels = c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis = axis.size2) 

axis(2, at = c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%", "18%", "20%"), 
        las=2, cex.axis=axis.size ) 



points(x=1 - SS_seq, y=RF_xIMM_meas_prev_quant[,1],
     type='l', lwd=2, col="aquamarine4" )

polygon( x = c( 1-SS_seq, rev(1-SS_seq) ), 
	   y = c( RF_xIMM_meas_prev_quant[,2], rev(RF_xIMM_meas_prev_quant[,3]) ), 
         col=rgb(69/256,139/256,116/256,0.2), border=NA)




#####################################
#####################################
##                                 ## 
##  PANEL 3                        ##
##  IMM adjusted prevalence        ##
##  Variation in specificity       ##
##                                 ##
#####################################
##################################### 

line_seq_x <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1)
line_seq_y <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,0.1), ylim=c(0.0,0.2),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Adjusted sero-prevalence (Paris, Apr 2020)",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3, 
cex=lab.size, 
text="sero-prevalence")

for(i in 1:length(line_seq_y))
{
	points(x=c(-1,0.1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

axis(1, at = c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels = c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis = axis.size2) 

axis(2, at = c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%", "18%", "20%"), 
        las=2, cex.axis=axis.size ) 



points(x=1 - SS_seq, y=RF_xIMM_adj_prev_quant[,1],
     type='l', lwd=2, col="aquamarine4" )

polygon( x = c( 1-SS_seq, rev(1-SS_seq) ), 
	   y = c( RF_xIMM_adj_prev_quant[,2], rev(RF_xIMM_adj_prev_quant[,3]) ), 
         col=rgb(69/256,139/256,116/256,0.2), border=NA)



dev.off()



















