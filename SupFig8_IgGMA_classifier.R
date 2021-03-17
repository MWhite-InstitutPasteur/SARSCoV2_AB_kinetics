library(binom)

AB_data <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Manuscript/Database/IDEA_db.csv")


##################################
## Exclude unknown and post-2020 negative

AB_data <- AB_data[ which(AB_data$covid_status %in% c("flow-positive", "pcr-positive", "pre-epidemic-negative")),]
AB_data <- AB_data[,-(38:39)]


AB_data$status <- rep( "pos", nrow(AB_data) )
AB_data$status[ which(AB_data$covid_status=="pre-epidemic-negative") ] <- "neg"



##################################
## Exclude individuals with missing days since symptom
## and with days <= 14

AB_data <- AB_data[-which( AB_data$status == "pos" & is.na(AB_data$days_pso)==TRUE ),]

AB_data <- AB_data[-which( AB_data$status == "pos" & AB_data$days_pso<=14 ),]



###################################################
###################################################
##                                               ##
##   ####  #   ## #####     ####  #   ## ######  ##
##  ##  ## ##  ## ##       ##  ## ##  ##   ##    ##
##  ##  ## ### ## ####     ###### ### ##   ##    ##
##  ##  ## ## ### ##       ##  ## ## ###   ##    ## 
##   ####  ##  ## #####    ##  ## ##  ##   ##    ##
##                                               ##
###################################################
###################################################


N_ant <- 9*3
N_SS  <- 50000

dil_cut <- exp(seq(from=log(1e-6), to=log(0.3), length=N_SS))

sens_mat <- matrix(NA, ncol=N_ant, nrow=N_SS)
spec_mat <- matrix(NA, ncol=N_ant, nrow=N_SS)

colnames(sens_mat) <- colnames(AB_data)[11:37]

colnames(spec_mat) <- colnames(AB_data)[11:37]

N_pos_vec <- rep(NA, N_ant)
N_neg_vec <- rep(NA, N_ant)

for(j in 1:(N_ant))
{
	ab_pos <- AB_data[which(AB_data$status=="pos"),10+j]
	ab_pos <- ab_pos[which(is.na(ab_pos)==FALSE)]

	ab_neg <- AB_data[which(AB_data$status=="neg"),10+j]
	ab_neg <- ab_neg[which(is.na(ab_neg)==FALSE)]

	for(i in 1:N_SS)
	{
		sens_mat[i,j] <- length(which( ab_pos > dil_cut[i] ))/length(ab_pos)
		spec_mat[i,j] <- length(which( ab_neg < dil_cut[i] ))/length(ab_neg)
	}

	N_pos_vec[j] <- length(ab_pos)
	N_neg_vec[j] <- length(ab_neg)
}


######################################
## Calculate Area Under Curve (AUC)


AUC_one_ant <- rep(NA, N_ant)

for(j in 1:(N_ant))
{
	AUC_one_ant[j] <- sum( (sens_mat[1:(N_SS-1),j] - sens_mat[2:N_SS,j])*
                             0.5*(spec_mat[1:(N_SS-1),j] + spec_mat[2:N_SS,j]) )
}





sens_target <- matrix(NA, nrow=3, ncol=N_ant)
colnames(sens_target) <- colnames(sens_mat)
rownames(sens_target) <- c("sens_99", "sens_spec", "spec_99")

spec_target <- matrix(NA, nrow=3, ncol=N_ant)
colnames(spec_target) <- colnames(spec_mat)
rownames(spec_target) <- c("sens_99", "sens_spec", "spec_99")

for(j in 1:(N_ant))
{
	sens_target[1,j] <- sens_mat[max(which(sens_mat[,j] > 0.99)),j]
	spec_target[1,j] <- spec_mat[max(which(sens_mat[,j] > 0.99)),j]

	sens_target[2,j] <- sens_mat[which.min(abs(sens_mat[,j] - spec_mat[,j])),j]
	spec_target[2,j] <- spec_mat[which.min(abs(sens_mat[,j] - spec_mat[,j])),j]

	sens_target[3,j] <- sens_mat[min(which(spec_mat[,j] > 0.99)),j]
	spec_target[3,j] <- spec_mat[min(which(spec_mat[,j] > 0.99)),j]
}




sens_target_lwr <- sens_target

for(i in 1:nrow(sens_target_lwr))
{
	for(j in 1:ncol(sens_target_lwr))
	{
		if( is.na(sens_target_lwr[i,j]) == FALSE )
		{
			sens_target_lwr[i,j] <- binom.confint( sens_target_lwr[i,j]*N_pos_vec[j], N_pos_vec[j], method="wilson")[1,5]
			sens_target_lwr[i,j] <- round( 100*sens_target_lwr[i,j], 1)  
		}
	}
}


sens_target_upr <- sens_target

for(i in 1:nrow(sens_target_upr))
{
	for(j in 1:ncol(sens_target_upr))
	{
		if( is.na(sens_target_upr[i,j]) == FALSE )
		{
			sens_target_upr[i,j] <- binom.confint( sens_target_upr[i,j]*N_pos_vec[j], N_pos_vec[j], method="wilson")[1,6]
			sens_target_upr[i,j] <- round( 100*sens_target_upr[i,j], 1)
		}  
	}
}


spec_target_lwr <- spec_target

for(i in 1:nrow(spec_target_lwr))
{
	for(j in 1:ncol(spec_target_lwr))
	{
		if( is.na(spec_target_lwr[i,j]) == FALSE )
		{
			spec_target_lwr[i,j] <- binom.confint( spec_target_lwr[i,j]*N_neg_vec[j], N_neg_vec[j], method="wilson")[1,5]
			spec_target_lwr[i,j] <- round( 100*spec_target_lwr[i,j], 1)  
		}
	}
}

spec_target_upr <- spec_target

for(i in 1:nrow(spec_target_upr))
{
	for(j in 1:ncol(spec_target_upr))
	{
		if( is.na(spec_target_upr[i,j]) == FALSE )
		{
			spec_target_upr[i,j] <- binom.confint( spec_target_upr[i,j]*N_neg_vec[j], N_neg_vec[j], method="wilson")[1,6]
			spec_target_upr[i,j] <- round( 100*spec_target_upr[i,j], 1)  
		}
	}
}

sens_target <- round( 100*sens_target, 1)



spec_target <- round( 100*spec_target, 1)





##############################################################
##############################################################
##                                                          ##
##  ##    # ##  ## ##   ###### ####    ####  #   ## ######  ##
##  ##   ## ##  ## ##     ##    ##    ##  ## ##  ##   ##    ##
##  ####### ##  ## ##     ##    ##    ###### ### ##   ##    ##
##  ## # ## ##  ## ##     ##    ##    ##  ## ## ###   ##    ## 
##  ##   ##  ##### #####  ##   ####   ##  ## ##  ##   ##    ##
##                                                          ##
##############################################################
##############################################################

set.seed(1234)

N_tree <- 10000


library(MASS)
library(ROCR)
library(randomForest)
library(pROC)

status <- as.factor(AB_data$status)

AB <- log(AB_data[,11:37])


####################################
## Variable Importance
## (all antigens)

RF_all = randomForest( status ~ ., data=AB, 
                       importance=TRUE, ntree=N_tree)

rf.roc_all <- roc(status, RF_all$votes[,2])

varImpPlot( RF_all )



RF_summary <- function( antigen_list )
{
	AB_trim <- AB[,antigen_list]
	status_trim <- status

	index_trim <- unique(which(is.na(AB_trim)==TRUE, arr.ind=TRUE)[,1])
	if( length(index_trim) > 0 )
	{
		AB_trim <- AB_trim[-index_trim,]
		status_trim <- status[-index_trim]
	}


	RF_ant = randomForest( status_trim ~ ., data=AB_trim, 
                                      importance=TRUE, ntree=N_tree)

	RF_ant_roc <- roc(status_trim, RF_ant$votes[,2])

	c( RF_ant_roc$sensitivities[min(which(RF_ant_roc$specificities > 0.99))],
         RF_ant_roc$auc,
         RF_ant_roc$specificities[max(which(RF_ant_roc$sensitivities > 0.99))] )
}


####################################
## 2 antigens: Spike IgG & other antigens

RF_2ant_summary <- matrix(NA, nrow=26, ncol=3)
rownames(RF_2ant_summary) <- colnames(AB)[2:27]
colnames(RF_2ant_summary) <- c("high_spec", "auc", "high_sens")

for(i in 1:26)
{
	RF_2ant_summary[i,] <- RF_summary( c(1,1+i) )
}


## Choose NP IgG (3)


RF_2ant = randomForest( status ~ ., data=AB[,c(1,3)], 
                        importance=TRUE, ntree=N_tree)

RF_2ant_roc <- roc(status, RF_2ant$votes[,2])


####################################
## 3 antigens: Spike IgG, NP IgG, & other antigens

RF_3ant_summary <- matrix(NA, nrow=25, ncol=3)
rownames(RF_3ant_summary) <- colnames( AB)[c(2,4:27)]
colnames(RF_3ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(2,4:27)

for(i in 1:25)
{
	RF_3ant_summary[i,] <- RF_summary( c(1,3,test_seq[i]) )
}

## Choose RBD IgG (2)


RF_3ant = randomForest( status ~ ., data=AB[,c(1,2,3)], 
                        importance=TRUE, ntree=N_tree)

RF_3ant_roc <- roc(status, RF_3ant$votes[,2])


####################################
## 4 antigens: Spike IgG, RBD IgG, NP IgG & other antigens

RF_4ant_summary <- matrix(NA, nrow=24, ncol=3)
rownames(RF_4ant_summary) <- colnames(AB)[c(4:27)]
colnames(RF_4ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(4:27)

for(i in 1:24)
{
	RF_4ant_summary[i,] <- RF_summary( c(1,2,3,test_seq[i]) )
}

## Choose spike IgA (19)

AB_trim <- AB[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1]),]
status_trim <- status[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])]

RF_4ant = randomForest( status_trim ~ ., data=AB_trim[,c(1,2,3,19)], 
                        importance=TRUE, ntree=N_tree)

RF_4ant_roc <- roc(status_trim, RF_4ant$votes[,2])




####################################
## 5 antigens: Spike IgG, RBD IgG, NP IgG 
##             Spike IgA & other antigens

RF_5ant_summary <- matrix(NA, nrow=23, ncol=3)
rownames(RF_5ant_summary) <- colnames(AB)[c(4:18,20:27)]
colnames(RF_5ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(4:18,20:27)

for(i in 1:23)
{
	RF_5ant_summary[i,] <- RF_summary( c(1,2,3,19,test_seq[i]) )
}

## Choose spike IgM (10)

AB_trim <- AB[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1]),]
status_trim <- status[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])]

RF_5ant = randomForest( status_trim ~ ., data=AB_trim[,c(1,2,3, 10, 19)], 
                        importance=TRUE, ntree=N_tree)

RF_5ant_roc <- roc(status_trim, RF_5ant$votes[,2])



####################################
## 6 antigens: Spike IgG, RBD IgG, NP IgG 
##             Spike IgA, Spike IgM & other antigens

RF_6ant_summary <- matrix(NA, nrow=22, ncol=3)
rownames(RF_6ant_summary) <- colnames(AB)[c(4:9,11:18,20:27)]
colnames(RF_6ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(4:9,11:18,20:27)

for(i in 1:22)
{
	RF_6ant_summary[i,] <- RF_summary( c(1,2,3,10,19,test_seq[i]) )
}

## Choose NP IgM (12)

AB_trim <- AB[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1]),]
status_trim <- status[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])]

RF_6ant = randomForest( status_trim ~ ., data=AB_trim[,c(1,2,3,10,12,19)], 
                        importance=TRUE, ntree=N_tree)

RF_6ant_roc <- roc(status_trim, RF_6ant$votes[,2])


####################################
## 7 antigens: Spike IgG, RBD IgG, NP IgG 
##             Spike IgA, Spike IgM, NP IgM & other antigens

RF_7ant_summary <- matrix(NA, nrow=21, ncol=3)
rownames(RF_7ant_summary) <- colnames(AB)[c(4:9,11,13:18,20:27)]
colnames(RF_7ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(4:9,11,13:18,20:27)

for(i in 1:21)
{
	RF_7ant_summary[i,] <- RF_summary( c(1,2,3,10,12,19,test_seq[i]) )
}

## Choose NP IgA (21)

AB_trim <- AB[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1]),]
status_trim <- status[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])]

RF_7ant = randomForest( status_trim ~ ., data=AB_trim[,c(1,2,3,10,12,19,21)], 
                        importance=TRUE, ntree=N_tree)

RF_7ant_roc <- roc(status_trim, RF_7ant$votes[,2])


####################################
## 8 antigens: Spike IgG, RBD IgG, NP IgG 
##             Spike IgA, Spike IgM, NP IgM, NP IgA & other antigens

RF_8ant_summary <- matrix(NA, nrow=20, ncol=3)
rownames(RF_8ant_summary) <- colnames(AB)[c(4:9,11,13:18,20,22:27)]
colnames(RF_8ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(4:9,11,13:18,20,22:27)

for(i in 1:20)
{
	RF_8ant_summary[i,] <- RF_summary( c(1,2,3,10,12,19,21,test_seq[i]) )
}

## Choose RBD IgA (20)

AB_trim <- AB[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1]),]
status_trim <- status[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])]

RF_8ant = randomForest( status_trim ~ ., data=AB_trim[,c(1,2,3,10,12,19,20,21)], 
                        importance=TRUE, ntree=N_tree)

RF_8ant_roc <- roc(status_trim, RF_8ant$votes[,2])



####################################
## 9 antigens: Spike IgG, RBD IgG, NP IgG 
##             Spike IgA, Spike IgM, NP IgM, NP IgA & other antigens

RF_9ant_summary <- matrix(NA, nrow=19, ncol=3)
rownames(RF_9ant_summary) <- colnames(AB)[c(4:9,11,13:18,22:27)]
colnames(RF_9ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(4:9,13:18,20,22:27)

for(i in 1:19)
{
	RF_9ant_summary[i,] <- RF_summary( c(1,2,3,10,12,19,20,21,test_seq[i]) )
}

## Choose RBD IgM (11)

AB_trim <- AB[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1]),]
status_trim <- status[-unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])]

RF_9ant = randomForest( status_trim ~ ., data=AB_trim[,c(1,2,3,10,11,12,19,20,21)], 
                        importance=TRUE, ntree=N_tree)

RF_9ant_roc <- roc(status_trim, RF_9ant$votes[,2])











####################################
##

Sens_Combo <- matrix(NA, nrow=9, ncol=3)

Sens_Combo[1,] <- c( sens_target[3,1], sens_target_lwr[3,1], sens_target_upr[3,1] )/100

Sens_Combo[2,1] <- RF_2ant_summary[2,1]

Sens_Combo[3,1] <- RF_3ant_summary[1,1]

Sens_Combo[4,1] <- RF_4ant_summary[16,1]

Sens_Combo[5,1] <- RF_5ant_summary[7,1]

Sens_Combo[6,1] <- RF_6ant_summary[8,1]

Sens_Combo[7,1] <- RF_7ant_summary[15,1]

Sens_Combo[8,1] <- RF_8ant_summary[14,1]

Sens_Combo[9,1] <- RF_9ant_summary[7,1]

N_pos <- length(which(AB_data$status=="pos"))
N_neg <- length(which(AB_data$status=="neg"))




for(k in 2:9)
{
	Sens_Combo[k,2] <- binom.confint( Sens_Combo[k,1]*N_pos, N_pos, method="wilson")[1,5]
	Sens_Combo[k,3] <- binom.confint( Sens_Combo[k,1]*N_pos, N_pos, method="wilson")[1,6]
}

rownames(Sens_Combo) <- c("spike_IgG", "NP_IgG", "RBD_IgG", 
                          "spike_IgA", "spike_IgM", "NP_IgM",
                          "NP_IgA", "RBD_IgA", "RBD_IgM")




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


N_tree <- 10000

N_rep <- 100

N_SS_cut <- 5000

SS_seq <- seq(from=0, to=1, length=N_SS_cut)




####################################
## Uncertainty in sensitivity and specificity 
## for algorithmvia cross-validation - randomly
## leaving 1/3 out

#############################################################################
## Random forests classifier for disjoint training and testing data

RF_multi_ant_SS = function( AB_train, status_train, AB_test, status_test, antigen_list )
{
	AB_train = as.matrix(AB_train[,antigen_list])

	if( length(which(is.na(AB_train))) > 0 )
	{
		status_train = status_train[-which(is.na(AB_train), arr.ind=TRUE)[,1]]

		AB_train = AB_train[-which(is.na(AB_train), arr.ind=TRUE)[,1],]
	}


	AB_test = as.matrix(AB_test[,antigen_list])

	if( length(which(is.na(AB_test))) > 0 )
	{
		status_test = status_test[-which(is.na(AB_test), arr.ind=TRUE)[,1]]

		AB_test = AB_test[-which(is.na(AB_test), arr.ind=TRUE)[,1],]
	}

	############################################
	## Fit Rforest and create prediction object

	tryCatch(
	{
		RF_multi_ant = randomForest( as.factor(status_train) ~ ., data=as.data.frame(AB_train), 
                                      importance=TRUE, ntree=N_tree )

		RF_multi_ant_pred_obj = predict( RF_multi_ant, newdata=AB_test, predict.all=TRUE)

		RF_multi_ant_votes = rowSums(RF_multi_ant_pred_obj$individual=="neg")/N_tree

		RF_multi_ant_roc <- roc(status_test, RF_multi_ant_votes)

	}, error=function(e){ NULL }
	)

	SS_mat = cbind( RF_multi_ant_roc$sensitivities, RF_multi_ant_roc$specificities  )
	colnames(SS_mat) = c("sens", "spec" )

	SS_mat
}


#############################################################################
## Cross-validated random forests classification algorithm

RF_multi_ant_SS_xval = function( AB_dd, status_data, antigen_list, train_prop, N_rep )
{
	AB_class <- list()	

	for(n in 1:N_rep)
	{
		############################################
		## Prepare testing and training data

		index_train = sample( nrow(AB_dd), train_prop*nrow(AB_dd) )
		index_test  = setdiff( 1:nrow(AB_dd), index_train )

		status_train = status_data[index_train]
		status_test  = status_data[index_test]

		AB_train = AB_dd[index_train,]
		AB_test  = AB_dd[index_test,]

		
		AB_class[[n]] = RF_multi_ant_SS( AB_train, status_train, AB_test, status_test, antigen_list )
	}

	AB_class
}



#########################################
##                                     ##
##   1 antigens                        ##
##                                     ##
#########################################


spike_IgG_pos <- AB_data$spike_IgG[which(AB_data$status == "pos")]
spike_IgG_neg <- AB_data$spike_IgG[which(AB_data$status == "neg")]


sens_1ant_xval_mat <- matrix(NA, nrow=N_SS, ncol=N_rep)
spec_1ant_xval_mat <- matrix(NA, nrow=N_SS, ncol=N_rep)

for(n in 1:N_rep)
{
	spike_IgG_pos_sub <- spike_IgG_pos[sample( length(spike_IgG_pos), (2/3)*length(spike_IgG_pos) )]		
	spike_IgG_neg_sub <- spike_IgG_neg[sample( length(spike_IgG_neg), (2/3)*length(spike_IgG_neg) )]		

	for(i in 1:N_SS)
	{
		sens_1ant_xval_mat[i,n] <- length(which( spike_IgG_pos_sub > dil_cut[i] ))/length(spike_IgG_pos_sub)
 		spec_1ant_xval_mat[i,n] <- length(which( spike_IgG_neg_sub < dil_cut[i] ))/length(spike_IgG_neg_sub)
	}
}


sens_1ant_xval_quant <- matrix(NA, nrow=N_SS, ncol=3)
spec_1ant_xval_quant <- matrix(NA, nrow=N_SS, ncol=3)

colnames(sens_1ant_xval_quant) <- c( "med", "lwr", "upr" )
colnames(spec_1ant_xval_quant) <- c( "med", "lwr", "upr" )

for(i in 1:N_SS)
{
	sens_1ant_xval_quant[i,] <- quantile( sens_1ant_xval_mat[i,], prob=c(0.5, 0.025, 0.975) )
	spec_1ant_xval_quant[i,] <- quantile( spec_1ant_xval_mat[i,], prob=c(0.5, 0.025, 0.975) )
}

index_1 <- which.min(abs(dil_cut - quantile( spike_IgG_neg, prob=0.99)))

sens_xval_1ant <- sens_1ant_xval_quant[index_1,]


#########################################
##                                     ##
##   2 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_2_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,3), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_2_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_2_ant_xval_se_vary_sp_fix[i,n] <- RF_2_ant_xval_rep[[n]][min(which(RF_2_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_2_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_2_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_2_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_2_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_2_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_2 <- which.min(abs(SS_seq - 0.99))


sens_xval_2ant <- RF_2_ant_xval_se_vary_sp_fix_quant[index_2,]



#########################################
##                                     ##
##   3 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_3_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,2,3), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_3_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_3_ant_xval_se_vary_sp_fix[i,n] <- RF_3_ant_xval_rep[[n]][min(which(RF_3_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_3_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_3_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_3_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_3_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_3_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}

index_3 <- which.min(abs(SS_seq - 0.99))


sens_xval_3ant <- RF_3_ant_xval_se_vary_sp_fix_quant[index_3,]



#########################################
##                                     ##
##   4 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_4_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,2,3,19), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_4_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_4_ant_xval_se_vary_sp_fix[i,n] <- RF_4_ant_xval_rep[[n]][min(which(RF_4_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_4_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_4_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_4_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_4_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_4_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_4 <- which.min(abs(SS_seq - 0.99))


sens_xval_4ant <- RF_4_ant_xval_se_vary_sp_fix_quant[index_4,]




#########################################
##                                     ##
##   5 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_5_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,2,3,10,19), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_5_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_5_ant_xval_se_vary_sp_fix[i,n] <- RF_5_ant_xval_rep[[n]][min(which(RF_5_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_5_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_5_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_5_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_5_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_5_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_5 <- which.min(abs(SS_seq - 0.99))


sens_xval_5ant <- RF_5_ant_xval_se_vary_sp_fix_quant[index_5,]




#########################################
##                                     ##
##   6 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_6_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,3,2,10,12,19), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_6_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_6_ant_xval_se_vary_sp_fix[i,n] <- RF_6_ant_xval_rep[[n]][min(which(RF_6_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_6_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_6_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_6_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_6_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_6_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_6 <- which.min(abs(SS_seq - 0.99))


sens_xval_6ant <- RF_6_ant_xval_se_vary_sp_fix_quant[index_6,]






#########################################
##                                     ##
##   7 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_7_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,2,3,10,12,19,21), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_7_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_7_ant_xval_se_vary_sp_fix[i,n] <- RF_7_ant_xval_rep[[n]][min(which(RF_7_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_7_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_7_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_7_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_7_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_7_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_7 <- which.min(abs(SS_seq - 0.99))


sens_xval_7ant <- RF_7_ant_xval_se_vary_sp_fix_quant[index_7,]






#########################################
##                                     ##
##   8 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_8_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,2,3,10,12,19,20,21), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_8_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_8_ant_xval_se_vary_sp_fix[i,n] <- RF_8_ant_xval_rep[[n]][min(which(RF_8_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_8_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_8_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_8_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_8_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_8_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_8 <- which.min(abs(SS_seq - 0.99))


sens_xval_8ant <- RF_8_ant_xval_se_vary_sp_fix_quant[index_8,]






#########################################
##                                     ##
##   9 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_9_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,2,3,10,11,12,19,20,21), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_9_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_9_ant_xval_se_vary_sp_fix[i,n] <- RF_9_ant_xval_rep[[n]][min(which(RF_9_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_9_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_9_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_9_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_9_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_9_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_9 <- which.min(abs(SS_seq - 0.99))


sens_xval_9ant <- RF_9_ant_xval_se_vary_sp_fix_quant[index_9,]











sens_target_combi <- matrix(NA, nrow=8, ncol=3)
rownames(sens_target_combi) <- c("NP_IgG", "RBD_IgG", 
                                 "spike_IgA", "spike_IgM", "NP_IgM",
                                 "NP_IgA", "RBD_IgA", "RBD_IgM")
colnames(sens_target_combi) <- c("spec_99", "sens_spec", "sens_99")

sens_target_combi_lwr <- sens_target_combi
sens_target_combi_upr <- sens_target_combi

spec_target_combi <- sens_target_combi
spec_target_combi_lwr <- sens_target_combi
spec_target_combi_upr <- sens_target_combi



for(i in 1:8)
{
	if( i == 1 )
	{
		RFX <- RF_2ant_roc
	}

	if( i == 2 )
	{
		RFX <- RF_3ant_roc
	}		

	if( i == 3 )
	{
		RFX <- RF_4ant_roc
	}		

	if( i == 4 )
	{
		RFX <- RF_5ant_roc
	}		

	if( i == 5 )
	{
		RFX <- RF_6ant_roc
	}

	if( i == 6 )
	{
		RFX <- RF_7ant_roc
	}

	if( i == 7 )
	{
		RFX <- RF_8ant_roc
	}

	if( i == 8 )
	{
		RFX <- RF_9ant_roc
	}


	################################
	## High specificity target

	sens_target_combi[i,1]     <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_pos, N_pos, method="wilson")[1,4]
	sens_target_combi_lwr[i,1] <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_pos, N_pos, method="wilson")[1,5]
	sens_target_combi_upr[i,1] <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_pos, N_pos, method="wilson")[1,6]

	spec_target_combi[i,1]     <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_neg, N_neg, method="wilson")[1,4]
	spec_target_combi_lwr[i,1] <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_neg, N_neg, method="wilson")[1,5]
	spec_target_combi_upr[i,1] <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_neg, N_neg, method="wilson")[1,6]


	################################
	## Balanced sensitivty and specificity target

	sens_target_combi[i,2]     <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_pos, N_pos, method="wilson")[1,4]
	sens_target_combi_lwr[i,2] <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_pos, N_pos, method="wilson")[1,5]
	sens_target_combi_upr[i,2] <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_pos, N_pos, method="wilson")[1,6]

	spec_target_combi[i,2]     <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_neg, N_neg, method="wilson")[1,4]
	spec_target_combi_lwr[i,2] <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_neg, N_neg, method="wilson")[1,5]
	spec_target_combi_upr[i,2] <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_neg, N_neg, method="wilson")[1,6]


	################################
	## High sensitivity target

	sens_target_combi[i,3]     <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_pos, N_pos, method="wilson")[1,4]
	sens_target_combi_lwr[i,3] <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_pos, N_pos, method="wilson")[1,5]
	sens_target_combi_upr[i,3] <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_pos, N_pos, method="wilson")[1,6]

	spec_target_combi[i,3]     <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_neg, N_neg, method="wilson")[1,4]
	spec_target_combi_lwr[i,3] <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_neg, N_neg, method="wilson")[1,5]
	spec_target_combi_upr[i,3] <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_neg, N_neg, method="wilson")[1,6]

}



spec_target_combi     <- round(100*spec_target_combi,1)
spec_target_combi_lwr <- round(100*spec_target_combi_lwr,1)
spec_target_combi_upr <- round(100*spec_target_combi_upr,1)

sens_target_combi     <- round(100*sens_target_combi,1)
sens_target_combi_lwr <- round(100*sens_target_combi_lwr,1)
sens_target_combi_upr <- round(100*sens_target_combi_upr,1)


################################################################################
################################################################################
##                                                                            ##
##   ####   ####  #####  #####  ##### ##     ####  ###### ####  ####  #   ##  ##
##  ##  ## ##  ## ##  ## ##  ## ##    ##    ##  ##   ##    ##  ##  ## ##  ##  ##
##  ##     ##  ## #####  #####  ####  ##    ######   ##    ##  ##  ## ### ##  ##
##  ##  ## ##  ## ## ##  ## ##  ##    ##    ##  ##   ##    ##  ##  ## ## ###  ##
##   ####   ####  ##  ## ##  ## ##### ##### ##  ##   ##   ####  ####  ##  ##  ##
##                                                                            ##
################################################################################
################################################################################

library(fields)

AB_cor <- matrix(NA, nrow=27, ncol=27)

cor_seq <- c(11:37)


for(i in 1:27)
{
	for(j in 1:27)
	{
		AB_i <- log(AB_data[,cor_seq[i]])
		AB_j <- log(AB_data[,cor_seq[j]])

		index_ij <- which( is.na(AB_i)==FALSE & is.na(AB_j)==FALSE )

		AB_i <- AB_i[index_ij]
		AB_j <- AB_j[index_ij]

		AB_cor[i,j] = cor( AB_i, AB_j, method="spearman")
	}
}


save.image("SupFig8_IgGMA_classifier.RData")





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

ant_cols <- c("red", "gold", "forestgreen", "lightgreen", "salmon",
              "lightslateblue", "lightskyblue", "navy", "purple",
              "red",  "gold", "forestgreen", "lightgreen", "salmon",
              "lightslateblue", "lightskyblue", "navy", "purple",
              "red", "gold", "forestgreen", "lightgreen", "salmon",
              "lightslateblue", "lightskyblue", "navy", "purple")

plot_names <- c("S IgG", "RBD IgG", "N IgG", "S2 IgG", "ME IgG", 
                "NL63 IgG", "229E IgG", "OC43 IgG", "HKU1 IgG",
                "S IgM", "RBD IgM", "N IgM", "S2 IgM", "ME IgM", 
                "NL63 IgM", "229E IgM", "OC43 IgM", "HKU1 IgM",
                "S IgA", "RBD IgA", "N IgA", "S2 IgA", "ME IgA", 
                "NL63 IgA", "229E IgA", "OC43 IgA", "HKU1 IgA" )


combi_cols <- c("red", "gold", "forestgreen", 
                "salmon", "magenta", "limegreen",
                "seagreen1", "yellow1", "orange3") 
                


lab.size   = 1.5
axis.size  = 1.25	
main.size  = 1.1
line.size  = 2
dash.line.size <- 1



tiff(file="SupFig8_IgGMA_classifier.tif", width=30, height=18, units="cm", res=500)

lay.mat <- rbind( c( 1, 2, 3, 7 ),
                  c( 4, 5, 6, 7 ),
                  c( 4, 5, 6, 8 ) )
layout(lay.mat, heights=c(10,8,2), widths=c(1,1,1,2))
layout.show(8)




########################
########################
##                    ##              
##  PANEL 1           ##            
##  IgG ROC analysis  ##
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
main="(A) IgG ROC analysis",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

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




for(k in 1:9)
{
	points( x=1-spec_mat[,k], y=sens_mat[,k], 
    		  type='S', lwd=2, col=ant_cols[k] )
}


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=1.0) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=1.5 ) 




########################
########################
##                    ##              
##  PANEL 2           ##            
##  IgM ROC analysis  ##
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
main="(B) IgM ROC analysis",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

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




for(k in 10:18)
{
	points( x=1-spec_mat[,k], y=sens_mat[,k], 
    		  type='S', lwd=2, col=ant_cols[k] )
}


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=1.0) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=1.5 ) 




########################
########################
##                    ##              
##  PANEL 3           ##            
##  IgA ROC analysis  ##
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
main="(C) IgA ROC analysis",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

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




for(k in 19:27)
{
	points( x=1-spec_mat[,k], y=sens_mat[,k], 
    		  type='S', lwd=2, col=ant_cols[k] )
}


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=1.0) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=1.5 ) 





########################
########################
##                    ##              
##  PANEL 4           ##            
##  AUC analysis      ##
##                    ##
########################
########################

par(mar=c(7,6,2,1))
par(mgp=c(2.5, 1.0,0))


gap <- 0.1

line_seq_y <- c(0.0, 0.25, 0.5, 0.75, 1.0)





plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,27*(1+gap)), ylim=c(0,1.002),
xaxs='i',
# yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(D) Area under ROC curve",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)




mtext(side = 2, line = 4, 
cex=1.0, 
text="area under ROC curve")




for(k in 1:27)
{
	polygon( x=gap*(k-1) + c(k-1,k,k,k-1), c(0,0,AUC_one_ant[k],AUC_one_ant[k]), col=ant_cols[k], border=NA)
}



for(i in 1:length(line_seq_y))
{
	points(x=c(-1000,365), y=rep(line_seq_y[i],2), 
             type='l', lwd=dash.line.size, col="grey", lty="dashed")
}


axis(2, at=c(0.0, 0.25, 0.5, 0.75, 1), 
        labels=c("0.0", "0.25", "0.5", "0.75", "1.0"), 
        las=2, cex.axis=1.5 ) 



axis(1, at=seq(from=0.5, by=1+gap, length=27), label=plot_names, las=2, cex.axis=0.5)


########################
########################
##                    ##              
##  PANEL 5           ##            
##  Multi-antigen ROC ##
##                    ##
########################
########################


par(mar=c(7,6,2,1))
par(mgp=c(2.5, 1.25,0))


line_seq_x <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1)
line_seq_y <- c(0.90, 0.92, 0.94, 0.96, 0.98, 1)



plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,0.1), ylim=c(0.9,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(E) Multi-antigen analysis",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


mtext(side = 1, line = 3, 
cex=1.0, 
text="1 - specificity")

mtext(side = 2, line = 4, 
cex=1.0, 
text="sensitivity")


for(i in 1:length(line_seq_y))
{
	points(x=c(-1,1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


points( x=1-spec_mat[,1], y=sens_mat[,1], 
    	  type='S', lwd=2, col=ant_cols[1] )


points( x=1-RF_2ant_roc$specificities, y=RF_2ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[2] )

points( x=1-RF_3ant_roc$specificities, y=RF_3ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[3] )

points( x=1-RF_4ant_roc$specificities, y=RF_4ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[4] )

points( x=1-RF_5ant_roc$specificities, y=RF_5ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[5] )

points( x=1-RF_6ant_roc$specificities, y=RF_6ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[6] )

points( x=1-RF_7ant_roc$specificities, y=RF_7ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[7] )

points( x=1-RF_8ant_roc$specificities, y=RF_8ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[8] )

points( x=1-RF_9ant_roc$specificities, y=RF_9ant_roc$sensitivities, 
    	  type='S', lwd=2, col=combi_cols[9] )


axis(1, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis=1.0) 

axis(2, at=c(0.90, 0.92, 0.94, 0.96, 0.98, 1), 
        labels=c("0%", "92%", "94%", "96%", "98%", "100%"), 
        las=2, cex.axis=1.5 ) 





########################
########################
##                    ##              
##  PANEL 6           ##            
##  Spex Xval         ##
##                    ##
########################
########################


par(mar=c(7,6,2,1))
par(mgp=c(2.5, 1.25,0))

gap <- 0.1

line_seq_y <- c(0.90, 0.92, 0.94, 0.96, 0.98, 1)



plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,9*(1+gap)), ylim=c(0.9,1.002),
xaxs='i',
# yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(F) High spec (99%) target",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



mtext(side = 2, line = 4, 
cex=1.0, 
text="sensitivity")


combi_cols2 <- combi_cols[c(1,3,2,5,4,6:9)]

for(k in 1:9)
{
	polygon( x=gap*(k-1) + c(k-1,k,k,k-1), c(0,0,Sens_Combo[k,1],Sens_Combo[k,1]), col=combi_cols2[k], border=NA)
}

arrows( x0=gap*(1-1) + 1 - 0.5, x1=gap*(1-1) + 1 - 0.5, y0=sens_xval_1ant[2], y1=sens_xval_1ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(2-1) + 2 - 0.5, x1=gap*(2-1) + 2 - 0.5, y0=sens_xval_2ant[2], y1=sens_xval_2ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(3-1) + 3 - 0.5, x1=gap*(3-1) + 3 - 0.5, y0=sens_xval_3ant[2], y1=sens_xval_3ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(4-1) + 4 - 0.5, x1=gap*(4-1) + 4 - 0.5, y0=sens_xval_4ant[2], y1=sens_xval_4ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(5-1) + 5 - 0.5, x1=gap*(5-1) + 5 - 0.5, y0=sens_xval_5ant[2], y1=sens_xval_5ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(6-1) + 6 - 0.5, x1=gap*(6-1) + 6 - 0.5, y0=sens_xval_6ant[2], y1=sens_xval_6ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(7-1) + 7 - 0.5, x1=gap*(7-1) + 7 - 0.5, y0=sens_xval_7ant[2], y1=sens_xval_7ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(8-1) + 8 - 0.5, x1=gap*(8-1) + 8 - 0.5, y0=sens_xval_8ant[2], y1=sens_xval_8ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  

arrows( x0=gap*(9-1) + 9 - 0.5, x1=gap*(9-1) + 9 - 0.5, y0=sens_xval_9ant[2], y1=sens_xval_9ant[3], 
        code=3, col="black", angle=90, lwd=2, length=0.05)  





points( x=gap*(1-1) + 1 - 0.5, y=sens_xval_1ant[1], pch=19, cex=1, col="black") 

points( x=gap*(2-1) + 2 - 0.5, y=sens_xval_2ant[1], pch=19, cex=1, col="black") 

points( x=gap*(3-1) + 3 - 0.5, y=sens_xval_3ant[1], pch=19, cex=1, col="black") 

points( x=gap*(4-1) + 4 - 0.5, y=sens_xval_4ant[1], pch=19, cex=1, col="black") 

points( x=gap*(5-1) + 5 - 0.5, y=sens_xval_5ant[1], pch=19, cex=1, col="black") 

points( x=gap*(6-1) + 6 - 0.5, y=sens_xval_6ant[1], pch=19, cex=1, col="black") 

points( x=gap*(7-1) + 7 - 0.5, y=sens_xval_7ant[1], pch=19, cex=1, col="black") 

points( x=gap*(8-1) + 8 - 0.5, y=sens_xval_8ant[1], pch=19, cex=1, col="black") 

points( x=gap*(9-1) + 9 - 0.5, y=sens_xval_9ant[1], pch=19, cex=1, col="black") 

for(i in 1:length(line_seq_y))
{
	points(x=c(0,100), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}



axis(2, at=c(0.9, 0.92, 0.94, 0.96, 0.98, 1), 
        labels=c("90%", "92%", "94%", "96%", "98%", "100%"), 
        las=2, cex.axis=1.5 ) 



axis(1, at=seq(from=0.5, by=1+gap, length=9), 
           label=c("S IgG", "+ N IgG", "+ RBD IgG", 
                   "+ S IgA", "+ S IgM", "+ N IgM",
                   "+ N IgA", "+ RBD IgA", "+ RBD IgM" ), 
            las=2, cex.axis=1)




########################
########################
##                    ##              
##  PANEL 7           ##            
##  IgG correlation   ##
##                    ##
########################
########################


par(mar=c(5,7,2,1))
par(mgp=c(2.5, 1.25,0))




cor_cols = c("springgreen4", "springgreen", "palegreen", "yellowgreen", "yellow", 
             "gold", "orange", "orangered", "firebrick1", "red3" )

N_cor_steps = length(cor_cols)

plot(x=100, y=100,
xlim=c(0,27), ylim=c(0,27),
xlab="", ylab="",
main="(G) Correlation between antibodies",
xaxt='n', yaxt='n', xaxs='i', yaxs='i',
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:27)
{
	for(j in 1:27)
	{
		polygon(x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                    border=NA, col="blue" )

		polygon(x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                    border=NA, col=cor_cols[ceiling(N_cor_steps*AB_cor[i,j])] )
	}
}


points(x=c(0,27), y=c(9,9), type='l', lwd=2)
points(x=c(0,27), y=c(18,18), type='l', lwd=2)

points(x=c(9,9), y=c(0,27), type='l', lwd=2)
points(x=c(18,18), y=c(0,27), type='l', lwd=2)


axis(1, at=seq(from=0.5, by=1, length=27), label=plot_names, las=2, cex.axis=1)
axis(2, at=seq(from=0.5, by=1, length=27), label=plot_names, las=2, cex.axis=1)


#############
## LEGEND

par(mar=c(3,6,1,2))


plot(x=100, y=100,
xlim=c(-10,100), ylim=c(0,1),
xlab="", ylab="",
main="",
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n')

polygon(y=c(0,1,1,0), x=c(-10,-10,0,0),
              border=NA, col="blue")	

for(i in 1:length(cor_cols))
{
	polygon(y=c(0,1,1,0), x=c(i-1,i-1,i,i)*100/N_cor_steps,
              border=NA, col=cor_cols[i])	
}


axis(1, at=100*c(-0.1,0,0.25,0.5,0.75,1), label=c("-100%", "0%", "25%", "50%", "75%", "100%"), cex.axis=1)




dev.off()








