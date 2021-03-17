

IMM_db <- read.csv("C:/U/CoronaVirus/PELLEAU_et_al/Databases/IDEA/IMM_sero_with_neut.csv")



index_neg <- which( IMM_db$IgGMA_status == "neg" )

index_pos_r1 <- which( IMM_db$IgGMA_status == "pos" & IMM_db$site == "IMM_r1" )

index_pos_r1 <- index_pos_r1[-which(IMM_db$days_pso[index_pos_r1] > 80)]

index_pos_r2 <- which( IMM_db$IgGMA_status == "pos" & IMM_db$site == "IMM_r2"  )

index_pos <- c(index_pos_r1, index_pos_r2)

IMM_db <- IMM_db[index_pos,]





plot( x = IMM_db$days_pso[index_pos_r1], 
      y = IMM_db$NT_T1_titre[index_pos_r1],
      xlim=c(-100,400), pch=19, log="y" )


points( x = IMM_db$days_pso[index_pos_r2], 
        y = IMM_db$NT_T2_titre[index_pos_r2] )




plot( x = IMM_db$days_pso[index_pos_r1], 
      y = IMM_db$FRNT_titre[index_pos_r1],
      xlim=c(-100,400), pch=19, log="y" )




plot( x = IMM_db$NT_T1_titre[index_pos_r1], 
      y = IMM_db$FRNT_titre[index_pos_r1],
      pch=19, log="xy" )






plot( x = IMM_db$NT_T1_titre[index_pos_r1], 
      y = IMM_db$RBD_IgG[index_pos_r1],
      pch=19, log="xy" )

points( x = IMM_db$NT_T2_titre[index_pos_r2], 
        y = IMM_db$RBD_IgG[index_pos_r2],
        log="xy" )




plot( x = IMM_db$NT_T1_titre[index_pos_r1], 
      y = IMM_db$RBD_IgA[index_pos_r1],
      pch=19, log="xy" )

points( x = IMM_db$NT_T2_titre[index_pos_r2], 
        y = IMM_db$RBD_IgA[index_pos_r2],
        log="xy" )




IMM_db$NT_T1_titre[which(IMM_db$NT_T1_titre == 0)] <- NA

IMM_db$NT_T2_titre[which(IMM_db$NT_T2_titre == 0)] <- NA

IMM_db$FRNT_titre[which(IMM_db$FRNT_titre == 0)] <- NA



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

line_seq_dil_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)
line_seq_dil_y <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)


tiff(file="SupFig4_neutralization.tif", width=40, height=30, units="cm", res=500)


lay.mat <- rbind( c( 1, 1, 2 ),
                  c( 3, 4, 5 ) )
layout(lay.mat, heights=c(1,1), widths=c(1,1,1))
layout.show(5)





par(mar=c(8,8,3,1))
par(mgp=c(4.5, 1, 0))

point.size = 1.5
lab.size   = 2
axis.size  = 1.5
main.size  = 2.5


###################################
## Panel 1: Spike_IPP IgG vs S1RBD_NA IgG

line_seq_dil_x <- 30*c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
line_seq_dil_y <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)



plot( x = IMM_db$days_pso[which(IMM_db$site=="IMM_r1")], 
      y = IMM_db$NT_T1_titre[which(IMM_db$site=="IMM_r1")],
      xlim=c(0, 360), ylim=c(1, 30000), log="y",
      pch=19, col="deeppink", cex=point.size,
	yaxt='n', xaxt='n', bty='n', 
	xlab="time post symptoms (months)", 
	ylab="reciprocal dilution (neutralization IC50)", 
	main="(A) S Fuse neutralization",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_y))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_x))
{
	points(x=rep(line_seq_dil_x[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points( x = IMM_db$days_pso[which(IMM_db$site=="IMM_r2")], 
        y = IMM_db$NT_T2_titre[which(IMM_db$site=="IMM_r2")],
        pch=19, col="grey11", cex=point.size )


axis(1, at=30*c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        label=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
        cex.axis=axis.size )

axis(2, at=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000), 
        label=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000),
        las=2, cex.axis=axis.size )



ID_draw <- unique(IMM_db$part_ID)

decay_slope <- rep(NA, length(ID_draw))

for(i in 1:length(ID_draw))
{
	index_plot <- which( IMM_db$part_ID == ID_draw[i] )

	if( is.na(IMM_db$NT_T1_titre[index_plot])==FALSE &
          is.na(IMM_db$NT_T2_titre[index_plot])==FALSE )
	{
		points( x = IMM_db$days_pso[index_plot],
                    y = c( IMM_db$NT_T1_titre[index_plot[1]], IMM_db$NT_T2_titre[index_plot[2]] ),
                    type='l', col="grey" )

		decay_slope[i] <- ( log(IMM_db$NT_T1_titre[index_plot[1]]) - log(IMM_db$NT_T2_titre[index_plot[2]]) )/
                              ( IMM_db$days_pso[index_plot[1]] - IMM_db$days_pso[index_plot[2]] ) 
	} 
}





###################################
## Panel 2: Spike_IPP IgG vs S1RBD_NA IgG

line_seq_dil_x <- 30*c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
line_seq_dil_y <- c(1, 3, 10, 30, 100, 300, 1000)



plot( x = IMM_db$days_pso[which(IMM_db$site=="IMM_r1")], 
      y = IMM_db$FRNT_titre[which(IMM_db$site=="IMM_r1")],
      xlim=c(0, 90), ylim=c(1, 1000), log="y",
      pch=19, col="deeppink", cex=point.size,
	yaxt='n', xaxt='n', bty='n', 
	xlab="time post symptoms (months)", 
	ylab="FRNT 50 titre", 
	main="(B) Foci Reduction Neutralization Test",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=0.75*main.size)

for(i in 1:length(line_seq_dil_y))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_x))
{
	points(x=rep(line_seq_dil_x[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}


axis(1, at=30*c(0, 1, 2, 3), 
        label=c(0, 1, 2, 3),
        cex.axis=axis.size )

axis(2, at=c(1, 3, 10, 30, 100, 300, 1000), 
        label=c(1, 3, 10, 30, 100, 300, 1000),
        las=2, cex.axis=axis.size )




###################################
## Panel 3: Spike_IPP IgG vs S1RBD_NA IgG

line_seq_dil_x <- c(1, 3, 10, 30, 100, 300, 1000)
line_seq_dil_y <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)



plot( x = IMM_db$FRNT_titre[which(IMM_db$site=="IMM_r1")], 
      y = IMM_db$NT_T1_titre[which(IMM_db$site=="IMM_r1")],
      xlim=c(1, 1000), ylim=c(1, 30000), log="xy",
      pch=19, col="deeppink", cex=point.size,
	yaxt='n', xaxt='n', bty='n', 
	xlab="FRNT 50 titre", 
	ylab="reciprocal dilution (neutralization IC50)", 
	main="(C) S Fuse vs FRNT",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_y))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_x))
{
	points(x=rep(line_seq_dil_x[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points( x = IMM_db$RBD_IgG[which(IMM_db$site=="IMM_r2")], 
        y = IMM_db$NT_T2_titre[which(IMM_db$site=="IMM_r2")],
        pch=19, col="grey11" )


axis(1, at=c(1, 3, 10, 30, 100, 300, 1000), 
        label=c(1, 3, 10, 30, 100, 300, 1000),
        cex.axis=axis.size )

axis(2, at=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000), 
        label=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000),
        las=2, cex.axis=axis.size )






###################################
## Panel 4: Spike_IPP IgG vs S1RBD_NA IgG

line_seq_dil_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)
line_seq_dil_y <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)



plot( x = IMM_db$RBD_IgG[which(IMM_db$site=="IMM_r1")], 
      y = IMM_db$NT_T1_titre[which(IMM_db$site=="IMM_r1")],
      xlim=c(0.00001, 0.01), ylim=c(1, 30000), log="xy",
      pch=19, col="deeppink", cex=point.size,
	yaxt='n', xaxt='n', bty='n', 
	xlab="anti-RBD IgG (RAU)", 
	ylab="reciprocal dilution (neutralization IC50)", 
	main="(D) S Fuse vs anti-RBD IgG",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_y))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_x))
{
	points(x=rep(line_seq_dil_x[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points( x = IMM_db$RBD_IgG[which(IMM_db$site=="IMM_r2")], 
        y = IMM_db$NT_T2_titre[which(IMM_db$site=="IMM_r2")],
        pch=19, col="grey11",  cex=point.size )


axis(1, at=c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01), 
        label=c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01),
        cex.axis=axis.size )

axis(2, at=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000), 
        label=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000),
        las=2, cex.axis=axis.size )






###################################
## Panel 5: Spike_IPP IgG vs S1RBD_NA IgG

line_seq_dil_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)
line_seq_dil_y <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)



plot( x = IMM_db$RBD_IgA[which(IMM_db$site=="IMM_r1")], 
      y = IMM_db$NT_T1_titre[which(IMM_db$site=="IMM_r1")],
      xlim=c(0.00001, 0.01), ylim=c(1, 30000), log="xy",
      pch=19, col="deeppink", cex=point.size,
	yaxt='n', xaxt='n', bty='n', 
	xlab="anti-RBD IgA (RAU)", 
	ylab="reciprocal dilution (neutralization IC50)", 
	main="(E) S Fuse vs anti-RBD IgA",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_y))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_x))
{
	points(x=rep(line_seq_dil_x[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points( x = IMM_db$RBD_IgA[which(IMM_db$site=="IMM_r2")], 
        y = IMM_db$NT_T2_titre[which(IMM_db$site=="IMM_r2")],
        pch=19, col="grey11", cex=point.size)


axis(1, at=c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01), 
        label=c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01),
        cex.axis=axis.size )

axis(2, at=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000), 
        label=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000),
        las=2, cex.axis=axis.size )



dev.off()



quantile( decay_slope, prob=c(0.5, 0.25, 0.75), na.rm=TRUE )*365

exp(quantile( decay_slope, prob=c(0.5, 0.25, 0.75), na.rm=TRUE )*365)







