# Load packages, set working directory
library("readxl")
library("stringr")

##--------------------------------------------------------------
# Load and prep respirometry RAW data (found in https://github.com/JPGibert/Microbial_Responses_Prots_Temp/tree/master/Raw%20Data/T12_data and T24_data)
##--------------------------------------------------------------

## These two loops both load and prep the data, and produce Figs S1 and S2 of the Appendix.
# Data needs to be separated in the two folders T12_data and T24_data for this code to run

setwd("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/Respirometry/T12_data")
names_12 <- dir()
respirometry_12 <- rep(0,length(names))
Jar_12 <- sub(".+_","", names_12)
par(mfrow=c(9,8), mar = c(0, 0, 0, 0), oma = c(3.5, 3.5, 1, 1), tck=0.025,cex.axis=0.8)
for(i in 1:length(names_12)){
	# Load data
	data <- read_excel(names_12[i], sheet=6, col_names=TRUE)
	# Keep data after initial 5 mins
	data <- data[which(data[,3]>5),]
	respirometry_12[i] <- abs(summary(lm(t(t(data[,7]))~t(t(data[,3]))))$coefficients[2])
	## Plotting Fig S1
	if(i %in% c(1,1+8*c(seq(1,7,1))) ){
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(25,300),col="grey",xaxt='n', yaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
		axis(2,at=c(50,150,250), padj=2)
	}else if(i %in% c(9*8-seq(0,6,1))){
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(25,300),col="grey",yaxt='n',xaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
		axis(1,at=c(5,15,25), padj=-2)
	}else if(i == 65){
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(25,300),col="grey",yaxt='n',xaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
		axis(2,at=c(50,150,250), padj=2)
		axis(1,at=c(5,15,25), padj=-2)
	}else{
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(25,300),col="grey",xaxt='n',yaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
	}	
}
mtext('Time (min.)',1, line=2, outer=TRUE)
mtext('Oxygen concentration (micro mol/L)',2, line=2, outer=TRUE)

## This preps the respirometry data for Day 24. 
setwd("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/Respirometry/T24_data")
names_24 <- dir()
respirometry_24 <- rep(0,length(names))
Jar_24 <- str_sub(names_24,start=-11)
par(mfrow=c(9,8), mar = c(0, 0, 0, 0), oma = c(3.5, 3.5, 1, 1), tck=0.025,cex.axis=0.8)
for(i in 1:length(names_24)){
	# Load data
	data <- read_excel(names_24[i], sheet=6, col_names=TRUE)
	# Keep data after initial 5 mins
	data <- data[which(data[,3]>5),]
	respirometry_24[i] <- abs(summary(lm(t(t(data[,7]))~t(t(data[,3]))))$coefficients[2])
	## Plotting Fig S2
	if(i %in% c(1,1+8*c(seq(1,7,1))) ){
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(0,300),col="grey",xaxt='n', yaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
		axis(2,at=c(50,150,250), padj=2)
	}else if(i %in% c(9*8-seq(0,6,1))){
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(0,300),col="grey",yaxt='n',xaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
		axis(1,at=c(5,15,25), padj=-2)
	}else if(i == 65){
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(0,300),col="grey",yaxt='n',xaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
		axis(2,at=c(50,150,250), padj=2)
		axis(1,at=c(5,15,25), padj=-2)
	}else{
		plot(t(t(data[,7]))~t(t(data[,3])),ylab="",xlab="",cex=0.5,ylim=c(0,300),col="grey",xaxt='n',yaxt='n')
		abline(lm(t(t(data[,7]))~t(t(data[,3]))), col="red", lwd=1)
	}	
}
mtext('Time (min.)',1, line=2, outer=TRUE)
mtext('Oxygen concentration (micro mol/L)',2, line=2, outer=TRUE)


##-------------------------------------------------------------------------------------
### This piece of code preps the respiration data by putting it in a data_frame that can be used to do stats.
##-------------------------------------------------------------------------------------

## Set treatments 
	# Colp, Tetra or Control?
Treat <- rep("na",72)
Treat[grep("Con",names_12, value=FALSE)] <- "Control"
Treat[grep("Temp",names_12, value=FALSE)] <- "Control"
Treat[grep("Tet",names_12, value=FALSE)] <- "Tetra"
Treat[grep("Col",names_12, value=FALSE)] <- "Colp"

Treat_24 <- rep("na",72)
Treat_24[grep("Con",names_24, value=FALSE)] <- "Control"
Treat_24[grep("Tet",names_24, value=FALSE)] <- "Tetra"
Treat_24[grep("Col",names_24, value=FALSE)] <- "Colp"
Treat <- rbind(t(t(Treat)),t(t(Treat_24)))

	# Temp or not? 
Temp <- rep("noT",72)
Temp[grep("Temp",names_12, value=FALSE)] <- "+3"

Temp_24 <- rep("noT",72)
Temp_24[grep("_T_",names_24, value=FALSE)] <- "+3"
Temp <- rbind(t(t(Temp)),t(t(Temp_24)))

	# Measurement Temp 
Temp_m <- rep("22",72)
Temp_m[grep("25",names_12, value=FALSE)] <- "25"

Temp_m_24 <- rep("22",72)
Temp_m_24[grep("25",names_12, value=FALSE)] <- "25"
Temp_m <- rbind(t(t(Temp_m)),t(t(Temp_m_24)))


	# Put together dataset
Time <- rep(c(12,24),each=72)

resp_dat <- data.frame(cbind(Time,rbind(t(t(respirometry_12)), t(t(respirometry_24)))))
resp_dat <- cbind(resp_dat,rbind(t(t(Jar_12)),t(t(Jar_24))))
resp_dat$Treat <- factor(Treat,levels = c("Control", "Tetra", "Colp"))
resp_dat$Temp <- factor(Temp,levels = c("noT","+3"))
resp_dat$Temp_m <- factor(Temp_m)
resp_dat$Time <- factor(resp_dat$Time)

colnames(resp_dat) <- c("Time","Resp", "Jar","Treat","Temp","Temp_m")

# Add Jar names to dataset
tojar <- c(116,124,17,27,7,70,116,124,17,27,7,70,120,3,41,57,72,75,120,3,41,57,72,75,10,33,50,61,64,69,10,33,50,61,64,69,102,31,102,31,128,2,8,81,128,2,8,81,103,109,123,34,87,99,103,109,123,34,87,99,18,35,63,80,82,97,18,35,63,80,82,97,04,38,04,55,65,77,04,38,04,55,65,77,27,30,48,66,78,84,27,30,48,66,78,84,35,5,51,60,08,35,5,51,60,08,58,58,29,39,29,39,11,76,85,90,11,76,85,90,15,47,15,47,12,20,12,20,01,7,56,95,1,7,56,95,18,12,26,83,18,13,26,83)

resp_dat$Jar2 <- resp_dat$Jar
resp_dat$Jar2 <- tojar

##-------------------------------------------------------------
## The code below does preliminary data analysis and model selection
##-------------------------------------------------------------
library("MuMIn")
full.model <- lm(Resp ~ Treat*Time*Temp*Temp_m,data=resp_dat, na.action = "na.fail")       
dd <- dredge(full.model)
# When doing Treat*Time*Temp*Temp_mm Temp_m isn't in the best model, so the temperature at which respirometry was measured can be dropped can get dropped from model.
imp_dd <- importance(dd)
## This line reproduces *Fig S4* of the appendix
barplot(t(imp_dd), main="Resp", ylab="Proportion od models containing each variable")

# Run again without Temp_m 
full.model <- lm(Resp ~ Treat*Time*Temp,data=resp_dat, na.action = "na.fail")       
dd <- dredge(full.model)

mods <- get.models(dd, subset=TRUE)

imp_dd <- importance(dd)
barplot(t(imp_dd), main="Resp", ylab="Proportion od models containing each variable")
## THREE double interactions, only the triple one is less important.

##----------------------------------------------------------------------
## Data Analysis
##----------------------------------------------------------------------

## Since Time_m has weak effect, we average respiration rates for the two measurements.
## By time point (Taking the average of the two resp measurements for each jar)

##Day 12
data_12 <- resp_dat[which(resp_dat$Time==12),]
uni_Jar_12 <- unique(Jar_12)
mean_resp <- rep(0,length(uni_Jar_12))
data_12_final <- cbind(data_12[!duplicated(data_12$Jar),],mean_resp) 
for(i in 1:length(uni_Jar_12)){
	#mean_resp_to[i] <- mean(data_12$Resp[which(data_12$Jar == uni_Jar_12[i])])
	data_12_final$mean_resp[which(data_12_final$Jar == uni_Jar_12[i])] <- mean(data_12$Resp[which(data_12$Jar == uni_Jar_12[i])])
}
data_12 <- data_12_final

## Day 24
data_24 <- resp_dat[which(resp_dat$Time==24),]
uni_Jar_24 <- unique(Jar_24)
mean_resp <- rep(0,length(uni_Jar_24))
data_24_final <- cbind(data_24[!duplicated(data_24$Jar),],mean_resp) 
for(i in 1:length(uni_Jar_24)){
	#mean_resp_to[i] <- mean(data_24$Resp[which(data_24$Jar == uni_Jar_24[i])])
	data_24_final$mean_resp[which(data_24_final$Jar == uni_Jar_24[i])] <- mean(data_24$Resp[which(data_24$Jar == uni_Jar_24[i])])
}
data_24 <- data_24_final


## All together
data_tot <- rbind(data_12,data_24)
# This model is as reported in the main text
summary(mod_TOT <- lm(mean_resp~Treat*Temp+Temp*Time+Treat*Time,data=data_tot))

yellow "#E69F00",
green "#009E73", 
blue "#0072B2",
pink "#CC79A7",
red "#D55E00"



##----------------------------------------------------------------------------
## Reproducing Figs 1b and 2b of the main text
##----------------------------------------------------------------------------

#pdf("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/1_Manuscript/Figures/3rd_version/Fig_2/Fig2_b.pdf",height=6,width=7.38)
par(mfrow=c(1,2))
boxplot(mean_resp ~ Treat+Temp,data=data_tot, outline=FALSE, col="white",medcol=c("#E69F00","#CC79A7","#009E73")[unique(data_tot$Treat)],whiskcol=c("#E69F00","#CC79A7","#009E73")[unique(data_tot$Treat)],boxcol=c("#E69F00","#CC79A7","#009E73")[unique(data_tot$Treat)], staplecol=c("#E69F00","#CC79A7","#009E73")[unique(data_tot$Treat)], boxwex=0.6,axes=FALSE, ylab="",xlab="",range=2) ## points twice interquartile distance are outliers
points(jitter(rep(0.85,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Control" & data_tot$Time=="12")])),2), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Control" & data_tot$Time=="12")],col="#009E73", pch=1)
points(jitter(rep(1.15,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Control" & data_tot$Time=="24")])),2), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Control" & data_tot$Time=="24")],col="#009E73", pch=2)
points(jitter(rep(1.85,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Tetra" & data_tot$Time=="12")])),1), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Tetra" & data_tot$Time=="12")],col="#E69F00", pch=1)
points(jitter(rep(2.15,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Tetra" & data_tot$Time=="24")])),1), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Tetra" & data_tot$Time=="24")],col="#E69F00", pch=2)
points(jitter(rep(2.85,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Colp" & data_tot$Time=="12")])),0.8), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Colp" & data_tot$Time=="12")],col="#CC79A7", pch=1)
points(jitter(rep(3.15,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Colp" & data_tot$Time=="24")])),0.8), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Treat=="Colp" & data_tot$Time=="24")],col="#CC79A7", pch=2)
points(jitter(rep(3.85,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Control" & data_tot$Time=="12")])),0.5), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Control" & data_tot$Time=="12")],col="#009E73", pch=16)
points(jitter(rep(4.15,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Control" & data_tot$Time=="24")])),0.5), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Control" & data_tot$Time=="24")],col="#009E73", pch=17)
points(jitter(rep(4.85,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Tetra" & data_tot$Time=="12")])),0.25), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Tetra" & data_tot$Time=="12")],col="#E69F00", pch=16)
points(jitter(rep(5.15,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Tetra" & data_tot$Time=="24")])),0.25), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Tetra" & data_tot$Time=="24")],col="#E69F00", pch=17)
points(jitter(rep(5.85,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Colp" & data_tot$Time=="12")])),0.25), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Colp" & data_tot$Time=="12")],col="#CC79A7", pch=16)
points(jitter(rep(6.15,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Colp" & data_tot$Time=="24")])),0.25), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Treat=="Colp" & data_tot$Time=="24")],col="#CC79A7", pch=17)
box(lwd=2, bty='l')
axis(1,at=seq(1,6,1), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(0,6,1), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('Oxygen consumption rate',2, line=2.5, cex=1.5)
#dev.off()

#pdf("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/1_Manuscript/Figures/3rd_version/Fig_1/Fig1_b.pdf",height=6,width=6.38)
boxplot(mean_resp ~ Time+Temp,data=data_tot, outline=FALSE, medcol=c("#0072B2","#D55E00")[unique(data_tot$Temp)],whiskcol=c("#0072B2","#D55E00")[unique(data_tot$Temp)],boxcol=c("#0072B2","#D55E00")[unique(data_tot$Temp)], staplecol=c("#0072B2","#D55E00")[unique(data_tot$Temp)], boxwex=0.6,axes=FALSE, ylab="",xlab="",range=2)
points(jitter(rep(1,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Time=="12")])),13), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Time=="12")],col="#0072B2", pch=1)
points(jitter(rep(2,length(data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Time=="24")])),6.5), data_tot$mean_resp[which(data_tot$Temp=="noT" & data_tot$Time=="24")],col="#D55E00", pch=16)
points(jitter(rep(3,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Time=="12")])),4.5), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Time=="12")],col="#0072B2", pch=2)
points(jitter(rep(4,length(data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Time=="24")])),3.2), data_tot$mean_resp[which(data_tot$Temp=="+3" & data_tot$Time=="24")],col="#D55E00", pch=17)
box(lwd=2, bty='l')
axis(1,at=seq(1,6,1), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(0,6,1), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('Oxygen consumption rate',2, line=2.5, cex=1.5)
#dev.off()







### THE END