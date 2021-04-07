
## Set working directory
setwd("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/Flow_Cam")

## Load Needed packages
library("factoextra")
library("car")
library("vegan")

##------------------------------------------------------------------------------
## Data prepping for protist traits and abundances
##------------------------------------------------------------------------------

dat <- read.csv("Master_dataset_pr_Awesome.csv")
head(dat)
dat$Class <- factor(dat$Class,levels = c("Tetrahymena", "Colpidium"))
dat$Temp <- factor(dat$Temp,levels = c("Original","22","25"))
dat$Time <- factor(dat$Time,levels = c("Original","12","24"))
dat$Geodesic.Width <- dat$Geodesic.Length*dat$Geodesic.Aspect.Ratio

dat_noIN <- dat[-which(dat$Jar=="INIT"),] 

# INIT densities need to be corrected because the dataset reports the densities in the original cultures. We used 2mls of the original cultures at the densities described in the dataset. 
dat$Density_inds.ml[which(dat$Jar=="INIT" & dat$Class=="Colpidium" )] <- (311*2)/20  
dat$Density_inds.ml[which(dat$Jar=="INIT" & dat$Class=="Tetrahymena" )] <- (1244*2)/20

jars <- unique(dat$Jar)
dat_den <- 0
for(i in 1:length(jars)){
	if(i == 1){
		dat_tobind <- rbind(dat[which(dat$Jar==jars[i]),][1,],
		dat[which(dat$Jar==jars[i]),][283,])
		dat_den <- rbind(dat_den,dat_tobind)
	}else{
	dat_tobind <- dat[which(dat$Jar==jars[i]),][1,]
	dat_den <- rbind(dat_den,dat_tobind)
	}
}
dat_den <- dat_den[-1,]
dat_den <- dat_den[-which(dat_den$Density_inds.ml>1800),] # Identified as massive outlier and dropped (distance to median larger than 1.5 interquartile distance)

plot(Density_inds.ml[which(Class=="Colpidium")]~Time[which(Class=="Colpidium")], data=dat_den)
boxplot(Density_inds.ml[which(Class=="Tetrahymena")]~Time[which(Class=="Tetrahymena")], data=dat_den)

plot(Density_inds.ml[which(Class=="Colpidium")]~Temp[which(Class=="Colpidium")], data=dat_den)
plot(Density_inds.ml[which(Class=="Tetrahymena")]~Temp[which(Class=="Tetrahymena")], data=dat_den)

dat_den2 <- dat_den[-c(1,2),]
dat_den2$Temp <- droplevels(dat_den2$Temp)
dat_den2$Time <- droplevels(dat_den2$Time)

##--------------------------------------------------------------------------
## Reproducing Figs 3a and 3b of the main text
##--------------------------------------------------------------------------
#pdf("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/1_Manuscript/Figures/3rd_version/Fig_3/Fig_SX.pdf", height=5, width=12)
par(mfrow=c(1,2),oma=c(0,4,0,0),mar=c(4,1,1,1.2)) 
boxplot(Density_inds.ml~Time+Class, data=dat_den2, outline=FALSE, medcol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"),whiskcol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"),boxcol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"), staplecol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"), boxwex=0.6, axes=FALSE, ylab="", xlab="")
points(jitter(rep(1,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Time=="12")])),13.5), dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Time=="12")],col="#E69F00", pch=c(1,16)[dat_den2$Temp[which(dat_den2$Class=="Tetrahymena" & dat_den2$Time=="12")]])
points(jitter(rep(2,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Time=="24")])),7.2), dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Time=="24")],col="#E69F00", pch=c(2,17)[dat_den2$Temp[which(dat_den2$Class=="Tetrahymena" & dat_den2$Time=="24")]])
points(jitter(rep(3,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Time=="12")])),5), dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Time=="12")],col="#CC79A7", pch=c(1,16)[dat_den2$Temp[which(dat_den2$Class=="Colpidium" & dat_den2$Time=="12")]])
points(jitter(rep(4,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Time=="24")])),3.5), dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Time=="24")],col="#CC79A7", pch=c(2,17)[dat_den2$Temp[which(dat_den2$Class=="Colpidium" & dat_den2$Time=="24")]])
box(lwd=2, bty='l')
axis(1,at=seq(1,4,1), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(0,800,200), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('Treatment',1, line=2.2, cex=1.5)
mtext('Population density (inds/ml)',2,line=2.5, cex=1.5)

boxplot(Density_inds.ml~Temp+Class, data=dat_den2, outline=FALSE, medcol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"),whiskcol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"),boxcol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"), staplecol=c("#E69F00","#E69F00","#CC79A7","#CC79A7"), boxwex=0.6, axes=FALSE, ylab="", xlab="")
points(jitter(rep(1,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Temp=="22")])),13.5), dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Temp=="22")],col="#E69F00", pch=c(1,2)[dat_den2$Time[which(dat_den2$Class=="Tetrahymena" & dat_den2$Temp=="22")]])
points(jitter(rep(2,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Temp=="25")])),7.8), dat_den2$Density_inds.ml[which(dat_den2$Class=="Tetrahymena" & dat_den2$Temp=="25")],col="#E69F00", pch=c(16,17)[dat_den2$Time[which(dat_den2$Class=="Tetrahymena" & dat_den2$Temp=="25")]])
points(jitter(rep(3,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Temp=="22")])),5), dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Temp=="22")],col="#CC79A7", pch=c(1,2)[dat_den2$Time[which(dat_den2$Class=="Colpidium" & dat_den2$Temp=="22")]])
points(jitter(rep(4,length(dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Temp=="25" & dat_den2$Density_inds.ml<400)])),3.5), dat_den2$Density_inds.ml[which(dat_den2$Class=="Colpidium" & dat_den2$Temp=="25" & dat_den2$Density_inds.ml<400)],col="#CC79A7", pch=c(16,17)[dat_den2$Time[which(dat_den2$Class=="Tetrahymena" & dat_den2$Temp=="25")]])
box(lwd=2, bty='l')
axis(1,at=seq(1,4,1), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(0,800,200), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('Treatment',1, line=2.2, cex=1.5)
#mtext('Population density (inds/ml)',2,line=2.5, cex=1.5)
dev.off()

##--------------------------------------------------------------------------
## Analyzing response of abundances over Time and with Temp
##--------------------------------------------------------------------------
dat_denColp <- dat_den[which(dat_den$Class=="Colpidium" & dat_den$Temp!="Original"),]
dat_denColp$Temp <- droplevels(dat_denColp$Temp)
dat_denTetra <- dat_den[which(dat_den$Class=="Tetrahymena" & dat_den$Temp!="Original"),]
dat_denTetra$Temp <- droplevels(dat_denTetra$Temp)

# These are the analyses reported in the main text for Figs 3a and 3b 
summary(lm(Density_inds.ml[which(Class=="Colpidium")]~Temp[which(Class=="Colpidium")], data=dat_denColp))
summary(lm(Density_inds.ml[which(Class=="Tetrahymena")]~Temp[which(Class=="Tetrahymena")], data=dat_denTetra))
## Neither species abundances respond to temperature

dat_den$Temp <- droplevels(dat_den$Temp)
summary(lm(Density_inds.ml[which(Class=="Colpidium")]~Time[which(Class=="Colpidium")], data=dat_den))
summary(lm(Density_inds.ml[which(Class=="Tetrahymena")]~Time[which(Class=="Tetrahymena")], data=dat_den))
## But both respond to Time



##--------------------------------------------------------------------------
## Trait analyses I: PCAs
##--------------------------------------------------------------------------

## Run PCA
 ## For colpidium
dat_colp <- dat[which(dat$Class=="Colpidium" & dat$Temp!="Original"),]
## Subset the traits to those reported in the main text
dat_sub_colp <- dat_colp[,c(1,2,3,4,5,6,11,13,18,20,21,22,23)]
head(dat_sub_colp)

# Reproduces Table S3 of Appendix 4
dat_pca_colp <- prcomp(dat_sub_colp[c(5:13)], center = TRUE, scale = TRUE)
summary(dat_pca_colp)

# Reproduces Figs S8 of Appendix 4
par(mfrow=c(1,2))
screeplot(dat_pca_colp, type = "l", npcs = 9, main = "Screeplot of the first 9 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(dat_pca_colp$sdev^2 / sum(dat_pca_colp$sdev^2))
plot(cumpro[0:9], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 3, col="blue", lty=5)
abline(h = 0.867, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC3"),
       col=c("blue"), lty=5, cex=0.6)

 ## For tetrahymena
dat_tetra <- dat[which(dat$Class=="Tetrahymena" & dat$Temp!="Original"),]
dat_sub_tetra <- dat_tetra[,c(1,2,3,4,5,6,11,13,18,20,21,22,23)]
head(dat_sub_tetra)

# Reproduces Table S3 of Appendix 4
dat_pca_tetra <- prcomp(dat_sub_tetra[c(5:13)], center = TRUE, scale = TRUE)
summary(dat_pca_tetra)

# Reproduces Figs S8 of Appendix 4
par(mfrow=c(1,2))
screeplot(dat_pca_tetra, type = "l", npcs = 9, main = "Screeplot of the first 9 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(dat_pca_tetra$sdev^2 / sum(dat_pca_colp$sdev^2))
plot(cumpro[0:9], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 3, col="blue", lty=5)
abline(h = 0.8766, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC3"),
       col=c("blue"), lty=5, cex=0.6)

  
 #### ---------------------------------------
 ## Reproduces Figs 3c-f of the main text 

#### Changes with Temp and Time for Tetra and Colpidium
post_pca_colp <- cbind(dat_colp,PCA1=dat_pca_colp$x[,1],PCA2=dat_pca_colp$x[,2])
post_pca_colp$Temp <- droplevels(post_pca_colp$Temp)
post_pca_colp$Time <- droplevels(post_pca_colp$Time)
post_pca_tetra <- cbind(dat_tetra,PCA1=dat_pca_tetra$x[,1],PCA2=dat_pca_tetra$x[,2])
post_pca_tetra$Temp <- droplevels(post_pca_tetra$Temp)
post_pca_tetra$Time <- droplevels(post_pca_tetra$Time)

PCA1Time_colp_12 <- mean(post_pca_colp$PCA1[which(post_pca_colp$Time=="12")])
PCA1Time_colp_24 <- mean(post_pca_colp$PCA1[which(post_pca_colp$Time=="24")])
PCA2Time_colp_12 <- mean(post_pca_colp$PCA2[which(post_pca_colp$Time=="12")])
PCA2Time_colp_24 <- mean(post_pca_colp$PCA2[which(post_pca_colp$Time=="24")])

PCA1Temp_colp_22 <- mean(post_pca_colp$PCA1[which(post_pca_colp$Temp=="22")])
PCA1Temp_colp_25 <- mean(post_pca_colp$PCA1[which(post_pca_colp$Temp=="25")])
PCA2Temp_colp_22 <- mean(post_pca_colp$PCA2[which(post_pca_colp$Temp=="22")])
PCA2Temp_colp_25 <- mean(post_pca_colp$PCA2[which(post_pca_colp$Temp=="25")])

PCA1Time_tetra_12 <- mean(post_pca_tetra$PCA1[which(post_pca_tetra$Time=="12")])
PCA1Time_tetra_24 <- mean(post_pca_tetra$PCA1[which(post_pca_tetra$Time=="24")])
PCA2Time_tetra_12 <- mean(post_pca_tetra$PCA2[which(post_pca_tetra$Time=="12")])
PCA2Time_tetra_24 <- mean(post_pca_tetra$PCA2[which(post_pca_tetra$Time=="24")])

PCA1Temp_tetra_22 <- mean(post_pca_tetra$PCA1[which(post_pca_tetra$Temp=="22")])
PCA1Temp_tetra_25 <- mean(post_pca_tetra$PCA1[which(post_pca_tetra$Temp=="25")])
PCA2Temp_tetra_22 <- mean(post_pca_tetra$PCA2[which(post_pca_tetra$Temp=="22")])
PCA2Temp_tetra_25 <- mean(post_pca_tetra$PCA2[which(post_pca_tetra$Temp=="25")])

par(mfrow=c(2,2), mar=c(3, 2.2, 0.5, 0.5),oma=c(1.5,2,1,1))
plot(PCA2~PCA1,col=c("#0072B2","#D55E00")[post_pca_colp$Temp],data=post_pca_colp, xlab="", ylab = "", main = "", pch=16, cex=0.55, axes=FALSE, ylim=c(-5,5))
dataEllipse(post_pca_colp$PCA1,post_pca_colp$PCA2,groups=post_pca_colp$Temp, plot.points=FALSE,add=TRUE, group.labels="", center.pch=FALSE, col=c("#0072B2","#D55E00"), levels=c(0.95))
arrows(0,0,dat_pca_colp$rotation[,1]*max(dat_pca_colp$x[,1]),
dat_pca_colp$rotation[,2]*max(dat_pca_colp$x[,2]),length=0.1,
angle=20, col="black", lwd=1.5)
text(dat_pca_colp$rotation[,1]*max(dat_pca_colp$x[,1])*1.1, dat_pca_colp$rotation[,2]*max(dat_pca_colp$x[,2])*1.1,labels=rownames(dat_pca_colp$rotation), col="black",cex=1)
points(PCA2Temp_colp_22~PCA1Temp_colp_22, cex=2, pch = 21, bg="#0072B2" , col = "black")
points(PCA2Temp_colp_25~PCA1Temp_colp_25, cex=2, pch = 21, bg="#D55E00" , col = "black")
box(lwd=2, bty='l')
axis(1,at=seq(-6,6,3), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(-6,6,2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('PC1 (32.5%)',1, line=2.2, cex=1.2)
mtext('PC2 (30.3%)',2, line=2.2, cex=1.2)

plot(PCA2~PCA1,col=c("#CC79A7","grey")[post_pca_colp$Time],data=post_pca_colp, xlab="", ylab = "", main = "", pch=16, cex=0.55, axes=FALSE, ylim=c(-5,5))
dataEllipse(post_pca_colp$PCA1,post_pca_colp$PCA2,groups=post_pca_colp$Time, plot.points=FALSE,add=TRUE, group.labels="", center.pch=FALSE, col=c("#CC79A7","grey"), levels=c(0.95))
arrows(0,0,dat_pca_colp$rotation[,1]*max(dat_pca_colp$x[,1]),
dat_pca_colp$rotation[,2]*max(dat_pca_colp$x[,2]),length=0.1,
angle=20, col="black", lwd=1.5)
text(dat_pca_colp$rotation[,1]*max(dat_pca_colp$x[,1])*1.1, dat_pca_colp$rotation[,2]*max(dat_pca_colp$x[,2])*1.1,labels=rownames(dat_pca_colp$rotation), col="black",cex=1)
points(PCA2Time_colp_12~PCA1Time_colp_12, cex=2, pch = 21, bg="#CC79A7" , col = "black")
points(PCA2Time_colp_24~PCA1Time_colp_24, cex=2, pch = 21, bg="grey" , col = "black")
box(lwd=2, bty='l')
axis(1,at=seq(-6,6,3), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(-6,6,2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('PC1 (32.5%)',1, line=2.2, cex=1.2)
#mtext('PC2 (30.3%)',2, line=2.2, cex=1.2)

plot(PCA2~PCA1,col=c("#0072B2","#D55E00")[post_pca_tetra$Temp],data=post_pca_tetra, xlab="", ylab = "", main = "", pch=16, cex=0.55, axes=FALSE)
dataEllipse(post_pca_tetra$PCA1,post_pca_tetra$PCA2,groups=post_pca_tetra$Temp, plot.points=FALSE,add=TRUE, group.labels="", center.pch=FALSE, col=c("#0072B2","#D55E00"), levels=c(0.95))
arrows(0,0,dat_pca_tetra$rotation[,1]*max(dat_pca_tetra$x[,1]),
dat_pca_tetra$rotation[,2]*max(dat_pca_tetra$x[,2]),length=0.1,
angle=20, col="black", lwd=1.5)
text(dat_pca_tetra$rotation[,1]*max(dat_pca_tetra$x[,1])*1.1, dat_pca_tetra$rotation[,2]*max(dat_pca_tetra$x[,2])*1.1,labels=rownames(dat_pca_tetra$rotation), col="black",cex=1)
points(PCA2Temp_tetra_22~PCA1Temp_tetra_22, cex=2, pch = 21, bg="#0072B2" , col = "black")
points(PCA2Temp_tetra_25~PCA1Temp_tetra_25, cex=2, pch = 21, bg="#D55E00" , col = "black")
box(lwd=2, bty='l')
axis(1,at=seq(-6,6,3), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(-6,6,2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('PC1 (38%)',1, line=2.2, cex=1.2)
mtext('PC2 (30%)',2, line=2.2, cex=1.2)

plot(PCA2~PCA1,col=c("#E69F00","grey")[post_pca_tetra$Time],data=post_pca_tetra, xlab="", ylab = "", main = "", pch=16, cex=0.55, axes=FALSE)
dataEllipse(post_pca_tetra$PCA1,post_pca_tetra$PCA2,groups=post_pca_tetra$Time, plot.points=FALSE,add=TRUE, group.labels="", center.pch=FALSE, col=c("#E69F00","grey"), levels=c(0.95))
arrows(0,0,dat_pca_tetra$rotation[,1]*max(dat_pca_tetra$x[,1]),
dat_pca_tetra$rotation[,2]*max(dat_pca_tetra$x[,2]),length=0.1,
angle=20, col="black", lwd=1.5)
text(dat_pca_tetra$rotation[,1]*max(dat_pca_tetra$x[,1])*1.1, dat_pca_tetra$rotation[,2]*max(dat_pca_tetra$x[,2])*1.1,labels=rownames(dat_pca_tetra$rotation), col="black",cex=1)
points(PCA2Time_tetra_12~PCA1Time_tetra_12, cex=2, pch = 21, bg="#E69F00" , col = "black")
points(PCA2Time_tetra_24~PCA1Time_tetra_24, cex=2, pch = 21, bg="grey" , col = "black")
box(lwd=2, bty='l')
axis(1,at=seq(-6,6,3), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(-6,6,2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('PC1 (38%)',1, line=2.2, cex=1.2)
#mtext('PC2 (30%)',2, line=2.2, cex=1.5)

##--------------------------------------------------------------------------
## Trait analyses II: perMANOVA
##--------------------------------------------------------------------------

## PERMANOVA results reproduce Tables S4-S7 of Appendix 4 and the results reported in the main text
PC1_dist_colp <- dist(post_pca_colp[,29])
PC2_dist_colp <- dist(post_pca_colp[,30])
# Everything super significant
mod_PC1_colp <- adonis(PC1_dist_colp~Time*Temp, data=post_pca_colp, permut=99)
mod_PC2_colp <- adonis(PC2_dist_colp~Time*Temp, data=post_pca_colp, permut=99)

PC1_dist_tetra <- dist(post_pca_tetra[,29])
PC2_dist_tetra <- dist(post_pca_tetra[,30])
mod_PC1_tetra <- adonis(PC1_dist_tetra~Time*Temp, data=post_pca_tetra, permut=99)
mod_PC2_tetra <- adonis(PC2_dist_tetra~Time*Temp, data=post_pca_tetra, permut=99)

#simple lms
modlm_PC1_colp <- lm(PCA1~Time*Temp, data=post_pca_colp)
modlm_PC2_colp <- lm(PCA2~Time*Temp, data=post_pca_colp)

summary(modlm_PC1_colp)
summary(modlm_PC2_colp) 

modlm_PC1_tetra <- lm(PCA1~Time*Temp, data=post_pca_tetra)
modlm_PC2_tetra <- lm(PCA2~Temp*Temp, data=post_pca_tetra)

summary(modlm_PC1_tetra)
summary(modlm_PC2_tetra) 
 
## THE END