# Set working directory
setwd("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/Spectrophotometry")
#Load data file
spec_dat <- read.csv("Spec.csv")

# Prep data for analysis
spec_dat$Spec <- spec_dat$Spec*3 # Because samples had to be diluted 
spec_dat$Treat <- factor(spec_dat$Treat,levels = c("Control", "Tetra", "Colp"))
spec_dat$Temp <- factor(spec_dat$Temp,levels = c("noT","3"))
spec_dat$Time <- factor(spec_dat$Time,levels = c("12","24"))

## Preliminary data analysis and model selection
#install.packages("MuMIn") # Install if missing
library("MuMIn")
# Start from full model
full.model <- lm(Spec ~ Treat*Time*Temp,data=spec_dat, na.action = "na.fail")       
dd <- dredge(full.model) # Run dd for Appendix Table S1
mods <- get.models(dd, subset=TRUE)
imp_dd <- importance(dd)

# Produce Appendix Fig S3
barplot(t(imp_dd), main="OD", ylab="Proportion od models containing each variable")

## Understanding the details
summary(mod <- lm(Spec~Treat+Temp*Time,data=spec_dat)) # Additive effect of Treatment and the Time*Temperature interaction (Can be plotted separately)

## Treat effects
yellow "#E69F00",
green "#009E73", 
blue "#0072B2",
pink "#CC79A7",
red "#D55E00"

##-----------------------------------------------------------------------------
#### Produce Figs 1a and 2a of the main text 
par(mfrow=c(1,2),oma=c(0,4,0,0),mar=c(4,1,1,1.2)) 

#pdf("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/1_Manuscript/Figures/3rd_version/Fig_2/Fig2_a.pdf",height=5.5,width=6)
boxplot(Spec~Treat,data=spec_dat, outline=FALSE, col="white", medcol=c("#009E73","#E69F00","#CC79A7")[unique(spec_dat$Treat)],whiskcol=c("#009E73","#E69F00","#CC79A7")[unique(spec_dat$Treat)],boxcol=c("#009E73","#E69F00","#CC79A7")[unique(spec_dat$Treat)], staplecol=c("#009E73","#E69F00","#CC79A7")[unique(spec_dat$Treat)], boxwex=0.5,axes=FALSE,range=2)
points(jitter(rep(0.85,length(spec_dat$Spec[which(spec_dat$Treat=="Control" & spec_dat$Time=="12")])),3), spec_dat$Spec[which(spec_dat$Treat=="Control"& spec_dat$Time=="12")],col="#009E73", pch=c(1,16)[spec_dat$Temp[spec_dat$Treat=="Control"]])
points(jitter(rep(1.15,length(spec_dat$Spec[which(spec_dat$Treat=="Control" & spec_dat$Time=="24")])),2), spec_dat$Spec[which(spec_dat$Treat=="Control"& spec_dat$Time=="24")],col="#009E73", pch=c(2,17)[spec_dat$Temp[spec_dat$Treat=="Control"]])
points(jitter(rep(1.85,length(spec_dat$Spec[which(spec_dat$Treat=="Tetra" & spec_dat$Time=="12")])),1), spec_dat$Spec[which(spec_dat$Treat=="Tetra"& spec_dat$Time=="12")],col="#E69F00", pch=c(1,16)[spec_dat$Temp[spec_dat$Treat=="Tetra"]])
points(jitter(rep(2.15,length(spec_dat$Spec[which(spec_dat$Treat=="Tetra" & spec_dat$Time=="24")])),0.9), spec_dat$Spec[which(spec_dat$Treat=="Tetra"& spec_dat$Time=="24")],col="#E69F00", pch=c(2,17)[spec_dat$Temp[spec_dat$Treat=="Tetra"]])
points(jitter(rep(2.85,length(spec_dat$Spec[which(spec_dat$Treat=="Colp" & spec_dat$Time=="12")])),0.5), spec_dat$Spec[which(spec_dat$Treat=="Colp"& spec_dat$Time=="12")],col="#CC79A7", pch=c(1,16)[spec_dat$Temp[spec_dat$Treat=="Colp"]])
points(jitter(rep(3.15,length(spec_dat$Spec[which(spec_dat$Treat=="Colp" & spec_dat$Time=="24")])),0.4), spec_dat$Spec[which(spec_dat$Treat=="Colp"& spec_dat$Time=="24")],col="#CC79A7", pch=c(2,17)[spec_dat$Temp[spec_dat$Treat=="Colp"]])
box(lwd=2, bty='l')
axis(1,at=seq(1,3,1), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(0,3,0.2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('Treatment',1, line=2.2, cex=1.5)
mtext('Bacterial Biomass (OD600)',2,line=2.5, cex=1.5)
#dev.off()


#pdf("~/Desktop/JP/Papers_in_progress/Jenny_Bact_Prot_Temp/1_Manuscript/Figures/3rd_version/Fig_1/Fig1_a.pdf",height=5.5,width=6)
boxplot(Spec~Time*Temp,data=spec_dat, outline=FALSE, medcol=c("#0072B2","#D55E00")[unique(spec_dat$Temp)],whiskcol=c("#0072B2","#D55E00")[unique(spec_dat$Temp)],boxcol=c("#0072B2","#D55E00")[unique(spec_dat$Temp)], staplecol=c("#0072B2","#D55E00")[unique(spec_dat$Temp)], boxwex=0.6,axes=FALSE, range=2)
points(jitter(rep(1,length(spec_dat$Spec[which(spec_dat$Temp=="noT" & spec_dat$Time=="12")])),13), spec_dat$Spec[which(spec_dat$Temp=="noT" & spec_dat$Time=="12")],col="#0072B2", pch=1)
points(jitter(rep(2,length(spec_dat$Spec[which(spec_dat$Temp=="3" & spec_dat$Time=="12")])),6.5), spec_dat$Spec[which(spec_dat$Temp=="3" & spec_dat$Time=="12")],col="#D55E00", pch=16)
points(jitter(rep(3,length(spec_dat$Spec[which(spec_dat$Temp=="noT" & spec_dat$Time=="24")])),4.5), spec_dat$Spec[which(spec_dat$Temp=="noT" & spec_dat$Time=="24")],col="#0072B2", pch=2)
points(jitter(rep(4,length(spec_dat$Spec[which(spec_dat$Temp=="3" & spec_dat$Time=="24")])),3.2), spec_dat$Spec[which(spec_dat$Temp=="3" & spec_dat$Time=="24")],col="#D55E00", pch=17)
box(lwd=2, bty='l')
axis(1,at=seq(1,4,1), tck=0.015, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(0,3,0.2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext('Treatment',1, line=2.2, cex=1.5)
#dev.off()

boxplot(Spec~Treat+Time,data=spec_dat)
boxplot(Spec~Time+Temp+Treat,data=spec_dat)

## Summary of results: 
#1) Colpidium significantly decreases total microbial biomass
#2) While Tetrahymena also decreases biomass, the effect is not significant and half the size of that of Colpidium
#3) On average, biomass increases over time, but the interactive effect with temperature suggests that at first, warmer treatments lead to more rapid increase in biomass, while longer term colder treatments have higher biomass. Whether this is a snapshot of a transient, or a longer term effect, we cannot tell.




## THE END