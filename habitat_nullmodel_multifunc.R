library(plyr)
library(multifunc)
library(reshape2)
library(ggplot2)
library(GGally)
library(gridExtra)


source("functions.R")



mono.habitat.func <- read.table("habitat_mono_function.txt", sep="\t", header=T)

# take average of replicates

mono.long <- melt(mono.habitat.func, id.vars=c("habitat", "season"))
mono.mean <- ddply(mono.long, .(habitat,season,variable), 
                   summarize, mean=mean(value))

#colnames to match the requirements for the funcitons below
colnames(mono.mean)[c(1,3,4)] <- c("Species", "Functions", "Funcval")


######### some graphs ######

ggplot(mono.mean, aes(x=Species, y=Funcval, colour=season, shape=season))+
  facet_wrap(~Functions, scales="free")+
  geom_point(size=3)+
  theme_bw(base_size=15)
  
         


################################################################################

#library(devtools)
#install_github("multifunc", "jebyrnes")




#### habitats = "species"###
spec <- unique(as.character(mono.mean$Species)) 
specnum <- length(spec)

## subset function per season ##
FuncMat_autumn <- mono.mean[mono.mean$season == "autumn",c(1,3,4)] 
FuncMat_summer <- mono.mean[mono.mean$season == "summer",c(1,3,4)] 
FuncMat_spring <-   mono.mean[mono.mean$season == "spring",c(1,3,4)] 

funcnum <- length(unique(mono.mean$Functions))

# set up species matrix with specnum species 

SpecMat <- SpeciesMatrix(spec = spec)

# calculate null.model ### choose your season ###
AvFunc <- AverageFunction(SpecMat, FuncMat_spring , method = "av")


 #################
# In the following I replicate the analysis steps described 
# In:
# Supplementary Information 1: Using the multifunc package for analysis of 
# Biodiversity Ecosystem Multifunctionality Relationships in R
# from: 
# Byrnes et al 2014; Investigating the relationship between biodiversity and 
# ecosystem multifunctionality: challenges and solutions

# The analysis above uses the german BIODEPTH data. 
# For comarision, I replicate all performed here also with the original data
# see Biodepth_example_JByrnes.R

# If I have to deviate from the analysis I specify it in the script
###################


######################
# THRESHOLD APPROACH #
######################

# extract function names
func.names <- colnames( AvFunc[ ( specnum + 3) : ncol( AvFunc)])

# add on the new (standardized) functions along with the averaged multifunctional index
AvFunc <- cbind(AvFunc, getStdAndMeanFunctions(AvFunc, func.names))


mixedThresh <- getFuncsMaxed(AvFunc, func.names, threshmin=0.05, threshmax=0.99, 
                             prepend=c("plot","Richness"), maxN=7)

gcPlot_mixed <- subset(mixedThresh, mixedThresh$thresholds %in% qw(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
gcPlot_mixed$percent <- paste(100*gcPlot_mixed$thresholds, "%", sep="")

# note: Jerret uses GLMs to specify the best fit lines. As they don't converge in all cases for this dataframe,
# I use just lm(). For the cases where the glm converged, the differences seem to be minimal.

# plot 4 selected threshold values
ggplot(gcPlot_mixed, aes(x=Richness, y=funcMaxed))+
  geom_point()+
  facet_wrap(~percent)+
  stat_smooth(method="glm",colour="red", lwd=1.2 )+
  labs(y = expression("Number of Functions" >= Threshold), x = ("Species Richness"))+
  theme_bw(base_size=15)


# plot N Function over Threshold ~ Richness Slopes for all Thresholds
mixedThresh$percent <- 100*mixedThresh$thresholds 

ggplot(data=mixedThresh, aes(x=Richness, y=funcMaxed, group=percent)) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("Species Richness") +
  stat_smooth(method="lm", lwd=0.8, fill=NA, aes(color=percent))+
  theme_bw(base_size=14) +
  scale_y_continuous(limits=c(0,max(mixedThresh$funcMaxed)))+
  geom_hline(aes(yintercept = max(mixedThresh$funcMaxed)),size=1)+
  geom_hline(aes(yintercept = 0), size=1)+
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red")

# plot the slopes of the relationship against Threshold values

# note that I take fun = lm as glm doesn't converge for all thresholds over 67%
# however, the parts that converge look almost identical so it's no big deal

mixedLinearSlopes<-getCoefTab(funcMaxed ~ Richness, fun = lm,  data=mixedThresh, 
                              coefVar="Richness")

colnames(mixedLinearSlopes) <- c("thresholds", "Estimate",  "Std. Error", "t value", "Pr(>|t|)")


ggplot(mixedLinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate- 1.96*mixedLinearSlopes[["Std. Error"]],
                                 ymax=Estimate+1.96*mixedLinearSlopes[["Std. Error"]])) + 
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Change in Number of Functions per Addition of 1 Species\n") + xlab("\nThreshold (%)") +
  stat_abline(intercept=0, slope=0, lwd=1, linetype=2) + 
  theme_bw(base_size=14)

# highlight critical slopes

# We can't calculate Tmin, Tmax and Tmde because the calculation rely on the fact that the slopes get non-significant at
# some point. Becasue these are simulated data with 100 replicates for each richness level, this never happens (all slopes
# are significantly differnt from 0 at all threshold levels) . 



