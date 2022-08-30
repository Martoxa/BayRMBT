# The following code provides an example of how to use the BayRMBT lakes
# calibration to calibrate MBT5Me data. The example data is from:
#
# Miller et al. (2018) A 900-year New England temperature reconstruction
# from in situ seasonally produced branched glycerol dialkyl glycerol
# tetraethers (brGDGTs), Climate of the Past 
# https://doi.org/10.5194/cp-14-1653-2018
#
# These data also appear in Figure 9 of the BayMBT paper:
#
# Martínez-Sosa, Pablo, Tierney, J. E., Stefanescu, I. C., Dearing
# Crampton-Flood, E., Shuman, B. N. & Routson, C (2021). BayMBT: A global
# Bayesian temperature calibration for lacustrine brGDGTs, Geochimica et
# Cosmochimica Acta, in revision.

library(ggplot2)
#load the bayRmbt_predict function
source("bayRmbt_predict.R") 

#read in data
lakeData<-read.csv(file('BasinPondMiller.csv'),head=TRUE,sep=',')
#pull out mbt5me for ease
mbt5me<-lakeData$MBT5Me
#set prior mean, standard deviation, and model type
prior_mean<-10
prior_std<-10
model<-"T0"

#calibrate both uncorrected and corrected mbt
calibratedData<-bayrmbt_predict(mbt5me,prior_mean,prior_std,Tmodel = model,Type = "lake")

#get all the necessary values in a single data frame for ease
Data_final<-as.data.frame(cbind(lakeData$Year,calibratedData$T))
colnames(Data_final)[1]<-c("Year")

#plot the median values and 1-sigma confidence intervals
ggplot(data=Data_final,aes(x=Year,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+coord_cartesian(xlim=c(1110,2015),ylim=c(2,20),expand = FALSE)+
  geom_ribbon(alpha=0.3,fill="#8da0cb")+
  geom_line(size=1.2,color="#fc8d62")+
  theme_bw()+
  theme(panel.border=element_rect(size=1,fill=NA,colour = "black"),axis.text=element_text(colour="black"),text=element_text(face="bold",size=12))+
  ggtitle("Basin Pond MAF reconstruction (1 sigma CI)")+labs(x="Year",y="Temperature (C)")
