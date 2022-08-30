# The following code provides an example of how to read in fractional
# brGDGT abundances and apply BayRMBT. The example data is from:
#
# Crampton-Flood, E. D., Peterse, F., Munsterman, D., & Damsté, J. S. S.
# (2018). Using tetraether lipids archived in North Sea Basin sediments to
# extract North Western European Pliocene continental air temperatures.
# Earth and Planetary Science Letters, 490, 193-205.
#
# These are Pliocene GDGT data from the Hank core. Crampton-Flood describe
# a way to correct brGDGT data for marine contributions, an example of that
# correction is also given here.
#
# These data also appear in Figure 9 of the BayMBT paper:
#
# Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
# & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
# for branched glycerol dialkyl glycerol tetraethers in soils and peats.
# Geochimica et Cosmochimica Acta, 268, 142-159.

#load functions to calculate GDGT relevant indices (https://github.com/Martoxa/GDGTindices)
source("Indices.R") 
#load the bayRmbt_predict function
source("bayRmbt_predict.R")

#just some packages to make pretty plots
library(ggplot2)
library(cowplot)

#read in data
hank<-read.csv(file("Hank_core.csv"),head=TRUE,sep=',')
#change names to run the MBT5 function
colnames(hank)<-c("Depth","Ia","Ib","Ic","IIa","IIa'","IIb","IIb'","IIc","IIc'","IIIa","IIIa'","IIIb","IIIb'","IIIc","IIIc'","percentTerrestrial")

#calculate MBT'5Me
Hmbt5me<-MBT5(hank)

#adjust for marine contribution according to Crampton-Flood et al., (2018)
bwt<-14.5
mbtmarine<-(bwt+23.7)/59.5
mbtterr<- (Hmbt5me - mbtmarine*(1-hank$percentTerrestrial/100))/(hank$percentTerrestrial/100)

#set prior mean, standard deviation, and model type
prior_mean<-10
prior_std<-10
model<-"T0"

#calibrate both uncorrected and corrected MBT'5Me
uncorrected<-bayrmbt_predict(Hmbt5me,prior_mean,prior_std,Tmodel=model,Type = "soil")
corrected<-bayrmbt_predict(mbtterr,prior_mean,prior_std,Tmodel=model,Type = "soil")

#arrange the data in data frames to plot
UnF<-as.data.frame(cbind(hank$Depth,uncorrected$T))
colnames(UnF)[1]<-"Depth"
CF<-as.data.frame(cbind(hank$Depth,corrected$T))
colnames(CF)[1]<-"Depth"
All<-as.data.frame(cbind(hank$Depth,uncorrected$T[,2],corrected$T[,2]))

#plot each reconstruction with their 1-sigma CI
p1<-ggplot(data=UnF,aes(x=Depth,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+coord_cartesian(xlim=c(130,410),ylim=c(-1,23),expand = FALSE)+
  geom_ribbon(alpha=0.3,fill="#a6dba0")+
  geom_line(size=1.2,color="#008837")+
  theme_bw()+
  theme(panel.border=element_rect(size=1,fill=NA,colour = "black"),axis.text=element_text(colour="black"),text=element_text(face="bold",size=12))+
  ggtitle("Uncorrected MAF reconstruction (1 sigma CI)")+labs(x="Depth (m)",y="Temperature (C)")

p2<-ggplot(data=CF,aes(x=Depth,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+coord_cartesian(xlim=c(130,410),ylim=c(-1,23),expand = FALSE)+
  geom_ribbon(alpha=0.3,fill="#c2a5cf")+
  geom_line(size=1.2,color="#7b3294")+
  theme_bw()+
  theme(panel.border=element_rect(size=1,fill=NA,colour = "black"),axis.text=element_text(colour="black"),text=element_text(face="bold",size=12))+
  ggtitle("Corrected MAF reconstruction (1 sigma CI)")+labs(x="Depth (m)",y="Temperature (C)")

plot_grid(p1,p2,ncol=1)