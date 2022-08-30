bayrmbt_predict<-function(mbt5me,prior_mean,prior_std,Tmodel,Type){

# BAYMBT prediction model for MBT5Me measured in soils and peats.
# Predicts Mean annual air temperature or mean temperatures above zero.
# ----- Inputs -----
# mbt5me: A scalar or vector of mbt5me values (1 x N) or (N x 1)
#
# prior_mean: A scalar prior mean value of T in degrees C.
#
# prior_std: A scalar prior standard deviation value of T in degrees C.
#
# Tmodel: A string corresponding the model you want to use. Options are:
# "T" = Calculate mean annual air temperature (BayMBT, for soils only)
# "T0" = Calculate mean annual temperatures above zero (BayMBT0)
#
# Type: A string corresponding to the data type. Options are:
# "soil" = use the soil calibration
# "lake" = use the lake calibration
#  
# Example: bayRmbt_predict(mbt5me,10,10,Tmodel="T0",Type="lake") 
#
# ----- Outputs -----
#
# output is a list containing:
# output$prior_mean: User choice of prior mean
# output$prior_std: User choice of prior standard deviation
# output$T: 2.5%, 50%, and 97.5% confidence intervals of posterior SST
# output$ens: full ensemble of posterior SST (N x 1000);
#
# ----- Citation -----
# Please cite the following papers when using the BayMBT calibrations:
#
# For soils:
# Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
# & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
# for branched glycerol dialkyl glycerol tetraethers in soils and peats.
# Geochimica et Cosmochimica Acta, 268, 142-159.
#
# For lakes:
# Martínez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Dearing
# Crampton-Flood, E., Shuman, B. N., Routson, C. (2021) A global Bayesian
# temperature calibration for lacustrine brGDGTs. 
#Geochimica et Cosmochimica Acta, 305, 87-105.
#  
  
  #ensure vector
  if(is.vector(mbt5me) == FALSE){mbt5me<-as.vector(mbt5me)}
  #load appropriate model
  if(Type == "soil" & Tmodel == "T"){
    params<-readRDS("baymbt_params_soil.RDS")
  } else if(Type == "soil" & Tmodel == "T0"){
    params<-readRDS("baymbt0_params_soil.RDS")
  } else if(Type == "lake" & Tmodel == "T"){
    params<-readRDS("baymbt_params_lake.RDS")
  } else if(Type == "lake" & Tmodel == "T0"){
    params<-readRDS("baymbt0_params_lake.RDS")
  } else print("Type or TModel not recognized")
  
  b_draws_final<-params[[1]]
  tau2_draws_final<-params[[2]]
  
  #get dimensions of time series and draws
  nd<-length(mbt5me)
  n_draws<-length(tau2_draws_final)
  
  #parse parameters
  alpha<-b_draws_final[,2]
  betaT<-b_draws_final[,1]
  sigmab<-sqrt(tau2_draws_final)
  
  #prior mean and inverse covariance matrix
  pmu<- replicate(n_draws,rep(1,nd)*prior_mean)
  pinv_cov<-(replicate(n_draws,rep(prior_std,nd)))^(-2)
  
  #posterior calculation
  post_mean_num<- pinv_cov*pmu+(t(replicate(nd,sigmab)))^(-2)*t(replicate(nd,betaT))*(mbt5me-t(replicate(nd,alpha)))
  post_mean_den<- pinv_cov+(t(replicate(nd,betaT)))^2*(t(replicate(nd,sigmab)))^(-2)
  post_mean<- post_mean_num/post_mean_den
  post_sig<-sqrt(post_mean_den^(-1))
  output.ens<-post_mean+matrix(rnorm(nd*n_draws),nrow = nd,ncol = n_draws)*post_sig
  output.prior_mean<-prior_mean
  output.prior_std<-prior_std
  
  #if using BayMBT0, truncate at T<0
  if(Tmodel == "T0"){
    output.ens[output.ens<0]<-0
  }
  
  output.T<-t(apply(output.ens,1,FUN = quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))
  output<-list(output.ens,output.prior_mean,output.prior_std,output.T)
  names(output)<-c("ens","prior_mean","prior_std","T")
  return(output)
  
}