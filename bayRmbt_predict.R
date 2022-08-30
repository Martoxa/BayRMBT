bayrmbt_predict<-function(mbt5me,prior_mean,prior_std,Tmodel,Type){
  
  if(is.vector(mbt5me) == FALSE){mbt5me<-as.vector(mbt5me)}
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
  
  nd<-length(mbt5me)
  n_draws<-length(tau2_draws_final)
  
  alpha<-b_draws_final[,2]
  betaT<-b_draws_final[,1]
  sigmab<-sqrt(tau2_draws_final)
  
  pmu<- replicate(n_draws,rep(1,nd)*prior_mean)
  pinv_cov<-(replicate(n_draws,rep(prior_std,nd)))^(-2)
  
##Posterior calculation
  post_mean_num<- pinv_cov*pmu+(t(replicate(nd,sigmab)))^(-2)*t(replicate(nd,betaT))*(mbt5me-t(replicate(nd,alpha)))
  post_mean_den<- pinv_cov+(t(replicate(nd,betaT)))^2*(t(replicate(nd,sigmab)))^(-2)
  post_mean<- post_mean_num/post_mean_den
  post_sig<-sqrt(post_mean_den^(-1))
  output.ens<-post_mean+matrix(rnorm(nd*n_draws),nrow = nd,ncol = n_draws)*post_sig
  output.prior_mean<-prior_mean
  output.prior_std<-prior_std
  
  if(Tmodel == "T0"){
    output.ens[output.ens<0]<-0
  }
  
  output.T<-t(apply(output.ens,1,FUN = quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))
  output<-list(output.ens,output.prior_mean,output.prior_std,output.T)
  names(output)<-c("ens","prior_mean","prior_std","T")
  return(output)
  
}