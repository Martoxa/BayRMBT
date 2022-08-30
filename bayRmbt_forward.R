bayrmbt_forward<-function(t,Tmodel,Type){
  
  if(is.vector(t) == FALSE){t<-as.vector(t)}
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
  
  Nobs<-length(t)
  
  alpha<-b_draws_final[,2]
  betaT<-b_draws_final[,1]
  sigmab<-sqrt(tau2_draws_final)
  
  mbt5me<-matrix(rnorm(Nobs*1000,mean=alpha+t%*%t(as.matrix(betaT)),sd=t(replicate(Nobs,sigmab))),ncol = 1000,byrow = TRUE)
  mbt5me[mbt5me<0]<-0
  mbt5me[mbt5me>1]<-1
  
  return(mbt5me)
}