bayrmbt_model<-function(Tmodel,Type){
  
  # BAYMBT calibration model for MBT5Me measured in soils and peats.
  # Computes slope, intercept, error variance terms using Bayesian linear
  # regression.
  #
  # ----- Inputs -----
  # Tmodel: A string corresponding the temperature model you want to calculate. Options are:
  # "T" = Calculate mean annual air temperature (BayMBT)
  # "T0" = Calculate mean annual temperatures above zero (BayMBT0)
  #
  # Type: A string corresponding to the data type. Options are:
  # "soil" = use the soil calibration data
  # "lake" = use the lake calibration data
  #
  # Example: bayRmbt_model(Tmodel="T",Type="lake")
  #
  # ----- Citation -----
  # This BayMBT calibration method is published here:
  #
  # Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
  # & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
  # for branched glycerol dialkyl glycerol tetraethers in soils and peats.
  # Geochimica et Cosmochimica Acta, 268, 142-159.
  
  
  ##Required packages
  if(!require(MASS)){cat("informative error message")}
  if(!require(Rlab)){cat("informative error message")}
  if(!require(ks)){cat("informative error message")}
  
  ##Select type of sample
  if(Type == "soil"){
    calibs<-read.csv(file("soil_calibration_data.csv"),head=TRUE,sep=',')
  } else if(Type == "lake"){
    calibs<-read.csv(file("lake_calibration_data.csv"),head=TRUE,sep=',')
    calibs<-calibs[calibs$Outliers==0,]
  } else
    print('Type not recognized - choices are "soil" and "lake"')
  
  
  ##Select temperature model and set up for Bayesian regression of form Y=X*B + e    
  if(Tmodel == "T"){
    X<-as.matrix(cbind(calibs$MAT,rep(1,dim(calibs)[1])))
  } else if(Tmodel == "T0"){
    X<-as.matrix(cbind(calibs$MAT_0,rep(1,dim(calibs)[1])))
  } else
    print('TModel not recognized - choices are "T" and "T0"')
  
  Y<-as.matrix(calibs$MBT5Me)
  
  iObs<-dim(X)[1]
  iR<-dim(X)[2]
  
  ##Least square approximation of coefficient B (Bhat)
  
  Bhat<-solve(t(X) %*% X,t(X) %*% Y)
  tau2hat<-(1/(iObs-iR))%*%t(Y-X%*%Bhat)%*%(Y-X%*%Bhat)
  
  ##Priors
  #We use a multivariate normal prior for B
  m_p<-Bhat #prior on the mean, let's use least square values
  
  #prior covariance matrix
  var_p<-as.vector((2*Bhat)^2) #variance, choose a large number, scale by Bhat
  cov_p<-diag(var_p) #variance times identity matrix
  
  #we use an inverse gamma prior for tau^2. For your own applications,
  #you might have to do some experimentation to find good starting values
  a<-0.05
  b<-0.006
  
  ## Gibbs sampler
  #set number of chains to do, and draws for each chain:
  Nc<-5 #number of chains
  Nd<-1000 #number of draws
  b_draws<-array(NaN,c(Nc,Nd,iR))
  tau2_draws<-array(NaN,c(Nc,Nd,1))
  N<-dim(X)[1]
  
  for (i in 1:Nc) {
    #start each chain by setting the initial tau^2 value.
    tau2_val<-as.numeric(tau2hat*10*runif(1))
    
    #now start the draws
    for (kk in 1:Nd) {
      #calculate conditional posterior distribution of b:
      b_covar<- solve(((t(X)%*%X)/tau2_val+solve(cov_p))) #covariance
      b_mn<- t(solve(cov_p,m_p)+(t(X)%*%Y)/tau2_val)%*%b_covar #mean
      #make a draw from the distribution
      b_val<- t(t(MASS::mvrnorm(n=1,b_mn,b_covar)))
      
      #now calculate marginal posterior distribution of tau^2:
      p_a<- a+N/2 #shape parameter
      p_b<- b+(1/2)*(t(Y-X%*%b_val))%*%(Y-X%*%b_val) #scale parameter
      #make a draw from the distribution
      tau2_val<-1/Rlab::rgamma(n=1,alpha = p_a, beta = (1/p_b))
      #put in a ceiling in case this blows up (it can happen sometimes)
      if(tau2_val > 20){tau2_val<-20}
      
      #fill in the vectors
      for(w in 1:2) {
        b_draws[i,kk,w]<-b_val[w]
      }
      tau2_draws[i,kk,]<-tau2_val
    }
    
  }
  
##discard first 200 samples of each chain as burn-in
  b_draws<-b_draws[,201:dim(b_draws)[2],]
  tau2_draws<-tau2_draws[,201:dim(tau2_draws)[2],]
  
  #calculate the Rhat statistics to assess convergence
  m<-dim(b_draws)[1]
  n<-dim(b_draws)[2]
  Bs_var<-n/(m-1)*sum((rowMeans(b_draws[,,1])-mean(rowMeans(b_draws[,,1])))^2)
  Ws_var<-mean(apply(b_draws[,,1],1,var))
  Rhat<-sqrt(((n-1)/n*Ws_var+1/n*Bs_var)/Ws_var)
  #Rhat should be close to 1 if convergence is achieved.
  
##Plot posteriors with priors
  
  #Slope
  slope<-b_draws[,,1]
  xts<-seq(min(slope)-0.01,max(slope)+0.01,by=0.0005)
  posts<-kde(c(slope),eval.points = xts) #x<-post$eval, y<-post$estimate
  priors<-dnorm(xts,mean = m_p[1],sd=sqrt(mean(cov_p[1,])))
  
  #Intercept
  intercept<-b_draws[,,2]
  xti<-seq(min(intercept)-0.01,max(intercept)+0.01,by=0.0005)
  posti<-kde(c(intercept),eval.points = xti)
  priori<-dnorm(xti,mean = m_p[1],sd=sqrt(mean(cov_p[1,])))
  
  #Error variance
  xte<-seq(0,0.02,by=0.0001)
  poste<-kde(c(tau2_draws),eval.points = xte)
  priore<-b^a/gamma(a)*xte^(-a-1)*exp(-b/xte)
  
  #Plot
  par(mfrow=c(1,3))
  plot(xts,posts$estimate,type='l',col="red",lwd=2, main = "Slope",ylab = "Prob. density",xlab="")
  lines(xts,priors,col="black",lwd=2)
  legend("topleft",legend = c("Prior","Posterior"),col = c("black","red"),lty=1,lwd = 2)
  
  plot(xti,posti$estimate,type='l',col="red",lwd=2, main = "Intercept",ylab="",xlab="")
  lines(xti,priori,col="black",lwd=2)
  legend("topleft",legend = c("Prior","Posterior"),col = c("black","red"),lty=1,lwd = 2)
  
  plot(xte,poste$estimate,type='l',col="red",lwd=2, main = "Error variance (t^2)",ylab="",xlab="")
  lines(xte,priore,col="black",lwd=2)
  legend("topleft",legend = c("Prior","Posterior"),col = c("black","red"),lty=1,lwd = 2)
  
##Reshape and thin draws, save parameters
  b_draws_r<-matrix(b_draws,nrow=dim(b_draws)[1]*dim(b_draws)[2],ncol=dim(b_draws)[3])
  tau2_draws_r<-matrix(tau2_draws,dim(tau2_draws)[1]*dim(tau2_draws)[2],1)
  b_draws_final<-b_draws_r[seq(1,dim(b_draws_r)[1],by=4),]
  tau2_draws_final<-tau2_draws_r[seq(1,dim(tau2_draws_r)[1],by=4),]
  param_final<-list(b_draws_final,tau2_draws_final)
  
  if(Type == "soil" & Tmodel == "T") {
    saveRDS(param_final,file = "baymbt_params_soil.RDS")
  } else if(Type == "soil" & Tmodel == "T0"){
    saveRDS(param_final,file = "baymbt0_params_soil.RDS")
  } else if(Type == "lake" & Tmodel == "T"){
    saveRDS(param_final,file = "baymbt_params_lake.RDS")
  } else if(Type == "lake" & Tmodel == "T0"){
    saveRDS(param_final,file = "baymbt0_params_lake.RDS")
  }
}  