bayrmbt_model<-function(Tmodel,Type){

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
  m_p<-Bhat
  
  var_p<-as.vector((2*Bhat)^2)
  cov_p<-diag(var_p)
  
  a<-0.05
  b<-0.006
  
## Gibbs sampler
  Nc<-5
  Nd<-1000
  b_draws<-array(NaN,c(Nc,Nd,iR))
  tau2_draws<-array(NaN,c(Nc,Nd,1))
  N<-dim(X)[1]
  
  for (i in 1:Nc) {
    tau2_val<-as.numeric(tau2hat*10*runif(1))
    #tau2_val<-0.109283810301172
    for (kk in 1:Nd) {
      b_covar<- solve(((t(X)%*%X)/tau2_val+solve(cov_p)))
      b_mn<- t(solve(cov_p,m_p)+(t(X)%*%Y)/tau2_val)%*%b_covar
      b_val<- t(t(MASS::mvrnorm(n=1,b_mn,b_covar)))
      #b_val<- matrix(c(0.0202,0.3477),nrow = 2,ncol = 1)
      
      p_a<- a+N/2
      p_b<- b+(1/2)*(t(Y-X%*%b_val))%*%(Y-X%*%b_val)
      tau2_val<-1/Rlab::rgamma(n=1,alpha = p_a, beta = (1/p_b))
      #tau2_val<-0.017570583757519
      if(tau2_val > 20){tau2_val<-20}
      
      for(w in 1:2) {
        b_draws[i,kk,w]<-b_val[w]
      }
      tau2_draws[i,kk,]<-tau2_val
    }
    
  }
##Burn in first 200 draws
  b_draws<-b_draws[,201:dim(b_draws)[2],]
  tau2_draws<-tau2_draws[,201:dim(tau2_draws)[2],]
  
  m<-dim(b_draws)[1]
  n<-dim(b_draws)[2]
  Bs_var<-n/(m-1)*sum((rowMeans(b_draws[,,1])-mean(rowMeans(b_draws[,,1])))^2)
  Ws_var<-mean(apply(b_draws[,,1],1,var))
  Rhat<-sqrt(((n-1)/n*Ws_var+1/n*Bs_var)/Ws_var) #Rhat should be close to 1 if convergence is achieved.
  
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
  
##Reshape and save
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