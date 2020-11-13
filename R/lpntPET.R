#' Generate data from lp-ntPET model
#'
#' This function generate the data from the lp-ntPET model as in Normandin (2012).
#' @param Ct vector of concentration values of the tracer in the target tissue.
#' @param Cr vector of concentration values in the reference regions. 
#' @param Ti vector of time points of observations.
#' @param R1 value of the R1 parameter.
#' @param K2 value of the K2 parameter.
#' @param K2a value of the K2a parameter.
#' @param gamma value of the activation parameter.
#' @param tD delay time at which the response starts.
#' @param tP peak time of maximal response.
#' @param alpha sharpeness parameters of the response function. 
#' @return M1 data under the model with no activation (gamma=0).
#' @return M2 data under the model with activation (gamma different from 0).
#' @keywords PETabc
#' @export
#generate the data from Normandin 2012 model
GenCurve=function(Ct, Cr, Ti, R1, K2, K2a, gamma, tD, tP, alpha){  
  col1=Cr
  col2=cumsum(Cr)
  col3=-cumsum(Ct)
  Ind=(Ti-tD)>0
  ht=pmax(0, (Ti-tD)/(tP-tD))^(alpha)*exp(alpha*(1- (Ti-tD)/(tP-tD)))*Ind
  col4=-cumsum(Ct*ht)
  BigMat=cbind(col1, col2, col3, col4)
  #for model M1 with gamma=0 and M2 with gamma\neq 0
  theta1=matrix(c(R1, K2, K2a, 0), ncol=1, nrow=4)
  theta2=matrix(c(R1, K2, K2a, gamma), ncol=1, nrow=4)
  
  return(list(M1=BigMat%*%theta1, M2=BigMat%*%theta2))
}



#' PET ABC function for the lp-nt model.
#'
#' This function peforms ABC inference for the lp-nt model. 
#' @param Ct vector of concentration values of the tracer in the target tissue.
#' @param Cr vector of concentration values in the reference regions. 
#' @param Ti vector of time points of observations.
#' @param S number of simulations. Default at 10^5.
#' @param R1a lower bound of the uniform prior distribution for R1. Default at 0.5.
#' @param R1a upper bound of the uniform prior distribution for R1. Default at 1.6.
#' @param K2alpha lower bound of the uniform prior distribution for K2. Default at 0.
#' @param K2beta upper bound of the uniform prior distribution for K2. Default at 1.
#' @param K2a.alpha lower bound of the uniform prior distribution for K2a. Default at 0.
#' @param K2a.beta upper bound of the uniform prior distribution for K2a. Default at 0.6.
#' @param gamma.a lower bound of the uniform prior distribution for gamma. Default at 0.
#' @param gamma.b upper bound of the uniform prior distribution for gamma. Default at 0.2.
#' @param tD.a lower bound of the uniform prior distribution for tD. Default at 18.
#' @param tD.b upper bound of the uniform prior distribution for tD. Default at 22.
#' @param tP.b upper bound of the uniform prior distribution for tP. Default at 40.
#' @param alpha.a lower bound of the uniform prior distribution for alpha. Default at 0.
#' @param alpha.b upper bound of the uniform prior distribution for alpha. Default at 3.
#' @return ABCout_act a matrix with values simulated from the posterior distribution of the parameters of the model with activation.
#' @return Smat_act a matrix with values of summary statistics of the simulated curves for the model with activation.
#' @return ABCout_accepted_act a matrix with values simulated from the posterior distribution of the parameters 
#' of the model with activation; it shows only the values accepted according to the automatically selected threshold.
#' @return error_act vector of computed squared differences among the observed and the simulated summary statistics for the model with activation. 
#' @return tol_act automatically selected tolerance level for the model with activation; this tolerance level is used to define the matrix ABCout_accepted.
#' @return ABCout_noact a matrix with values simulated from the posterior distribution of the parameters of the model with no activation.
#' @return Smat_noact a matrix with values of summary statistics of the simulated curves for the model with no activation.
#' @return ABCout_accepted_noact a matrix with values simulated from the posterior distribution of the parameters 
#' of the model with no activation; it shows only the values accepted according to the automatically selected threshold.
#' @return error_noact vector of computed squared differences among the observed and the simulated summary statistics for the model with no activation. 
#' @return tol_noact automatically selected tolerance level for the model with no activation; this tolerance level is used to define the matrix ABCout_accepted.
#' @keywords PETabc
#' @export
lp_ntPETabc <- function(Ct,Cr,Ti,S=10^5,R1a=0.5,R1b=1.6,K2alpha=0,K2beta=1,K2a.alpha=0,K2a.beta=0.6,
                        gamma.a=0,gamma.b=0.2,tD.a=18,tD.b=22,tP.b=40,alpha.a=0,alpha.b=3)
{

  # Observed ummary statistics
  nT=length(Ti)
  Sobs=smooth.spline(c(1:nT), Ct, df=15)$y
  
  errorM1=errorM2=c()
  parMat_act=matrix(0, nrow=S, ncol=7)
  parMat_noact=matrix(0, nrow=S, ncol=3)
  Smat_act <- matrix(NA, nrow=S,ncol=length(Ti))
  Smat_noact <- matrix(NA, nrow=S,ncol=length(Ti))
  
  for(iter in 1:S){
    # Simulation of the parameters from their prior distributions
    R1=runif(1, R1a, R1b)
    K2=runif(1, K2alpha, K2beta)
    K2a=runif(1, K2.alpha, K2.beta)
    gamma=runif(1, gamma.a, gamma.b)
    tD=runif(1, tD.a, tD.b)
    tP=runif(1, tD+1, tP.b)
    alpha=runif(1, alpha.a, alpha.b)
    
    data_sim=GenCurve(Ct, Cr, Ti, R1, K2, K2a, gamma, tD, tP, alpha)
    Ssim1=data_sim$M1
    Ssim2=data_sim$M2
    Smat_act[iter,] <- Ssim2
    Smat_noact[iter,] <- Ssim1
    
    errorM1[iter]=sum(abs(Sobs-Ssim1))
    errorM2[iter]=sum(abs(Sobs-Ssims2))
    parMat_act[iter,] <- c(R1,K2,K2a,gamma,tD,tP,alpha)
    parMat_noact[iter,] <- c(R1,K2,K2a)
    
    if(round(iter/10000)==iter/10000){
      
      cat(iter, "..\n")
      
      parMat_temp <- parMat_noact[1:iter,]
      error_temp <- errorM1[1:iter,]
      h1=quantile(error_temp, probs=0.05)
      out1=parMat_temp[(error_temp<h1[1])==1,]
      
      parMat_temp <- parMat_act[1:iter,]
      error_temp <- errorM2[1:iter,]
      h2=quantile(error_temp, probs=0.05)
      out2=parMat_temp[(error_temp<h2[1])==1,]
      
      # par(mfrow=c(1,2))
      # plot(density(out1[,1]),main="K1",
      #      xlab="K1")
      # abline(v=mean(out1[,1]))
      # plot(density(out1[,2]),main="k2",
      #      xlab="k2")
      # abline(v=mean(out1[,2]))
      
    } #intermediate plot

  } 

  # Chosen threshold for no activation
  h1=apply(errorM1, 2, quantile, probs=0.05)
  # Select the values respecting the threshold
  out1=parMat_noact[(errorM1[,1]<h1[1])==1,]
  
  # Chosen threshold for no activation
  h2=apply(errorM2, 2, quantile, probs=0.05)
  # Select the values respecting the threshold
  out2=parMat_act[(errorM2[,1]<h2[1])==1,]
  
  return(list(ABCout_act=parMat_act,Smat_act=Smat_act,ABCout_act_accepted=out2,error_act=error2,tol_act=h1,
              ABCout_noact=parMat_noact,Smat_noact=Smat_noact,ABCout_noact_accepted=out1,error_noact=error1,tol_noact=h1) )
  
}

