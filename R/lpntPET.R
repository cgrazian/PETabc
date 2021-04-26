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
#' @import pracma
#' @keywords PETabc
#' @export
#generate the data from Normandin 2012 model
GenCurve=function(Ct, Cr, Ti, R1, K2, K2a, gamma, tD, tP, alpha){
  #Cr=TACtrueC
  #Ct=TACtrueS
  col1=Cr
  #col2=cumsum(Cr)
  col2=cumtrapz(Ti, Cr)
  #col3=-cumsum(Ct)
  col3=-cumtrapz(Ti, Ct)
  Ind=(Ti-tD)>0
  ht=pmax(0, (Ti-tD)/(tP-tD))^(alpha)*exp(alpha*(1- (Ti-tD)/(tP-tD)))*Ind
  #col4=-cumsum(Ct*ht)
  col4=-cumtrapz(Ti, Ct*ht)
  #yy=smooth.spline(c(1:60), Ct, cv=T)$y
  #col4=-cumtrapz(Ti, yy*ht)
  #col4=smooth.spline(c(1:60), col4, cv=T)$y
  BigMat=cbind(col1, col2, col3, col4)
  theta=matrix(c(R1, K2, K2a, gamma), ncol=1, nrow=4)
  M2=BigMat%*%theta
  theta0=matrix(c(R1, K2, K2a, 0), ncol=1, nrow=4)
  M1=BigMat%*%theta0
  return(list(M1=M1,M2=M2))
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
#' @param lsq type of least squares estimation: 1 for non-negative least squares, 2 for weighted non-negative least squares. Default at 1.
#' @param PLOT If plots have to be produced. Default at TRUE.
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
#' @import pracma
#' @import nnls
#' @keywords PETabc
#' @export
# Similar to rat [11C]raclopride studies:
# R1 = 1;K2 = 0.43;K2a = 0.23; = 0.07; tD = 20; tP = 35; alpha= 1
# Choose flag for lsq=1,2 or 3 depending on the least square method to apply

lp_ntPETabc <- function(Ct,Cr,Ti,S=10^5,R1a=0.8,R1b=1.1,K2alpha=0.2,K2beta=0.6,K2a.alpha=0.1,K2a.beta=0.4,
                        gamma.a=0,gamma.b=0.2,tD.a=15,tD.b=25,tP.b=40,alpha.a=0,alpha.b=3,lsq=1, PLOT=T) {
  #### nnls - option 1
  # model with no activation: MRTM
  col1=Cr
  if (lsq==1) {
    col2=cumtrapz(Ti, Cr)
    col3=-cumtrapz(Ti, Ct)
    BigMat=cbind(col1, col2, col3)
    obj_nnls_noa <- nnls(BigMat,Ct)

    # model with activation: lp-ntPET
    # Keep the same range for the sampling of td, tp and alpha
    count=0
    res_vec <- c()
    for (tD in seq(tD.a,tD.b,0.25)) {
      for (tP in seq (tD+1,tP.b,0.25)) {
        for (alpha in seq(alpha.a,alpha.b,0.25)) {
          count <- count + 1
          Ind=(Ti-tD)>0
          ht=pmax(0, (Ti-tD)/(tP-tD))^(alpha)*
            exp(alpha*(1- (Ti-tD)/(tP-tD)))*Ind
          col4=-cumtrapz(Ti, Ct*ht)
          BigMat=cbind(col1, col2, col3, col4)
          obj_nnls_Ma <- nnls(BigMat,Ct)
          if (count==1){
            nnls_mat <- matrix(c(obj_nnls_Ma$x,tD,tP,alpha),ncol=7,nrow=1)
          } else {
            nnls_mat <- rbind(nnls_mat,c(obj_nnls_Ma$x,tD,tP,alpha))
          }
          res_vec <- c(res_vec,sum(obj_nnls_Ma$residuals))
        }
      }
    }
    obj_nnls_a <- nnls_mat[res_vec==min(res_vec),]
  }

  #### wnnls - option 2
  if (lsq==2) {
    # model with no activation
    col1=Cr
    col2=cumtrapz(Ti, Cr)
    col3=-cumtrapz(Ti, Ct)
    weights=diag(1.0/(sqrt(Ct)+eps))
    wCt=weights%*%Ct
    BigMat=cbind(col1, col2, col3)
    wBigMat=weights%*%BigMat
    obj_wnnls_noa <- nnls(wBigMat,wCt)
    sum(obj_wnnls_noa$residuals^2)
    weightvec=1.0/(sqrt(Ct)+eps)
    sum(weightvec*(Ct-obj_wnnls_noa$fitted)^2)

    # model with activation
    count=0
    res_vec <- c()
    for (tD in seq(td.a,td.b,0.25))
    {
      for (tP in seq (tD+1,tP.b,0.25))
      {
        for (alpha in seq(alpha.a,alpha.b,0.25))
        {
          count <- count + 1
          Ind=(Ti-tD)>0
          ht=pmax(0, (Ti-tD)/(tP-tD))^(alpha)*
            exp(alpha*(1- (Ti-tD)/(tP-tD)))*Ind
          col4=-cumtrapz(Ti, Ct*ht)
          BigMat=cbind(col1, col2, col3, col4)
          wBigMat=weights%*%BigMat
          obj_wnnls_Ma <- nnls(wBigMat,wCt)
          if (count==1){
            wnnls_mat <- matrix(c(obj_wnnls_Ma$x,tD,tP,alpha),ncol=7,nrow=1)
          } else {
            wnnls_mat <- rbind(wnnls_mat,c(obj_wnnls_Ma$x,tD,tP,alpha))
          }
          res_vec <- c(res_vec,sum(obj_wnnls_Ma$residuals^2))
        }
      }
    }

    #the residuals are calculated for all combinations of implicit parameters, td, tp, alpha
    #RSS or res_vec=sum(weights*(Ct-FittedCt)*(Ct-FittedCt))
    obj_wnnls_a <- wnnls_mat[res_vec==min(res_vec),] # estimated parameters selected from minimum res_vec
  }

  # Observed summary statistics
  nT=length(Ti)
  Sobs=smooth.spline(c(1:nT), Ct, df=15)$y

  errorM1=errorM2=c()
  parMat_act=matrix(0, nrow=S, ncol=7)
  parMat_noact=matrix(0, nrow=S, ncol=3)
  Smat_act <- matrix(NA, nrow=S,ncol=length(Ti))
  Smat_noact <- matrix(NA, nrow=S,ncol=length(Ti))

  for (iter in 1:S){
    # Simulation of the parameters from their prior distributions
    R1=runif(1, R1a, R1b)
    K2=runif(1, K2alpha, K2beta)
    K2a=runif(1, K2a.alpha, K2a.beta)
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
    errorM2[iter]=sum(abs(Sobs-Ssim2))
    parMat_act[iter,] <- c(R1,K2,K2a,gamma,tD,tP,alpha)
    parMat_noact[iter,] <- c(R1,K2,K2a)

    if(round(iter/10000)==iter/10000){

      cat(iter, "..\n")

      parMat_noact_temp <- parMat_noact[1:iter,]
      error_temp <- errorM1[1:iter]
      h1=quantile(error_temp, probs=0.05)
      out1=parMat_noact_temp[(error_temp<h1[1])==1,]

      parMat_act_temp <- parMat_act[1:iter,]
      error_temp <- errorM2[1:iter]
      h2=quantile(error_temp, probs=0.05)
      out2=parMat_act_temp[(error_temp<h2[1])==1,]

      par(mfrow=c(2,2))

      plot(density(out2[,5]),main="lp-ntPET with act",xlab=expression(t[D]))
      abline(v=mean(out2[,5]))

      plot(density(out2[,6]),main="lp-ntPET with act",xlab=expression(t[P]))
      abline(v=mean(out2[,6]))

      plot(density(out2[,4]),main="lp-ntPET with act",xlab=expression(gamma))
      abline(v=mean(out2[,4]))

      plot(density(out2[,3]),main="lp-ntPET with act",xlab=expression(k[2][a]))
      abline(v=mean(out2[,3]))

    } #intermediate plot

  }

  # Chosen threshold for no activation
  h1=quantile(errorM1, probs=0.05)
  # Select the values respecting the threshold
  out1=parMat_noact[(errorM1<h1[1])==1,]

  # Chosen threshold for activation
  h2=quantile(errorM2, probs=0.05)
  # Select the values respecting the threshold
  out2=parMat_act[(errorM2<h2[1])==1,]

  ### Saving the output

  ### Posterior plots
  if(PLOT==T){

      pdf("posteriors_noactivation.pdf")
      par(mfrow=c(2,2))
      plot(density(out1[,1]),main=expression(R[1]), xlab=expression(R[1]))
      abline(v=mean(out1[,1]))

      plot(density(out1[,2]),main=expression(k[2]), xlab=expression(k[2]))
      abline(v=mean(out1[,2]))

      plot(density(out1[,3]),main=expression(k[2][a]), xlab=expression(k[2][a]))
      abline(v=mean(out1[,3]))

      dev.off()

      pdf("posteriors_activation_explicit.pdf")
      par(mfrow=c(2,2))
      plot(density(out2[,1]),main=expression(R[1]), xlab=expression(R[1]))
      abline(v=mean(out2[,1]))

      plot(density(out2[,2]),main=expression(k[2]), xlab=expression(k[2]))
      abline(v=mean(out2[,2]))

      plot(density(out2[,3]),main=expression(k[2][a]), xlab=expression(k[2][a]))
      abline(v=mean(out2[,3]))

      plot(density(out2[,4]),main=expression(gamma), xlab=expression(gamma))
      abline(v=mean(out2[,4]))

      dev.off()

      pdf("posteriors_activation_implicit.pdf")
      par(mfrow=c(2,2))
      plot(density(out2[,5]),main=expression(t[D]), xlab=expression(t[D]))
      abline(v=mean(out2[,5]))

      plot(density(out2[,6]),main=expression(t[P]), xlab=expression(t[P]))
      abline(v=mean(out2[,6]))

      plot(density(out2[,7]),main=expression(alpha), xlab=expression(alpha))
      abline(v=mean(out2[,7]))

      dev.off()
  } # end of plots

  ### Output: files
    write(t(parMat_act), file="parMat_act.out", ncol=4)
    write(t(Smat_act), file="SMat_act.out", ncol=60)
    write(t(errorM2), file="error_act.out", ncol=3)

    write(t(parMat_noact), file="parMat_noact.out", ncol=4)
    write(t(Smat_noact), file="SMat_noact.out", ncol=60)
    write(t(errorM1), file="error_noact.out", ncol=3)

  #modify returned values depending on the least square method used
  return(list(ABCout_act=parMat_act,Smat_act=Smat_act,ABCout_act_accepted=out2,error_act=errorM2,tol_act=h2,
              ABCout_noact=parMat_noact,Smat_noact=Smat_noact,ABCout_noact_accepted=out1,error_noact=errorM1,tol_noact=h1,
              nnls_noa=obj_nnls_noa$x,nnls_a=obj_nnls_a) )#modify depending on the least square method used

}
