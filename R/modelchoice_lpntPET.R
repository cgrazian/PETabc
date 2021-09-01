#' Mode function.
#'
#' This function computes the modal value of a kernel density estimate of a vector of
#' (possibly) continuous values.
#' @param xvec vector of values for which to compute the kernel density estimate.
#' @return The modal value of the kernel density estimate of the vector used as input.
#' @keywords PETabc
#' @export
mode_value <- function(xvec){
  density_estimate <- density(xvec)
  mode_value <- density_estimate$x[which.max(density_estimate$y)]
  return(mode_value)
}





#' Model choice function for the lp-ntPET model.
#'
#' This function peforms model choice between the lp-ntPET model with activation and without activation.
#' @param Ct vector of concentration values of the tracer in the target tissue.
#' @param Cr vector of concentration values in the reference regions.
#' @param Ti vector of time points of observations.
#' @param abc_out output of the lp_ntPETabc() function for the lp-ntPET model (both with and without activation).
#' @param tol1 tolerance level for the lp-ntPET model with no activation. The tolerance is given in terms
#' of percentage of posterior values to keep.
#' If tol1=NULL, the tolerance level automatically chosen by the lp_ntPETabc()
#' function will be used. Default NULL.
#' @param tol2 tolerance level for the lp-ntPET model with activation. The tolerance is given in terms
#' of percentage of posterior values to keep.
#' If tol2=NULL, the tolerance level automatically chosen by the lp_ntPETabc()
#' function will be used. Default NULL.
#' @return postProb posterior probabilities of the lp-ntPET model with no activation.
#' @return RMSEnoact root mean squared error of the lp-ntPET model with no activation (on a grid of tolerance levels).
#' @return RMSEact root mean squared error of the lp-ntPET model with activation (on a grid of tolerance levels).
#' @param PLOT If plots have to be produced. Default at FALSE.
#' @keywords PETabc
#' @export
MC_abc_lpntPET <- function(Ct,Cr,Ti,abc_out,tol1=NULL,tol2=NULL,PLOT=F)
{
  # abc_single: output of the PETabc() function for the single-tissue single-tissue compartment model
  # abc_two   : output of the PETabc() function for the two-tissue single-tissue compartment model
  # PLOT      : TRUE/FALSE if plots have to be produced

  parM1_nnls <- abc_out$nnls_noa
  parM2_nnls <- abc_out$nnls_a

  error1 <- abc_out$error_noact
  error2 <- abc_out$error_act

  parMat1 <- abc_out$ABCout_noact
  parMat2 <- abc_out$ABCout_act

  if(is.null(tol1)==T){
    parM1_mean <- apply(abc_out$ABCout_noact_accepted,2,mean)
    parM1_mode <- apply(abc_out$ABCout_noact_accepted,2,mode_value)
    parM1_median <- apply(abc_out$ABCout_noact_accepted,2,median)
  } else {
    h1 <- quantile(error1,probs=tol1)
    out1=parMat1[(error1<h1)==1,]
    parM1_mean <- apply(out1,2,mean)
    parM1_mode <- apply(out1,2,mode_value)
    parM1_median <- apply(out1,2,median)
  }

  if(is.null(tol2)==T){
    parM2_mean <- apply(abc_out$ABCout_act_accepted,2,mean)
    parM2_mode <- apply(abc_out$ABCout_act_accepted,2,mode_value)
    parM2_median <- apply(abc_out$ABCout_act_accepted,2,median)
  } else {
    h2 <- quantile(error2,probs=tol2)
    out2=parMat2[(error2<tol2)==1,]
    parM2_mean <- apply(out2,2,mean)
    parM2_mode <- apply(out2,2,mode_value)
    parM2_median <- apply(out2,2,median)
  }

  lprob_grid <- seq(log(0.00001),log(0.05),length=20)
  hgrid1=quantile(error1,prob=exp(lprob_grid))
  hgrid2=quantile(error2,prob=exp(lprob_grid))
  RMSE1=RMSE2=matrix(NA, nrow=1, ncol=length(lprob_grid))
  mprob1=mprob2=matrix(NA, ncol=length(lprob_grid), nrow=1)

  for (j in c(1:length(hgrid1))) {
    ind=error1<hgrid1[j]
    RMSE1[1, j]=sqrt(mean((apply(parMat1[ind==1,], 2, mean)-parM1_mean)^2))
    mprob1[1, j]=sum(error1 < min(hgrid1[j],hgrid2[j]))
  }

  for (j in c(1:length(hgrid2))) {
    ind=error2<hgrid2[j]
    RMSE2[1, j]=sqrt(mean((apply(parMat2[ind==1,], 2,mean)-parM2_mean)^2))
    mprob2[1, j]=sum(error2 < min(hgrid1[j],hgrid2[j]))
  }
  prob=mprob1/(mprob1+mprob2)

  if(PLOT==T){

    pdf("modelprob_TAC.pdf")
    data1sim_mean=GenCurve(Ct, Cr, Ti, R1=parM1_mean[1], K2=parM1_mean[2],
                      K2a=parM1_mean[3], gamma=0, tD=parM2_mean[5],
                      tP=parM2_mean[6], alpha=parM2_mean[7])$M1
    data2sim_mean=GenCurve(Ct, Cr, Ti, R1=parM2_mean[1], K2=parM2_mean[2], K2a=parM2_mean[3],
                      gamma=parM2_mean[4], tD=parM2_mean[5], tP=parM2_mean[6],
                      alpha=parM2_mean[7])$M2

    data1sim_mode=GenCurve(Ct, Cr, Ti, R1=parM1_mode[1], K2=parM1_mode[2],
                           K2a=parM1_mode[3],gamma=0, tD=parM2_mode[5],
                           tP=parM2_mode[6], alpha=parM2_mode[7])$M1
    data2sim_mode=GenCurve(Ct, Cr, Ti, R1=parM2_mode[1], K2=parM2_mode[2],
                           K2a=parM2_mode[3],gamma=parM2_mode[4],
                           tD=parM2_mode[5], tP=parM2_mode[6], alpha=parM2_mode[7])$M2

    data1sim_median=GenCurve(Ct, Cr, Ti, R1=parM1_median[1], K2=parM1_median[2],
                           K2a=parM1_median[3],gamma=0, tD=parM2_median[5],
                           tP=parM2_median[6], alpha=parM2_median[7])$M1
    data2sim_median=GenCurve(Ct, Cr, Ti, R1=parM2_median[1], K2=parM2_median[2],
                           K2a=parM2_median[3],gamma=parM2_median[4],
                           tD=parM2_median[5], tP=parM2_median[6], alpha=parM2_median[7])$M2

    data1sim_nnls=GenCurve(Ct, Cr, Ti, R1=parM1_nnls[1], K2=parM1_nnls[2],
                          K2a=parM1_nnls[3],gamma=0, tD=parM2_nnls[5],
                          tP=parM2_nnls[6], alpha=parM2_nnls[7])$M1
    data2sim_nnls=GenCurve(Ct, Cr, Ti, R1=parM2_nnls[1], K2=parM2_nnls[2],
                          K2a=parM2_nnls[3],gamma=parM2_nnls[4],
                          tD=parM2_nnls[5], tP=parM2_nnls[6],
                          alpha=parM2_nnls[7])$M2

    par(mfrow=c(1,1))
    plot(Ti, Ct, type="p", lwd=3, xlab="time", ylab="Observations",
      ylim=c(0,max(data1sim_mean,data2sim_mean,Ct)),cex.lab=1, cex.axis=0.8)
    #time interval and sampling
    lines(Ti, data1sim_mean, lty=1, lwd=2, col=1)
    lines(Ti, data2sim_mean, lty=2, lwd=2, col=1)
    lines(Ti, data1sim_median, lty=1, lwd=2, col=2)
    lines(Ti, data2sim_median, lty=2, lwd=2, col=2)
    lines(Ti, data1sim_nnls, lty=1, lwd=2, col=3)
    lines(Ti, data2sim_nnls, lty=2, lwd=2, col=3)
    legend("topright",c("no activ.-mean","activ.-mean","no activ.-median",
                        "activ.-median","no activ.-nnls","activ.-nnls"),
           lty=c(1,2,1,2,1,2),
           col=c(1,1,2,2,3,3),lwd=2,cex=0.8)
    dev.off()

    pdf("modelprobWW_M1.pdf")
    par(mfrow=c(1,1))
    plot(exp(lprob_grid), prob[1,], type="l", lwd=3, col=1, lty=1,
         xlab="tolerance", ylab="model prob", ylim=c(0,1),
         cex.lab=1.6,cex.axis=0.8)
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    indgrid2=which(RMSE2[1,]==min(RMSE2[1,], na.rm=T))
    hh=(exp(lprob_grid)[indgrid1]+ exp(lprob_grid)[indgrid2])/2
    abline(v=hh)
    abline(h=0.5)
    dev.off()

    pdf("modelprobWW_RMSE.pdf") #WW=1
    par(mfrow=c(1,1))
    plot(exp(lprob_grid), RMSE1[1,], type="l", lty=1, lwd=3, col=1, xlab="tolerance",
         ylab="RMSE", cex.lab=1.6, cex.axis=0.8,ylim=c(0,max(RMSE1,RMSE2,na.rm=T)))
    lines(exp(lprob_grid), RMSE2[1,], lty=2, lwd=3, col=1)
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    #half way betwwen the two minimums
    indgrid2=which(RMSE2[1, ]==min(RMSE2[1, ], na.rm=T))
    hh=(exp(lprob_grid)[indgrid1]+ exp(lprob_grid)[indgrid2])/2
    abline(v=hh)
    legend("topright",c("no activ.","activ."),lty=c(1,2),lwd=2)
    dev.off()

    pdf("modelprobWW_RMSE_two.pdf") #WW=1
    par(mfrow=c(1,2))
    plot(exp(lprob_grid), RMSE1[1,], type="l", lty=1, lwd=3, col=1, xlab="tolerance",main="no act.",
         ylab="RMSE", cex.lab=1.6, cex.axis=0.8,ylim=c(0,max(RMSE1,RMSE2,na.rm=T)))
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    abline(v=exp(lprob_grid)[indgrid1])

    plot(exp(lprob_grid), RMSE2[1,], type="l", lty=2, lwd=3, col=1, xlab="tolerance",main="act.",
         ylab="RMSE", cex.lab=1.6, cex.axis=0.8,ylim=c(0,max(RMSE1,RMSE2,na.rm=T)))
    indgrid2=which(RMSE2[1,]==min(RMSE2[1,], na.rm=T))
    abline(v=exp(lprob_grid)[indgrid2])

    dev.off()

    }

    hdf1 <- data.frame(perc=paste(round(as.numeric(exp(lprob_grid))*100,2) , "%" ),
                        err=round(as.numeric(hgrid1,3)) )
    hdf2 <- data.frame(perc=paste(round(as.numeric(exp(lprob_grid))*100,2) , "%" ),
                       err=round(as.numeric(hgrid2,3)) )

    return(list(postProb=prob,RMSEnoact=RMSE1,RMSEact=RMSE2,
                err_noact=hdf1, err_act=hdf2,
                parMnoact_mean=parM1_mean,parMact_mean=parM2_mean,
                parMnoact_mode=parM1_mode,parMact_mode=parM2_mode,
                parMnoact_med=parM1_median,parMact_med=parM2_median))
}
