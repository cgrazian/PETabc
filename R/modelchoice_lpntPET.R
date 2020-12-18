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
#' @param tol1 tolerance level for the lp-ntPET model with no activation.
#' If tol1=NULL, the tolerance level automatically chosen by the lp_ntPETabc()
#' function will be used. Default NULL.
#' @param tol1 tolerance level for the lp-ntPET model with activation.
#' If tol1=NULL, the tolerance level automatically chosen by the lp_ntPETabc()
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

  error1 <- abc_out$error_noact
  error2 <- abc_out$error_act

  parMat1 <- abc_out$ABCout_noact
  parMat2 <- abc_out$ABCout_act

  parM1_min=parMat1[(error1==min(error1,na.rm=T))==1,]
  parM2_min=parMat2[(error2==min(error2,na.rm=T))==1,]

  if(is.null(tol1)==T){
    parM1 <- apply(abc_out$ABCout_noact_accepted,2,mean)
    parM1_mode <- apply(abc_out$ABCout_noact_accepted,2,mode_value)
  } else {
    out1=parMat1[(error1<tol1)==1,]
    parM1 <- apply(out1,2,mean)
    parM1_mode <- apply(out1,2,Mode)
  }

    if(is.null(tol2)==T){
    parM2 <- apply(abc_out$ABCout_act_accepted,2,mean)
    parM2_mode <- apply(abc_out$ABCout_act_accepted,2,mode_value)
  } else {
    out2=parMat2[(error2<tol2)==1,]
    parM2 <- apply(out2,2,mean)
    parM2_mode <- apply(out2,2,Mode)
  }

#  hgrid=seq(0.05, 0.4, length=10)
  hgrid=seq(quantile(error1,prob=0.001),quantile(error1,prob=0.08),length=10)
  RMSE1=RMSE2=matrix(NA, nrow=1, ncol=length(hgrid))
  mprob1=mprob2=matrix(NA, ncol=length(hgrid), nrow=1)

  for (j in c(1:length(hgrid))) {
    ind=error1<hgrid[j]
    RMSE1[1, j]=sqrt(mean((apply(parMat1[ind==1,], 2, mean)-parM1)^2))
    mprob1[1, j]=sum(ind)
  }

  for (j in c(1:length(hgrid))) {
    ind=error2<hgrid[j]
    RMSE2[1, j]=sqrt(mean((apply(parMat2[ind==1,], 2,mean)-parM2)^2))
    mprob2[1, j]=sum(ind)
  }
  prob=mprob1/(mprob1+mprob2)

  if(PLOT==T){

    pdf("modelprob_TAC.pdf")
    data1sim=GenCurve(Ct, Cr, Ti, R1=parM1[1], K2=parM1[2], K2a=parM1[3], gamma=0, tD=parM2[5], tP=parM2[6], alpha=parM2[7])$M1
    data2sim=GenCurve(Ct, Cr, Ti, R1=parM2[1], K2=parM2[2], K2a=parM2[3], gamma=parM2[4], tD=parM2[5], tP=parM2[6], alpha=parM2[7])$M2
    data1sim_mode=GenCurve(Ct, Cr, Ti, R1=parM1_mode[1], K2=parM1_mode[2], K2a=parM1_mode[3],
                           gamma=0, tD=parM2_mode[5], tP=parM2_mode[6], alpha=parM2_mode[7])$M1
    data2sim_mode=GenCurve(Ct, Cr, Ti, R1=parM2_mode[1], K2=parM2_mode[2], K2a=parM2_mode[3],
                           gamma=parM2_mode[4], tD=parM2_mode[5], tP=parM2_mode[6], alpha=parM2_mode[7])$M2
    data1sim_min=GenCurve(Ct, Cr, Ti, R1=parM1_min[1], K2=parM1_min[2],
                          K2a=parM1_min[3],gamma=0, tD=parM2_min[5],
                          tP=parM2_min[6], alpha=parM2_min[7])$M1
    data2sim_min=GenCurve(Ct, Cr, Ti, R1=parM2_min[1], K2=parM2_min[2],
                          K2a=parM2_min[3],gamma=parM2_min[4],
                          tD=parM2_min[5], tP=parM2_min[6],
                          alpha=parM2_min[7])$M2

    par(mfrow=c(1,1))
    plot(Ti, Ct, type="p", lwd=3, xlab="time", ylab="Observations",
      ylim=c(0,max(data1sim,data2sim,Ct)),cex.lab=1, cex.axis=0.8)
    #time interval and sampling
    lines(Ti, data1sim, lty=1, lwd=3, col=1)
    lines(Ti, data2sim, lty=2, lwd=3, col=1)
    lines(Ti, data1sim_min, lty=1, lwd=3, col=2)
    lines(Ti, data2sim_min, lty=2, lwd=3, col=2)
    legend("topright",c("no activ.-mean","activ.-mean","no activ.-min err","activ.-min err"),lty=c(1,2,1,2),
           col=c(1,1,2,2),lwd=2,cex=0.8)
    dev.off()

    pdf("modelprobWW_M1.pdf")
    par(mfrow=c(1,1))
    plot(hgrid, prob[1,], type="l", lwd=3, col=1, lty=1,
         xlab="tolerance", ylab="model prob", ylim=c(0,1),
         cex.lab=1.6,cex.axis=0.8)
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    indgrid2=which(RMSE2[1,]==min(RMSE2[1,], na.rm=T))
    hh=(hgrid[indgrid1]+ hgrid[indgrid2])/2
    abline(v=hh)
    abline(h=0.5)
    dev.off()

    pdf("modelprobWW_RMSE.pdf") #WW=1
    par(mfrow=c(1,1))
    plot(hgrid, RMSE1[1,], type="l", lty=1, lwd=3, col=1, xlab="tolerance",
         ylab="RMSE", cex.lab=1.6, cex.axis=0.8,ylim=c(0,max(RMSE1,RMSE2,na.rm=T)))
    lines(hgrid, RMSE2[1,], lty=2, lwd=3, col=1)
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    #half way betwwen the two minimums
    indgrid2=which(RMSE2[1, ]==min(RMSE2[1, ], na.rm=T))
    hh=(hgrid[indgrid1]+ hgrid[indgrid2])/2
    abline(v=hh)
    legend("topright",c("no activ.","activ."),lty=c(1,2),lwd=2)
    dev.off()

    pdf("modelprobWW_RMSE_two.pdf") #WW=1
    par(mfrow=c(1,2))
    plot(hgrid, RMSE1[1,], type="l", lty=1, lwd=3, col=1, xlab="tolerance",main="no act.",
         ylab="RMSE", cex.lab=1.6, cex.axis=0.8,ylim=c(0,max(RMSE1,RMSE2,na.rm=T)))
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    abline(v=hgrid[indgrid1])

    plot(hgrid, RMSE2[1,], type="l", lty=2, lwd=3, col=1, xlab="tolerance",main="act.",
         ylab="RMSE", cex.lab=1.6, cex.axis=0.8,ylim=c(0,max(RMSE1,RMSE2,na.rm=T)))
    indgrid2=which(RMSE2[1,]==min(RMSE2[1,], na.rm=T))
    abline(v=hgrid[indgrid2])

    dev.off()

    }

    return(list(postProb=prob,RMSEnoact=RMSE1,RMSEact=RMSE2,parMnoact_mean=parM1,parMact_mean=parM2,
                parMnoact_mode=parM1_mode,parMact_mode=parM2_mode,
                parMnoact_min=parM1_min,parMact_min=parM2_min))
}
