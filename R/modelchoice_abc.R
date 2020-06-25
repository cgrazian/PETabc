#' Model choice function.
#'
#' This function peforms model choice between the single-tissue and the 
#' two-tissue models based on the outputs of the PETabc() function.
#' @param y vector of measured radioactivity concentrations in a voxel for 
#' each time-point. This can be in terms of observations or directly in terms 
#' of voxel time activity curve.
#' @param tspan vector of time points of observations.
#' @param abc_single output of the PETabc() function for the single-tissue 
#' single-tissue compartment model
#' @param abc_two output of the PETabc() function for the two-tissue 
#' single-tissue compartment model
#' @param tol1 tolerance level for the single-tissue model. 
#' If tol1=NULL, the tolerance level automatically chosen by the PETabc() 
#' function will be used. Default NULL.
#' @param tol2 tolerance level for the two-tissue model. 
#' If tol2=NULL, the tolerance level automatically chosen by the PETabc() 
#' function will be used. Default NULL.
#' @param inputfunction A function describing the input function. 
#' Default inputfunction().
#' @param type The type of model in the input function.
#' @param PLOT If plots have to be produced. Default at FALSE.
#' @return ABCres_single an object of class abc for the single-tissue compartment model.
#' @return ABCres_two an object of class abc for the two-tissue compartment model.
#' @keywords PETabc
#' @export
MC_abc <- function(y,tspan,abc_single, abc_two,tol1=NULL,tol2=NULL,
                   inputfunction.=inputfunction, 
                   type=2,PLOT=F)
{
  # abc_single: output of the PETabc() function for the single-tissue single-tissue compartment model
  # abc_two   : output of the PETabc() function for the two-tissue single-tissue compartment model
  # PLOT      : TRUE/FALSE if plots have to be produced   

  error1 <- abc_single$error
  error2 <- abc_two$error
  
  parMat1 <- abc_single$ABCout
  parMat2 <- abc_two$ABCout

  if(is.null(tol1)==T){
    parM1 <- apply(abc_single$ABCout_accepted,2,mean)
  } else {
    out1=parMat1[(error1[,1]<tol1)==1,]
    parM1 <- apply(out1,2,mean)
  }
  
  if(is.null(tol2)==T){
    parM2 <- apply(abc_two$ABCout_accepted,2,mean)
  } else {
    out2=parMat2[(error2[,1]<tol2)==1,]
    parM2 <- apply(out2,2,mean)
  }

  hgrid=seq(0.05, 0.4, length=10)
  RMSE1=RMSE2=matrix(NA, nrow=1, ncol=length(hgrid))
  mprob1=mprob2=matrix(NA, ncol=length(hgrid), nrow=1)

  for (j in c(1:length(hgrid))) {
    ind=error1[,1]<hgrid[j]
    RMSE1[1, j]=sqrt(mean((apply(parMat1[ind==1,], 2, mean)-parM1)^2))
    mprob1[1, j]=sum(ind)
  }

  for (j in c(1:length(hgrid))) {
    ind=error2[,1]<hgrid[j]
    RMSE2[1, j]=sqrt(mean((apply(parMat2[ind==1,], 2,mean)-parM2)^2))
    mprob2[1, j]=sum(ind)
  }
  prob=mprob1/(mprob1+mprob2)

  singtissue= function(t, y, parms){
    K1c=parms[1]
    k2c=parms[2]
    type=parms[3]
    ipt=inputfunction.(t, type)
    #differential equation describing the concentration in the reference tissue
    #input function taken from...
    dydt= K1c*ipt-k2c*y
    return(list(c(dydt)))
  } # end of single-tissue function

  twotissue = function(t, y, parms) {
    #differential equation describing the concentration in the target tissue
    #based on a 2 tissue compartment model
    K1=parms[1]
    k2=parms[2]
    k3=parms[3]
    k4=parms[4]
    type=parms[5]
    ipt=inputfunction.(t, type)
    dydt=c()
    dydt[1]= ((K1*ipt)-(k2*y[1])-(k3*y[1])+(k4*y[2]))
    dydt[2]=(k3*y[1]-k4*y[2])
    
    return(list(c(dydt)))
  } # end of two-tissue function
  
  pdf("modelprob_TAC.pdf")
  par(mfrow=c(1,1))
  plot(y, type="p", lwd=3, xlab="time", ylab="Observations",
     ylim=c(0, 0.3), cex.lab=1, cex.axis=0.8)
  parms1=c(parM1,  type)
  #time interval and sampling
  out=ode(0, tspan, singtissue, parms1, method="ode45")
  data1sim=as.vector(t(out[,2]))
  lines(data1sim, lty=1, lwd=3, col=1)

  parms2=c(parM2,  type)
  out=ode(c(0,0), tspan, twotissue, parms2, method="ode45")
  data2sim=as.vector(t(out[,2])+t(out[,3]))
  lines(data2sim, lty=2, lwd=3, col=1)

  legend("topright",c("single","two"),lty=c(1,2),lwd=2)
  dev.off()

#    par(mfrow=c(1,1))
#    plot(hgrid, prob[1,], ylim=c(0,1),type="l")
#    lines(hgrid, RMSE1[1,], col=2)
#    lines(hgrid, RMSE2[1,], col=3)

  pdf("modelprobWW_M1.pdf")
    par(mfrow=c(1,1))
    plot(hgrid, prob[1,], type="l", lwd=3, col=1, lty=1,
       xlab="tolerance", ylab="model prob", ylim=c(0,1),
       cex.lab=1.6,cex.axis=0.8)
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    #half way betwwen the two minimums
    indgrid2=which(RMSE2[1,]==min(RMSE2[1,], na.rm=T))
    hh=(hgrid[indgrid1]+ hgrid[indgrid2])/2
    abline(v=hh)
  dev.off()

  pdf("modelprobWW_RMSE.pdf") #WW=1
    par(mfrow=c(1,1))
    plot(hgrid, RMSE1[1,], type="l", lty=1, lwd=3, col=1, xlab="tolerance",
       ylab="RMSE", cex.lab=1.6, cex.axis=0.8)
    lines(hgrid, RMSE2[1,], lty=2, lwd=3, col=1)
    indgrid1=which(RMSE1[1,]==min(RMSE1[1,], na.rm=T))
    #half way betwwen the two minimums
    indgrid2=which(RMSE2[1, ]==min(RMSE2[1, ], na.rm=T))
    hh=(hgrid[indgrid1]+ hgrid[indgrid2])/2
    abline(v=hh)
  dev.off()

}
