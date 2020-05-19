#' PET ABC function
#'
#' This function peforms ABC inference for the compartimental models.
#' @param y vector of measured radioactivity concentrations in a voxel for each time-point. This can be in terms of observations or directly in terms of voxel time activity curve.
#' @param t vector of time points of observations.
#' @param N Number of simulations. Default to 10^5.
#' @param inputfunction A function describing the input function. Default inputfunction().
#' @param type The type of model in the input function.
#' @param a1 Minimum of the uniform prior for K1. Default to 0 .
#' @param b1 Maximum of the uniform prior for K1. Default to 0.2 .
#' @param a2 Minimum of the uniform prior for K2. Default to 0.3 .
#' @param b2 Maximum of the uniform prior for K2. Default to 0.5 .
#' @param a3 Minimum of the uniform prior for K3. Default to 0 .
#' @param b3 Maximum of the uniform prior for K3. Default to 0.1 .
#' @param a4 Minimum of the uniform prior for K4. Default to 0 .
#' @param b4 Maximum of the uniform prior for K4. Default to 0.2 .
#' @param tol Tolerance level of the abc function.
#' @return ABCres_single an object of class abc for the single-tissue compartment model.
#' @return ABCres_two an object of class abc for the two-tissue compartment model.
#' @keywords PETabc
#' @export
PETabc <- function(y,t,N=100000, inputfunction.=inputfunction, type=2,
                   a1=0,b1=0.2,a2=0.3,b2=0.5,a3=0,b3=0.1,a4=0,b4=0.2,
                   tol)
{
  # y: vector of measured radioactivity concentrations in a voxel for each time-point. This can be in terms of observations or directly in terms of voxel time activity curve.
  # t: vector of time points of observations
  # N: number of simulations
  # Tmax: maximum time range
  # type: type for the singtissue and twotissue functions
  # other parameters: bounds of the prior distributions
  #

  # Storing matrices
  parMat1=matrix(NA, ncol=2, nrow=N)
  colnames(parMat1) <- c("K1","K2")
  parMat2=matrix(NA, ncol=4, nrow=N)
  colnames(parMat2) <- c("K1","K2","K3","K4")
  errorMat1=errorMat2=c()

  # singtissue2 <- function(t,y,parms){
  #   K1=parms[1]
  #   k2=parms[2]
  #   dydt= K1*ipt[1]-k2*y[1]
  #   return(list(c(dydt)))
  # }

  singtissue= function(t, y, parms){
    K1c=parms[1]
    k2c=parms[2]
    type=parms[3]
    ipt=inputfunction.(t, type)
    #differential equation describing the concentration in the reference tissue
    #input function taken from...
    dydt= K1c*ipt-k2c*y
    return(list(c(dydt)))
  }

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
  }

  Sobs=smooth.spline(time_vec,y, cv=T)$y
  #Sobs= ksmooth(c(1:61), y2, bandwidth=2)$y
  #Smat1=Smat2=matrix(NA, ncol=length(Sobs), nrow=N)
  tspan=time_vec #time interval and sampling

  Smat1 <- matrix(NA, nrow=N,ncol=length(time_vec))
  Smat2 <- matrix(NA, nrow=N,ncol=length(time_vec))
  names.Smat <- c()
    for(j in 1:length(time_vec)){
      names.Smat[j] <- paste("y",j,sep="")
  }
  colnames(Smat1) <- names.Smat
  colnames(Smat2) <- names.Smat

  for (i in c(1:N)) {

    if(round(i/10000)==i/10000){
      cat("i", i, "\n")
    }
    # Initialisation
      K1=runif(1, a1, b1)
      K2=runif(1, a2, b2)
      K3=runif(1, a3, b3)
      K4=runif(1, a4, b4)

      parMat1[i,] <- c(K1,K2)
      parMat2[i,] <- c(K1,K2,K3,K4)

      parms1=c(K1, K2, type)
      parms2=c(K1, K2, K3, K4, type)
    # Ct0=0 #initial condition : cncentration=0 at t=0

    #out=ode(0, tspan, singtissue, parms1, method="ode45")
    out=ode(0, tspan, singtissue, parms1, method="ode45")
#    out=ode(0, tspan, singtissue2, parms1[1:2], method="ode45")

    data1=t(out[,2])
    Smat1[i,]=as.vector(data1)
    #Smat1=smooth.spline(time_vec, data1, cv=T)$y

    out=ode(c(0,0), tspan, twotissue, parms2, method="ode45")
    data2=t(out[,2])+t(out[,3])
    #Smat2=as.vector(data2)
    Smat2[i,]=smooth.spline(time_vec, data2, cv=T)$y

#    write(c(K1, K2, K3, K4), file="parMat.out", ncol=4, append=T)
#    write(t(Smat1), file="SMat1.out", ncol=length(time_vec), append=T)
#    write(t(Smat2), file="SMat2.out", ncol=length(time_vec), append=T)
#    write(t(smooth.spline(time_vec, Smat1, cv=T)$y), ncol=length(time_vec),
#          file="Smat1S.out", append=T)
#    write(t(smooth.spline(time_vec, Smat2, cv=T)$y), ncol=length(time_vec),
#          file="Smat2S.out", append=T)

  }
  out1=abc(target=Sobs, param=parMat1, sumstat=Smat1, tol=tol,
           method="neuralnet")
  out2=abc(target=Sobs, param=parMat2, sumstat=Smat2, tol=tol,
           method="neuralnet")
  #  adj=out$adj.values
  #  unadj=out$unadj.values

  return(list(ABCres_single=out1,ABCres_two=out2))

}

#
#
#
# #simulate N=200000 parameters and store
# N=1000
# parMat1=matrix(NA, ncol=2, nrow=N)
# parMat2=matrix(NA, ncol=4, nrow=N)
# errorMat1=errorMat2=c()
# Sobs=smooth.spline(c(1:61), y1, cv=T)$y
# #Sobs= ksmooth(c(1:61), y2, bandwidth=2)$y
# #Smat1=Smat2=matrix(NA, ncol=length(Sobs), nrow=N)
#
# ntries=0
#
# for (i in c(1:N)) {
#   cat("i", i, ".. ntries", ntries, "\n")
#   error=999
#   ntries=0
#
#   while (error > 0.7) {
#     ntries=ntries+1
#     K1=runif(1, 0, 0.2)
#     K2=runif(1, 0.3, 0.5)
#     # K3=runif(1, 0, 0.1)
#     # K4=runif(1, 0, 0.2)
#
#
#     parms1=c(K1, K2, 2)
#     #parms2=c(K1, K2, K3, K4, 2)
#     tspan=c(0:60) #time interval and sampling
#     # Ct0=0 #initial condition : cncentration=0 at t=0
#
#     out=ode(0, tspan, singtissue, parms1, method="ode45")
#     data1=t(out[,2])
#     Smat1[i,]=smooth.spline(c(1:61), data1, cv=T)$y
#     error=sum(abs(Smat1[i,]-Sobs))
#     errorMat1[i]=error
#     parMat1[i,]=c(K1, K2)
#   }
#
#   # out=ode(c(0,0), tspan, twotissue, parms2, method="ode45")
#   # data2=t(out[,2])+t(out[,3])
#   #Smat2=smooth.spline(c(1:61), data2, cv=T)$y
#   # Smat2=data2
#
#   Ind=errorMat1< thresh   #0.685
#   par(mfrow=c(2,2))
#   hist(parMat1[Ind==1,1])
# #  abline(v=0.0918, lwd=3, col=2)
# #  abline(v=mean(parMat1[Ind==1,1]), lwd=3, col=3)
#   hist(parMat1[Ind==1,2])
# #  abline(v=0.4484, lwd=3, col=2)
# #  abline(v=mean(parMat1[Ind==1,2]), lwd=3, col=3)
#   plot( parMat1[Ind==1,1], errorMat1[Ind==1])
#   plot( parMat1[Ind==1,2], errorMat1[Ind==1])
#
#
#   library(abc)
#   out=abc(Sobs, parMat1, Smat1, tol=1, method="loclinear")
#   adj=out$adj.values
#   unadj=out$unadj.values
#
#   return(ABCres=out)
# }
