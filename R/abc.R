#' PET ABC function
#'
#' This function peforms ABC inference for the compartimental models.
#' @param y vector of measured radioactivity concentrations in a voxel for each time-point. This can be in terms of observations or directly in terms of voxel time activity curve.
#' @param tspan vector of time points of observations.
#' @param N Number of simulations. Default to 10^5.
#' @param inputfunction A function describing the input function. Default inputfunction().
#' @param type The type of model in the input function.
#' @param model The type of kinetic model: "single" performs the analysis for the single-tisse
#' model, "two" performs the analysis for the two-tissue model. Default "single".
#' @param a1 Minimum of the uniform prior for K1. Default to 0 .
#' @param b1 Maximum of the uniform prior for K1. Default to 0.2 .
#' @param a2 Minimum of the uniform prior for K2. Default to 0.3 .
#' @param b2 Maximum of the uniform prior for K2. Default to 0.5 .
#' @param a3 Minimum of the uniform prior for K3. Default to 0 .
#' @param b3 Maximum of the uniform prior for K3. Default to 0.1 .
#' @param a4 Minimum of the uniform prior for K4. Default to 0 .
#' @param b4 Maximum of the uniform prior for K4. Default to 0.2 .
#' @param PLOT If plots have to be produced. Default at TRUE.
#' @return ABCout a matrix with values simulated from the posterior distribution of the parameters of the selected model.
#' @return Smat a matrix with values of summary statistics of the simulated curves.
#' @return ABCout_accepted a matrix with values simulated from the posterior distribution of the parameters
#' of the selected model; it shows only the values accepted according to the automatically selected threshold.
#' @return error vector of computed squared differences among the observed and the simulated summary statistics.
#' @return tol automatically selected tolerance level; this tolerance level is used to define the matrix ABCout_accepted.
#' @keywords PETabc
#' @export
PETabc <- function(y,tspan,N=100000, inputfunction.=inputfunction, type=2,model="single",
                   a1=0,b1=0.2,a2=0.3,b2=0.5,a3=0,b3=0.1,a4=0,b4=0.2, PLOT=T)
{
  # y: vector of measured radioactivity concentrations in a voxel for each time-point. This can be in terms of observations or directly in terms of voxel time activity curve.
  # t: vector of time points of observations
  # N: number of simulations
  # Tmax: maximum time range
  # type: type for the singtissue and twotissue functions
  # other parameters: bounds of the prior distributions
  #

  # observed summary statistics
  Sobs=smooth.spline(tspan,y, cv=T)$y

  names.Smat <- c()
  for(j in 1:length(tspan)){
    names.Smat[j] <- paste("y",j,sep="")
  }

  if(model=="single"){

    # Storing matrices
    parMat1=matrix(NA, ncol=2, nrow=N)
    colnames(parMat1) <- c("K1","K2")
    error1=matrix(NA, ncol=1, nrow=N)
    Smat1 <- matrix(NA, nrow=N,ncol=length(tspan))
    colnames(Smat1) <- names.Smat

    # single-tissue function
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

  } else {

    # Storing matrices
    parMat2=matrix(NA, ncol=4, nrow=N)
    colnames(parMat2) <- c("K1","K2","K3","K4")
    error2=matrix(NA, ncol=1, nrow=N)
    Smat2 <- matrix(NA, nrow=N,ncol=length(tspan))
    colnames(Smat2) <- names.Smat

    # two-tissue function
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


  } # end of if

  #Sobs= ksmooth(c(1:61), y2, bandwidth=2)$y
  #Smat1=Smat2=matrix(NA, ncol=length(Sobs), nrow=N)
  #time interval and sampling

  ### ABC simulations

  for (i in 1:N) {

    # printin every 10,000 simulations
    if(round(i/10000)==i/10000){
      cat("i", i, "\n")
    } #end printing

    # Simulation from the prior
    K1=runif(1, a1, b1)
    K2=runif(1, a2, b2)

    if(model=="single"){
      parms1=c(K1, K2, type)

      out=ode(0, tspan, singtissue, parms1, method="ode45")
      data1=t(out[,2])
      #Smat1[i,]=as.vector(data1)
      Smat1[i,]=smooth.spline(tspan, data1, cv=T)$y
      error1[i,]=sum((y-Smat1[i,])^2)
      parMat1[i,] <- c(K1,K2)

      if(round(i/10000)==i/10000){
        parMat_temp <- parMat1[1:i,]
        error_temp <- error1[1:i,]
        h1=quantile(error_temp, probs=0.05)
        out1=parMat_temp[(error_temp<h1[1])==1,]

        par(mfrow=c(1,2))
        plot(density(out1[,1]),main="K1",
             xlab="K1")
        abline(v=mean(out1[,1]))
        plot(density(out1[,2]),main="k2",
             xlab="k2")
        abline(v=mean(out1[,2]))

      } #intermediate plot

    } else { # end of single-tissue

      K3=runif(1, a3, b3)
      K4=runif(1, a4, b4)
      parms2=c(K1, K2, K3, K4, type)

      out=ode(c(0,0), tspan, twotissue, parms2, method="ode45")
      data2=t(out[,2])+t(out[,3])
      #Smat2=as.vector(data2)
      Smat2[i,]=smooth.spline(tspan, data2, cv=T)$y

      error2[i,]=sum((y-Smat2[i,])^2)
      parMat2[i,] <- c(K1,K2,K3,K4)

      if(round(i/10000)==i/10000){
        parMat_temp <- parMat2[1:i,]
        error_temp <- error2[1:i,]
        h2=quantile(error_temp, probs=0.05)
        out2=parMat_temp[(error_temp<h2[1])==1,]

        par(mfrow=c(2,2),oma=c(3,3,0,0),mar=c(3,3,2,2))
        plot(density(out2[,1]),main="K1",
             xlab="K1")
        abline(v=mean(out2[,1]))
        plot(density(out2[,2]),main="k2",
             xlab="k2")
        abline(v=mean(out2[,2]))
        plot(density(out2[,3]),main="k3",
             xlab="K1")
        abline(v=mean(out2[,3]))
        plot(density(out2[,4]),main="k4",
             xlab="K1")
        abline(v=mean(out2[,4]))

      } #intermediate plot

    } # end of if (two-tissue)

    # Ct0=0 #initial condition : cncentration=0 at t=0

    #out=ode(0, tspan, singtissue, parms1, method="ode45")
#    out=ode(0, tspan, singtissue2, parms1[1:2], method="ode45")



#    write(c(K1, K2, K3, K4), file="parMat.out", ncol=4, append=T)
#    write(t(Smat1), file="SMat1.out", ncol=length(time_vec), append=T)
#    write(t(Smat2), file="SMat2.out", ncol=length(time_vec), append=T)
#    write(t(smooth.spline(time_vec, Smat1, cv=T)$y), ncol=length(time_vec),
#          file="Smat1S.out", append=T)
#    write(t(smooth.spline(time_vec, Smat2, cv=T)$y), ncol=length(time_vec),
#          file="Smat2S.out", append=T)

  } # end of for loop

  ### Saving the output

  if(model=="single"){

    # Chosen threshold
    h1=apply(error1, 2, quantile, probs=0.05)
    # Select the values respecting the threshold
    out1=parMat1[(error1[,1]<h1[1])==1,]

  } else{ # end of single-tissue

    # Chosen threshold
    h2=apply(error2, 2, quantile, probs=0.05)
    # Select the values respecting the threshold
    out2=parMat2[(error2[,1]<h2[1])==1,]

  } # end of two-tissue

  ### Posterior plots
  if(PLOT==T){

    if(model=="single"){
      pdf("posteriors_singletissue.pdf")
      par(mfrow=c(1,2))
      plot(density(out1[,1]),main="K1",
           xlab="K1")
      abline(v=mean(out1[,1]))
      plot(density(out1[,2]),main="k2",
           xlab="k2")
      abline(v=mean(out1[,2]))
      dev.off()
    } else { # end of plots for single-tissue
      pdf("posteriors_twotissue.pdf")
      par(mfrow=c(2,2))
      plot(density(out2[,1]),main="K1",
           xlab="K1")
      abline(v=mean(out2[,1]))
      plot(density(out2[,2]),main="k2",
           xlab="k2")
      abline(v=mean(out2[,2]))
      plot(density(out2[,3]),main="k3",
           xlab="K1")
      abline(v=mean(out2[,3]))
      plot(density(out2[,4]),main="k4",
           xlab="K1")
      abline(v=mean(out2[,4]))
      dev.off()
    } # end of plots for two-tissue
  } # end of plots


  ### Output: files
  if(model=="single"){

    write(t(parMat1), file="parMat1.out", ncol=4)
    write(t(Smat1), file="SMat1.out", ncol=60)
    write(t(error1), file="error1.out", ncol=3)

    return(list(ABCout=parMat1,Smat=Smat1,ABCout_accepted=out1,error=error1,tol=h1))

  } else { # end of single-tissue

    write(t(parMat2), file="parMat2.out", ncol=4)
    write(t(Smat2), file="SMat2.out", ncol=60)
    write(t(error2), file="error2.out", ncol=3)

    return(list(ABCout=parMat2,Smat=Smat2,ABCout_accepted=out2,error=error2,tol=h2) )
  } # end of two-tissue

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
