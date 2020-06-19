#' Input function
#'
#' Input function, values from Feng et al. 1993, Models for computer simulation
#' studies of input functions for tracer kinetic modeling with positron
#' emission tomography.
#'
#' @param t time
  #' @param type type of the model. If type=1 is Model 1 in Feng et al. Type=2 represents Model 2 in Feng et al.
#' @return A vector of concentrations for each time point.
#' @keywords BayesPET
#' @export
inputfunction= function(t,type) {

  #Input function, values from Feng et al. 1993, Models for computer simulation
  #studies of input functions for tracer kinetic modeling with positron
  #emission tomography

  if (type==1) {
    Ct=c(851.1, 21.9, 20.8)
    bt=c(-4.13, -0.12, -0.01)
  }

  if  (type==2) {
    Ct = c(12, 1.8, 0.45)
    bt = c(-4, -0.5, -0.008)
  }


  if (all((type==c(1,2))==FALSE)) {
    cat("Incorrect type for input function parameters", "\n")
    cat("Default values", "\n")
    Ct=c(851.1, 21.9, 20.8)
    bt=c(-4.13, -0.12, -0.01)
  }

  #cp=(A*t - B - C)*exp(a*t) + B*exp(b*t) + C*exp(c*t);

  y=exp(bt[1]*t)* (Ct[1]*t - Ct[2] - Ct[3])+ Ct[2]*exp(bt[2]*t) + Ct[3]*exp(bt[3]*t)

  return(y)
}

#' Single-tissue compartment model
#'
#' This function produce the output of the differential equation describing the
#' concentration in the reference tissue input function. The model considers a
#' bidirectional flux between tissue and blood.
#'
#' @param t time
#' @param y a matrix or dataframe of measured radioactivity concentration in a voxel
#' for each time-point. This can be in terms of observations or
#' directly in terms of voxel time activity curve. The first column refers
#' to the observations from the tissue.
#' @param parms vector of two elements: K1, constant for the flux from blood
#' to tissue, k2, constant for the flux from tissue to blood.
#' @param inputfunction input function. Default inputfunction()
#' @param type type of the input function.
#' @return A list with the output of the differential equation.
#' @keywords BayesPET
#' @export
singtissue= function(t, yobs, parms, inputfunction.=inputfunction,
                     type=2){

  ipt <- inputfunction(t, type)
  y   <- yobs[,1]
  K1c <- parms[1]
  k2c <- parms[2]
#  type=parms[3]


  #differential equation describing the concentration in the reference tissue
  #input function taken from...

  dydt= K1c*ipt-k2c*y

  return(list(c(dydt)))
}

#' Two-tissue compartment model
#'
#' This function produce the output of the differential equation describing the
#' concentration in a two-tissue compartment model. First compartment is for the input function,
#' the concentration of tracer in plasma or blood curve as a function of time C0(t).
#' The following two compartments are for the two distinct kinetic compartments
#' in the tissue (C1(t) and C2(t)), representing, for instance, the free
#' and receptor-bound tracer concentrations.
#' Compartments are connected with four rate constants, K1, k2, k3, and k4.
#'
#' @param t time
#' @param yobs matrix or dataframe of measured radioactivity concentration in a
#' voxel for each time-point. This can be in terms of observations or
#' directly in terms of voxel time activity curve. The
#' first column refers to the observations from the
#' first tissue, the second column refers to the observations from the
#' second tissue.
#' @param parms vector of three elements: K1, constant for the flux from blood to
#' tissue, k2, constant for the flux from tissue to blood, k3, constant
#' for the flux from the first tissue to the second, k4, constant for the flux
#' from the second to the first tissue.
#' @param inputfunction input function. Default inputfunction()
#' @param type type of the input function.
#' @return A list with the output of the differential equation
#' @keywords BayesPET
#' @export
twotissue = function(t, yobs, parms, inputfunction.=inputfunction,
                     type=2) {
  #differential equation describing the concentration in the target tissue
  #based on a 2 tissue compartment model
  ipt <- inputfunction(t, type)
  y   <- yobs[,1:2]
  K1  <- parms[1]
  k2  <- parms[2]
  k3  <- parms[3]
  k4  <- parms[4]

  dydt=matrix(NA,ncol=2,nrow=length(ipt))
  dydt[,1]= ((K1*ipt)-(k2*y[,1])-(k3*y[,1])+(k4*y[,2]))
  dydt[,2]=(k3*y[,1]-k4*y[,2])

  return(list(dydt))
}

#' Simulation for single-tissue compartment model
#'
#' Time Activity Curves Simulations with Normal distributed and Poisson noise for
#' the single-tissue compartment model.
#'
#' @param K1 constant for the flux from the blood to the tissue. Default to 0.0918.
#' @param k2 constant for the flux from the tissue to the blood. Default to 0.4484.
#' @param noise type of error noise. Default Poisson (Pois), alternative
#' Gaussian ("norm") and Poisson-normal ("mixed")
#' @param PLOT if a plot has to be produced. Default FALSE.
#' @param l2 noise level for Poisson model. Default decreasing 1:20.
#' @param l1 noise level for the normal model. Default increasing 1:10.
#' @param tspan vector of time points. Default 1:60.
#' @param inputfunction input function. Default inputfunction()
#' @param type type of the input function.
#' @return A list with two vectors with different type of noise.
#' @keywords BayesPET
#' @export
TAC_One_Compartment=function(K1=0.0918, k2=0.4484, noise="Pois", PLOT=F,
                             l2=1, l1=3,tspan=1:60,
                             inputfunction.=inputfunction,
                             type=2) {
#Time Activity Curves Simulations with Normal distributed and Poisson noise
#One compartment model l2 (noise level for Pois, decreasing 1:20),
#  l1 (noise level for Normal, increasing 1:10)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Parameters Values from Yoder et al. 2004
  #K1=0.0918
  #k2=0.4484

  Ct0=0 #initial condition : concentration=0 at t=0

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

  parms=c(K1=K1, k2=k2,type=type)
  out=ode(Ct0, tspan, singtissue, parms, method="ode45")
#  parms=c(K1, k2, type)
#  out=ode(Ct0, tspan, singtissue, parms, inputfunction, method="ode45")
  Ct=t(out[,2])
  #solving the equation over the time interval tspan [0:1:50] using ode45
  #input function described in the function, second argument is the type

  if (noise=="norm") {
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Adding Normally distributed noise
  #Level of noise from Gunn et al. 2005
  nmin=0.001
  nmax=0.1
  n=10^(seq(log10(nmin),log10(nmax),length=10) )

  noiselevel=n[l1] #choosinglevel of noise for normally distributed noise
  x=rnorm(length(Ct))
  ynoise1=as.vector(Ct+noiselevel*x) #adding noise, randn normally distributed random numbers
  #cat("norm", ynoise1, "\n")

  if (PLOT==T) {
    par(mfrow=c(1,2))
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(tspan,ynoise1, type="l", col=2, xlab="time", ylab="Concentration of the radiotracer")
    lines(tspan, Ct, lty=2, col=3)
    title("TAC Target tissue")
    } #plot
  return(ynoise1)
  }

  if (noise=="Pois") {
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Adding poisson noise

  ynoise1=matrix(NA, ncol=1, nrow=length(Ct))
  ll=length(Ct)

  counts=10^(seq(3,6,length=20))

  calib=counts[l2]/sum(Ct) #choosing level of noise

  ynoise1=rpois(ll, calib*Ct)/calib #adding poisson noise

  if (PLOT==T) {
    par(mfrow=c(1,2))
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(tspan,ynoise1, type="l", col=2, xlab="time", ylab="Concentration of the radiotracer")
    lines(tspan, Ct, lty=2, col=3)
    title("TAC Target tissue")
    } #plot
  return(ynoise1)
  }  #pois

  if(noise=="mixed"){
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  #Adding normally distributed noise to poisson noise
  ynoise1=ynoise2=matrix(NA, ncol=1, nrow=length(Ct))
  ll=length(Ct)

  poissonN=10^(seq(3, 6, length=20))

  calib=poissonN[l2]/sum(Ct)  #choosing level of noise

  ynoise1=rpois(ll,calib*Ct)/calib  #adding poisson noise

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Adding Normally distributed noise
  #Level of noise from Gunn et al. 2005
  nmin=0.001
  nmax=0.1
  n=10^(seq(log10(nmin),log10(nmax),length=10) )

  X=rnorm(length(ynoise1))
  #ynoise2=Ct+noiselevel*X   #adding noise, poisson and normally distributed
  # ynoise2=Ct+n[nL]*X
  ynoise2=ynoise1+n[l1]*X

  if (PLOT==T) {
    par(mfrow=c(1,2))
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(tspan,ynoise1, type="l", col=2, xlab="time", ylab="Concentration of the radiotracer")
    lines(tspan, Ct, lty=2, col=3)
    title("TAC Target tissue")

    plot(tspan,ynoise1, type="l", col=2, xlab="time", ylab="Concentration of the radiotracer")
    lines(tspan,ynoise2, lty=2, col=3 )
    lines(tspan, Ct, lty=3, col=4 )
    title("TAC Target tissue")
  } #plot

  return(list(ynoise1=ynoise1, ynoise2=as.vector(ynoise2)))
  } #mixed

}

#' Simulation for two-tissue compartment model
#'
#' Time Activity Curves Simulations with Normal distributed and Poisson noise for
#' the two-tissue compartment model.
#'
#' @param K1 constant for the flux from the blood to the tissue. Default to 0.0918.
#' @param k2 constant for the flux from the tissue to the blood. Default to 0.4484.
#' @param k3 constant for the flux from the first to the second tissue. Default to 0.0482.
#' @param k4 constant for the flux from the second to the first tissue. Default to 0.1363.
#' @param noise type of error noise. Default Poisson (Pois), altnernative Gaussian ("norm"), or Poisson-normal ("mixed")
#' @param PLOT if a plot has to be produced. Default FALSE.
#' @param l2 noise level for Poisson model. Default decreasing 1:20.
#' @param l1 noise level for the normal model. Default increasing 1:10.
#' @param tspan vector of time points. Default 1:60
#' @param inputfunction input function. Default inputfunction()
#' @param type type of the input function.
#' @return A list with two vectors with different type of noise.
#' @keywords BayesPET
#' @export
TAC_2_Compartment=function(K1=0.0918, k2=0.4484, k3=0.0482, k4=0.1363,
                           noise="Pois", PLOT=F, l2=1, l1=3,
                           tspan=1:60,inputfunction.=inputfunction,
                           type=2) {
  #Time Activity Curves Simulations with Normal distributed and Poisson noise
  #Two compartment model l2 (noise level for Pois, decreasing 1:20), l1 (noise level for Normal, increasing 1:10)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Parameters Values from Yoder et al. 2004
  #K1=0.0918
  #k2=0.4484
  #k3=0.0482
  #k4=0.1363

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

  y0=c(0,0) #+initial condition : cncentration=0 at t=0
  parms=c(K1,k2,k3,k4,type)

  y= ode(y0, tspan, twotissue, parms,method="ode45")
  #solving the equation over the time interval tspan [0:1:60] using ode45
  #input function described in the function, second argument is the type

  #total concentration measured in the region of interest: Ct
  Ct= y[,2] + y[,3]

  if (noise=="norm") {
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #Adding Normally distributed noise
    #Level of noise from Gunn et al. 2005
    nmin=0.001
    nmax=0.1
    n=10^seq(log10(nmin),log10(nmax),length=10)

    noiselevel=n[l1]#choosinglevel of noise
    x=rnorm(length(Ct))
    ynoise1=Ct+noiselevel*x #adding noise, randn normally distributed random numbers

    if (PLOT==T) {
      par(mfrow=c(1,2))
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      plot(tspan,ynoise1, xlab="time", yla="Concentration of the radiotracer", type="l", col=2)
      lines(tspan, Ct, lty=3, col=3)
      #plot(t,ynoise1(:,j),'r',t,Ct,'g');
      title("TAC Target tissue Two-compartment model")

    } #plot

    return(ynoise1)
  }

  if (noise=="Pois") {
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #Adding poisson noise - 20 different levels of noise from 0.001 to 1.0
    ynoise1=ynoise2=matrix(NA, ncol=1, nrow=length(Ct))
    ll=length(Ct)

    counts=10^(seq(3,6,length=20))

    calib=counts[l2]/sum(Ct) #choosing level of noise

    ynoise1=rpois(ll, calib*Ct)/calib #adding poisson noise

    if (PLOT==T) {
      par(mfrow=c(1,2))
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      plot(tspan,ynoise1, xlab="time", yla="Concentration of the radiotracer", type="l", col=2)
      lines(tspan, Ct, lty=3, col=3)
      #plot(t,ynoise1(:,j),'r',t,Ct,'g');
      title("TAC Target tissue Two-compartment model")

    } #plot

    return(ynoise1)

  }

  if (noise=="mixed"){
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #Adding normally distributed noise to poisson noise
    #Level of noise from Gunn et al. 2005
    ynoise1=ynoise2=matrix(NA, ncol=1, nrow=length(Ct))
    ll=length(Ct)

    counts=10^(seq(3,6,length=20))

    calib=counts[l2]/sum(Ct) #choosing level of noise

    ynoise1=rpois(ll, calib*Ct)/calib #adding poisson noise

    nmin=0.001
    nmax=0.1
    n=10^(seq(log10(nmin),log10(nmax),length=10) )

    X=rnorm(ll)
    #ynoise2=Ct+noiselevel*X  #adding noise, poisson and normally distributed
    ynoise2=ynoise1+n[l1]*X

    if (PLOT==T) {
      par(mfrow=c(1,2))
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      plot(tspan,ynoise1, xlab="time", yla="Concentration of the radiotracer", type="l", col=2)
      lines(tspan, Ct, lty=3, col=3)
      #plot(t,ynoise1(:,j),'r',t,Ct,'g');
      title("TAC Target tissue Two-compartment model")

      plot(tspan,ynoise1, type="l", col=4,xlab="time", yla="Concentration of the radiotracer")
      lines(tspan,ynoise2,lty=2, col=5)
      lines(tspan,Ct,lty=3, col=6)
      title('TAC Target tissue Two-compartment model')
    } #plot

    return(list(ynoise1=ynoise1, ynoise2=as.vector(ynoise2)) )
  }

}
