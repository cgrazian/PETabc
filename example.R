rm(list=ls())
library(PETabc)

data <- read.csv("ToTest-Real-1.csv",header=T)
head(data)

obj <- lp_ntPETabc(Ct=data$Ctreal,Cr=data$Cr,Ti=data$Time.min.,S=10^5)

obj2 <- lp_ntPETabc(Ct=data$Ctreal,Cr=data$Cr,Ti=data$Time.min.,S=10^5,
                    R1a=0.8,R1b=1.2,K2alpha=0.1,K2beta=0.4,K2a.alpha=0,
                    K2a.beta=0.1,gamma.a=0,gamma.b=0.2,tD.a=18,tD.b=22,
                    tP.b=40,alpha.a=0,alpha.b=3)

str(obj)

obj_MC <- MC_abc_lpntPET(Ct=data$Ctreal,Cr=data$Cr,Ti=data$Time.min.,
                         abc_out=obj,tol1=NULL,tol2=NULL,PLOT=T)

obj2_MC <- MC_abc_lpntPET(Ct=data$Ctreal,Cr=data$Cr,Ti=data$Time.min.,
                          abc_out=obj2,tol1=NULL,tol2=NULL,PLOT=T)


# In case you want to change tolerance level (by default, it uses 5% of the 
# distances)

err_act_new <- quantile(obj2$error_act,probs=0.01)
ABCout_act_new <- obj2$ABCout_act[obj2$error_act < err_act_new,]

par(oma=c(3,3,3,3)) # all sides have 3 lines of space
par(mfrow=c(3,2))

plot(density(ABCout_act_new[,1]),main="", xlab=expression(R[1]))
abline(v=mean(ABCout_act_new[,1]))

plot(density(ABCout_act_new[,2]),main="",xlab=expression(k[2]))
abline(v=mean(ABCout_act_new[,2]))

plot(density(ABCout_act_new[,3]),main="",xlab=expression(k[2][a]))
abline(v=mean(ABCout_act_new[,3]))

plot(density(ABCout_act_new[,4]),main="",xlab=expression(gamma))
abline(v=mean(ABCout_act_new[,4]))

plot(density(ABCout_act_new[,5]),main="",xlab=expression(t[D]))
abline(v=mean(ABCout_act_new[,5]))

plot(density(ABCout_act_new[,6]),main="", xlab=expression(t[P]))
abline(v=mean(ABCout_act_new[,6]))

mtext(paste("Posterior distributions, tol=",round(err_act_new,2)), outer = TRUE, cex = 1)