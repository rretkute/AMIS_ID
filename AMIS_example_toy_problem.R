###############################################
#
#  Retkute et al. 2020
#  "Integrating geostatistical maps and transmission models using 
# multiple impotance sampling
#  Code for a toy model with banana shape target
#
###############################################



library(tmvtnorm)
library(mnormt)
library("mclust")
library(ggplot2)
library(ggpubr)
require(gridExtra)

source("AMIS_source.r")

###################################################################	
#          Model and prevalences
####################################################################	

# Banana shape
# Models takes parameter values and returns prevalence value in [0,1].
model<-function(x){
  sig<-100; b<-0.1
  if (is.vector(x)) x <- t(as.matrix(x))
  p <- ncol(x)
  x[,2] <- x[,2]+b*x[,1]^2-b*sig
  
  Mu <- rep(0,p)
  Sigma <- diag(p)
  Sigma[1,1] <- sig
  
  dmnorm(x,mean=Mu,varcov=Sigma,log=FALSE)/0.01591549
}

# Test model output
# model(c(0.1,0.1))
n.param<-2 #  number of parameters

n.pixels<-100  # Number of pixels
n.map.sampl<-1000 # Number of sammples for each pixel
ESS.R<-200 # Desired effective sample 
delta<-5 # delta value (width for the Radon-Nikodym derivative)

# Make synthetic map ~ bimodal distribution
# Alternatively load prevalences as a matrix and set n.pixels accordingly
prev<-matrix(NA, nrow=n.pixels, ncol=n.map.sampl)
for(i in 1:n.pixels){
	 if(i<=round(n.pixels/2)){
	    prev[i,]<-rtmvnorm(n.map.sampl, mean=25, sigma=100*runif(1), lower=1, upper=100)
	 } else {
	   	prev[i,]<-rtmvnorm(n.map.sampl, mean=75, sigma=100*runif(1), lower=1, upper=100)
	 }
}

hist(prev, 100)
mean.prev<-sapply(1:n.pixels, function(a) mean(prev[a,])) # Mean prevalence at each pixel


###################################################################	
#          AMIS setup
####################################################################

# Set distribution for proposal: Student's t distribution
proposal=mvtComp(df=3); mixture=mclustMix(); 
dprop <- proposal$d
rprop <- proposal$r

# Set prior distribution: uniform
dprop0<-function(a,b){
  return(dunif(a, min=-30, max=30)* dunif(b, min=-50, max=30))
}
rprop0<-function(n){
  return(list(runif(n, min=-30, max=30), runif(n, min=-50, max=30)))
}


T<-200; # max number of iterations
NN<-1000  # Number of parameter sets in each iteration
N<-rep(NN,T)  # This allows to have different number of parameters sampled each iteration. Here it's the same

param<-matrix(NA, ncol=n.param+1, nrow=T*NN)  # Matrix for parameter values + corresponding prevalence 
Sigma <- list(NA, 10*T)
Mean<-list(NA, 10*T)
PP<-list(NA,T)
GG<-list(NA,T)

###################################################################	
#          Iteration 1. 
####################################################################	


t<-1  # Iteration 
tmp<-rprop0(N[t]) # Sample from prior
x<-tmp[[1]]
y<-tmp[[2]]


# Run the model
sim<-c()
for(i in 1:N[t]){
	sim<-c(sim, model(c(x[i], y[i])))
}

# Make prevalence from into percents
sim<-100*sim

param[1:N[1],1]<-x
param[1:N[1],2]<-y
param[1:N[1],3]<-sim

# Calculate weights
# For iteration 1, proposal==prior, so weight ratio  is 1
w1<-rep(1, length(sim))

tmp<-get.WW.and.ESS(prev, sim, w1)

WW<-tmp[[1]]
ess<-tmp[[2]]

ESS<-matrix(ess, nrow=1, ncol=n.pixels)

# Visualise results
pp<-data.frame(x=param[1:sum(N[1:(t)]),1], y=param[1:sum(N[1:(t)]),2],  prevalence=param[1:sum(N[1:(t)]),3])
pp<-pp[order(pp$prevalence),]
f1<-ggplot(pp, aes(x,y, colour = prevalence))+   geom_point()  +scale_color_gradientn(colours = rainbow(5))
f2<- ggplot(pp, aes(prevalence)) + geom_histogram()
xx<-data.frame(mean.prevalence=mean.prev, ESS=ESS[nrow(ESS),])
f3<-qplot(mean.prevalence, ESS, data = xx)
grid.arrange(f1, f2, f3, ncol=3, widths=c(1.25,1,1))


###################################################################	
#          Iteration 2+
####################################################################	


stop<-0
while(stop==0){

	t<-t+1
	cat(c("Iteration: ", t,", min(ESS): ", min(ess),"\n"))

	wh<-which(ess>=ESS.R)
	W1<-WW; W1[wh,]<-0

	w1<- c(colSums(W1))


	J<-sample(1:sum(N[1:(t-1)]), NN, prob= w1, replace=T)
	xx<-param[J,1:2]
	clustMix <- mixture(xx)

	G <- clustMix$G
	cluster <- clustMix$cluster

	### Components of the mixture
	ppt <- clustMix$alpha
	muHatt <- clustMix$muHat
	varHatt <- clustMix$SigmaHat
	GG[[t-1]]<-G
	G1<-0; G2<-G
	if(t>2) {
		G1<-sum(sapply(1:(t-2), function(a) GG[[a]]))
		G2<-sum(sapply(1:(t-1), function(a) GG[[a]]))
	}
	for(i in 1:G){
		Sigma[[i+G1]] <- varHatt[,,i]
		Mean[[i+G1]] <- muHatt[i,]	
		PP[[i+G1]]<-ppt[i]
	}
			
	### Sample new from the mixture...
ans<-c(); x<-c(); y<-c()
while(length(ans)<N[t+1]){
  compo <- sample(1:G,1,prob=ppt) 
  x1 <- t(rprop(1,muHatt[compo,], varHatt[,,compo]))
  new.param<-as.numeric(x1)
  if(dprop0(new.param[1],new.param[2])>0){	
    ans<-c(ans, model(new.param))
		x<-c(x, new.param[1])
		y<-c(y, new.param[2])
	}
	i<-i+1
}

param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),1]<-x
param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),2]<-y
param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),3]<-100*ans

sim<-param[1:sum(N[1:(t)]),3]

sim<-param[1:sum(N[1:t]),3]
w1 <- sapply(1:sum(N[1:t]), function(b)  dprop0(param[b,1], param[b,2]))/(sapply(1:sum(N[1:t]), function(b)  dprop0(param[b,1], param[b,2]) + sum(sapply(1:G2, function(a) PP[[a]] * dprop(param[b,1:2],mu= Mean[[a]], Sig=Sigma[[a]])))))

tmp<-get.WW.and.ESS(prev, sim, w1)

WW<-tmp[[1]]
ess<-tmp[[2]]

cat( c("min(ESS)=", min(ess),  ", max(ESS)=", max(ess), "\n"))

ESS<-rbind(ESS, as.numeric(ess))

if(min(ess)>=ESS.R) stop<-1
if(t>= T) stop<-1
# Visualise resylts
pp<-data.frame(x=param[1:sum(N[1:(t)]),1], y=param[1:sum(N[1:(t)]),2],  prevalence=param[1:sum(N[1:(t)]),3])
pp<-pp[order(pp$prevalence),]
f1<-ggplot(pp, aes(x,y, colour = prevalence))+   geom_point()  +scale_color_gradientn(colours = rainbow(5))
f2<- ggplot(pp, aes(prevalence)) + geom_histogram()
xx<-data.frame(mean.prevalence=mean.prev, ESS=ESS[nrow(ESS),])
f3<-qplot(mean.prevalence, ESS, data = xx)
grid.arrange(f1, f2, f3, ncol=3, widths=c(1.25,1,1))
}


