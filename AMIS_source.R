
H <-function(x) as.numeric(x>0)


get.WW.and.ESS<-function(prev, sim, w1){
  n.sim<-length(sim)
  n.pixels<-nrow(prev)
  n.samples<-ncol(prev)
  WW.cnt<-matrix(0, nrow=n.pixels, ncol=n.sim)
  g<-rep(0, n.sim)
  for(j in 1:n.sim){
    x1<-H(delta/2-abs(prev-sim[j]))
    x2<-matrix(x1, ncol=n.samples, nrow=n.pixels)
    x3<-rowSums(x2)
    WW.cnt[,j]<-x3
    g[j]<-sum(w1[which(abs(sim-sim[j])<=delta/2)])/sum(w1)
  } 
  
  ess<-c()
  WW<-matrix(0, nrow=n.pixels, ncol=n.sim)
  for(i in 1:n.pixels){
    f<-WW.cnt[i,]
    ww<-w1*(f/g)  
    if(sum(ww) >0)
      ww<-ww/sum(ww)
    WW[i,]<-ww
    if(sum(ww)>0) {
      www<-(sum((ww)^2))^(-1)
    } else {
      www<-0
    }
    ess[i]<- www
  }
  
  return(list(WW, ess))
}


mvtComp<-function(df=3){
	list("d"=function(xx,mu=rep(0,ncol(xx)),Sig=diag(1,ncol(xx),ncol(xx)),log=FALSE){
		dmt(xx,mean=mu,S=Sig,df=df,log=log)
		},"r"=function(n=1,mu=0,Sig=1){
			rmt(n,mean=mu,S=Sig,df=df)
			})
}

mclustMix<-function(G=1:10){
	
	if (any(as.numeric(G)<=0)) stop("G must be positive")
	if(!require("mclust")) stop("Package mclust is not installed...")
	
	function(xx){
		
		clustering <- fitMclust(xx,modelName="VVV",G= G)
		
		G <- clustering$G
		
		if(G==1) clustering$parameters$pro <- 1

		return(list(alpha=clustering$parameters$pro, muHat=t(clustering$parameters$mean), SigmaHat=clustering$parameters$variance$sigma,G=G,cluster=clustering$classification))
		}		
}



fitMclust<-function(xx,modelName="VVV",G= G){
	
	options(warn=-1)
	
	control <- emControl(eps=sqrt(.Machine$double.eps))
	
	n <- nrow(xx)
	p <- ncol(xx)
	
	clustering <-Gout <- BIC <- NA
	
	if (G[1] == 1) {
		clustering <- mvn(modelName = modelName, data = xx)
		BIC <- bic(modelName=modelName,loglik=clustering$loglik,n=n,d=p,G=1)
		Gout <- 1
		G <- G[-1]
		}
	
	if (p != 1) {
		if (n > p) {
			hcPairs <- hc(modelName="VVV",data=xx)
			}else {
				hcPairs <- hc(modelName="EII",data=xx)
				}
		}else hcPairs <- NULL
	if (p > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
	
	for (g in G) {
		
		if (p > 1 || !is.null(hcPairs)) {
			cl <- clss[, as.character(g)]
			}else {
				cl <- .qclass(data[subset], as.numeric(g))
				}
		z <- unmap(cl, groups = 1:max(cl))

		new <- me(modelName=modelName,data=xx,z=z,control=control)
		
		if(!is.na(new$loglik)){
			
			BICnew <- bic(modelName=modelName,loglik=new$loglik,n=n,d=p,G=g,equalPro=control$equalPro)
			
			if(is.na(BIC)){
				clustering <- new
				BIC <- BICnew
				Gout <- g
				}else{
					if(BICnew>BIC){
						clustering <- new
						BIC <- BICnew
						Gout <- g
						}
					}
			}
		}
	
	options(warn=0)
	
	return(c(clustering,G=Gout))
		
	}


