rm(list=ls())

source(file="script/functions-for-chemogenomic-project.R");
source(file="script/pdfTools.R");


M1=read.csv(file="data/dataTable3-afterRescaling.csv",sep=",",row.names=1)

mu1<- myLog10LinearHybrid(0)
mu2<- myLog10LinearHybrid(1/nrow(M1))
sigma=(mu2-mu1)/5;



treshold=10^-4.5
value.at.minus.treshold <- (log(treshold)-2)/ log(10)
value.at.plus.treshold <-  log10(treshold)
	


#####################################################################################################
#####################################################################################################
#####################################################################################################


transform.param2<-function( param)
{
	
	p<-list();
	
	p$mean1 <- param[1]
	p$mean2 <- param[2]
	
	p$alpha1 <- 1/(1+exp(-param[3]))
	left3<-p$alpha1*p$mean1+(1-p$alpha1)*p$mean2
	p$alpha2 <- 1/(1+exp(-param[9]))
	right3<-p$alpha2*p$mean1+(1-p$alpha2)*p$mean2
	p$mean3 <- (right3+left3)/2
	p$pl <- abs(right3-left3)
	
	p$sigma1 <- 1/(1+exp(-param[4]))
	p$sigma2 <- 1/(1+exp(-param[5]))
	p$sigma3 <- min(c(p$sigma1,p$sigma2)) + 1/(1+exp(-param[6]))
	
	
	q4 <- 1/(1+exp(-param[7]))  # sigmoid function. target range:[0,+1]
	p$q3 <- 1-q4
	q.tmp<- 1/(1+exp(-param[8]))
	p$q1 <- q4*q.tmp  
	p$q2 <- q4*(1-q.tmp)
	
	
	if(p$mean2 < p$mean1){
		p2<-p;
		p2$mean1 <- p$mean2
		p2$mean2 <- p$mean1
		p2$sigma2 <- p$sigma1
		p2$sigma1 <- p$sigma2
		p2$q1<-p$q2
		p2$q2<-p$q1
		
		p<-p2;
	} 
	
	return(p)
}


f2<-function(y, param)
{
	p <- transform.param2( param)
	

	logLiklehood<-mean(log(
					p$q1*dnorm(y, mean = p$mean1, sd = p$sigma1, log = FALSE)+
							p$q2*dnorm(y, mean = p$mean2, sd = p$sigma2, log = FALSE)+
							p$q3*mixDistributipnPdf(y, mu = p$mean3, sigma = p$sigma3, pl=p$pl)
			), na.rm=TRUE)
	
	return(-logLiklehood);
	
}


fit.mixed.normal.distribution<- function(y)
{

	start.parameter<-c(mu1,mu2,-0.1,0,0,0,2,0,0.1)
	opt1<-optim(par=start.parameter, fn=function(param) f2(y,param),  method="BFGS", control =list(maxit=4000))
	
	q1<-quantile(y,probs=c(0.25,0.75))
	start.parameter<-c(q1[1],q1[2],-0.1,-0.7,-0.7,0,2,0,0.1)
	opt2<-optim(par=start.parameter, fn=function(param) f2(y,param),  method="BFGS", control =list(maxit=4000))
	
	if(opt1$value<opt2$value)
	{ 
		opt<-opt1
	}	else
	{
		opt<-opt2
	}
	
	
	p<-transform.param2(opt$par)
	
	
	
	# purge outliers
	tr1<-qnorm(1/(p$q1*length(y)) , mean = p$mean1, sd = p$sigma1, log = FALSE)
	tr2<-qnorm(1-1/(p$q2*length(y)) , mean = p$mean2, sd = p$sigma2, log = FALSE)
	# abline(v=c(tr1,tr2))
	yp <- y[y>tr1 & y<tr2] ; 
	
	start.parameter<-opt$par
	start.parameter[3] <- 0.1 #  make mean3  the mean of mean1 and mean2 
	start.parameter[9] <- -0.1 # plateau width = 0
	start.parameter[6] <- -10 # make sigma3 approx. mean of sigma1 and sigma2 
	start.parameter[7] <- 1.5 # make q3 approx 20%
	#transform.param2(start.parameter)
	
	opt3<-optim(par=start.parameter, fn=function(param) f2(yp,param),  method="BFGS", control =list(maxit=4000))
	
	
#	Sys.time()
	p<-transform.param2(opt3$par)
#	print(unlist(p))
#	abline(v=c(p$mean1,p$mean2), col="blue",lwd=3, lty=4)
	
	return(p)
}




#####################################################################################################
#####################################################################################################
#####################################################################################################

transform.param1<-function( param)
{
	
	p<-list();
	
	
	p$c <- param[[1]]
	
	p$d <- param[[2]]
	# We need to solve equalition for (a,b) , the  (mu1, mu2 ,d) are known
	#	mu1 = a+b*mu1
	#	mu2 = a+b*(mu2 - d)
	#	
	#	mu2-mu1 = b*(mu2 - mu1 -d)
	#	b=(mu2-mu1)/(mu2 - mu1 -d)
	p$b <- (mu2-mu1)/(mu2 - mu1 -p$d)
	p$a <- mu1 - p$b*mu1
	
	return(p)
}



transform.data<-function(x,p)
{
	y<-p$a+p$b*myLog10LinearHybrid(x+p$c )
	return(y)
}





######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################


M2<-data.frame(row.names=rownames(M1))


graphics.off()
out.dir1<-"out/cg05-normalisation/"
if(!dir.exists(out.dir1)){dir.create(out.dir1	, recursive = TRUE)}
unlink(paste0(out.dir1,"*"),recursive = TRUE, force = FALSE)

for(col1 in names(M1)){
	#col1<- "PXB_1"# names(M1)[7]
	print(col1)
	
	#pdf(file=paste0(out.dir1,col1,".pdf"), width=8,height=4.5)
#	png(filename=paste0(out.dir1,col1,".png"), width=1920,height=1080)
	
	
	x<-M1[,col1]	
	
	y<-myLog10LinearHybrid(x)
	
#	par(mfrow=c(5,2))
	
	#  parameter c determines the shift of the first mode
	#  parameter d determines the shift of the second mode
	
	pc<-0
	pd<-0
	for(i in 1:8)
	{
		
		png(filename=paste0(out.dir1,col1,"-A",i,".png"), width=1920,height=1080)
			
		breaks=seq(min(y)-0.05,max(y)+0.05,by=0.05)
		plot(c(), xlim=c(-8,0),ylim=c(0,500), axes=FALSE, main=col1)
		axis(1, at=-8:0)
		axis(2, las=2)
		abline(v=c(mu1,mu2), col="red")
		abline(v=c(value.at.minus.treshold,value.at.plus.treshold), col="gray")
		h1<-hist(y,breaks=breaks,  axes=FALSE, add=TRUE, col="gray")
		
		
		fit1<-fit.mixed.normal.distribution(y)
		fit1<-fit1[order(names(fit1))]
		
		x1<-seq(from=-10, to=-2, by=0.01)
		y1<-dnorm(x1, mean = fit1$mean1, sd = fit1$sigma1, log = FALSE)
		y2<-dnorm(x1, mean = fit1$mean2, sd = fit1$sigma2, log = FALSE)
		y3<-mixDistributipnPdf(x1, mu = fit1$mean3, sigma = fit1$sigma3,pl=fit1$pl )
		w<-length(y)*(max(h1$mids)-min(h1$mids))/length(h1$mids) # area covered by histogram
#	lines(x1,0.333*y1*w,lwd=2, col="pink", lty=2)
#	lines(x1,0.333*y2*w,lwd=2, col="pink", lty=2)
#	lines(x1,0.333*y3*w,lwd=2, col="pink", lty=2)
		lines(x1,fit1$q3*y3*w,lwd=4, col="cyan")
		lines(x1,fit1$q1*y1*w,lwd=2, col="blue")
		lines(x1,fit1$q2*y2*w,lwd=2, col="blue")
		
		lines(x1,(fit1$q1*y1+fit1$q2*y2+fit1$q3*y3)*w, lwd=1.5, lty=2)
		
		axis(1,at= c(fit1$mean1,fit1$mean2),col="blue", lwd.tick=3, labels=rep("",2))
		axis(1,at= fit1$mean3,col="cyan", lwd.tick=3, labels="")
		
		usr=par("usr")
		text(x=usr[1], y=usr[4],labels=paste0(sprintf("\n %s:%0.2f",names(fit1),fit1),collapse = ""),adj=c(0,1), col="blue",cex=1)
		

		is.ok <- abs(mu1-fit1$mean1)< 1e-4
		
		ap<-transform.param1(c(c= -pc,d=pd))
		lab <- sprintf("correcting LEFT peak\n step: A%i\n\nmu1-peak1 = %0.5f\n",i,mu1-fit1$mean1)
		lab <- paste0(lab,"\nactual transfom: \ny=a+b*log(x+c)  ",collapse = "")
		lab <- paste0(lab,paste0(sprintf("\n%s:%f",names(ap),ap),collapse = ""),collapse = "")
		lab <- paste0(lab,"\nis OK? :", is.ok ,collapse = "")
		text(x=usr[2], y=usr[4],labels=lab ,adj=c(1,1), col="blue",cex=1)
		
		
		
		
		dev.off()
		
		
		cat("delta:",mu1-fit1$mean1,"\n")
		if(is.ok)
		{
			break;
		}else if(abs(mu1-fit1$mean1)< 0.3 )
		{
			pc <-pc + myExp10LinearHybrid(fit1$mean1)
		}else 
		{
			pc <-pc + 0.5*myExp10LinearHybrid(fit1$mean1)
		}		
		p<-transform.param1(c(c= -pc,d=pd))
		cat("transform parameters:",sprintf("%s:%f",names(p),p), "\n")
		
#		text(x=-0.1, y=390,labels=paste0("calculated transfom\n",paste0(sprintf("\n%s:%f",names(p),p),collapse = ""),collapse = ""),adj=c(1,1), col="green",cex=1.5)
		
		y<-p$a+p$b*myLog10LinearHybrid(x+p$c )
	}
	
	
	for(i in 1:8)
	{
		
		png(filename=paste0(out.dir1,col1,"-B",i,".png"), width=1920,height=1080)
		
		breaks=seq(min(y)-0.05,max(y)+0.05,by=0.05)
		plot(c(), xlim=c(-8,0),ylim=c(0,500), axes=FALSE, main=col1)
		axis(1, at=-8:0)
		axis(2, las=2)
		abline(v=c(mu1,mu2), col="red")
		abline(v=c(value.at.minus.treshold,value.at.plus.treshold), col="gray")
		h1<-hist(y,breaks=breaks,  axes=FALSE, add=TRUE, col="gray")
		
		fit1<-fit.mixed.normal.distribution(y)
		fit1<-fit1[order(names(fit1))]
		
		x1<-seq(from=-10, to=-2, by=0.01)
		y1<-dnorm(x1, mean = fit1$mean1, sd = fit1$sigma1, log = FALSE)
		y2<-dnorm(x1, mean = fit1$mean2, sd = fit1$sigma2, log = FALSE)
		y3<-mixDistributipnPdf(x1, mu = fit1$mean3, sigma = fit1$sigma3,pl=fit1$pl )
		w<-length(y)*(max(h1$mids)-min(h1$mids))/length(h1$mids) # area covered by histogram
		lines(x1,fit1$q3*y3*w,lwd=4, col="cyan")
		lines(x1,fit1$q1*y1*w,lwd=2, col="blue")
		lines(x1,fit1$q2*y2*w,lwd=2, col="blue")
		
		lines(x1,(fit1$q1*y1+fit1$q2*y2+fit1$q3*y3)*w, lwd=1.5, lty=2)
		
		axis(1,at= c(fit1$mean1,fit1$mean2),col="blue", lwd.tick=3, labels=rep("",2))
		axis(1,at= fit1$mean3,col="cyan", lwd.tick=3, labels="")
		
		usr=par("usr")
		text(x=usr[1], y=usr[4],labels=paste0(sprintf("\n %s:%0.2f",names(fit1),fit1),collapse = ""),adj=c(0,1), col="blue",cex=1)
		
		
		is.ok <- abs(mu2-fit1$mean2)< 5e-4
		
		ap<-transform.param1(c(c= -pc,d=pd))
		lab <- sprintf("correcting RIGHT peak\n step: B%i\n\nmu2-peak2 = %0.5f\n",i,mu2-fit1$mean2)
		lab <- paste0(lab,"\nactual transfom: \ny=a+b*log(x+c)  ",collapse = "")
		lab <- paste0(lab,paste0(sprintf("\n%s:%f",names(ap),ap),collapse = ""),collapse = "")
		lab <- paste0(lab,"\nis OK? :", is.ok ,collapse = "")
		text(x=usr[2], y=usr[4],labels=lab ,adj=c(1,1), col="blue",cex=1)
		
		
		dev.off()
		
		
		cat("delta:",mu2-fit1$mean2,"\n")
		if(is.ok)
		{
			break;
		}
		else if(abs(mu2-fit1$mean2)< 0.2 )
		{
			pd <- pd + (mu2-fit1$mean2)
			
		}
		else 
		{
			pd <- pd + 0.5*(mu2-fit1$mean2)
		}		
		p<-transform.param1(c(c= -pc,d=pd))
		cat("transform parameters:",sprintf("%s:%f",names(p),p), "\n")
		
		#text(x=-0.1, y=390,labels=paste0(sprintf("\n%s:%f",names(p),p),collapse = ""),adj=c(1,1), col="blue",cex=1.5)
		
		y<-p$a+p$b*myLog10LinearHybrid(x+p$c )
	}
	
	
	M2[,col1]<-y
	
	
	
	
}

write.table(M2, file="data/dataTable5-normalisedLogReadNumbers.csv",sep=",",quote = FALSE,row.names=TRUE, col.names=NA)
