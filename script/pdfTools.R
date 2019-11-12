#  Helper functions for the normalisation
#


# it is the pdf (probability distribution function) of a special distribution:
# a combination of uniform distribution and a normal distribution.
# there is a plateau of the pdf. Th pl parameter is the length of the plateau.
mixDistributipnPdf <- function(x,mu, sigma,pl)
{

# Norm PDF		exp(- ((x-mu)^2)/(2*sigma^2)    )/(sigma*sqrt(2*pi))  

	idx1 <- x < mu-pl/2
	idx3 <- x > mu+pl/2
	idx2 <- ! (idx1 | idx3)
	
	
	s <- 1 + pl/(sigma*sqrt(2*pi))
	
	k<-(sigma*sqrt(2*pi)*s)
	
	p<-x
	
	p[idx1] <- exp(- ((x[idx1]-(mu-pl/2))^2)/(2*sigma^2)    )/k
	p[idx3] <- exp(- ((x[idx3]-(mu+pl/2))^2)/(2*sigma^2)    )/k
	p[idx2] <- 1/k
 
	return(p)
}



# We have two bell curve and we want to calculate the intersectionpoint
calculate.intersections.of.gaus.curves <- function(mu1,mu2,sigma1,sigma2,q1,q2)
{
#	 We need to solve the next equalisation for x. All of the other parameters are known.
#	
#	This is the formula of the pdf of normal distribution
#	
#	q1*exp( - (x-mu1)^2/(2*sigma1^2) ) / sigma1 = q2*exp( - (x-mu2)^2/(2*sigma2^2) ) / sigma2
#	
#	( - (x-mu1)^2/(2*sigma1^2) ) - log(sigma1) + log(q1)  = ( - (x-mu2)^2/(2*sigma2^2) ) - log(sigma2) + log(q2)
#	
#	(x-mu1)^2/(2*sigma1^2)  + log(sigma1) - log(q1) =   (x-mu2)^2/(2*sigma2^2)  + log(sigma2) - log(q2)
#	
#	(x-mu1)^2/(2*sigma1^2)  - (x-mu2)^2/(2*sigma2^2) =  log(sigma2) - log(sigma1) + log(q1) - log(q2)
#	
#	( (x-mu1)^2*sigma2^2 -(x-mu2)^2*sigma1^2 ) /(2*sigma1^2*sigma2^2) =  log(sigma2) - log(sigma1) + log(q1) - log(q2)
#	
#	(x-mu1)^2*sigma2^2 -(x-mu2)^2*sigma1^2  =  (log(sigma2) - log(sigma1) + log(q1) - log(q2) ) * (2*sigma1^2*sigma2^2)
#	
#	
#	(sigma2^2 - sigma1^2) * x^2
#	+ (-2*mu1*sigma2^2 + 2*mu2*sigma1^2) * x
#	+ mu1^2*sigma2^2 - mu2^2*sigma1^2 - (log(sigma2) - log(sigma1) + log(q1) - log(q2) ) * (2*sigma1^2*sigma2^2) = 0 
#	
#   Now it is a  quadratic equation. I can use the formula of the roots	
	
	
	A <- (sigma2^2 - sigma1^2)
	B <- (-2*mu1*sigma2^2 + 2*mu2*sigma1^2)
	C <- mu1^2*sigma2^2 - mu2^2*sigma1^2 - (log(sigma2) - log(sigma1)+ log(q1) - log(q2)) * (2*sigma1^2*sigma2^2)
	
	if(A !=0){
	
		D2 <- B*B - 4*A*C  #discriminant
		if(D2 > 0)
		{
			D <- sqrt(D2);
			r1 <- (-B + D ) / (2*A)
			r2 <- (-B - D ) / (2*A)
			
			return(c(r1,r2))
			
		}else if(D2 == 0)
		{
			r1 <- -B / (2*A)
			return(c(r1))
		}else
		{
			return(c())
		}
		
	r1<-(-B +sqrt(B*B - 4*A*C)) /(2*A)
	r2<-(-B -sqrt(B*B - 4*A*C)) /(2*A)
	
	return(c(r1,r2))
	
	} else if (B!=0){
		 # B*x + C =0
		r1 <- -C/B
		return(c(r1))
	}else{
		return(c())
	}
}	


calculate_intersections_of_Gauss_curve_and_a_line <- function(mu,sigma,q,h)
{
	# q*exp(- ((x-mu)^2)/(2*sigma^2)    )/(sigma*sqrt(2*pi)) = h
	#
	# exp(- ((x-mu)^2)/(2*sigma^2)    ) = h*sigma*sqrt(2*pi) /q
	#
	# - ((x-mu)^2)/(2*sigma^2)     =  log(h*sigma*sqrt(2*pi) /q )
	#
	#  (x-mu)^2     =  -log(h*sigma*sqrt(2*pi) /q ) * (2*sigma^2)
	#
	#  x     = mu +- sqrt( -log(h*sigma*sqrt(2*pi) /q ) * (2*sigma^2))
	
	d2 <- -log(h*sigma*sqrt(2*pi) /q ) * (2*sigma^2)
	if(d2>0)
	{
		d<- sqrt( -log(h*sigma*sqrt(2*pi) /q ) * (2*sigma^2))
		x1     = mu + d
		x2     = mu - d
		return(c(x1,x2))
	}else
	{
		return(c())
	}
	
}


calculate_intersections_of_gaus_and_mixed_curves <- function(mu1,mu2,sigma1,sigma2,q1,q2, pl2)
{
	s <- 1 + pl2/(sigma2*sqrt(2*pi))
	k <- (sigma2*sqrt(2*pi)*s)
	
#	
#	x <- seq(-8,0,by=0.001)
#	y1 <- q1*dnorm(x, mean=mu1,sd=sigma1)
#	y2L <- q2*dnorm(x, mean=mu2-pl2/2,sd=sigma2)/s
#	y2R <- q2*dnorm(x, mean=mu2+pl2/2,sd=sigma2)/s
#	points(x,y1,cex=0.1,col="blue")
#	points(x,y2L,cex=0.1,col="green")
#	points(x,y2R,cex=0.1,col="green")
	
	r1 <- calculate.intersections.of.gaus.curves(mu1,mu2-pl2/2 ,sigma1,sigma2,q1,q2/s)
	r2 <- calculate.intersections.of.gaus.curves(mu1,mu2+pl2/2 ,sigma1,sigma2,q1,q2/s)
	
	r3 <-calculate_intersections_of_Gauss_curve_and_a_line(mu1,sigma1,q1,h=q2/k)
	
	r1<-r1[r1 < mu2-pl2/2]
	r2<-r2[r2 > mu2+pl2/2]
	r3<-r3[mu2-pl2/2 <= r3 & r3 <= mu2+pl2/2]
	
	return(c(r1,r2,r3));
}


##Example:
##
##mu1=0
##mu2=2
##sigma1=2
##sigma2=2
##q1=0.5
##q2=0.7
##pl2=4
##
#
#mu1=fit1$mean2
#mu2=fit1$mean3
#sigma1=fit1$sigma2
#sigma2=fit1$sigma3
#q1=fit1$q2
#q2=fit1$q3
#pl2=fit1$pl
#
#x <- seq(-8,0,by=0.001)
#y1 <- q2*mixDistributipnPdf(x, mu=mu2,sigma=sigma2,pl=pl2)
#y2 <- q1*dnorm(x, mean=mu1,sd=sigma1)
##y3 <- dnorm(x, mean=2,sd=2)
#plot(x,y1, cex=0.1, ylim=c(0,0.7))
#points(x,y2, cex=0.1, col="red")
##points(x,y3, cex=0.1, col="purple")
#
##abline(v=calculate.intersections.of.gaus.curves(mu=1.2,mu2=2,sigma1=2,sigma2=2,q1=1,q2=1))
#
##abline(h=0.19)
##abline(v=calculate_intersections_of_Gauss_curve_and_a_line(mu=1.2,sigma=2,q=1,h=0.19))
#abline(v=calculate_intersections_of_gaus_and_mixed_curves(mu1,mu2,sigma1,sigma2,q1,q2, pl2))
#
