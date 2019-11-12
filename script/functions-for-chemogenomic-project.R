## file of helper functions



estimate_mode <- function(x) {
	x<-x[is.finite(x)]
	d <- density(x)
	d$x[which.max(d$y)]
}


# this function is designed for reorder data columns. The key.list gives regexp pattarns, and 
# the return value will be ordered that the firsts match to the regexp
my.order<-function(orig.names, key.list)
{

	reordered.names=c();
	for(pattern in key.list){
		idx=grepl(pattern, orig.names)
		reordered.names=c(reordered.names,orig.names[idx])
		orig.names=orig.names[!idx]
	}
	
	reordered.names=c(reordered.names,orig.names)
	return(reordered.names)
}





# this function can write a text upon a plot aligned to a line,
# the line is givver by the lm1 parameter, what is the result of the lm() function. 
# The x is in parameter. The y and the angle is calculated automatically. The rest of the parameters are forwarded to the text() function
#
# example:
# tex.on.line(lm1=lm1,x=5.3,labels="svalamiss")
tex.on.line<-function(lm1, x, ...){

	#x=4.5
	uy <- diff(grconvertY(1:2,"user","inches"))
	ux <- diff(grconvertX(1:2,"user","inches"))
	asp<-	uy/ux
	angle=180/pi*atan(lm1$coefficients[2]*asp)
	text(x=x, y=lm1$coefficients[1]+lm1$coefficients[2]*x, srt=angle, adj=c(0.5,0) , ...)
}


# writes text with a little shadow. It is easier to read on the plot.
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
		theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
	
	xy <- xy.coords(x,y)
	xo <- r*strwidth('x')
	yo <- r*strheight('x')
	
	for (i in theta) {
		text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
	}
	text(xy$x, xy$y, labels, col=col, ... )
}





# mixture of log and linear function.
# if x> treshold then it is a log,
# if x<= treshold the it is a linear function. 
# The slope and y parameters of the linear parts are determined to be smooth at the treshold.
log10LinearHybrid<-function(x, treshold=1e-4){

	idx<-(x>treshold)
	
	x[idx]<-log10(x[idx])
	
	x[!idx]<-( x[!idx]/treshold + log(treshold)-1)/ log(10)
	
	return(x)
}

#inverse oflog10LinearHybrid()-
exp10LinearHybrid<-function(y, treshold=1e-4){
	treshold2<-log10(treshold)
	
	idx<-(y>treshold2)
	
	y[idx]<-10^(y[idx])
	
	# y= a*x+b
	# y= 1/(treshold*log(10)) * x + (log(treshold)-1)/ log(10)
	#
	#inverse:
	# x= (y-b)/a
	# x= (y - (log(treshold)-1)/ log(10))*treshold*log(10)
	# x= ( y*log(10)  - (log(treshold)-1) ) * treshold
	
	y[!idx]<- ( y[!idx]*log(10) - (log(treshold)-1) ) * treshold
	
	
	return(y)
}

# 
# if x > treshold then it is a log10(x)
# if x in [-treshold,+treshold] then it is linear
# if x < -treshold then it is a shifted version log10(-x)
# the function s smooth at both  x=treshold and x=-treshold position.
#
# It does not result huge negative values at negative numbers like the log10LinearHybrid()
myLog10LinearHybrid<-function(x, treshold=10^-4.5)
{
	#value.at.0<-(log(treshold)-1)/ log(10)
	value.at.minus.treshold <- (log(treshold)-2)/ log(10)
	value.at.plus.treshold <-  log10(treshold)

	idx1 <- (x> treshold)
	idx2 <- (x< -treshold)
	idx3 <- !(idx1 | idx2)

	x[idx1]<- log10(x[idx1])
	x[idx2]<- value.at.minus.treshold+value.at.plus.treshold- log10(-x[idx2])

	x[idx3]<-( x[idx3]/treshold + log(treshold)-1)/ log(10)

	return(x)
}


#############

myExp10LinearHybrid<-function(y, treshold=10^-4.5){

	value.at.minus.treshold <- (log(treshold)-2)/ log(10)
	value.at.plus.treshold <-  log10(treshold)


	idx1 <- (y> value.at.plus.treshold)
	idx2 <- (y< value.at.minus.treshold)
	idx3 <- !(idx1 | idx2)

	y[idx1]<- 10^(y[idx1])
	y[idx2]<- -10^((value.at.minus.treshold+value.at.plus.treshold)-y[idx2])


	# y= a*x+b
	# y= 1/(treshold*log(10)) * x + (log(treshold)-1)/ log(10)
	#
	#inverse:
	# x= (y-b)/a
	# x= (y - (log(treshold)-1)/ log(10))*treshold*log(10)
	# x= ( y*log(10)  - (log(treshold)-1) ) * treshold

	y[idx3]<- ( y[idx3]*log(10) - (log(treshold)-1) ) * treshold


	return(y)
}
