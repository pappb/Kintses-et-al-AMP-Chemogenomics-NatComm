#
# This program generates a lot of figures to survey the data, but has no effect to the rest of the data analysis
#
####################################
# This program generates one barplot for each of the 4441 genes
#



rm(list=ls())

source("script/functions-for-chemogenomic-project.R");

batch.df <- read.csv(file="data/batch-structure.csv",sep=",")

M0=read.csv(file="data/dataTable1-read_numbers_of_sequencing.csv",sep="," )
rownames(M0)<-M0$Feature_ID
names(M0)<-sapply(names(M0), FUN=function(x)  sprintf("%s_B%i",x, batch.df$batch.num[ batch.df$cname==x]))


M1=read.csv(file="data/dataTable3-afterRescaling.csv",sep=",",row.names=1)
names(M1)<-sapply(names(M1), FUN=function(x)  sprintf("%s_B%i",x, batch.df$batch.num[ batch.df$cname==x]))

names1<-my.order(sort(names(M1)), c( "NT_"))
M1<-M1[,names1]
M0<-M0[,names1]

M0<-M0[rownames(M0)!="lacI" & rownames(M0)!="lacZ",]
M0<-as.data.frame(as.matrix(M0) %*% diag(1/colSums(M0)))
names(M0)<-names1

stopifnot(all(row.names(M0)== row.names(M1)))



tr1<-gsub("^([^_]+)_\\d_B\\d","\\1",names(M1))

walls1<-which(c("xx",tr1)!= c(tr1,"xx"))-1 # positions of separator vertical lines on the barplots

tmp<-M1

graphics.off()
out.dir1<-"out/cg04-raw-barplots"
if(!file.exists(out.dir1)){dir.create(out.dir1	, recursive = TRUE)}
unlink(paste0(out.dir1,"/*"),recursive = TRUE, force = FALSE)


for(i in 1:nrow(tmp)){
	name1<- row.names(tmp)[i]
	png(filename=paste(out.dir1,"/",name1,".png", sep=""), width=1280, height=720, pointsize = 24)
	par(mar=c( 5.1, 4.1, 2.1, 2.1))
	x<-as.numeric(tmp[i,])
	x[x<=0] <- 1e-10

	color1<-rep(rgb(0.85,0.85,0.85,alpha = 0.5),ncol(M1))  #"gray85"
	
	y<-as.numeric(M0[i,])
	y[y<=0] <- 1e-10
	
	barplot(y, names.arg=NULL ,space=c(0.2, rep(0.2,length(x)-1)) , horiz=FALSE,  log="y", ylim=c(5e-6, 5e-2),  main=name1,  ylab="log of proportion of read-numbers", col="gray30", axes=FALSE)
	barplot(x, names.arg=names(M1) ,space=c(0.3, rep(0.2,length(x)-1)), horiz=FALSE, las=2 ,cex.names=0.8, log="y", ylim=c(5e-6, 5e-2),   col=color1 , add=TRUE) 
	abline(v=0.1+walls1*1.2, col="lightgray")
	dev.off()
	
	
}

#############################
