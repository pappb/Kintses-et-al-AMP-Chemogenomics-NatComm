#
# This program generates a lot of figures to survey the data, but has no effect to the rest of the data analysis
#

rm(list=ls())

source("script/functions-for-chemogenomic-project.R");

M1=read.csv(file="data/dataTable3-afterRescaling.csv",sep=",",row.names=1)

batch.df <- read.csv(file="data/batch-structure.csv",sep=",")


########################################################################################
########################################################################################
########################################################################################


graphics.off()
out.dir1="out/cg03-pie_charts-afterBgRemoveAndNormToOne"
unlink(out.dir1, recursive = TRUE)
dir.create(out.dir1, recursive = TRUE)

name.list1<-names(M1)
row.idx<-rep(TRUE,nrow(M1))
for(column.name1 in name.list1)
{
	batch.label<-batch.df$batch.num[ batch.df$cname==column.name1]
	
	png(filename=paste(out.dir1,"/",column.name1,".png", sep=""), width=800, height=750, pointsize=24)
	par(mar=c(0,0,2,0))
	x<-M1[row.idx,column.name1]
	x[x<0]<-0
	names1<-as.character( rownames(M1)[row.idx])
	names1[x<2*quantile(x,0.99)]<-""
	pie( x, labels=names1, main = c(column.name1,paste0("(batch: ",batch.label,")")))
	dev.off()
}

rm(column.name1,name.list1,names1,out.dir1, row.idx, x,batch.label)    



########################################################################################
########################################################################################
########################################################################################

color1<-function(x) {
	if(x!=0) {return ("blue");}else{return ("gray");}
}

out.dir1="out/cg03-hist-logScale-afterBgRemoveAndNormToOne"
unlink(out.dir1, recursive = TRUE)
dir.create(out.dir1, recursive = TRUE)

names1<-names(M1) 
range1<-range(unlist(log10(M1)), finite=TRUE)
for(t1.col in names1 )
{
	batch.label<-batch.df$batch.num[ batch.df$cname==t1.col]
	
	
	pdf(file=paste(out.dir1,"/",t1.col,".pdf", sep=""), width=6, height=3.5)
	
	t1<-M1[,t1.col]
	names(t1)<-row.names(M1)
	t1<-log10(t1)
	cnt.inf=sum(t1==Inf, na.rm=TRUE)
	cnt.minus.inf=sum(t1==-Inf, na.rm=TRUE)
	cnt.nan=sum(is.nan(t1))
	cnt.na=sum(is.na(t1))-cnt.nan
	
	t1<-t1[is.finite(t1)]
	
	
	
	cnt.of.nonfinties<-nrow(M1)-length(t1)
	
	breaks=seq(range1[1]-0.1,range1[2]+0.1,by=0.05)
	hist(t1, main=sprintf("%s (Batch: %i)",t1.col,batch.label), breaks=breaks, xlim=range1 , xlab="log10(read number)")
	
	N=6
	t2<-head(sort(t1, decreasing = TRUE),N)
	text(x=t2,y=jitter(rep(25,N),amount=25) ,labels=names(t2),srt=90, cex=0.7,col=1+(1:N),adj=c(0,0.5))
	
	
	
	usr1<-par("usr")
	x=usr1[1]
	y=mean(usr1[c(3,4)]);
	cex=0.6
	text(x=x, y=y, label=sprintf(" -Inf cnt: %i",cnt.minus.inf ) , adj=c(0,0.5) , col=color1(cnt.minus.inf), cex=cex )
	text(x=x, y=y, label=sprintf("\n +Inf cnt: %i",cnt.inf ) , adj=c(0,0.5) , col=color1(cnt.inf), cex=cex )
	text(x=x, y=y, label=sprintf("\n\n\n  NA  cnt: %i",cnt.na ) , adj=c(0,0.5) , col=color1(cnt.na), cex=cex )
	text(x=x, y=y, label=sprintf("\n\n\n\n  NaN cnt: %i",cnt.nan ) , adj=c(0,0.5) , col=color1(cnt.nan), cex=cex )
	
	dev.off()
}

rm(t1.col,t1, breaks, out.dir1)  









#########################


color1<-function(x) {
	if(x!=0) {return ("blue");}else{return ("gray");}
}

out.dir1="out/cg03-hist-logLinearScale-afterBgRemoveAndRescale"
unlink(out.dir1, recursive = TRUE)
dir.create(out.dir1, recursive = TRUE)

names1<- names(M1) 
range1<-range(unlist(myLog10LinearHybrid(M1)), finite=TRUE)
for(t1.col in names1 )
{
	batch.label<-batch.df$batch.num[ batch.df$cname==t1.col]
	
	
	pdf(file=paste(out.dir1,"/",t1.col,".pdf", sep=""), width=6, height=3.5)
	
	t1<-M1[,t1.col]
	names(t1)<-row.names(M1)
	t1<-myLog10LinearHybrid(t1)
	cnt.inf=sum(t1==Inf)
	cnt.minus.inf=sum(t1==-Inf)
	cnt.na=sum(is.na(t1))
	cnt.nan=sum(is.nan(t1))
	
	t1<-t1[is.finite(t1)]
	
	
	
	cnt.of.nonfinties<-nrow(M1)-length(t1)
	
	breaks=seq(range1[1]-0.1,range1[2]+0.1,by=0.05)
	hist(t1, main=sprintf("%s (Batch: %i)",t1.col,batch.label), breaks=breaks, xlim=range1 , xlab="log10(read number)")
	
	abline(v=myLog10LinearHybrid(c(0,1/length(t1)) ) , col="red" )
	

	N=6
	t2<-head(sort(t1, decreasing = TRUE),N)
	text(x=t2,y=jitter(rep(25,N),amount=25) ,labels=names(t2),srt=90, cex=0.7,col=1+(1:N),adj=c(0,0.5))
	
	
	
	usr1<-par("usr")
	x=usr1[1]
	y=mean(usr1[c(3,4)]);
	cex=0.6
	text(x=x, y=y, label=sprintf(" -Inf cnt: %i",cnt.minus.inf ) , adj=c(0,0.5) , col=color1(cnt.minus.inf), cex=cex )
	text(x=x, y=y, label=sprintf("\n +Inf cnt: %i",cnt.inf ) , adj=c(0,0.5) , col=color1(cnt.inf), cex=cex )
	text(x=x, y=y, label=sprintf("\n\n\n  NA  cnt: %i",cnt.na ) , adj=c(0,0.5) , col=color1(cnt.na), cex=cex )
	text(x=x, y=y, label=sprintf("\n\n\n\n  NaN cnt: %i",cnt.nan ) , adj=c(0,0.5) , col=color1(cnt.nan), cex=cex )
	
	dev.off()
}

rm(t1.col,t1, breaks, out.dir1)  



###############################################################################
###############################################################################
###############################################################################
###############################################################################
##  Plot correlation between replicates

out.dir1="out/cg03-correlationOfReplicatesAfterNormalization"
graphics.off()
Sys.sleep(1) # 1 sec pause
unlink(out.dir1, recursive = TRUE)
Sys.sleep(1) # 1 sec pause
dir.create(out.dir1, recursive = TRUE)

range1<-range(M1)

treatments <- unique(gsub("([^_]+)_\\d","\\1",names(M1)))
treatments0 <- setdiff(treatments,"NT")

for(tr in treatments){
	columns1<-sort(grep(paste0(tr,"_\\d"), names(M1), value=TRUE))
	
	if(length(columns1)>=2){
		for(i in 1:(length(columns1)-1)) {
			for(j in (i+1):length(columns1)) {
				
				filename<-paste0(columns1[i],"-",gsub("[^_]+_(\\d)","\\1",columns1[j]))
				pdf(file=paste0(out.dir1,"/",filename,".pdf"), width=4, height=4)
				
				x <- M1[,columns1[i]]
				y <- M1[,columns1[j]]
				
				plot(x,y,xlab=columns1[i], ylab=columns1[j], cex=0.2, asp=1, xlim=range1,ylim=range1)
				abline(a=0,b=1, lwd=3,col="gray")
				abline(v=myLog10LinearHybrid(c(0,1/nrow(M1)) ), col="red",lwd=1)
				abline(h=myLog10LinearHybrid(c(0,1/nrow(M1)) ), col="red",lwd=1)
				
				usr1<-par("usr")
				x1=usr1[1]
				y1=sum(usr1[c(3,4)] * c(0.1,0.9))
				text(x=x1, y=y1, label=sprintf("  Pearson rho: %0.3f", cor(x,y) ) , adj=c(0,0.5) , col="black", cex=1 )
				
				dev.off()				
			}
		}
	}
}

