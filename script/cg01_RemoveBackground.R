###############################
## We estimate the amount of background noise coming from genomic DNA ,
# and then we will remove it.

rm(list=ls())
source("script/functions-for-chemogenomic-project.R");

#####################################
# loading read numbers
M1=read.csv(file="data/dataTable1-read_numbers_of_sequencing.csv",sep=",")
rownames(M1)<-M1$Feature_ID

# batch structure describes witch data coluns are coming from the same assay
batch.df <- read.csv(file="data/batch-structure.csv",sep=",")
#########################


# names of columns of treatments
t.column.idx<-names(M1)
t.column.idx<-t.column.idx[t.column.idx!="Feature_ID" &t.column.idx!=  "b_numbers"]


####################################################################################################
####################################################################################################
####################################################################################################
#################################  Generating plots  ############################################### 

graphics.off()
out.dir1="out/cg01-RemoveBackground"
unlink(out.dir1, recursive = TRUE)
dir.create(out.dir1, recursive = TRUE)


cs1 <- colSums(M1[,t.column.idx]) # sum over columns
lacZ.values <- unlist(M1[rownames(M1)=="lacZ",t.column.idx])
cs1<- cs1-lacZ.values

lacI.values <- unlist(M1[rownames(M1)=="lacI",t.column.idx])
cs2 <- cs1-lacI.values
mx<-rbind(lacI.values,cs2)
#mx<-rbind(lacI.values,lacZ.values,cs2)

batch.num<-sapply(colnames(mx),FUN=function(x) batch.df$batch.num[batch.df$cname==x])
stopifnot(length(batch.num)==ncol(mx))
walls1<-which( batch.num[-1]!=batch.num[-length(batch.num)]) # position of separators on the figures

png(filename=paste(out.dir1,"/barplot-absoluteRedNumbers.png", sep=""), width=1920, height=1080, pointsize=24)
par( mar= c(7.5, 4.1, 4.1,2.1))
mp<-barplot(mx, las=2 , main="lacI vs. other ORFs", cex.names=0.6)
abline(v=(mp[walls1]+ mp[1+walls1])/2,col="blue", lwd=3)
text.positions<-c(1,1+walls1);
mtext(text=sprintf("Batch: %i",batch.num[text.positions]), side=3, at=mp[text.positions], adj=c(0,0), cex=0.7, col="blue")

dev.off()




png(filename=paste(out.dir1,"/barplot-relativeRedNumbers.png", sep=""), width=1920, height=1080, pointsize=24)
par( mar= c(7.5, 4.1, 4.1,2.1))
mp<-barplot(t(t(mx)/cs1), las=2 , main="lacI vs. other ORFs", axes=FALSE,cex.names=0.6)
abline(v=(mp[walls1]+ mp[1+walls1])/2,col="blue", lwd=3)
text.positions<-c(1,1+walls1);
mtext(text=sprintf("Batch: %i",batch.num[text.positions]), side=3, at=mp[text.positions], adj=c(0,0), cex=0.7, col="blue")
ticks1<-c(10,20,30,40,45,50,60,70,80,90,100)
axis(side=2,at=ticks1/100, labels=sprintf("%i%%",ticks1),las=2, cex.axis=0.7)
abline(h=0.5, col="blue", lwd=1)
abline(h=0.45, col="blue", lwd=1)
abline(h=0.4, col="blue", lwd=1)
dev.off()

rm(cs1,cs2,lacI.values,lacZ.values,mx, walls1)


####################################################################################################
####################################################################################################
####################################################################################################


##############################################################
# compose an aggregated  background column


stopifnot(1==sum(M1$Feature_ID=="lacI")) # check if we have lacI and lacZ data rows
stopifnot(1==sum(M1$Feature_ID=="lacZ"))


bg.df<-M1[,grep("^background_\\d",names(M1), value=TRUE)]


# I assign weights to the background columns. Each batch has the same sum of weight
tmp.df1<-data.frame(cname=names(bg.df), order1=1:ncol(bg.df))
tmp.df1<-merge(batch.df, tmp.df1, by="cname")
tmp.df2<-aggregate(tmp.df1$batch.num,by=list(batch.num=tmp.df1$batch.num),FUN=function(x) 1/length(x) )
names(tmp.df2)[2]<-"weight"
tmp.df1<-merge(tmp.df2, tmp.df1, by="batch.num")
tmp.df1<-tmp.df1[order(tmp.df1$order1), ]
tmp.df1$weight<-tmp.df1$weight/sum(tmp.df1$weight)

weight1<-tmp.df1$weight


bg.aggregated<-apply(bg.df,1,	function(x) sum(weight1*x))


# a lacI and lacZ will have the median of the others
idx<-names(bg.aggregated)=="lacI" | names(bg.aggregated)=="lacZ"
bg.aggregated[idx]<-median(bg.aggregated[!idx])


bg.df$bg.aggregated<-bg.aggregated

M1$aggregated_background<-bg.aggregated;

rm(weight1,tmp.df1, tmp.df2,  bg.df,bg.aggregated, idx)


########################################################################################


sum.bg<-sum(M1$aggregated_background)

sum1<-apply(M1[,t.column.idx],2,sum)
value.lacI<-unlist(M1[M1$Feature_ID=="lacI",t.column.idx])
value.lacZ<-unlist(M1[M1$Feature_ID=="lacZ",t.column.idx])

sum2<- sum1-value.lacZ
sum3<- sum2-value.lacI

stopifnot(all(names(value.lacI)==names(sum1)))

lacI.in.background <- M1["lacI","aggregated_background"]

#
# In the ideal case the amount of lacI would be equal with the sum of others.
# The alpha value shows the surplus. It is needed to remove for balance.
#
# I am looking for alpha with property:
# sum3-alpha*sum.bg=value.lacI-alpha*lacI.in.background
#
# -alpha*sum.bg+alpha*lacI.in.background=value.lacI-sum3
#  alpha*(-sum.bg+lacI.in.background)=value.lacI-sum3
alpha <- (value.lacI-sum3)/(lacI.in.background-sum.bg)

rm(value.lacI, value.lacZ, sum1, sum2, sum3,sum.bg)

alpha<-alpha[!grepl("^background_\\d$",names(alpha))]

png(filename=paste0(out.dir1,"/alpha-barplot.png"), width=1920, height=1080, pointsize=24)
batch.num<-sapply(names(alpha),FUN=function(x) batch.df$batch.num[batch.df$cname==x])
stopifnot(length(batch.num)==ncol(alpha))
walls1<-which( batch.num[-1]!=batch.num[-length(batch.num)])
mp<-barplot(alpha, las=2, ylim=c(-1,2.5), cex.names=0.6, ylab="alpha", main=c("Background Correction Alpha", " y = x-alpha*bg" ))
grid()
abline(v=(mp[walls1]+ mp[1+walls1])/2,col="blue", lwd=3)
text.positions<-c(1,1+walls1);
mtext(text=sprintf("Batch: %i",batch.num[text.positions]), side=3, at=mp[text.positions], adj=c(0,0), cex=0.7, col="blue")
dev.off()



alpha[alpha<0]<-0

# M2 will be the matrix of corrected data
M2<-data.frame(row.names=rownames(M1))
for( i in names(alpha)){
	alpha1<-alpha[i]
	M2[,i]<-M1[,i]-alpha1*M1$aggregated_background
}
rm(alpha1, lacI.in.background,i, walls1,ticks1)
rm(M1)

############################################################################################
# I am generating figures of the after correction status

graphics.off()

cs1 <- colSums(M2)
lacZ.values <- unlist(M2[rownames(M2)=="lacZ",])
cs1<- cs1-lacZ.values

lacI.values <- unlist(M2[rownames(M2)=="lacI",])
cs2 <- cs1-lacI.values
mx<-rbind(lacI.values,cs2)



batch.num<-sapply(colnames(mx),FUN=function(x) batch.df$batch.num[batch.df$cname==x])
stopifnot(length(batch.num)==ncol(mx))
walls1<-which( batch.num[-1]!=batch.num[-length(batch.num)])


png(filename=paste(out.dir1,"/barplot-absoluteRedNumbers-afterLacICorrection.png", sep=""), width=1920, height=1080, pointsize=24)
par( mar= c(7.5, 4.1, 4.1,2.1))
mp<-barplot(mx, las=2 , main="lacI vs. other ORFs - After lacI Correction" , cex.names=0.6)
abline(v=(mp[walls1]+ mp[1+walls1])/2,col="blue", lwd=3)
text.positions<-c(1,1+walls1);
mtext(text=sprintf("Batch: %i",batch.num[text.positions]), side=3, at=mp[text.positions], adj=c(0,0), cex=0.7, col="blue")
dev.off()



png(filename=paste(out.dir1,"/barplot-relativeRedNumbers-afterLacICorrection.png", sep=""), width=1920, height=1080, pointsize=24)
par( mar= c(7.5, 4.1, 4.1,2.1))
mp<-barplot(t(t(mx)/cs1), las=2 , main="lacI vs. other ORFs-After lacI Correction" , axes=FALSE,cex.names=0.7 )
ticks1<-c(10,20,30,40,50,60,70,80,90,100)
axis(side=2,at=ticks1/100, labels=sprintf("%i%%",ticks1),las=2, cex.axis=0.7)

abline(v=(mp[walls1]+ mp[1+walls1])/2,col="blue", lwd=3)
text.positions<-c(1,1+walls1);
mtext(text=sprintf("Batch: %i",batch.num[text.positions]), side=3, at=mp[text.positions], adj=c(0,0), cex=0.7, col="blue")
dev.off()

rm(cs1,cs2,lacI.values,lacZ.values,mx,out.dir1, walls1)    

############################################################################################


write.table(M2, file="data/dataTable2-afterBackgroundRemove.csv",sep=",",quote = FALSE,row.names=TRUE, col.names=NA)
