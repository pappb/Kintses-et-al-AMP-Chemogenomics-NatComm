rm(list=ls())
M1=read.csv(file="data/dataTable2-afterBackgroundRemove.csv",sep=",",row.names=1)



# remove rows of lacI and LacZ  (These are parts of induction casette.  Each plasmid contains them.)
M2<-M1[rownames(M1)!="lacI" & rownames(M1)!="lacZ",] 

# scale each column to sum=1

colsums1<-colSums(M2)
for( i in names(M2)){
	M2[,i]=M2[,i]/colsums1[i] 
}      


write.table(M2, file="data/dataTable3-afterRescaling.csv",sep=",",quote = FALSE,row.names=TRUE, col.names=NA)