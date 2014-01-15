library(ggm)
THRESH<-0.051

data<-read.table("lab/lab_data_high_qual.txt",header=T,row.names=1)
# tidy up
data<-data[rowMeans(is.na(data))<0.8,]
data<-data[,colMeans(is.na(data))<0.8]

# scale
z_standard<-function(x){
    mu<-mean(x,na.rm=TRUE)
    sigma<-sd(x,na.rm=TRUE)
    return((x-mu)/sigma)
}
data<-apply(data,2,z_standard)

# scaling mucks this up...
S<-cov(data,use="pairwise.complete.obs")
n<-ncol(data)

edge_mat<-matrix(1,n,n)
rownames(edge_mat)<-colnames(data)
colnames(edge_mat)<-colnames(data)
amat<-adjMatrix(edge_mat)
cg<-fitConGraph(amat, S, n)

congraph<-(cg$Shat>THRESH)*1

