library(ggm)
library(igraph)

#THRESH<-0.051*5
THRESH<-as.numeric(readline("Threshold? "))

data_type<-readline("Data type? (lab/enviro/both) ")
if(data_type=="lab"){
    data<-read.table("lab/lab_data_high_qual.txt",header=T,row.names=1)
    # tidy up
    data<-data[rowMeans(is.na(data))<0.8,]
    data<-data[,colMeans(is.na(data))<0.8]
    vertex_cols<-"dodgerblue3"
} else if(data_type=="enviro"){
    data<-read.table("enviro/enviro_non_cat.txt",header=T,row.names=1)
    # tidy up
    data<-data[rowMeans(is.na(data))<0.8,]
    data<-data[,colMeans(is.na(data))<0.8]
    vertex_cols<-"green3"
} else if(data_type=="both"){
    lab_data<-read.table("lab/lab_data_high_qual.txt",header=T,row.names=1)
    enviro_data<-read.table("enviro/enviro_non_cat.txt",header=T,row.names=1)
    enviro_indiv<-rownames(enviro_data)
    # lab contains enviro
    lab_data<-lab_data[enviro_indiv,]
    data<-data.frame(lab_data,enviro_data)
    # tidy up
    data<-data[rowMeans(is.na(data))<0.8,]
    data<-data[,colMeans(is.na(data))<0.8]
    vertex_cols<-ifelse(colnames(data)%in%colnames(lab_data),"dodgerblue3","green3")
} else{
    cat("Data type must be one of 'lab', 'enviro', or 'both'!\n")
    quit("no")
}


# scale
z_standard<-function(x){
    mu<-mean(x,na.rm=TRUE)
    sigma<-sd(x,na.rm=TRUE)
    return((x-mu)/sigma)
}
data<-apply(data,2,z_standard)

# scaling mucks this up...
S<-cov(data,use="pairwise.complete.obs")
# this is exceedingly dodgy... NAs arise when there are no shared values between the covariates... assume the correlation is 0
S[is.na(S)]<-0
n<-ncol(data)

edge_mat<-matrix(1,n,n)
rownames(edge_mat)<-colnames(data)
colnames(edge_mat)<-colnames(data)
amat<-adjMatrix(edge_mat)
cg<-fitConGraph(amat, S, n)
#covg<-fitCovGraph(amat, S, n)

congraph<-(cg$Shat>THRESH)*1
adjgraph<-congraph-diag(1,n)
#covgraph<-(covg$Shat>THRESH)*1
ticks<-seq(nrow(edge_mat))/nrow(edge_mat)
labels<-rownames(edge_mat)

# this isn't as pretty
#pdf(paste(data_type,"_ggm_network.pdf",sep=""))
#par(las=2,cex=0.5,mar=c(5, 4, 4, 2) + 3.5)
#rownames(adjgraph)<-colnames(adjgraph)<-labels
#image(congraph,col=c("white","black"),axes=FALSE)
#axis(1,labels=labels,at=ticks)
#axis(2,labels=labels,at=ticks)
#dev.off()

g<-graph.adjacency(adjgraph,add.colnames="name",mode="undirected")
pdf(paste(data_type,"_ggm_network_",THRESH,".pdf",sep=""))
cat("Drawing network!\n")
plot(g,vertex.size=8,vertex.color="gray90",vertex.label.family="Helvetica",vertex.label.cex=0.5,vertex.label.color=vertex_cols,vertex.frame.color=NA,main=paste("GGM for ",data_type," measurements (threshold=",THRESH,")",sep=""))
dev.off()

#pdf('covgraph_network.pdf')
#par(las=2,cex=0.5,mar=c(5, 4, 4, 2) + 3.5)
#image(covgraph,col=c("white","black"),axes=FALSE)
#axis(1,labels=labels,at=ticks)
#axis(2,labels=labels,at=ticks)
#dev.off()
