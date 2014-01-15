# Take a list of XPT files of molecular measurements,
# produce a combined table of all measurements for each individual.
# 
# Can then be combined with cancer status.

library(foreign)

xpt_files<-list.files(pattern='*.XPT')
n_files<-length(xpt_files)
file_list<-list()
r<-NULL
vars<-NULL
# just get the range
for (xpt in xpt_files){
    file_list[[xpt]]<-read.xport(xpt)
    rownames(file_list[[xpt]])<-file_list[[xpt]]$"SEQN"
    r<-range(c(r,range(file_list[[xpt]]$"SEQN")))
    vars<-c(vars,names(file_list[[xpt]]))
}

n_indiv<-max(r)-min(r)+1
vars<-unique(vars)
print(r)
cat(vars,"\n",sep="\",\"")

# create a giant table
master_list<-matrix(rep(NA,length(vars)*n_indiv),nrow=n_indiv,ncol=length(vars))
colnames(master_list)<-vars
rownames(master_list)<-seq(min(r),max(r))

# now to populate it (this is horrible and slow)
for (indiv in seq(min(r),max(r))){
    print(indiv)
    for (xpt_file in file_list){
        if (indiv %in% xpt_file$"SEQN"){
            vars<-names(xpt_file)
            for (var in vars){
                if (!var=="SEQN"){
                    master_list[toString(indiv),var]<-xpt_file[toString(indiv),var]
                }
            }
        }
    }
}
master_list[,1]<-rownames(master_list)

# delete all NAs...
master_list<-master_list[!rowMeans(is.na(master_list))==1,]
master_list<-master_list[,!colMeans(is.na(master_list))==1]

# delete duplicates
del_list<-NULL
for(i in 3:ncol(master_list)){
    curr<-as.numeric(master_list[,i])
    prev<-as.numeric(master_list[,(i-1)])
    indx<-which(!is.na(curr))[1]
    ratio<-curr[indx]/prev[indx]
    # checking if proportional..
    diff<-curr-ratio*prev
    mean_diff<-mean(diff<0.001,na.rm=TRUE)
    if(!is.na(mean_diff)&(mean_diff==1)){
#    if(!is.na(mean(curr==ratio*prev,na.rm=TRUE))&(mean(curr==ratio*prev,na.rm=TRUE)==1)){
        cat("Duplicated information (",colnames(master_list)[i],") - deleting.\n")
        del_list<-c(del_list,i)
    }
}

master_list<-master_list[,-del_list]
write.table(master_list,file="enviro_data.txt",col.names=T,row.names=F,quote=FALSE)
