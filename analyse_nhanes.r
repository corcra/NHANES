# Do hierarchical/k-means clustering, then test for enrichment.
# Also do a load of other stuff...

do_cluster<-FALSE
do_kmeans<-FALSE
do_heatmap<-FALSE
do_cluster_enrich<-FALSE
only_cancer<-FALSE
do_logistic<-FALSE
do_supervised<-TRUE

library(cluster)
library(gplots)

data_type<-readline("Data type? (lab/enviro) ")
if(data_type=="lab"){
    data<-read.table("lab/lab_data_high_qual.txt",header=T,row.names=1)
    # the high qual data is all the continuous measurements (hopefully)
    non_categorical<-colnames(data)
} else if(data_type=="enviro"){
#    data<-read.table("enviro/enviro_data.txt",header=T,row.names=1)
    data<-read.table("enviro/enviro_non_cat.txt",header=T,row.names=1)
    # for enviro data, assume categorical unless in this list...
#    non_categorical<-c("ALQ120Q","ALQ130U","ALQ130","ALQ141Q","ALQ155","BPD035","BPD058","CKQ070Q","DEQ038G","DEQ038Q","DED120","DED125","DID040","DID060","DID250","DID260","DIQ280","DIQ300S","DIQ300D","DID310S","DID310D","DID320","DID330","DID341","DID350","DUQ210","DUQ213","DUQ215Q","DUQ220Q","DUQ230","DUQ260","DUQ270Q","DUQ280","DUQ300","DUQ310Q","DUQ320","DUQ340","DUQ350Q","DUQ360","DUQ390","DUQ400Q","ECD010","ECD070A","ECD070B","HUD080","HOD050","IMQ090","PAQ706","PAQ610","PAD615","PAQ625","PAD630","PAQ640","PAD645","PAQ655","PAD660","PAQ670","PAD675","PAD680","RXQ525Q","RXD530","RDD040","RDD060","RDQ080","RDD120","SXD031","SXD171","SXD510","SXQ824","SXQ827","SXD633","SXQ636","SXQ639","SXD642","SXQ410","SXQ550","SXQ836","SXQ841","SXD621","SXQ624","SXQ627","SXD630","SXQ590","SXQ600","SXD101","SXD450","SXQ724","SXQ727","SXQ130","SXQ490","SLD010H","SMD030","SMQ050Q","SMD055","SMD057","SMD641","SMD650","SMD100TR","SMD100NI","SMD100CO","SMD630","SMD430","SMQ710","SMQ720","SMQ740","SMQ750","SMQ770","SMQ780","SMQ800","SMQ817","SMQ830")
    # in no particular order...
#    cat<-c("SMD100BR","SMDUPCA")
#    for now, using 'high qual' data
    non_categorical<-colnames(data)
} else {
    cat("Data type must be one of 'lab' or 'enviro'!\n")
    quit("no")
}

cat("Getting disease and demographic information!\n")
# get disease state
source("disease_state.r")
# get demographic info
demo_full<-read.xport("demo/DEMO_G.XPT")
# this is just gender, age, race
demographics<-demo_full[,c("RIAGENDR","RIDAGEYR","RIDRETH3")]
rownames(demographics)<-demo_full[,"SEQN"]

# restrict to cancer-diagnosed individuals
if(only_cancer){
    cat("Restricting to individuals with history of cancer.\n")
    # note this produces a lot of NAs... they are removed in the pruning, but still
    data<-data[disease_state[rownames(data),"cancer"]==1,]
    disease_state<-cancer_state
    }

cat("Cleaning up data!\n")
# remove rows which are all NAs, or close to it...
data<-data[rowMeans(is.na(data))<0.8,]
# now remove cols...
data<-data[,colMeans(is.na(data))<0.8]

# now scale data (also turn to categorical...)
# seems to be not a great idea right now
z_standard<-function(x){
    col<-data[,x]
    if(x %in% non_categorical){
        # even still, this isn't exactly perfect...
        mu<-mean(col,na.rm=TRUE)
        sigma<-sd(col,na.rm=TRUE)
        return((col-mu)/sigma)
    } else{
        return(as.factor(col))
    }
}
cat("Normalising data!\n")
data_new<-lapply(names(data),z_standard)
data_temp<-data.frame(data_new)
colnames(data_temp)<-colnames(data)
rownames(data_temp)<-rownames(data)
data<-data_temp
#data<-apply(data,2,z_standard)
#data2<-data.matrix(scale(data))

# restrict to my list of individuals
my_indiv<-rownames(data)
disease_state<-disease_state[my_indiv,]
bg_success<-colSums(disease_state==1,na.rm=TRUE)
if(only_cancer){
    bg_failure<-nrow(disease_state)-bg_success
} else{
    bg_failure<-colSums(disease_state==2,na.rm=TRUE)
}
demographics<-data.matrix(demographics[my_indiv,])

# get cancer status
if(!only_cancer){
    cat("Getting cancer status!\n")
    cancer_status<-1*(disease_state[my_indiv,"cancer"]==1)
    cancer_status[is.na(cancer_status)]<-0
} else{
    cancer_status<-rep(1,nrow(data))
}

if(do_cluster){
    cat("Doing hierarchical clustering!\n")
    n_clust<-as.integer(readline("number of clusters: "))
    d<-dist(data)
    h<-hclust(d,method="complete")
    h_clusters<-cutree(h,k=n_clust)
}

if(do_kmeans){
    cat("Doing k-means clustering!\n")
    n_clust<-as.integer(readline("number of clusters: "))
    k_means<-pam(data,k=n_clust,keep.diss=TRUE)
    k_clusters<-k_means$clustering
#    plot(k_means)
}

if(do_heatmap){
    cat("Making a heatmap!\n")
    heatmap.2(data,RowSideColors=ifelse(cancer_status==1,"seagreen3","grey37"),dendrogram="row",na.color="white",trace="none",Colv=FALSE,key=FALSE,ylab="Individuals",xlab="Molecular measurements",labRow="",cexCol=0.3,lhei=c(1,9))
}

if(do_cluster_enrich){
    cat("Searching for enrichments...\n")
    min_pvals<-NULL
    groups<-vector("list")
    for (i in 1:n_clust){
        this_group<-list()
        # h clusters
    #    g<-subset(data,h_clusters==i)
        # k-means clusters
        g<-subset(data,k_clusters==i)
        this_group[["number"]]<-i
        this_group[["data"]]<-g
        this_group[["size"]]<-nrow(g)
        this_group[["freq"]]<-colMeans(disease_state[rownames(g),]==1,na.rm=TRUE)
        this_group[["success"]]<-colSums(disease_state[rownames(g),]==1,na.rm=TRUE)
        this_group[["failure"]]<-colSums(disease_state[rownames(g),]==2,na.rm=TRUE)
        this_group[["disease_size"]]<-this_group[["success"]]+this_group[["failure"]]

        # check enrichment using hypergeometric
        enrich_p_vals<-phyper(this_group[["success"]],bg_success,bg_failure,this_group[["disease_size"]])

        # do logistic regression on the disease
        browser()
        #
        # coming soon
        this_group[["enrich_p_vals"]]<-enrich_p_vals

        groups[[i]]<-this_group
    } 

    # Get most enriched disease/cancer type for each group
    which_disease<-function(x){
        group_number<-x$"number"
        min_pval<-min(x$"enrich_p_vals")
        names<-colnames(disease_state)[which(x$"enrich_p_vals"==min_pval)]
        if(length(names)>1){
            cat("Multiple diseases:",paste(names,collapse=', '),"(choosing first)\n")
        }
        name<-names[1]
        group_success<-x$"success"[[name]]
        group_total<-x$"disease_size"[[name]]
        background_success<-bg_success[[name]]
        background_total<-bg_success[[name]]+bg_failure[[name]]
        cat(group_number,group_total,group_success,min_pval,background_total,background_success,"",name,"\n",sep="\t")
    }
    cat("Most enriched (p-val) disease for each cluster:\n")
    cat("cluster","size","cases","p-val","","total","all cases","disease","\n",sep="\t")
    lapply(groups,which_disease)
}

get_logistic<-function(x,cancer_status,demo){
    cat("get_logistic has been called\n")
    if((mean(is.na(x))>0.9|length(levels(x))==1)){
        return(1)
    }
    if(demo==1){
        lfit<-glm(cancer_status ~ demographics + x, family=binomial)
        index<-3
    } else{
        lfit<-glm(cancer_status ~ x, family=binomial)
        index<-2
    }
    p_val<-anova(lfit,test="Chisq")$"Pr(>Chi)"[index]
    #cat(mean(is.na(x)),p_val,"\n")
    return(p_val)
}

logit_p_vals<-NULL
logit_p_vals_nodemo<-NULL
if(do_logistic){
    cat("Checking each measurement for cancer association (logistic)!\n")
    for (i in 1:ncol(data)){
        if(i%%5==0){
            cat(100*i/ncol(data),"%\n",sep="")
        }
        col<-data[,i]
        logit_p_vals<-c(logit_p_vals,get_logistic(col,cancer_status,1))
        logit_p_vals_nodemo<-c(logit_p_vals,get_logistic(col,cancer_status,0))
        }
    names(logit_p_vals)<-colnames(data)
    names(logit_p_vals_nodemo)<-colnames(data)

    bonf<-0.05/ncol(data)
    which_sig<-which(logit_p_vals<=bonf)
    sig<-logit_p_vals[which_sig]
    #cat("Significant hits:\n")
    #print(sig)

    which_sig_nodemo<-which(logit_p_vals_nodemo<=bonf)
    sig_nodemo<-logit_p_vals_nodemo[which_sig_nodemo]
}

# --- supervised analysis --- #
# right now, this doesn't seem to work for the enviro option, and I'm not sure why... probably too many NAs?
if(do_supervised){
    cat("Performing supervised analysis!\n")
    # first, take a training set...
    train_frac<-0.5
    cat("Fraction of dataset used for training:",train_frac,"\n")
    training<-sample(seq(nrow(data)),floor(train_frac*nrow(data)))
    # not sure how well this takes categorical data
    training_data<-as.matrix(data[training,])
    training_cancer<-1*(disease_state[rownames(training_data),"cancer"]==1)
    training_cancer[is.na(training_cancer)]<-0
    training_demo<-demographics[rownames(training_data),]
    testing_data<-as.matrix(data[-training,])
    testing_cancer<-1*(disease_state[rownames(testing_data),"cancer"]==1)
    testing_cancer[is.na(testing_cancer)]<-0
    testing_demo<-demographics[rownames(testing_data),]
    # now fit a model based on the training data
   
    fit_data<-data.frame(Cancer=training_cancer,Demo=training_demo,Data=training_data)
    browser()
    fit_model<-glm(Cancer~.,family=binomial,data=fit_data)

    vali_data<-data.frame(Demo=testing_demo,Data=testing_data)
    test_predictions<-predict(fit_model,newdata=vali_data,type="response",se.fit=TRUE,na.action="na.omit")
    pred_vals<-1*(test_predictions$fit>=0.5)
    true_vals<-1*(disease_state[names(test_predictions$fit),"cancer"]==1)
    
    tp<-sum(pred_vals==1&true_vals==1,na.rm=TRUE)
    fp<-sum(pred_vals==1&true_vals==0,na.rm=TRUE)
    tn<-sum(pred_vals==0&true_vals==0,na.rm=TRUE)
    fn<-sum(pred_vals==0&true_vals==1,na.rm=TRUE)

    sens<-tp/(tp+fn)*100
    spec<-tn/(fp+tn)*100
    cat("Senstivity: ",sens,"%\nSpecificity: ",spec,"%\n",sep="")
}
