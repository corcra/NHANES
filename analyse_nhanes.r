# Do hierarchical/k-means clustering, then test for enrichment.
# Also do a load of other stuff...

do_cluster<-FALSE
do_kmeans<-FALSE
do_heatmap<-FALSE
do_cluster_enrich<-FALSE
only_cancer<-FALSE
do_logistic<-TRUE
do_supervised<-FALSE

library(cluster)
library(gplots)

data_type<-readline("Data type? (lab/enviro) ")
if(data_type=="lab"){
    data<-read.table("lab/lab_data.txt",header=T,row.names=1)
    continuous<-read.table("lab/lab_cont.txt",as.is=TRUE)[,1]
} else if(data_type=="enviro"){
    data<-read.table("enviro/enviro_data.txt",header=T,row.names=1)
    continuous<-read.table("enviro/enviro_cont.txt",as.is=TRUE)[,1]
} else {
    cat("Data type must be one of 'lab' or 'enviro'!\n")
    quit("no")
}
# restricting to continuous data, for now...
keep_cont<-intersect(colnames(data),continuous)
data<-subset(data,select=keep_cont)
# this is kinda redundant until i deal with categorical data properly
non_categorical<-colnames(data)

cat("Getting disease and demographic information!\n")
# get disease state
source("disease_state.r")
# get demographic info
demo_full<-read.xport("demo/DEMO_G.XPT")
# this is just gender, age, race
demographics<-demo_full[,c("RIAGENDR","RIDAGEYR","RIDRETH3")]
rownames(demographics)<-demo_full[,"SEQN"]
d<-apply(demographics,2,as.factor)
demographics<-data.frame(d[,1],demographics[,2],d[,3])
colnames(demographics)<-c("sex","age","ethnicity")

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
        # first: get rid of refused/missing... this will hit some true values most likely, but it shouldn't be the worst...
        missing<-NULL
        # what we're looking for is the max value being 9 or 99 or 999 or 9999 or 99999 (hopefully nothing higher than that)
        max_val<-max(col,na.rm=TRUE)
        # weird weird formatting here
        if(max_val==(9)|(max_val==7)){
            missing<-which(col==9|col==7)
        } else if(max_val==(99)|(max_val==77)){
            missing<-which(col==99|col==77)
        } else if(max_val==(999)|(max_val==777)){
            missing<-which(col==999|col==777)
        } else if(max_val==(9999)|(max_val==7777)){
            missing<-which(col==9999|col==7777)
        } else if(max_val==(99999)|(max_val==77777)){
            missing<-which(col==99999|col==77777)
        }

        col[missing]<-NA

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
rm(data_temp)
rm(data_new)

# restrict to my list of individuals
my_indiv<-rownames(data)
disease_state<-disease_state[my_indiv,]
bg_success<-colSums(disease_state==1,na.rm=TRUE)
if(only_cancer){
    bg_failure<-nrow(disease_state)-bg_success
} else{
    bg_failure<-colSums(disease_state==2,na.rm=TRUE)
}
demographics<-demographics[my_indiv,]

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
    if((mean(is.na(x))>0.9|length(levels(x))==1)){
        return(1)
    }
    if(demo==1){
        lfit<-glm(cancer_status ~ sex + age + ethnicity + x, family=binomial,data=demographics)
        sex_assoc<-coefficients(summary(lfit))["sex2",4]
        age_assoc<-coefficients(summary(lfit))["age",4]
    } else{
        lfit<-glm(cancer_status ~ x, family=binomial)
        sex_assoc<-1
        age_assoc<-1
    }
#   the anova way is testing if including the covariate is significant...
#    p_val<-anova(lfit,test="Chisq")$"Pr(>Chi)"[index]
    p_val<-coefficients(summary(lfit))["x",4]
    return(list("p_var"=p_val,"p_sex"=sex_assoc,"p_age"=age_assoc))
}

var_p_vals<-1
var_p_vals_nodemo<-1
sex_p_vals<-1
age_p_vals<-1

if(do_logistic){
    cat("Checking each measurement for cancer association (logistic)!\n")
    for (i in 2:ncol(data)){
        if(i%%5==0){
            cat(100*i/ncol(data),"%\n",sep="")
        }
        col<-data[,i]
        p_vals<-get_logistic(col,cancer_status,1)
        var_p_vals<-c(var_p_vals,p_vals$"p_var")
        sex_p_vals<-c(sex_p_vals,p_vals$"p_sex")
        age_p_vals<-c(age_p_vals,p_vals$"p_age")
        p_vals_nodemo<-get_logistic(col,cancer_status,0)
        var_p_vals_nodemo<-c(var_p_vals_nodemo,p_vals_nodemo$"p_var")
    }
    names(var_p_vals)<-colnames(data)
    names(sex_p_vals)<-colnames(data)
    names(age_p_vals)<-colnames(data)
    names(var_p_vals_nodemo)<-colnames(data)

    bonf<-0.05/ncol(data)
    which_sig<-which(var_p_vals<=bonf)
    sig<-var_p_vals[which_sig]
    cat("Significant hits:\n")
    print(sig)

    which_sig_nodemo<-which(var_p_vals_nodemo<=bonf)
    sig_nodemo<-var_p_vals_nodemo[which_sig_nodemo]
}

shist<-function(x,name){
    hist(x,col="deepskyblue2",breaks=50,main=name,xlab="p value")
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
