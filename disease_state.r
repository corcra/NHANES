# Obtain disease state information...
cat("(disease_state.r) Getting disease information!\n")

library(foreign)

meds<-read.xport("disease/MCQ_G.XPT")
individuals<-meds$"SEQN"

# note: for all of these, 1 means yes and 2 means no
# these are also 'ever been told you had...'

asthma<-meds$"MCQ010"
psoriasis<-meds$"MCQ070"
overweight<-meds$"MCQ080"
coeliac<-meds$"MCQ082"
arthritis<-meds$"MCQ180A"
gout<-meds$"MCQ160N"
congestiveheart<-meds$"MCQ160B"
coronaryheart<-meds$"MCQ180C"
angina<-meds$"MCQ160D"
heartattack<-meds$"MCQ160E"
stroke<-meds$"MCQ160F"
emphysema<-meds$"MCQ160G"
thyroid<-meds$"MCQ160M"
bronchitis<-meds$"MCQ160K"
liver<-meds$"MCQ160L"
cancer<-meds$"MCQ220"

disease_state<-data.frame(asthma,psoriasis,overweight,coeliac,arthritis,gout,congestiveheart,coronaryheart,heartattack,stroke,emphysema,thyroid,bronchitis,liver,cancer)
rownames(disease_state)<-individuals

# types of cancers
bladder<-1*(meds$"MCQ230A"==10|meds$"MCQ230B"==10|meds$"MCQ230C"==10|meds$"MCQ230D"==10)
blood<-1*(meds$"MCQ230A"==11|meds$"MCQ230B"==11|meds$"MCQ230C"==11|meds$"MCQ230D"==11)
bone<-1*(meds$"MCQ230A"==12|meds$"MCQ230B"==12|meds$"MCQ230C"==12|meds$"MCQ230D"==12)
brain<-1*(meds$"MCQ230A"==13|meds$"MCQ230B"==13|meds$"MCQ230C"==13|meds$"MCQ230D"==13)
breast<-1*(meds$"MCQ230A"==14|meds$"MCQ230B"==14|meds$"MCQ230C"==14|meds$"MCQ230D"==14)
cervix<-1*(meds$"MCQ230A"==15|meds$"MCQ230B"==15|meds$"MCQ230C"==15|meds$"MCQ230D"==15)
colon<-1*(meds$"MCQ230A"==16|meds$"MCQ230B"==16|meds$"MCQ230C"==16|meds$"MCQ230D"==16)
esophagus<-1*(meds$"MCQ230A"==17|meds$"MCQ230B"==17|meds$"MCQ230C"==17|meds$"MCQ230D"==17)
gallbladder<-1*(meds$"MCQ230A"==18|meds$"MCQ230B"==18|meds$"MCQ230C"==18|meds$"MCQ230D"==18)
kidney<-1*(meds$"MCQ230A"==19|meds$"MCQ230B"==19|meds$"MCQ230C"==19|meds$"MCQ230D"==19)
larynx<-1*(meds$"MCQ230A"==20|meds$"MCQ230B"==20|meds$"MCQ230C"==20|meds$"MCQ230D"==20)
leukemia<-1*(meds$"MCQ230A"==21|meds$"MCQ230B"==21|meds$"MCQ230C"==21|meds$"MCQ230D"==21)
liver<-1*(meds$"MCQ230A"==22|meds$"MCQ230B"==22|meds$"MCQ230C"==22|meds$"MCQ230D"==22)
lung<-1*(meds$"MCQ230A"==23|meds$"MCQ230B"==23|meds$"MCQ230C"==23|meds$"MCQ230D"==23)
lymphoma<-1*(meds$"MCQ230A"==24|meds$"MCQ230B"==24|meds$"MCQ230C"==24|meds$"MCQ230D"==24)
melanoma<-1*(meds$"MCQ230A"==25|meds$"MCQ230B"==25|meds$"MCQ230C"==25|meds$"MCQ230D"==25)
mouth<-1*(meds$"MCQ230A"==26|meds$"MCQ230B"==26|meds$"MCQ230C"==26|meds$"MCQ230D"==26)
nervous<-1*(meds$"MCQ230A"==27|meds$"MCQ230B"==27|meds$"MCQ230C"==27|meds$"MCQ230D"==27)
ovary<-1*(meds$"MCQ230A"==28|meds$"MCQ230B"==28|meds$"MCQ230C"==28|meds$"MCQ230D"==28)
pancreas<-1*(meds$"MCQ230A"==29|meds$"MCQ230B"==29|meds$"MCQ230C"==29|meds$"MCQ230D"==29)
prostate<-1*(meds$"MCQ230A"==30|meds$"MCQ230B"==30|meds$"MCQ230C"==30|meds$"MCQ230D"==30)
rectum<-1*(meds$"MCQ230A"==31|meds$"MCQ230B"==31|meds$"MCQ230C"==31|meds$"MCQ230D"==31)
skin_nonmelanoma<-1*(meds$"MCQ230A"==32|meds$"MCQ230B"==32|meds$"MCQ230C"==32|meds$"MCQ230D"==32)
skin_unknown<-1*(meds$"MCQ230A"==33|meds$"MCQ230B"==33|meds$"MCQ230C"==33|meds$"MCQ230D"==33)
soft_tissue<-1*(meds$"MCQ230A"==34|meds$"MCQ230B"==34|meds$"MCQ230C"==34|meds$"MCQ230D"==34)
stomach<-1*(meds$"MCQ230A"==35|meds$"MCQ230B"==35|meds$"MCQ230C"==35|meds$"MCQ230D"==35)
testis<-1*(meds$"MCQ230A"==36|meds$"MCQ230B"==36|meds$"MCQ230C"==36|meds$"MCQ230D"==36)
thyroid<-1*(meds$"MCQ230A"==37|meds$"MCQ230B"==37|meds$"MCQ230C"==37|meds$"MCQ230D"==37)
uterus<-1*(meds$"MCQ230A"==38|meds$"MCQ230B"==38|meds$"MCQ230C"==38|meds$"MCQ230D"==38)
other<-1*(meds$"MCQ230A"==39|meds$"MCQ230B"==39|meds$"MCQ230C"==39|meds$"MCQ230D"==39)
unknown<-1*(meds$"MCQ230A"==99|meds$"MCQ230B"==99|meds$"MCQ230C"==99|meds$"MCQ230D"==99)

cancer_state<-data.frame(bladder,blood,bone,brain,breast,cervix,colon,esophagus,gallbladder,kidney,larynx,leukemia,liver,lung,lymphoma,melanoma,mouth,nervous,ovary,pancreas,prostate,rectum,skin_nonmelanoma,skin_unknown,soft_tissue,stomach,testis,thyroid,uterus,other,unknown)
rownames(cancer_state)<-individuals
# saying 'not this one' otherwise...
cancer_state[is.na(cancer_state)]<-2
