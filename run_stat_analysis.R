usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
packages= c("numDeriv","RCurl","ggplot2","reshape2","plyr","randomForest","lme4","nlme","car","lmtest","coefplot2","GGally","scales",'Rmisc',"tidyr","effsize","boot",'grid')
lapply(packages,usePackage)

checkConv = function(mc) {
  #check gradient, has converged if tol<0.001
  derivs1 <- mc@optinfo$derivs
  sc_grad1 <- with(derivs1,solve(Hessian,gradient))
  g1 = max(abs(sc_grad1))
  
  #parallel 
  g2 = max(pmin(abs(sc_grad1),abs(derivs1$gradient)))
  #not singular
  m1_sc = mc
  tt <- getME(m1_sc,"theta")
  ll <- getME(m1_sc,"lower")
  g3 = min(tt[ll==0])
  c(g1<0.001,g2<0.001,g3>10e-6)
}
summarySE = function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

################ READ DATA ################
##READ IN ALL .CSV IN FOLDER, COMBINE UNDER ONE HEADING and common columns
csvdir = '../Patient_Output_ja_spatial_Szonset/Brodman/L_R/New/'
temp = list.files(path=csvdir,pattern="*.csv")
tables = lapply(paste(csvdir,temp,sep=''),read.csv,header=TRUE,as.is=TRUE)
dat = Reduce(function(x, y) merge(x, y, all=TRUE),tables)
dat[dat$sex==FALSE,]$sex = "F"

str(dat)
dat$patient = factor(dat$patient)
dat$sex = as.factor(dat$sex)
dat$session = factor(dat$session)
dat$task = factor(dat$task)
dat$trial = factor(dat$trial)
dat$word = factor(dat$word)
dat$recalled = factor(dat$recalled)

str(dat)
names(dat)
summary(dat)

# only use the encoding period
recallDat= dat[dat$task=="WORD",]
recallDat = recallDat[recallDat$recalled!=2,]
recallDat$recalled = factor(recallDat$recalled)

#Get list of patients
uniquePt = unique(recallDat$patient) 

##### ADDING CATEGORIES #####


## LATERALITY
leftSubj = NULL
k = 0
for (i in c(4,11,13,16)) {
  k = k + 1
  leftSubj[k] = sprintf('I022_P0%.2d',i)
}
for (i in c(3,6,7,9,10,11,12,14,21,22,26,32,35,40,43,44,45,51,54)) {
  k = k + 1
  leftSubj[k] = sprintf('I027_P0%.2d',i)
}

rightSubj = NULL
k = 0
for (i in c(1,5,7,14,15,17)) {
  k = k + 1
  rightSubj[k] = sprintf('I022_P0%.2d',i)
}
for (i in c(2,4,13,15,18,19,23,24,25,30,34,36,37,39,41,46,47,48,49,53)) {
  k = k + 1
  rightSubj[k] = sprintf('I027_P0%.2d',i)
}

bilateralSubj = NULL
k = 0
for (i in c(2,6,8,12)) {
  k = k + 1
  bilateralSubj[k] = sprintf('I022_P0%.2d',i)
}
for (i in c(8,17,28,29,33,42)) {
  k = k + 1
  bilateralSubj[k] = sprintf('I027_P0%.2d',i)
}

recallDat$onset = 'UNLOCALIZED'
recallDat[recallDat$patient %in% leftSubj,]$onset = 'LEFT'
recallDat[recallDat$patient %in% rightSubj,]$onset = 'RIGHT'
recallDat[(recallDat$patient %in% bilateralSubj),]$onset = 'BILATERAL'

onlysz = recallDat[recallDat$onset !='UNLOCALIZED',]
recallDat$onset = as.factor(recallDat$onset)
onlysz$onset = as.factor(onlysz$onset)

onlysz$LeftOnset = 0
onlysz$RightOnset = 0
onlysz[onlysz$patient %in% leftSubj,]$LeftOnset = 1
onlysz[onlysz$patient %in% rightSubj,]$RightOnset = 1
onlysz[onlysz$patient %in% bilateralSubj,]$LeftOnset = 1
onlysz[onlysz$patient %in% bilateralSubj,]$RightOnset = 1

length(unique(recallDat$patient))

length(unique(onlysz$patient)) 

## TEMPORAL LOBE PATIENTS
temporalSubj = NULL
k = 0
for (i in c(1,2,4,5,6,7,8,11,12,14,18)) {
  k = k + 1
  temporalSubj[k] = sprintf('I022_P0%.2d',i)
}
for (i in c(2,4,6,8,10,11,13,14,15,17,19,21,22,23,24,25,26,28,29,30,32,33,34,36,37,39,40,41,42,43,45,46,54)) {
  k = k + 1
  temporalSubj[k] = sprintf('I027_P0%.2d',i)
}

length(temporalSubj)
onlysz$TemporalOnset = 0
onlysz[onlysz$patient %in% temporalSubj,]$TemporalOnset = 1
onlysz$TemporalOnset = as.factor(onlysz$TemporalOnset)

#### SZ ONSET IN BA 20
szOnsetInLBA20 = NULL
k = 0
for (i in c(2,6,7,12)) {
  k = k + 1
  szOnsetInLBA20[k] = sprintf('I022_P0%.2d',i)
}
for (i in c(2,4,14,17,22,26,28,29,32,40,42,45,54)) {
  k = k + 1
  szOnsetInLBA20[k] = sprintf('I027_P0%.2d',i)
}

length(szOnsetInLBA20)
onlysz$szOnsetInLBA20 = 0
onlysz[onlysz$patient %in% szOnsetInLBA20,]$szOnsetInLBA20 = 1
onlysz$szOnsetInLBA20 = as.factor(onlysz$szOnsetInLBA20)

#### SZ ONSET IN BA 21
szOnsetInLBA21 = NULL
k = 0
for (i in c(4,6,12)) {
  k = k + 1
  szOnsetInLBA21[k] = sprintf('I022_P0%.2d',i)
}
for (i in c(4,21,22,28,29,32,34,40,43,46,54)) {
  k = k + 1
  szOnsetInLBA21[k] = sprintf('I027_P0%.2d',i)
}

length(szOnsetInLBA21)
onlysz$szOnsetInLBA21 = 0
onlysz[onlysz$patient %in% szOnsetInLBA21,]$szOnsetInLBA21 = 1
onlysz$szOnsetInLBA21 = as.factor(onlysz$szOnsetInLBA21)

length(unique(onlysz[onlysz$LeftOnset==1 & onlysz$RightOnset==0,]$patient))
length(unique(onlysz[onlysz$LeftOnset==0 & onlysz$RightOnset==1,]$patient))
length(unique(onlysz[onlysz$LeftOnset==1 & onlysz$RightOnset==1,]$patient))


#### FIGURE 1 ### SPIKE AND SZ LATERALIZATION #### 
tmp = melt(onlysz[,c('patient','onset','LeftCerebrum','RightCerebrum')])
tmpavg1 = ddply(tmp,.(patient,onset),summarize,totalSpikes=sum(value))
tmpavg2 = ddply(tmp,.(patient,onset,variable),summarize,totalSpikes=sum(value))
for (i in 1:length(tmpavg1$patient)) {
  pt = tmpavg1$patient[i]
  tmpspk = tmpavg2[tmpavg2$patient == pt,'totalSpikes'] 
  tmpspk = tmpspk/sum(tmpspk,na.rm=TRUE)
  tmpavg2[tmpavg2$patient == pt,'percentSpikes'] = tmpspk 
}

rlch = read.csv('rlch.csv')
rlch = as.data.frame(rlch)
tmpavg2 = cbind(tmpavg2,rlch$NumCh)
colnames(tmpavg2)[6] = 'NumCh'
tmpavg2$spikesperelec = tmpavg2$totalSpikes/tmpavg2$NumCh

#normalize by electrode
for (i in 1:length(tmpavg2$patient)) {
  pt = tmpavg2$patient[i]
  tmpspk = tmpavg2[tmpavg2$patient == pt,'spikesperelec'] 
  tmpspk = tmpspk/sum(tmpspk,na.rm=TRUE)
  tmpavg2[tmpavg2$patient == pt,'percentspikesperelec'] = tmpspk 
}

g = ggplot(tmpavg2,aes(x=onset,y=percentspikesperelec,colour=variable)) + geom_boxplot()
g = g + xlab('Seizure Onset') + ylab ('Percent Spikes/Electrode') +ggtitle('Spike and Seizure lateralization')
g = g + scale_color_manual("Spike lateralization",labels = c("Left Hemisphere", "Right Hemisphere"), values = c("blue", "red"))
g = g + theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6), legend.text=element_text(size=14),legend.title=element_text(size=14),legend.position='top') 
g
ggsave('spikelatvsonset.pdf',dpi=600,width=6,height=6)




t.test(tmpavg2[tmpavg2$onset=='LEFT' & tmpavg2$variable=='LeftCerebrum','percentspikesperelec'],tmpavg2[tmpavg2$onset=='LEFT' & tmpavg2$variable=='RightCerebrum','percentspikesperelec'])
cohen.d(tmpavg2[tmpavg2$onset=='LEFT' & tmpavg2$variable=='LeftCerebrum','percentspikesperelec'],tmpavg2[tmpavg2$onset=='LEFT' & tmpavg2$variable=='RightCerebrum','percentspikesperelec'],na.rm=TRUE)

t.test(tmpavg2[tmpavg2$onset=='RIGHT' & tmpavg2$variable=='LeftCerebrum','percentspikesperelec'],tmpavg2[tmpavg2$onset=='RIGHT' & tmpavg2$variable=='RightCerebrum','percentspikesperelec'])
cohen.d(tmpavg2[tmpavg2$onset=='RIGHT' & tmpavg2$variable=='LeftCerebrum','percentspikesperelec'],tmpavg2[tmpavg2$onset=='RIGHT' & tmpavg2$variable=='RightCerebrum','percentspikesperelec'],na.rm=TRUE)

################# SUMMARY STATISTICS AND MODELING #################
#calculate percent recall for each subject
avgdat<-ddply(recallDat, .(patient,age,sex,onset), summarize,numSessions = length(unique(session)), spikeMean=mean(spikes_ja_spatial), spikeTotal = sum(spikes_ja_spatial),recallMean=mean(as.numeric(recalled)-1))
mean(avgdat$age,na.rm=TRUE)
sd(avgdat$age,na.rm=TRUE)
#recall mean: 
mean(avgdat$recallMean)
sd(avgdat$recallMean)
min(avgdat$recallMean)
max(avgdat$recallMean)

#### SUPP FIG 1 #####

### age + sex + onset
avgdat<-ddply(onlysz, .(patient,age,sex,onset), summarize,numSessions = length(unique(session)), spikeMean=mean(spikes_ja_spatial), recallMean=mean(as.numeric(recalled)-1))
g1 = ggplot(data=avgdat,aes(x=age,y=recallMean)) + geom_point(size=2,aes(colour=sex)) + geom_smooth(method="lm") + xlab('Age') + ylab('Mean Recall Percentage')
linagemod = lm(dat=avgdat,recallMean ~ age + sex)
summary(linagemod)
g1

# SEX
sum(avgdat[!is.na(avgdat$sex),]$sex=='M')
sum(avgdat[!is.na(avgdat$sex),]$sex=='F')

#recall vs sex
g2 = ggplot(data=avgdat,aes(x=sex,y=recallMean)) + geom_boxplot() + xlab('Sex') + ylab('Mean Recall Percentage') 
t.test(avgdat[avgdat$sex=='M','recallMean'],avgdat[avgdat$sex=='F','recallMean'])
cohen.d(avgdat[avgdat$sex=='M','recallMean'],avgdat[avgdat$sex=='F','recallMean'])

#spikes vs sex
t.test(avgdat[avgdat$sex=='M','spikeMean'],avgdat[avgdat$sex=='F','spikeMean'])

## PLOT OF WORD ORDER
tmp<-ddply(recallDat, .(patient,age,sex,word_order), summarize, spikeMean=mean(spikes_ja_spatial), recallMean=mean(as.numeric(recalled)-1))
g3 = ggplot(tmp,aes(x=as.factor(word_order),y=recallMean)) + geom_boxplot() + xlab('Serial Word Position') + ylab('Mean Recall percentage')

## multiplot
pdf(file = 'supplementary_summ.pdf',width=8,height=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
print(g1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(g2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(g3,vp=viewport(layout.pos.row=2,layout.pos.col=1:2))
dev.off()


#### BUILDING THE GLMER ####
# Word order
logit1 = glmer(recalled ~(1|patient),data = recallDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit1b = glmer(recalled ~ word_order  + (1|patient),data = recallDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit1b,logit1)

# Sex
ageDat = recallDat
logit2 = glmer(recalled ~  word_order  + (1|patient),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit3 = glmer(recalled ~ word_order  + sex  + (1|patient),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit3,logit2) 

# Age
logit2 = glmer(recalled ~  word_order + (1|patient),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit3 = glmer(recalled ~ word_order  + age  +(1|patient),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit3,logit2) 

# Session
logit2 = glmer(recalled ~  word_order + age + (1|patient),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit3 = glmer(recalled ~ word_order  + age  + (1|patient/session),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit3,logit2) 

# Spikes
logit2 = glmer(recalled ~ spikes_ja_spatial+  word_order  + age  + (1|patient/session),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit3 = glmer(recalled ~ word_order  + age  + (1|patient/session),data = ageDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit3,logit2) 

#### SESSION ####
### number of sessions per patient
tmp = ddply(recallDat,.(patient,onset),summarize,numSessions= length(unique(session)))
ggplot(tmp,aes(x=patient,y=numSessions)) + geom_bar(stat='identity') + scale_x_discrete(labels=1:67) + xlab('Patient') + ylab('Number of Sessions') + ggtitle('Number of Sessions per patient')


############ Predicted probabilities vs spikes ###########
bootstrapCI = function(df,reg,mod,bootstrapN,bootstrapR) {
  #### confidence interval
  spikerange = 0:max(df[,reg],na.rm=T)
  my.bootstrap.predictions.f <- function(data, indices){
    return(mean(predict(mod, newdata = data[indices, ], type = "response", allow.new.levels=TRUE), na.rm=TRUE))
  }
  my.results <- matrix(nrow=length(spikerange), ncol = 4)
  for(x in 1:length(spikerange)){
    my.results[x, 1] =spikerange[x]
    new.df[,reg] = spikerange[x]
    my.boot.obj <- boot(data = new.df[sample(nrow(new.df), bootstrapN, replace = TRUE), ], statistic = my.bootstrap.predictions.f, R = bootstrapR)
    my.results[x, 2] <- my.boot.obj[[1]]
    my.results[x, 3:4] <- quantile(my.boot.obj[[2]], c(0.025, 0.975))
  }
  colnames(my.results) <- c("SpikeCount", "pp", "lower.ci", "upper.ci")
  my.results = data.frame(my.results)
  return(my.results)
}

#### BOOTSTRAP OUTSIDE SOZ #########
logit6= glmer(data=onlysz, recalled ~ word_order + age + non_SzOnset_spikes + (1|patient/session), family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
reg = 'non_SzOnset_spikes'
my.results = bootstrapCI(onlysz,reg,logit6,5000,1000)
g1 = ggplot(data=my.results,aes(x=SpikeCount,y=pp)) + geom_line() + geom_errorbar(aes(ymin=lower.ci,ymax=upper.ci)) + ylim(.2,.3) + xlab('Spike Count') + ylab('Probability of Recall') + ggtitle('Outside the SOZ') 
g1 = g1 + theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6))+ ylim(0.05,.3)
g1 = g1 + scale_x_continuous(breaks=0:maxSpikes)
ggsave(g1,file='pp.pdf',dpi=600)

######## BOOTSTRAP L BA 37 #############
tmpmod = glmer(data=recallDat, recalled ~ word_order + age + L_Brodmannarea37 + (1|patient/session), family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
reg = 'L_Brodmannarea37'
my.results = bootstrapCI(recallDat,reg,tmpmod,5000,1000)
g2 = ggplot(data=my.results,aes(x=SpikeCount,y=pp)) + geom_line() + geom_errorbar(aes(ymin=lower.ci,ymax=upper.ci)) + xlab('Spike Count') + ylab('Probability of Recall') + ggtitle('L BA 37') 
g2 = g2 + theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6))+ ylim(0.05,.3)
g2 = g2 + scale_x_continuous(breaks=0:maxSpikes)
ggsave(g2,file='ppba37.pdf',dpi=600)

######## L BA 20 #######
tmpmod = glmer(data=recallDat, recalled ~ word_order + age + L_Brodmannarea20 + (1|patient/session), family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
reg = 'L_Brodmannarea20'
my.results = bootstrapCI(recallDat,reg,tmpmod,5000,1000)
g3 = ggplot(data=my.results,aes(x=SpikeCount,y=pp)) + geom_line() + geom_errorbar(aes(ymin=lower.ci,ymax=upper.ci)) + xlab('Spike Count') + ylab('Probability of Recall') + ggtitle('L BA 20') 
g3 = g3 + theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6))+ ylim(0.05,.3)
g3 = g3 + scale_x_continuous(breaks=0:maxSpikes)
ggsave(g3,file='ppba20.pdf',dpi=600)

######## L BA 21 #########
tmpmod = glmer(data=recallDat, recalled ~ word_order + age + L_Brodmannarea21 + (1|patient/session), family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
reg = 'L_Brodmannarea21'
my.results = bootstrapCI(recallDat,reg,tmpmod,5000,1000)
g4 = ggplot(data=my.results,aes(x=SpikeCount,y=pp)) + geom_line() + geom_errorbar(aes(ymin=lower.ci,ymax=upper.ci)) + xlab('Spike Count') + ylab('Probability of Recall') + ggtitle('L BA 21') 
g4 = g4 + theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6)) + ylim(0.05,.3)
g4 = g4 + scale_x_continuous(breaks=0:maxSpikes)
ggsave(g4,file='ppba21.pdf',dpi=600)

pdf(file = 'ppcombined.pdf')
multiplot(g1,g2,g3,g4,cols=2)
dev.off()

############ GLMER MODELING ##########
# do you need session? yes
l1 = glmer(recalled ~ spikes_ja_spatial + word_order  + (1|patient),data = recallDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit1 = glmer(recalled ~ spikes_ja_spatial + word_order  + (1|patient),data = recallDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit1 = glmer(recalled ~ spikes_ja_spatial + word_order  + (1|patient),data = recallDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)

logit1 = glmer(recalled ~ spikes_ja_spatial + word_order  + (1|patient),data = recallDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit1b = glmer(recalled ~ spikes_ja_spatial + word_order  + (1|patient/session),data = recallDat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit1b,logit1)# (need session nest?) #yes

## SPIKES JA SPATIAL
logit1a = glmer(recalled ~ spikes_ja_spatial + word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(logit1a)
logit1b = glmer(recalled ~ word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(logit1b)
anova(logit1a,logit1b)

#sz onset spikes do not give more info than SzOnset_spikes
logitnull = glmer(recalled ~word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit = glmer(recalled ~ SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit,logitnull)

#once you take into account nonszonset, sz onset spikes do give more info
logit2 = glmer(recalled ~ SzOnset_spikes + non_SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit3 = glmer(recalled ~ non_SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(logit,logit2) #need non sz with sz? yes
anova(logit3,logit2) #need sz spikes with nonsz? yes

logit4 = glmer(recalled ~ word_order + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
logit5 = glmer(recalled ~ SzOnset_spikes + non_SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(logit5)

logit6 = glmer(recalled ~ non_SzOnset_spikes +word_order + age  +(1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(logit6)
dispersion_glmer(logit6) #shouldn't be over 1.4
cc <- confint(logit6,parm="beta_") 
ctab = cbind(est=fixef(logit6),cc)
print(ctab)

logit7 = glmer(recalled ~ SzOnset_spikes +word_order + age  +(1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(logit7)

anova(logit3,logit4) #need nonsz spikes? yes

summary(logit6)
exp(-coefficients(logit6)) 

modnoage = glmer(recalled ~ non_SzOnset_spikes +word_order + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
modage = glmer(recalled ~ non_SzOnset_spikes +word_order + age  +(1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
anova(modnoage,modage) #KEEP AGE? yes

############### REGIONAL ANALYSIS #####################
########## LATERALITY #########
#Difference in mean recall rates between left vs right sided onset?
tmp = ddply(onlysz,.(patient,age,sex,onset),summarize,meanrecall = mean(as.numeric(recalled)-1))

#recall worse with left onset

##### FIGURE LATERALITY #####
g = ggplot(tmp,aes(onset,meanrecall)) + geom_boxplot() + ylab('Mean Recall Percentage') + xlab('Seizure Onset') + ggtitle('Mean Recall by Seizure Onset Laterality')
g = g + theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6), legend.text=element_text(size=14),legend.title=element_text(size=14),legend.position='top') 
ggsave('lateralityvsmeanrecall.pdf',dpi=600,width=6,height=6)
t.test(tmp[tmp$onset=='RIGHT',]$meanrecall,tmp[tmp$onset=='LEFT',]$meanrecall)
cohen.d(tmp[tmp$onset=='RIGHT',]$meanrecall,tmp[tmp$onset=='LEFT',]$meanrecall)

#DOES ONSET PREDICT RECALL
logit1a = glmer(recalled ~ spikes_ja_spatial + word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(logit1a)
logit1b = glmer(recalled ~   spikes_ja_spatial  + onset + word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(logit1b)
anova(logit1a,logit1b)

## RIGHT ONSET
#in those with right sided onset, do spikes affect recall?
rnull = glmer(recalled ~ word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='RIGHT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(rnull)

rspikesall = glmer(recalled ~ spikes_ja_spatial +word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='RIGHT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(rspikesall)
rlr1 = anova(rnull,rspikesall)

#in those with right sided onset, do SZ onset spikes affect recall?
rszspikes = glmer(recalled ~ SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='RIGHT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(rnull)
checkConv(rszspikes)
rlr2 = anova(rnull,rszspikes)

#in those with right sided onset, do nonSZ onset spikes affect recall?
rnonszspikes = glmer(recalled ~ non_SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='RIGHT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(rnonszspikes)
rlr3 = anova(rnull,rnonszspikes)


#in those with right sided onset, do sz and nonSZ onset spikes affect recall?
rspikesboth = glmer(recalled ~ SzOnset_spikes + non_SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='RIGHT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(rspikesboth)
rlrb = anova(rnull,rspikesboth)
cc <- confint(rspikesboth,parm="beta_") 


## LEFT ONSET

#in those with left sided onset, do spikes affect recall?
lnull = glmer(recalled ~ word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='LEFT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
lspikesall = glmer(recalled ~ spikes_ja_spatial +word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='LEFT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(lnull)
checkConv(lspikesall)
anova(lnull,lspikesall)

#in those with left sided onset, do SZ onset spikes affect recall?
lszspikes = glmer(recalled ~ SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='LEFT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(lszspikes)
anova(lnull,lszspikes)

#in those with left sided onset, do nonSZ onset spikes affect recall?
lnonszspikes = glmer(recalled ~ non_SzOnset_spikes +word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='LEFT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(lnonszspikes)
anova(lnull,lnonszspikes)

#if already having nonsz spikes, do sz spikes give me information?
lspikesboth = glmer(recalled ~ non_SzOnset_spikes + SzOnset_spikes + word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='LEFT',], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(lspikesboth)
anova(lnull,lspikesboth)
cc2 <- confint(lspikesboth,parm="beta_") 

ltspikesboth = glmer(recalled ~ non_SzOnset_spikes + SzOnset_spikes + word_order + age + (1|patient/session),data = onlysz[onlysz$onset=='LEFT' & onlysz$TemporalOnset==1,], family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
checkConv(ltspikesboth)
anova(ltnull,ltspikesboth)

#if already having sz spikes, do nonsz spikes give me information? %YES
anova(logit3,logit2a)

########## LEFT VS RIGHT INTERACTION ############
#all pts
L_all = glmer(recalled ~ non_SzOnset_spikes*LeftOnset + SzOnset_spikes*LeftOnset + word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
cc3 <- confint(onsetinteraction,parm="beta_") 

R_all = glmer(recalled ~ non_SzOnset_spikes*RightOnset + SzOnset_spikes*RightOnset + word_order + age + (1|patient/session),data = onlysz, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
cc4 <- confint(LRinteraction,parm="beta_") 

#only unilateral
noBilat = onlysz[!(onlysz$onset=='BILATERAL'),]
L_uni = glmer(recalled ~ non_SzOnset_spikes*RightOnset + SzOnset_spikes*RightOnset + word_order + age + (1|patient/session),data = noBilat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
cc5 <- confint(lrin,parm="beta_") 

R_uni = glmer(recalled ~ non_SzOnset_spikes*LeftOnset + SzOnset_spikes*LeftOnset + word_order + age + (1|patient/session),data = noBilat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
cc6 <- confint(rin,parm="beta_") 

#only bilateral, too few patients
bilat = onlysz[(onlysz$onset=='BILATERAL'),]
L_b = glmer(recalled ~ non_SzOnset_spikes*LeftOnset + SzOnset_spikes*LeftOnset + word_order + age + (1|patient/session),data = bilat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
R_b = glmer(recalled ~ non_SzOnset_spikes*RightOnset + SzOnset_spikes*RightOnset + word_order + age + (1|patient/session),data = bilat, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)

####### REGIONAL ANALYSIS ########

#all patients included, for each analysis use if there are electrodes in the region
regions = names(onlysz[18:length(names(onlysz))])
idx = 1:170
specreg = NULL
mod = NULL
chip = NULL
chisq = NULL
numPt = NULL
b = NULL
k = 1
se = NULL
lowercc=NULL
uppercc=NULL
for (i in idx) {
  tmp = recallDat[!is.na(recallDat[,17+i]),]
  tmpnumPt = length(unique(tmp$patient))
  if (tmpnumPt > 20) {
    numPt[k] = tmpnumPt
    tmp = tmp[!is.na(tmp$age),]
    original = glmer(data=tmp,recalled~word_order + age + (1|patient/session),family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
    frmla = as.formula(paste('recalled ~ word_order + age + ',regions[i],'+ (1 + ',regions[i],'|patient/session)',sep=''))
    specreg[k] = regions[i]
    mod = glmer(data=tmp, frmla, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
    tmpsum = summary(mod)
    tmpcc = confint(mod,parm="beta_",method="Wald")
    lowercc[k] = tmpcc[4];
    uppercc[k] = tmpcc[8];
    summ = anova(original,mod)
    chisq[k] = summ$Chisq[2]
    chip[k] = summ$`Pr(>Chisq)`[2]
    b[k] = fixef(mod)[4]
    se[k] = tmpsum$coefficients[8]
    k = k +1
  }
}
idx = numPt>20 & !is.na(chip)
chisq2=chisq[idx]
p = chip[idx]
reg = specreg[idx]
numPt2 = numPt[idx]
b2 = b[idx]
se2 = se[idx]
ci = cbind(lowercc[idx],uppercc[idx])

tal = read.csv('TalairachRegions.csv')
## LEVEL5
regidx = NULL
for (i in 1:length(reg)) {
  if (length(grep(substr(reg[i],3,nchar(reg[i])),tal$Level5))>0) {
    regidx = append(regidx, i)
  }
}

lvl5idx = c(grep("Brodmann", reg),regidx)
newp = p[lvl5idx]
tmpreg = reg[lvl5idx]
pidx = p.adjust(newp,'holm')<0.05
tmpreg[pidx]
newp[pidx]
newchisq = chisq2[lvl5idx]

a = data.frame(NumPatients = numPt2[lvl5idx], Region = tmpreg, odds = exp(b2[lvl5idx]), ci_odds = exp(ci[lvl5idx,]), chisq=newchisq,P = newp, P.holm = p.adjust(newp,'holm'))
a = a[order(a$P),]
regionsAll = a
write.csv(regionsAll,file='regionsorig_rs.csv')


### Create figure with odds for each region ####
factorlist = levels(a$Region)
for (i in 1:length(factorlist)) {
  tmpstr = sub('_',' ',factorlist[i])
  tmpstr = gsub('_',',',tmpstr)
  tmpstr = sub('Brodmannarea','BA',tmpstr)
  tmpstr = sub('([[:alpha:]])([[:digit:]])','\\1 \\2',tmpstr)
  factorlist[i] = tmpstr
}

## LEVEL3
regidx = 1:length(reg)
regidx = regidx[reg %in% tal$Level3]
newp = p[regidx]
tmpreg = reg[regidx]
pidx = p.adjust(newp,'holm')<0.05
tmpreg[pidx]
newp[pidx]


## LEVEL2
regidx = 1:length(reg)
regidx = regidx[reg %in% tal$Level2]
newp = p[regidx]
tmpreg = reg[regidx]
pidx = p.adjust(newp,'holm')<0.05
tmpreg[pidx]
newp[pidx]

## LEVEL1
regidx = 1:length(reg)
regidx = regidx[reg %in% tal$Level1]
newp = p[regidx]
tmpreg = reg[regidx]
pidx = p.adjust(newp,'holm')<0.05
tmpreg[pidx]
newp[pidx]

pidx = p.adjust(p,'holm')<0.05
reg[pidx]
a = data.frame(NumPatients = numPt, Region = reg, log_odds = b2, se = se2, odds = exp(b2), P = p, P.holm = p.adjust(p,'holm'))
a = a[order(a$P),]

########################## FIGURE - ODDS FOR ALL REGIONS   ######################
shading <- data.frame(min = seq(from = 0.5, to = max(as.numeric(a$Region)), by = 1),
                      max = seq(from = 1.5, to = max(as.numeric(a$Region)) + 0.5, by = 1),
                      col = c(0,1))
g = ggplot(dat=a,aes(x=Region,y=exp(log_odds),colour=P.holm<0.05)) 
g = g + geom_point() + geom_errorbar(aes(ymin=exp(log_odds-se),ymax=exp(log_odds+se))) 
g = g + ylab('Odds') + xlab('Talairach Region') + ggtitle('Odds of Successful Recall Per Spike')
g = g + scale_color_manual("Corrected P\n",labels = c("P>0.05", "P<0.05"), values = c("blue", "red"))
g = g + theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6), legend.text=element_text(size=14)) 
g = g + geom_hline(yintercept=1)
g = g + coord_flip()
ggsave(g,file='regionodds.pdf',dpi=600)

## Create figure with odds for each region 
factorlist = levels(a$Region)
for (i in 1:length(factorlist)) {
  tmpstr = sub('_',' ',factorlist[i])
  tmpstr = gsub('_',',',tmpstr)
  tmpstr = sub('Brodmannarea','BA',tmpstr)
  tmpstr = sub('([[:alpha:]])([[:digit:]])','\\1 \\2',tmpstr)
  factorlist[i] = tmpstr
}
shading <- data.frame(min = seq(from = 0.5, to = max(as.numeric(a$Region)), by = 1),
                      max = seq(from = 1.5, to = max(as.numeric(a$Region)) + 0.5, by = 1),
                      col = c(0,1))
g = ggplot(dat=a,aes(x=Region,y=odds,colour=P.holm<0.05)) 
g = g + geom_point() + geom_errorbar(aes(ymin=ci_odds.1,ymax=ci_odds.2)) 
g = g + ylab('Odds') + xlab('Region') + ggtitle('Odds of Successful Recall Per Spike')
g = g + scale_color_manual("Corrected P\n",labels = c("P>0.05", "P<0.05"), values = c("blue", "red"))
g = g + theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=18, face="bold", margin = margin(10, 0, 10, 0),vjust=1, lineheight=0.6), legend.text=element_text(size=14),legend.title=element_text(size=14)) 
g = g + geom_hline(yintercept=1) + scale_x_discrete(labels=factorlist)
g = g + coord_flip()
ggsave(g,file='regionodds_BA.pdf',dpi=600)

########## RIGHT ONLY ########
regions = names(onlysz[18:length(names(onlysz))])
idx = 1:170
specreg = NULL
mod = NULL
chip = NULL
chisq = NULL
numPt = NULL
b = NULL
k = 1
se = NULL
lowercc=NULL
uppercc=NULL
for (i in idx) {
  #tmp = recallDat[!is.na(recallDat[,17+i]),]
  tmp = onlysz[onlysz$RightOnset==1,]
  tmp = tmp[!is.na(tmp[,17+i]),]
  #tmp = tmpDat[!is.na(tmpDat[,17+i]),]
  tmpnumPt = length(unique(tmp$patient))
  if (tmpnumPt > 20) {
    numPt[k] = tmpnumPt
    tmp = tmp[!is.na(tmp$age),]
    original = glmer(data=tmp,recalled~word_order + age + (1|patient/session),family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
    frmla = as.formula(paste('recalled ~ word_order + age + ',regions[i],'+ (1 + ',regions[i],'|patient/session)',sep=''))
    specreg[k] = regions[i]
    mod = glmer(data=tmp, frmla, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
    tmpsum = summary(mod)
    tmpcc = confint(mod,parm="beta_",method="Wald")
    lowercc[k] = tmpcc[4];
    uppercc[k] = tmpcc[8];
    summ = anova(original,mod)
    chisq[k] = summ$Chisq[2]
    chip[k] = summ$`Pr(>Chisq)`[2]
    b[k] = fixef(mod)[4]
    se[k] = tmpsum$coefficients[8]
    k = k +1
  }
}
idx = numPt>20 & !is.na(chip)
chisq2=chisq[idx]
p = chip[idx]
reg = specreg[idx]
numPt2 = numPt[idx]
b2 = b[idx]
se2 = se[idx]
ci = cbind(lowercc[idx],uppercc[idx])

tal = read.csv('TalairachRegions.csv')
## LEVEL5
regidx = NULL
for (i in 1:length(reg)) {
  if (length(grep(substr(reg[i],3,nchar(reg[i])),tal$Level5))>0) {
    regidx = append(regidx, i)
  }
}

lvl5idx = c(grep("Brodmann", reg),regidx)
newp = p[lvl5idx]
tmpreg = reg[lvl5idx]
pidx = p.adjust(newp,'holm')<0.05
tmpreg[pidx]
newp[pidx]
newchisq = chisq2[lvl5idx]

a = data.frame(NumPatients = numPt2[lvl5idx], Region = tmpreg, odds = exp(b2[lvl5idx]), ci_odds = exp(ci[lvl5idx,]), chisq=newchisq,P = newp, P.holm = p.adjust(newp,'holm'))
a = a[order(a$P),]
regionsRight = a
write.csv(regionsRight,file='regionsRight.csv')

########## LEFT ONLY ########
tmpDat = onlysz[onlysz$onset=='LEFT',]
regions = names(tmpDat[18:length(names(tmpDat))])
idx = 1:length(regions)
specreg = NULL
mod = NULL
chip = NULL
numPt = NULL
b = NULL
k = 1
se = NULL
lowercc=NULL
uppercc=NULL
for (i in idx) {
  #tmp = recallDat[!is.na(recallDat[,17+i]),]
  #tmp = onlysz[onlysz$LeftOnset==1 & onlysz$RightOnset==0,]
  #tmp = tmp[!is.na(tmp[,17+i]),]
  tmp = tmpDat[!is.na(tmpDat[,17+i]),]
  tmpnumPt = length(unique(tmp$patient))
  if (tmpnumPt > 20) {
    numPt[k] = tmpnumPt
    original = glmer(data=tmp,recalled~word_order + age + (1|patient/session),family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
    frmla = as.formula(paste('recalled ~ word_order + age + ',regions[i],'+ (1|patient/session)',sep=''))
    specreg[k] = regions[i]
    mod = glmer(data=tmp, frmla, family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=1)
    tmpsum = summary(mod)
    tmpcc = confint(mod,parm="beta_",method="Wald")
    lowercc[k] = tmpcc[4];
    uppercc[k] = tmpcc[8];
    summ = anova(original,mod)
    chip[k] = summ$`Pr(>Chisq)`[2]
    b[k] = fixef(mod)[4]
    se[k] = tmpsum$coefficients[8]
    k = k +1
  }
}
idx = numPt>20 & !is.na(chip)
p = chip[idx]
reg = specreg[idx]
numPt2 = numPt[idx]
b2 = b[idx]
se2 = se[idx]
ci = cbind(lowercc[idx],uppercc[idx])

tal = read.csv('TalairachRegions.csv')
## LEVEL5
regidx = NULL
for (i in 1:length(reg)) {
  if (length(grep(substr(reg[i],3,nchar(reg[i])),tal$Level5))>0) {
    regidx = append(regidx, i)
  }
}

lvl5idx = c(grep("Brodmann", reg),regidx)
newp = p[lvl5idx]
tmpreg = reg[lvl5idx]
pidx = p.adjust(newp,'holm')<0.05
tmpreg[pidx]
newp[pidx]

a = data.frame(NumPatients = numPt2[lvl5idx], Region = tmpreg, odds = exp(b2[lvl5idx]), ci_odds = exp(ci[lvl5idx,]), P = newp, P.holm = p.adjust(newp,'holm'))
a = a[order(a$P),]
regionsLeft = a
write.csv(regionsLeft,file='regionsleft.csv')


######## PLOT OF RAW DATA ACROSS REGIONS ########

#difference in spikes in various regions of the brain between recall and not recalled by onset
tmp = recallDat[recallDat$word_order>5,]
a = ddply(tmp, .(patient,onset,recalled),summarize,avgSpikes = mean(spikes_ja_spatial),
          avgSpikesNonOnset = mean(non_SzOnset_spikes),avgSpikesOnset = mean(SzOnset_spikes),
          lc = mean(LeftCerebrum),lt = mean(L_TemporalLobe), l37 = mean(L_Brodmannarea37), l21 = mean(L_Brodmannarea21),l20=mean(L_Brodmannarea20))
#all
b = ddply(a, .(patient,onset),summarize,dall=diff(avgSpikes)/max(abs(avgSpikes)),dASO=diff(avgSpikesOnset)/max(abs(avgSpikesOnset)),dASNO=diff(avgSpikesNonOnset)/max(abs(avgSpikesNonOnset)),dc=diff(lc)/max(abs(lc)),dt=diff(lt)/max(abs(lt)),
          dlc=diff(lc)/max(abs(lc)),d37=diff(l37)/max(abs(l37)),d21=diff(l21)/max(abs(l21)),d20=diff(l20)/max(abs(l20)))

c = melt(b)
c$value = c$value*-1
c = c[c$variable %in% c('dt','d37','d20','d21'),]
ggplot(dat=c,aes(x=variable,y=value,colour=onset)) + geom_boxplot() + geom_hline(yintercept=0) + xlab('Region') +
  ylab('Percent change in spike count') + scale_x_discrete(labels=c('L Temporal','L BA 20','L BA 21','L BA 37')) +
  scale_color_discrete(name="Seizure Onset") + ggtitle('Percent change in spike count during failed recall') + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold"),  plot.title = element_text(size=16, face="bold",vjust=1, lineheight=0.6), legend.text=element_text(size=12),legend.title=element_text(size=12),legend.position='right') 
ggsave('percentspikechange.pdf',dpi=600,height=8,width=10)





