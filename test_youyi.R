source("myhelper.R")
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')

#data( 'dat.505', package = 'HVTN505' )
#suppressWarnings( data( 'var.super', package = 'HVTN505' ) )
dat.505 <- read.table( 'dat_505.txt', header = T, sep = ',' )
dat.505 <- subset(dat.505, trt==1) # only treatment
var.super <- read.table( 'var_super.txt', header = T, sep = ',' )

# Scaling each variables: mean 0 and sd 1 #
dat.505[,c('age', 'BMI', 'bhvrisk', var.super$varname)]=scale(dat.505[,c('age', 'BMI', 'bhvrisk', var.super$varname)])

Y_vaccine <- dat.505$case
X_vaccine <- dat.505 %>% select(age, BMI, bhvrisk, var.super$varname) # without _bin variables
#rownames(X_vaccine) <- 1:150 # to assign 1:150 for treatment
weights_vaccine <- dat.505$wt
strata_vaccine <- dat.505$hg_strata


### 1. Setting process arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are need to a batch script
  # sim.setting: rv144ph2_n10000_1 for RV144 experiments in Section 3
  # fit.setting: 5fold for 5-fold cross validation
  # proj: you can choose one of tests listed in Section 3
  # ipw: uw for unweighted; sw for semi-weights; w for weighted
  Args=c(batch.size="10",batch.number="1",can.set='RF:tcell+glm:antibody',method="ST")
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; can.set=Args[i]
i=i+1; method=Args[i]



# Specifying screened variables index #
# Base domain for stacking #
#dat <- screen.dat.varidx(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set='all', screen.method='no')$dat

# For all markers #
screen.var.all <- screen.dat.varidx(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set='all', screen.method='lasso')$screen.index.var

# For t cell markers #
screen.var.tcell.temp <- screen.dat.varidx(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set='tcell', screen.method='lasso')$screen.index.var
screen.var.tcell <- rep(FALSE,ncol(X_vaccine))
for(i in 1:length(which(screen.var.tcell.temp))){
  temp <- which(names(X_vaccine) == names(which(screen.var.tcell.temp)[i]))
  screen.var.tcell[temp] <- TRUE
}
screen.var.tcell

# For antibody markers #
screen.var.antibody.temp <- screen.dat.varidx(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set='antibody', screen.method='lasso')$screen.index.var
screen.var.antibody <- rep(FALSE,ncol(X_vaccine))
for(i in 1:length(which(screen.var.antibody.temp))){
  temp <- which(names(X_vaccine) == names(which(screen.var.antibody.temp)[i]))
  screen.var.antibody[temp] <- TRUE
}
screen.var.antibody

var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.tcell)) # GLM: antidoby + RF: tcell

## Stacking #
#pred.vec = sapply(1:10, function(seed){
#    myprint(seed)
#    
##    Y=Y_vaccine; X=X_vaccine; obsWeights=weights_vaccine; strata=strata_vaccine; var.index=var.index.set; method='method.NNloglik'; mc.cores=1;
#
#    res.temp <- get.st.cvauc(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, strata=strata_vaccine, var.index=var.index.set, method='method.NNloglik', seed=seed, mc.cores=1)
#    mean(do.call(rbind,res.temp)[,'est.cvauc'])
#})
#mean(pred.vec)


Y=Y_vaccine; X=X_vaccine; obsWeights=weights_vaccine; strata=strata_vaccine; var.index=var.index.set; mc.cores=1;

var.index=var.index.set; mc.cores=1;



res=sapply(0:10 ,function(seed) {
    set.seed(seed)
    
    replace=seed>0
    
    idx.stratum1.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[18.4, 25)"), replace=replace)
    idx.stratum1.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[18.4, 25)"), replace=replace)
    idx.stratum2.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[25, 29.8)"), replace=replace)
    idx.stratum2.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[25, 29.8)"), replace=replace)
    #idx.stratum3.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[29.8, 40)"), replace=replace)
    idx.stratum3.case <- 134
    idx.stratum3.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[29.8, 40)"), replace=replace)
    idx.stratum4.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[18.4, 25)"), replace=replace)
    idx.stratum4.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[18.4, 25)"), replace=replace)
    idx.stratum5.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[25, 29.8)"), replace=replace)
    idx.stratum5.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[25, 29.8)"), replace=replace)
    idx.stratum6.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[29.8, 40)"), replace=replace)
    idx.stratum6.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[29.8, 40)"), replace=replace)
    
    idx.case <- c(idx.stratum1.case,idx.stratum2.case,idx.stratum3.case,idx.stratum4.case,idx.stratum5.case,idx.stratum6.case)
    idx.control <- c(idx.stratum1.control,idx.stratum2.control,idx.stratum3.control,idx.stratum4.control,idx.stratum5.control,idx.stratum6.control)
    
    Y_vaccine_boot <- Y_vaccine[c(idx.case, idx.control)]
    X_vaccine_boot <- X_vaccine[c(idx.case, idx.control),]
    weights_vaccine_boot <- weights_vaccine[c(idx.case, idx.control)]
    strata_vaccine_boot <- strata_vaccine[c(idx.case, idx.control)]
    
    Y=Y_vaccine_boot; X=X_vaccine_boot; obsWeights=weights_vaccine_boot; strata=strata_vaccine_boot; 
    
    split <- random.strata.splits(Y, X, strata, seed)[[1]]  
  
    # Target variable must be named as 'Y' #
    dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
    dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
    # Y=1,0 must be converted to factor that has two levels of 'case' (for Y=1) and 'control' (for Y=0) #
    dat.train$Y <- as.factor(dat.train$Y) ; levels(dat.train$Y) <- c('control','case')
    weights.train <- obsWeights[c(split$training$case,split$training$control)] ; 
    weights.test <- obsWeights[c(split$test$case,split$test$control)]
    
    # RF
    set.seed(123)
    fit.rf <- ranger( factor(Y)~., data = dat.train[,var.index$rf], probability = TRUE, min.node.size = 1, case.weights=weights.train ) # Weights
    fit.rf
    WeightedAUC(WeightedROC(guess=predict(fit.rf, dat.train)$predictions[,"case"] , label=dat.train$Y, weight=weights.train))    
    pred.rf = predict(fit.rf, dat.test)$predictions[,"case"]    
  
    
    # glm
    fit.glm <- glm( as.numeric(Y)~., data = dat.train[,var.index$rf], weights=weights.train ) # Weights
    pred.glm = predict(fit.glm, dat.test)
    
    cvauc.rf  <- WeightedAUC(WeightedROC(guess=pred.rf , label=dat.test$Y, weight=weights.test))    
    cvauc.glm <- WeightedAUC(WeightedROC(guess=pred.glm, label=dat.test$Y, weight=weights.test))    
    c(rf=cvauc.rf, glm=cvauc.glm)
    
})
res




#
#
#### 2. Experiments ###
#res=sapply(seeds, simplify="array", function (seed) {
#  
#  myprint(seed)
#  
#  table(strata_vaccine, Y_vaccine)
#  
#  # Bootstrap #
#  set.seed(seed)
#  idx.stratum1.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[18.4, 25)"), replace=TRUE)
#  idx.stratum1.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[18.4, 25)"), replace=TRUE)
#  idx.stratum2.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[25, 29.8)"), replace=TRUE)
#  idx.stratum2.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[25, 29.8)"), replace=TRUE)
#  #idx.stratum3.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[29.8, 40)"), replace=TRUE)
#  idx.stratum3.case <- 134
#  idx.stratum3.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[29.8, 40)"), replace=TRUE)
#  idx.stratum4.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[18.4, 25)"), replace=TRUE)
#  idx.stratum4.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[18.4, 25)"), replace=TRUE)
#  idx.stratum5.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[25, 29.8)"), replace=TRUE)
#  idx.stratum5.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[25, 29.8)"), replace=TRUE)
#  idx.stratum6.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[29.8, 40)"), replace=TRUE)
#  idx.stratum6.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[29.8, 40)"), replace=TRUE)
#  
#  idx.case <- c(idx.stratum1.case,idx.stratum2.case,idx.stratum3.case,idx.stratum4.case,idx.stratum5.case,idx.stratum6.case)
#  idx.control <- c(idx.stratum1.control,idx.stratum2.control,idx.stratum3.control,idx.stratum4.control,idx.stratum5.control,idx.stratum6.control)
#  Y_vaccine_boot <- Y_vaccine[c(idx.case, idx.control)]
#  X_vaccine_boot <- X_vaccine[c(idx.case, idx.control),]
#  weights_vaccine_boot <- weights_vaccine[c(idx.case, idx.control)]
#  strata_vaccine_boot <- strata_vaccine[c(idx.case, idx.control)]
#  
#  # Generating data where stacking is fitted #
#  dat.X <- screen.dat.index(Y=Y_vaccine_boot, X=X_vaccine_boot, obsWeights=weights_vaccine_boot, fit.set='all', screen.method='no')$dat
#  
#  dat.tmp=(data.frame(Y_vaccine_boot, X_vaccine_boot[,var.index.set$rf]))
#  
#    set.seed(123)
#    fit.rf <- ranger( factor(Y_vaccine_boot)~., data = dat.tmp, probability = TRUE, min.node.size = 1, case.weights=weights_vaccine_boot ) # Weights
#
#
#  if(can.set=='RF:tcell+glm:antibody'){
#    # RF:tcell + glm:antibody #
#    var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.tcell))
#  } else if(can.set=='RF:tcell+glm:all'){
#    # RF:tcell + glm:all #
#    var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.tcell))
#  } else if(can.set=='RF:all+glm:antibody'){
#    # RF:all + glm:antibody #
#    var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.all))
#  } else if(can.set=='RF:all+glm:all'){
#    # RF:all + glm:all #
#    var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.all))
#  }
#
#
##
#
#    
#
#
#  
#  # Stacking (CV-AUC) #
#  pred.vec <- c()
#  for( i in 1:500 ){
#    print(i)
#    seed <- i
#    # Stacking (CV-AUC) #
#    res.temp <- get.st.cvauc(dat=dat.X, obsWeights=weights_vaccine, strata=strata_vaccine, var.index=var.index.set, method='method.NNloglik', seed=seed)
#    pred.vec[i] <- mean(do.call(rbind,res.temp)[,'est.cvauc'])
#  }
#  mean(pred.vec)
#})
#res
#
#
#
## save res
#foldername="res_"%.%method%.%"_"%.%can.set%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
#save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")


library(kyotil)
library(aucm)
load("D:/downloads/RF_tcell+GLM_antibody_res1.Rdata")
load("D:/downloads/RF_tcell+GLM_antibody_res2.Rdata")

colMeans(res1.total[1:5,])

# take the first seed results
tmp=data.frame(res2.total[1:150,])
fast.auc(tmp[,"glm"], tmp[,"case"])
fast.auc(tmp[,"ranger"], tmp[,"case"])
fast.auc(tmp[,"stacking"], tmp[,"case"])

plot(glm~ranger, tmp, col=tmp[,"case"]+1)

myfigure(mfrow=c(1,3))
    myboxplot(ranger~case, tmp, main="ranger", names=c("controls","cases"), col=1+tmp$case); title(sub="CV-AUC = "%.%formatDouble(fast.auc(tmp[,"ranger"], tmp[,"case"]), 3), line=3)
    myboxplot(glm~case, tmp, main="glm", names=c("controls","cases"), col=1+tmp$case); title(sub="CV-AUC = "%.%formatDouble(fast.auc(tmp[,"glm"], tmp[,"case"]), 3), line=3)
    myboxplot(stacking~case, tmp, main="stacking", names=c("controls","cases"), col=1+tmp$case); title(sub="CV-AUC = "%.%formatDouble(fast.auc(tmp[,"stacking"], tmp[,"case"]), 3), line=3)
mydev.off(file="predscores_boxplot")


tmp[with(tmp, order(stacking, decreasing=T)[1:4]),]

myfigure(mfrow=c(1,2))
    tmp$ranger2=1*tmp$ranger
    plot(glm~ranger, tmp, col=tmp[,"case"]+1, xlim=c(0,1), ylim=c(0,1), cex=.7, pch=ifelse(tmp$ptid %in% c(180,183), 2, 1))
    abline(v=c(.45,.58), lty=2)
    abline(0,1)
    
    plot(stacking~ranger, tmp, col=tmp[,"case"]+1, xlim=c(0,1), ylim=c(0,1), cex=.7, pch=ifelse(tmp$ptid %in% c(180,183), 2, 1))
    abline(0,1)    
    abline(v=c(.45,.58), lty=2)
    abline(h=c(.4), lty=2)
mydev.off(file="predscores_scatterplot")
