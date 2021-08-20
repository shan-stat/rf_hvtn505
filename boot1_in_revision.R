##### Bootstrap method 1 #####
# Description #
# (1) bootstrap is for the variability of CV-AUC estimates
# (2) bootstrap is performed to this original data, and the bootstrapped sample is used to split training subset and validation subset in random stratified 4:1 cv#



# Required functions (in addition to all function in the "myhelper.R" file) #
# Random stratified 4:1 #
source("myhelper.R")
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')

random.strata.splits <- function (Y, X, strata, seed){
  set.seed(seed)
  case.val.idx <- sample(which(Y==1), size=sum(Y==1)/5, replace=FALSE) # randomly sample 1/5 cases 
  case.val.freq <- table(strata[case.val.idx]) # strata information for the sampled cases
  case.total.freq <- table(strata[which(Y==1)]) # strata information for the entire cases (for 25 cases, not the sampled cases)
  
  # compute the ratio of that how many cases are sampled from each stratum #
  # e.g. Suppose strata1 - 4cases/20controls. If 2cases are randomly sampled, the ratio will be 2/4=0.5 #
  ratio <- c()
  for(i in 1:length(case.val.freq)){
    freq.temp <- case.val.freq[i]
    ratio[i] <- freq.temp/case.total.freq[names(freq.temp)]
  }
  names(ratio) <- names(case.val.freq)
  
  # randomly sample controls depending on the ratio of that how many cases are sampled in each stratum #
  # e.g. if the ratio is 0.5, 20*0.5 controls will be sampled in strata 1 #
  control.val.idx <- c()
  for(i in 1:length(ratio)){
    control.stratum <- which(strata==names(ratio[i]) & Y==0)
    control.temp <- sample(control.stratum, size=round(length(control.stratum)*ratio[i]), replace=FALSE)
    control.val.idx <- c(control.val.idx, control.temp)
  }
  
  splits=list()
  case.train.idx <- setdiff(which(Y==1), case.val.idx)
  control.train.idx <- setdiff(which(Y==0), control.val.idx)
  splits[[1]] <- list(training=list(case=case.train.idx, control=control.train.idx),
                      test=list(case=case.val.idx, control=control.val.idx))
  splits
}

# Stacking with random stratified cv #
get.st.cvauc.random = function(Y, X, obsWeights, strata, var.index, method, seed=1, mc.cores=mc.cores){
  
  # Outer layer: Random stratified 4:1 cv #
  splits <- random.strata.splits(Y, X, strata, seed)
  
  # Inner layer: 10-fold CV (fitting learners and obtaining out-of-sample predictions) #
  my_control <- trainControl(
    method="cv",
    number=10,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=twoClassSummary
  )
  
  # Method of combining predictions from different learners #
  method <- get(method, mode = 'function')()
  if(!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), character.only = TRUE))
  }
  
  # Fitting candidate learners #
  cv.aucs <-  mclapply( splits, function(split){
    # Target variable must be named as 'Y' #
    dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
    dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
    # Y=1,0 must be converted to factor that has two levels of 'case' (for Y=1) and 'control' (for Y=0) #
    dat.train$Y <- as.factor(dat.train$Y) ; levels(dat.train$Y) <- c('control','case')
    weights.train <- obsWeights[c(split$training$case,split$training$control)] ; weights.test <- obsWeights[c(split$test$case,split$test$control)]
    
    # Setting three RF hyperparameters to their default values (without this codes, hyperparameter tuning is automatically done) #
    rf_grid <- expand.grid(mtry = floor(sqrt(sum(var.index$rf)-1)), splitrule = 'gini', min.node.size = 1)
    set.seed( 123 )
    suppressWarnings(
      model_list <- caretList(
        list(Y~., data=dat.train[,var.index$glm], weights = weights.train), # glm
        list(Y~., data=dat.train[,var.index$rf], tuneGrid = rf_grid, weights = weights.train), # RF
        trControl=my_control,
        methodList=c('glm', 'ranger')
      )
    )
    # Results of caretList() #
    # 1. model_list$ranger$pred: contains out-of-sample prediction scores (from 10-fold inner layer CV), and they are used for fitting meta-learner 
    # 2. model_list$ranger$finalModel: contains final RF model on the entire training data. At the final stage of stacking, test data is applied to this final model to generate test prediction scores 
    
    # For RF learners, final model should be set manually because of a random seed issue from caret.train() (for glm, it is not happned) #
    # caretList() fits RF using ranger() function, so if random seed is the same in caretList() and ranger() functions, the two results should be same #
    # but, the results are a bit different (note that their hyperparameter setting is the same) #
    # This issue is not appeared in GLM. The results from caretList() and glm() functions are the same under the same random seed #
    set.seed(123)
    fit.rf <- ranger( factor(Y)~., data = dat.train[,var.index$rf], probability = TRUE, min.node.size = 1, case.weights=weights.train ) # Weights
    model_list$ranger$finalModel <- fit.rf
    
    # stacking #
    set.seed( 123 )
    pred.cv <- makePredObsMatrix(model_list) # pred.cv is out-of-sample prediction scores. They are the same as "case" column in model_list$glm$pred and model_list$ranger$pred.
    res.st.fit <- method$computeCoef(Z = pred.cv$preds, Y =(as.numeric(pred.cv$obs)-1), obsWeights = weights.train, libraryNames = names(var.index), verbose = FALSE, control=list(trimLogit=0.001))
    pred.test <- predict.caretList(model_list, newdata=dat.test) # test data is applied to the final model from model_list
    pred.st <- method$computePred(predY = pred.test, coef = res.st.fit$coef, control=list(trimLogit=0.001))
    
    # Calculating CV-AUC #
    est.cvauc <- WeightedAUC(WeightedROC(guess=pred.st, label=dat.test$Y, weight=weights.test))
    # Calculating Pearson correlation coefficient #
    est.corr <- cor(pred.cv$preds[,1], pred.cv$preds[,2])
    c(est.cvauc=est.cvauc, est.corr=est.corr)
  }, mc.cores = mc.cores )
  cv.aucs
}



### 0. Pre-processing HVTN 505 dataset ###
data( 'dat.505', package = 'HVTN505' )
dat.505 <- subset(dat.505, trt==1) # only treatment
suppressWarnings( data( 'var.super', package = 'HVTN505' ) )

# Scaling each variables: mean 0 and sd 1 #
#dat.505[,c('age', 'BMI', 'bhvrisk', var.super$varname)]=scale(dat.505[,c('age', 'BMI', 'bhvrisk', var.super$varname)])
for( i in var.super$varname ){
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
  dat.505[[i %.% '_bin']] <- scale( dat.505[[i %.% '_bin']], center = mean(dat.505[[i %.% '_bin']][dat.505$trt == 1]), scale = sd(dat.505[[i %.% '_bin']][dat.505$trt == 1]))
}
for( i in c( 'age', 'BMI', 'bhvrisk' ) ) {
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
}

Y_vaccine <- dat.505$case
X_vaccine <- dat.505 %>% select(age, BMI, bhvrisk, var.super$varname, paste0(var.super$varname, '_bin'))
weights_vaccine <- dat.505$wt
strata_vaccine <- dat.505$hg_strata



### 1. Setting process arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are need to a batch script
  # can.set.st: one of c('RF:tcell+glm:antibody','RF:tcell+glm:all','RF:all+glm:antibody','RF:tcell+glm:all')
  # can.set.rf: one of c('all','tcell','antibody','no')
  # proj: ST or RF
  Args=c(batch.size="10",batch.number="1",can.set.st='RF:tcell+glm:antibody',can.set.rf='all',proj="ST")
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; can.set.st=Args[i]
i=i+1; can.set.rf=Args[i]
i=i+1; proj=Args[i]

### 2. Specifying screened variables index ###
# This process can be used for stacking and random forest #
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

### 3. Experiments ###
res=sapply(seeds, simplify="array", function (seed) {
  
  myprint(seed)
  replace=seed>0 # seed=0 for bootstrap without replacement; seed>0 for bootstrap with replacement
  
  # Bootstrap based on strata information #
  set.seed(seed)
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
  
  # Stacking #
  if(proj=='ST'){
    # specify candidate learners #
    if(can.set.st=='RF:tcell+glm:antibody'){
      # RF:tcell + glm:antibody #
      var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.tcell))
    } else if(can.set.st=='RF:tcell+glm:all'){
      # RF:tcell + glm:all #
      var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.tcell))
    } else if(can.set.st=='RF:all+glm:antibody'){
      # RF:all + glm:antibody #
      var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.all))
    } else if(can.set.st=='RF:all+glm:all'){
      # RF:all + glm:all #
      var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.all))
    }
    
    # cv-auc #
    pred.vec <- c()
    for( i in 1:500 ){
      print(i)
      seed <- i
      res.temp <- get.st.cvauc.random(Y=Y_vaccine_boot, X=X_vaccine_boot, obsWeights=weights_vaccine_boot, strata=strata_vaccine,
                                      var.index=var.index.set, method='method.NNloglik', seed=seed, mc.cores=1)
      pred.vec[i] <- mean(do.call(rbind,res.temp)[,'est.cvauc'])
    }
    mean(pred.vec)
    
  } else if (proj=='RF'){
    # specify candidate set #
    if(can.set.rf=='all'){
      var.index.set <- c(TRUE,screen.var.all)
    } else if(can.set.rf=='tcell'){
      var.index.set <- c(TRUE,screen.var.tcell)
    } else if(can.set.rf=='antibody'){
      var.index.set <- c(TRUE,screen.var.antibody)
    } else if(can.set.rf=='no'){
      var.index.set <- c(TRUE,rep(TRUE,3),rep(FALSE,420))
    }
    
    # random stratified 4:1 cv #
    split <- random.strata.splits(Y=Y_vaccine_boot, X=X_vaccine_boot, strata=strata_vaccine_boot, seed)[[1]]
    
    # Target variable must be named as 'Y' #
    dat.train <- cbind(Y=Y_vaccine_boot[c(split$training$case,split$training$control)], X_vaccine_boot[c(split$training$case,split$training$control),,drop=F])
    dat.test <- cbind(Y=Y_vaccine_boot[c(split$test$case,split$test$control)], X_vaccine_boot[c(split$test$case,split$test$control),,drop=F])
    # Y=1,0 must be converted to factor that has two levels of 'case' (for Y=1) and 'control' (for Y=0) #
    dat.train$Y <- as.factor(dat.train$Y) ; levels(dat.train$Y) <- c('control','case')
    weights.train <- weights_vaccine_boot[c(split$training$case,split$training$control)]
    weights.test <- weights_vaccine_boot[c(split$test$case,split$test$control)]
    
    # RF
    set.seed(123)
    fit.rf <- ranger( factor(Y)~., data = dat.train[,var.index.set], probability = TRUE, min.node.size = 1, case.weights=weights.train ) # Weights
    pred.rf = predict(fit.rf, dat.test)$predictions[,"case"]
    
    # glm
    fit.glm <- glm( as.numeric(Y)~., data = dat.train[,var.index.set], weights=weights.train ) # Weights
    pred.glm = predict(fit.glm, dat.test)
    
    cvauc.rf  <- WeightedAUC(WeightedROC(guess=pred.rf , label=dat.test$Y, weight=weights.test))    
    cvauc.glm <- WeightedAUC(WeightedROC(guess=pred.glm, label=dat.test$Y, weight=weights.test))    
    c(rf=cvauc.rf, glm=cvauc.glm)
  }
})
res

# save res
foldername="res_"%.%method%.%"_"%.%can.set%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
