source("myhelper.R")
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')



##### Bootstrap method 2 #####
# Description #
# (1) bootstrap is for the variability of CV-AUC estimates
# (2) 5-fold cv is performed ahead. And then, bootstrap is performed to the training subset and validation subset 


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


### 1. Specifying screened variables index ###
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


### 2. Setting process arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are need to a batch script
  # can.set.st: one of c('RF:tcell+glm:antibody','RF:tcell+glm:all','RF:all+glm:antibody','RF:tcell+glm:all')
  # can.set.rf: one of c('all','tcell','antibody','no')
  # proj: ST or RF
  # outer.feed: seed for 5-fold cv
  # fold: if fold=1, it is the first stage of the 5-fold cv
  Args=c(batch.size="10",batch.number="1",can.set.st='RF:tcell+glm:antibody',can.set.rf='all',proj="ST",outer.seed=1,fold=1)
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
i=i+1; outer.seed=as.numeric(Args[i])
i=i+1; fold=as.numeric(Args[i])


### 3. Experiments ###
if(proj=='ST'){
  
  # specify candidate set #
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
  
  # Outer layer: 5-fold CV #
  splits <- get.kfold.splits(Y=Y_vaccine, X=X_vaccine, k=5, seed=outer.seed)
  split <- splits[[fold]]
  
  # Inner layer: 10-fold CV (fitting learners and obtaining out-of-sample predictions) #
  my_control <- trainControl(
    method="cv",
    number=10,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=twoClassSummary
  )
  
  # Method of combining predictions from different learners #
  method='method.NNloglik'
  method <- get(method, mode = 'function')()
  if(!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), character.only = TRUE))
  }
  
  # Fitting candidate learners #
  # Target variable must be named as 'Y' #
  dat.train <- cbind(Y=Y_vaccine[c(split$training$case,split$training$control)], X_vaccine[c(split$training$case,split$training$control),,drop=F])
  dat.test <- cbind(Y=Y_vaccine[c(split$test$case,split$test$control)], X_vaccine[c(split$test$case,split$test$control),,drop=F])
  # Y=1,0 must be converted to factor that has two levels of 'case' (for Y=1) and 'control' (for Y=0) #
  dat.train$Y <- as.factor(dat.train$Y) ; levels(dat.train$Y) <- c('control','case')
  weights.train <- weights_vaccine[c(split$training$case,split$training$control)]
  weights.test <- weights_vaccine[c(split$test$case,split$test$control)]
  
  ### 2. Experiments ###
  res=sapply(seeds, simplify="array", function (seed) {
    
    myprint(seed)
    
    # Bootstrapping #
    set.seed(seed)
    case.train.idx <- sample(which(dat.train$Y=='case'),replace=TRUE); control.train.idx <- sample(which(dat.train$Y=='control'),replace=TRUE)
    case.test.idx <- sample(which(dat.test$Y==1),replace=TRUE); control.test.idx <- sample(which(dat.test$Y==0),replace=TRUE)
    dat.train.boot <- dat.train[c(case.train.idx,control.train.idx),,drop=F]
    dat.test.boot <- dat.test[c(case.test.idx,control.test.idx),,drop=F]
    weights.train.boot <- weights.train[c(case.train.idx,control.train.idx)]; weights.test.boot <- weights.test[c(case.test.idx,control.test.idx)]
    
    # Setting three RF hyperparameters to their default values (without the following code, hyperparameter tuning is automatically done) #
    rf_grid <- expand.grid(mtry = floor(sqrt(sum(var.index.set$rf)-1)), splitrule = 'gini', min.node.size = 1)
    set.seed( 123 )
    suppressWarnings(
      model_list <- caretList(
        list(Y~., data=dat.train.boot[,var.index.set$glm], weights = weights.train.boot), # glm
        list(Y~., data=dat.train.boot[,var.index.set$rf], tuneGrid = rf_grid, weights = weights.train.boot), # RF
        trControl=my_control,
        methodList=c('glm', 'ranger')
      )
    )
    
    # For RF learners, final model should be set manually because of a random seed issue from caret.train() (for glm, it is not happned) #
    set.seed(123)
    fit.rf <- ranger( factor(Y)~., data = dat.train.boot[,var.index.set$rf], probability = TRUE, min.node.size = 1, case.weights=weights.train.boot ) # Weights
    model_list$ranger$finalModel <- fit.rf
    
    # stacking #
    set.seed( 123 )
    pred.cv <- makePredObsMatrix(model_list)
    res.st.fit <- method$computeCoef(Z = pred.cv$preds, Y =(as.numeric(pred.cv$obs)-1), obsWeights = weights.train.boot, libraryNames = names(var.index.set), verbose = FALSE, control=list(trimLogit=0.001))
    pred.test <- predict.caretList(model_list, newdata=dat.test.boot)
    pred.st <- method$computePred(predY = pred.test, coef = res.st.fit$coef, control=list(trimLogit=0.001))
    
    # Calculating CV-AUC #
    est.cvauc <- WeightedAUC(WeightedROC(guess=pred.st, label=dat.test.boot$Y, weight=weights.test.boot))
    est.cvauc
  })
  res
  
} else if(proj=='RF'){
  
  # specify candidate set #
  if(can.set.rf=='all'){
    var.index.set <- screen.var.all
  } else if(can.set.rf=='tcell'){
    var.index.set <- screen.var.tcell
  } else if(can.set.rf=='antibody'){
    var.index.set <- screen.var.antibody
  }
  
  # Outer layer: 5-fold CV #
  splits <- get.kfold.splits(Y=Y_vaccine, X=X_vaccine, k=5, seed=outer.seed)
  
  # Fitting RF #
  split <- splits[[fold]]
  dat.train <- cbind(Y=Y_vaccine[c(split$training$case,split$training$control)], X_vaccine[c(split$training$case,split$training$control),,drop=F])
  dat.test <- cbind(Y=Y_vaccine[c(split$test$case,split$test$control)], X_vaccine[c(split$test$case,split$test$control),,drop=F])
  weights.train <- weights_vaccine[c(split$training$case,split$training$control)]
  weights.test <- weights_vaccine[c(split$test$case,split$test$control)]
  
  ### 2. Experiments ###
  res=sapply(seeds, simplify="array", function (seed) {
    
    myprint(seed)
    
    # Bootstrapping #
    set.seed(seed)
    case.train.idx <- sample(which(dat.train$Y==1),replace=TRUE); control.train.idx <- sample(which(dat.train$Y==0),replace=TRUE)
    case.test.idx <- sample(which(dat.test$Y==1),replace=TRUE); control.test.idx <- sample(which(dat.test$Y==0),replace=TRUE)
    dat.train.boot <- dat.train[c(case.train.idx,control.train.idx),,drop=F]
    dat.test.boot <- dat.test[c(case.test.idx,control.test.idx),,drop=F]
    weights.train.boot <- weights.train[c(case.train.idx,control.train.idx)]; weights.test.boot <- weights.test[c(case.test.idx,control.test.idx)]
    
    set.seed( 123 )
    fit.rf <- ranger( factor(Y)~., data = dat.train.boot[,c(TRUE,var.index.set)], case.weights = weights.train.boot, probability = TRUE, min.node.size = 1 ) # weights
    pred.rf <- predict( fit.rf, data = dat.test.boot )
    
    # Calculating CV-AUC #
    est.cvauc <- WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test.boot$Y, weight=weights.test.boot))
    est.cvauc
  })
  res
}

# save res
foldername="res_"%.%can.set%.%"_"%.%outer.seed%.%"_"%.%fold%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
