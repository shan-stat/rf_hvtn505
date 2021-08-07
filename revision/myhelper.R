# Required packages #
library(parallel); library(foreach); library(doParallel)
library(aucm); library(cvAUC); library(MASS)
library(kyotil); library(dplyr); library(vimp); library(WeightedROC)
library(glmnet); library(tuneRanger); library(caret); library(SuperLearner)

get.kfold.splits <- function (dat, k, seed){
  save.seed <- try(get(".Random.seed", .GlobalEnv), silent = TRUE)
  if (class(save.seed) == "try-error") {
    set.seed(1)
    save.seed <- get(".Random.seed", .GlobalEnv)
  }
  set.seed(seed)
  n0 = nrow(dat$control)
  n1 = nrow(dat$case)
  training.subsets = list()
  test.subsets = list()
  tmp1 = sample(1:n1)
  tmp0 = sample(1:n0)
  splits = list()
  for (ki in 1:k) {
    splits[[ki]] = list(training = list(case = tmp1[(1:n1)%%k != ki - 1], control = tmp0[(1:n0)%%k != ki - 1]), test = list(case = tmp1[(1:n1)%%k == ki - 1], control = tmp0[(1:n0)%%k == ki - 1]))
  }
  splits
  assign(".Random.seed", save.seed, .GlobalEnv)
  splits
}

random.strata.splits <- function (dat, strata, seed){
  set.seed(seed)
  case.idx <- rownames(dat$case)
  case.val.idx <- sample(case.idx, size=5, replace=FALSE) # 1/5 cases (total cases are 25)
  case.val.freq <- table(strata[floor(as.numeric(case.val.idx))])
  case.total.freq <- table(strata[floor(as.numeric(case.idx))])
  ratio <- c()
  for(i in 1:length(case.val.freq)){
    freq.temp <- case.val.freq[i]
    ratio[i] <- freq.temp/case.total.freq[names(freq.temp)]
  }
  names(ratio) <- names(case.val.freq)
  
  control.idx <- floor(as.numeric(rownames(dat$control)))
  control.val.idx <- c()
  for(i in 1:length(ratio)){
    stratum.idx <- which(strata==names(ratio[i]))
    control.stratum <- c()
    for(j in 1:nrow(dat$control)){
      control.stratum[j] <- sum(control.idx[j] == stratum.idx)!=0
    }
    control.temp <- sample(rownames(dat$control)[control.stratum], size=round(sum(control.stratum)*ratio[i]), replace=FALSE)
    control.val.idx <- c(control.val.idx, control.temp)
  }
  
  splits=list()
  case.train.idx <- setdiff(case.idx, case.val.idx)
  control.train.idx <- setdiff(rownames(dat$control), control.val.idx)
  splits[[1]] <- list(training=list(case=case.train.idx, control=control.train.idx),
                      test=list(case=case.val.idx, control=control.val.idx))
  splits
}


# Lasso variable screening #
screen_lasso <- function(Y, X, family, obsWeights=rep(1, nrow(X)), alpha = 1) {
  set.seed(123)
  res.ls <- cv.glmnet( x = as.matrix(X), y =  as.matrix(Y), weights = obsWeights, family = "binomial" , type.measure = 'auc', nfolds = 5, alpha = alpha ) # Lasso penalty
  vars.ls <- (coef( res.ls, s = res.ls$lambda.min ) != 0)[-1]
  vars <- vars.ls
  names(vars) <- colnames(X)
  # always keep clinical covariates (age, BMI, behavior risk score)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# Importing screened data and screened variable index #
screen.dat.index <- function(Y, X, fit.set = c('no','antibody','tcell','all'),
                             screen.dat.method = c('no','lasso'), obsWeights){
  
  p=ncol(X)-3
  
  # candidate set #
  if( fit.set == 'no' ){
    # 1. no markers (only clinical covariates : age, BMI, bhrisk) #
    var_set_none <- c( rep(TRUE, 3), rep(F, p) )
    dat <- cbind(Y = Y, X[,var_set_none])
  } else if( fit.set == 'antibody' ){
    # 2. antibody markers (IgG + IgA + IgG3 + phago + fcrR2a + fcrR3a) #
    var_set_igg_iga_igg3_fxab <- get_nms_group_all_antigens(X, assays = c("IgG", "IgA", "IgG3", "phago", "fcrR2a", "fcrR3a"))
    var_set_igg_iga_igg3_fxab[1:3] =T
    dat <- cbind(Y = Y, X[,var_set_igg_iga_igg3_fxab])
  } else if( fit.set == 'tcell' ){
    # 3. T cell markers (CD4 and CD8) #
    var_set_tcells <- get_nms_group_all_antigens(X, assays = c("CD4", "CD8"))
    var_set_tcells[1:3] =T
    dat <- cbind(Y = Y, X[,var_set_tcells])
  } else if( fit.set == 'all' ){
    # 4. all markers #
    var_set_all <- rep(TRUE, 3+p)
    dat <- cbind(Y = Y, X[,var_set_all])
  }
  
  # Screened data #
  if( screen.dat.method == 'no' ){
    # no screening #
    screen.dat.var <- rep(TRUE, (ncol(dat)-1))
  } else if( screen.dat.method == 'lasso' ){
    # lasso screening #
    screen.dat.var <- screen_lasso( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } 
  
  dat.X <- list( case = dat[dat$Y == 1, c(FALSE, screen.dat.var)], control = dat[dat$Y == 0, c(FALSE, screen.dat.var)])
  return(list(screen.index.var = screen.dat.var, dat = dat.X))
}

# Variable names #
get_nms_group_all_antigens <- function(X, assays, assays_to_exclude = "") {
  # set all variables to be false
  vars <- rep(FALSE, ncol(X))
  # set variables with assay in name to be true
  # may be more than one
  for (i in 1:length(assays)) {
    if (assays_to_exclude != "") {
      vars[grepl(assays[i], names(X)) & !grepl(assays_to_exclude, names(X))] <- TRUE
    } else {
      vars[grepl(assays[i], names(X))] <- TRUE
      if (assays[i] == "phago") {
        vars[grepl("ADCP1", names(X))] <- TRUE
      }
    }
  }
  return(vars)
}

# Random forests #
get.rf.cvauc = function(dat, obsWeights, method=c('RF','RF_under','RF_over','tRF'), ipw, seed=1){
  
  # 5-fold CV #
  splits <- get.kfold.splits(dat, k=5, seed)
  
  # Standard random forest (RF) #
  if( method == 'RF' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[floor(as.numeric(rownames(dat.train)))] ; weights.test <- obsWeights[floor(as.numeric(rownames(dat.test)))]
      
      set.seed( 123 )
      if(ipw == TRUE){
        fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE, min.node.size = 1 ) # weights
      } else if(ipw == FALSE) {
        fit.rf <- ranger( factor(Y)~., data = dat.train, probability = TRUE, min.node.size = 1 ) # no weights
      }
      pred.rf <- predict( fit.rf, data = dat.test )
      WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=weights.test))
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # RF with under-sampling (RF_under) #
  if( method == 'RF_under' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      weights.sampling <- rep(NA, nrow(dat.train))
      weights.sampling[dat.train$Y == 1] <- 5 ; weights.sampling[dat.train$Y == 0] <- 1
      
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, probability = TRUE, replace = TRUE, case.weights = weights.sampling, 
                        sample.fraction = 40/120, keep.inbag = TRUE, min.node.size = 1 )
      pred.rf <- predict( fit.rf, data = dat.test)
      WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=weights.test))
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # RF with over-sampling (RF_over) #
  if( method == 'RF_over' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed( 123 )
      idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE ) # Oversampling minority class
      dat.temp <- list( case = subset(dat.train, Y==1, select=-Y), control = subset(dat.train, Y==0, select=-Y) )
      dat.train.mod <- rbind( data.frame(Y=1, dat.temp$case[idx1,,drop=F]),   data.frame(Y=0, dat.temp$control) )
      weights.train <- obsWeights[floor(as.numeric(rownames(dat.train.mod)))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train.mod, probability = TRUE, replace = TRUE, case.weights = weights.train, min.node.size = 1 )
      pred.rf <- predict( fit.rf, data = dat.test)
      WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=weights.test))
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # Tuned random forest (tRF) #
  if( method == 'tRF' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train$Y <- as.factor(dat.train$Y) ; dat.test$Y <- as.factor(dat.test$Y)
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      set.seed( 123 )
      rf.task <- makeClassifTask(data = dat.train, target = 'Y')
      res.tunerf <- tuneRanger( rf.task, measure = list(auc), num.trees = 500, iters.warmup = 50, iters = 100, save.file.path = NULL, 
                                tune.parameters = c("mtry", "min.node.size","sample.fraction"), parameters = list(replace = TRUE) )
      pred.tunerf <- predict( res.tunerf$model, newdata = dat.test )
      WeightedAUC(WeightedROC(guess=pred.tunerf$data$prob.1, label=dat.test$Y, weight=weights.test))
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  cv.aucs
}

# Generalized linear models (GLM) #
get.glm.cvauc = function(dat, obsWeights, ipw, seed=1){
  
  # 5-fold CV #
  splits <- get.kfold.splits(dat, k=5, seed)
  
  # glm #
  cv.aucs <-  mclapply( splits, function(split){
    dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
    dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
    weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
    
    set.seed(123)
    if(ipw == TRUE){
      fit.mlr <- glm( factor(Y)~., dat.train, family=binomial, weights=weights.train ) # weights
    } else {
      fit.mlr <- glm( factor(Y)~., dat.train, family=binomial ) # no weights
    }
    pred.mlr <- predict( fit.mlr, newdata=dat.test, type='response' )
    WeightedAUC(WeightedROC(guess=pred.mlr, label=dat.test$Y, weight=weights.test))
  }, mc.cores = 4 )
  cv.aucs <- unlist( cv.aucs )
  
  cv.aucs
}

# Stacking #
get.st.cvauc = function(dat, obsWeights, strata, var.index, method, seed=1){
  
  # Outer layer: random stratified cv #
  splits <- random.strata.splits(dat, strata, seed)
  
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
    dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
    dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
    # Y=1,0 must be converted to factor that has two levels of 'case' (for Y=1) and 'control' (for Y=0) #
    dat.train$Y <- as.factor(dat.train$Y) ; levels(dat.train$Y) <- c('control','case')
    weights.train <- obsWeights[floor(as.numeric(rownames(dat.train)))] ; weights.test <- obsWeights[floor(as.numeric(rownames(dat.test)))]
    
    # Setting three RF hyperparameters to their default values (without the following code, hyperparameter tuning is automatically done) #
    rf_grid <- expand.grid(mtry = floor(sqrt(sum(var.index$rf)-1)), splitrule = 'gini', min.node.size = 1)
    set.seed( 123 )
    model_list <- caretList(
      list(Y~., data=dat.train[,var.index$glm], weights = weights.train), # glm
      list(Y~., data=dat.train[,var.index$rf], tuneGrid = rf_grid, weights = weights.train), # RF
      trControl=my_control,
      methodList=c('glm', 'ranger')
    )
    
    # For RF learners, final model should be set manually because of a random seed issue from caret.train() (for glm, it is not happned) #
    set.seed(123)
    fit.rf <- ranger( factor(Y)~., data = dat.train[,var.index$rf], probability = TRUE, min.node.size = 1, case.weights=weights.train ) # Weights
    model_list$ranger$finalModel <- fit.rf
    
    # stacking #
    set.seed( 123 )
    pred.cv <- makePredObsMatrix(model_list)
    res.st.fit <- method$computeCoef(Z = pred.cv$preds, Y =(as.numeric(pred.cv$obs)-1), obsWeights = weights.train, libraryNames = names(var.index), verbose = FALSE, control=list(trimLogit=0.001))
    pred.test <- predict.caretList(model_list, newdata=dat.test)
    pred.st <- method$computePred(predY = pred.test, coef = res.st.fit$coef, control=list(trimLogit=0.001))
  
    # Calculating CV-AUC #
    est.cvauc <- WeightedAUC(WeightedROC(guess=pred.st, label=dat.test$Y, weight=weights.test))
    # Calculating Pearson correlation coefficient #
    est.corr <- cor(pred.cv$preds[,1], pred.cv$preds[,2])
    c(est.cvauc=est.cvauc, est.corr=est.corr)
  }, mc.cores = 4 )
  cv.aucs
}
