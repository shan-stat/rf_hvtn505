# Required packages #
library(parallel); library(foreach); library(doParallel)
library(aucm); library(cvAUC); library(MASS)
library(HVTN505); library(kyotil); library(dplyr); library(vimp); library(WeightedROC)
library(glmnet); library(tuneRanger); library(caret); library(SuperLearner)

# K-fold cross validation #
get.kfold.splits <- function (Y, X, k, seed){
  set.seed(seed)
  n0 = sum(Y==0)
  n1 = sum(Y==1)
  training.subsets = list()
  test.subsets = list()
  tmp1 = sample(which(Y==1), replace=FALSE)
  tmp0 = sample(which(Y==0), replace=FALSE)
  splits = list()
  for (ki in 1:k) {
    splits[[ki]] = list(training = list(case = tmp1[(1:n1)%%k != ki - 1], control = tmp0[(1:n0)%%k != ki - 1]), test = list(case = tmp1[(1:n1)%%k == ki - 1], control = tmp0[(1:n0)%%k == ki - 1]))
  }
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
screen.dat.varidx <- function(Y, X, fit.set = c('no','antibody','tcell','all'), screen.method = c('no','lasso'), obsWeights){
  
  # Pre-processing #
  X_markers <- select(X,-c('age','BMI','bhvrisk')) # without clinical covariates
  
  # candidate set #
  if( fit.set == 'no' ){
    # 1. no markers (only clinical covariates : age, BMI, bhrisk) #
    var_set_none <- rep(FALSE, ncol(X_markers))
    var_set_none <- c( rep(TRUE, 3), var_set_none )
    dat <- cbind(Y = Y, X[,var_set_none])
  } else if( fit.set == 'antibody' ){
    # 2. antibody markers (IgG + IgA + IgG3 + phago + fcrR2a + fcrR3a) #
    var_set_igg_iga_igg3_fxab <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3", "phago", "fcrR2a", "fcrR3a"))
    var_set_igg_iga_igg3_fxab <- c( rep(TRUE, 3), var_set_igg_iga_igg3_fxab )
    dat <- cbind(Y = Y, X[,var_set_igg_iga_igg3_fxab])
  } else if( fit.set == 'tcell' ){
    # 3. T cell markers (CD4 and CD8) #
    var_set_tcells <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8"))
    var_set_tcells <- c( rep(TRUE, 3), var_set_tcells )
    dat <- cbind(Y = Y, X[,var_set_tcells])
  } else if( fit.set == 'all' ){
    # 4. all markers #
    var_set_all <- rep(TRUE, ncol(X_markers))
    var_set_all <- c( rep(TRUE, 3), var_set_all )
    dat <- cbind(Y = Y, X[,var_set_all])
  }
  
  # Screened data #
  if( screen.method == 'no' ){
    # no screening #
    screen.index <- rep(TRUE, (ncol(dat)-1))
  } else if( screen.method == 'lasso' ){
    # lasso screening #
    screen.index <- screen_lasso( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } 
  
  dat.X <- dat[,c(TRUE, screen.index)] # TRUE for Y
  return(list(screen.index.var = screen.index, dat = dat.X))
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
get.rf.cvauc = function(Y, X, obsWeights, method=c('RF','RF_under','RF_over','tRF'), ipw, seed=1){
  
  # 5-fold CV #
  splits <- get.kfold.splits(Y=Y, X=X, k=5, seed=seed)
  
  # Standard random forest (RF) #
  if( method == 'RF' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
      dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
      weights.train <- obsWeights[c(split$training$case,split$training$control)]; weights.test <- obsWeights[c(split$test$case,split$test$control)]
      
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
      dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
      dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
      weights.train <- obsWeights[c(split$training$case,split$training$control)]; weights.test <- obsWeights[c(split$test$case,split$test$control)]
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
      dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
      dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
      weights.train <- obsWeights[c(split$training$case,split$training$control)]; weights.test <- obsWeights[c(split$test$case,split$test$control)]
      set.seed( 123 )
      idx1 <- sample( which(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE ) # Oversampling minority class
      dat.train.mod <- dat.train[c(idx1,which(dat.train$Y==0)),,drop=F]
      weights.train <- weights.train[c(idx1,which(dat.train$Y==0))]
      
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
      dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
      dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
      dat.train$Y <- as.factor(dat.train$Y) ; dat.test$Y <- as.factor(dat.test$Y)
      weights.train <- obsWeights[c(split$training$case,split$training$control)]; weights.test <- obsWeights[c(split$test$case,split$test$control)]
      
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
get.glm.cvauc = function(Y, X, obsWeights, ipw, seed=1){
  
  # 5-fold CV #
  splits <- get.kfold.splits(Y=Y, X=X, k=5, seed=seed)
  
  # glm #
  cv.aucs <-  mclapply( splits, function(split){
    dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
    dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
    weights.train <- obsWeights[c(split$training$case,split$training$control)]; weights.test <- obsWeights[c(split$test$case,split$test$control)]
    
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
get.st.cvauc = function(Y, X, obsWeights, var.index, method, seed=1, mc.cores=4, detail=FALSE){
  
  # Outer layer: 5-fold CV #
  splits <- get.kfold.splits(Y=Y, X=X, k=5, seed=seed)
  
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
    
    # Setting three RF hyperparameters to their default values (without the following code, hyperparameter tuning is automatically done) #
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
    
    # For RF learners, final model should be set manually because of a random seed issue from caret.train() (for glm, it is not happened) #
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
    st.cvauc <- WeightedAUC(WeightedROC(guess=pred.st, label=dat.test$Y, weight=weights.test))
    st.corr <- cor(pred.cv$preds[,1], pred.cv$preds[,2]) # Pearson correlation coefficient
    glm.cvauc <- WeightedAUC(WeightedROC(guess=pred.test[,'glm'], label=dat.test$Y, weight=weights.test))
    rf.cvauc <- WeightedAUC(WeightedROC(guess=pred.test[,'ranger'], label=dat.test$Y, weight=weights.test))
    
    if(detail==TRUE){
      res <- cbind(as.numeric(rownames(dat.test)), dat.test$Y, pred.test, pred.st); colnames(res) <- c('ptid','case','glm','ranger','stacking')
      res
    } else{
      c(est.cvauc=st.cvauc, est.corr=st.corr, coef=res.st.fit$coef)
    }
  }, mc.cores = mc.cores )
  cv.aucs
}
