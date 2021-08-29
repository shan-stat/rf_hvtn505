##### Import helper functions #####
setwd("set your local directory where myhelper.R file is located")
source('myhelper.R')



##### 1. Import HVTN505 dataset ######
# 1-1) download "HVTN505_2019-4-25.tar.gz" file at https://atlas.scharp.org/cpas/project/HVTN%20Public%20Data/HVTN%20505/begin.view
# 1-2) install R package HVTN505 using "devtools::install_local("HVTN505_2019-4-25.tar.gz")"
data( 'dat.505', package = 'HVTN505' )
suppressWarnings( data( 'var.super', package = 'HVTN505' ) )

# Scaling each variables: mean 0 and sd 1 #
for( i in var.super$varname ){
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
  dat.505[[i %.% '_bin']] <- scale( dat.505[[i %.% '_bin']], center = mean(dat.505[[i %.% '_bin']][dat.505$trt == 1]), scale = sd(dat.505[[i %.% '_bin']][dat.505$trt == 1]))
}
for( i in c( 'age', 'BMI', 'bhvrisk' ) ) {
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
}

dat.505 <- subset(dat.505, trt==1) # only treatment
Y_vaccine <- dat.505$case
X_vaccine <- dat.505 %>% select(age, BMI, bhvrisk, var.super$varname, paste0(var.super$varname, '_bin'))
weights_vaccine <- dat.505$wt



##### 2. Experiments ######
##### 2-1. Random forest (RF) #####
can.set <- c('no','antibody','tcell','all') # (1) no markers, (2) antibody markers, (3) t cell markers, (4) all markers
screening <- 'no'; screening <- 'lasso' # no: no screening, lasso: lasso screening
covariates <- 'no'; covariates <- 'cvr' # no: without covariates, cvr: with covariates

# Empty matrix for containing results #
pred.mat <- matrix(NA, nrow=100, ncol=4 ); colnames(pred.mat) <- can.set

for(j in 1:length(can.set)){
  print(j)
  # Generating data #
  res.screen <- screen.dat.varidx(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set=can.set[j], screen.method=screening)
  if(covariates == 'no'){
    # Without covariates #
    dat <- res.screen$dat
    dat <- select(dat, select=-c('age','BMI','bhvrisk'))
  } else if(covariates == 'cvr'){
    # With covariates #
    dat <- res.screen$dat
  }
  
  if(can.set[j]=='no' & covariates=='no'){
    pred.mat[,can.set[j]] <- 0.5 # CV-AUC is 0.5 in no markers and no covariates #
  } else {
    # One-hundred replicates #
    for( i in 1:100 ){
      print(i)
      seed <- i
      # RF #
      res.temp <- get.rf.cvauc(Y=dat[,'Y'], X=select(dat,-"Y"), obsWeights=weights_vaccine, method='RF', ipw=FALSE, seed=seed)
      pred.mat[i, can.set[j]] <- mean( res.temp )
    }
  }
}
pred.mat



##### 2-2. Generalized linear models (GLM) #####
can.set <- c('no','antibody','tcell','all') # (1) no markers, (2) antibody markers, (3) t cell markers, (4) all markers
screening <- 'no'; screening <- 'lasso' # no: no screening, lasso: lasso screening
covariates <- 'no'; covariates <- 'cvr' # no: without covariates, cvr: with covariates

# Empty matrix for containing results #
pred.mat <- matrix(NA, nrow=100, ncol=4 ); colnames(pred.mat) <- can.set

for(j in 1:length(can.set)){
  print(j)
  # Generating data #
  res.screen <- screen.dat.varidx(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set=can.set[j], screen.method=screening)
  if(covariates == 'no'){
    # Without covariates #
    dat <- res.screen$dat
    dat <- select(dat, select=-c('age','BMI','bhvrisk'))
  } else if(covariates == 'cvr'){
    # With covariates #
    dat <- res.screen$dat
  }
  
  # One-hundred replicates #
  for( i in 1:100 ){
    print(i)
    seed <- i
    # GLM #
    res.temp <- get.glm.cvauc(Y=dat[,'Y'], X=select(dat,-"Y"), obsWeights=weights_vaccine, ipw=FALSE, seed=seed)
    pred.mat[i, can.set[j]] <- mean( res.temp )
  }
}
pred.mat



##### 2-3. Stacking (RF + GLM) #####
# The following four files are required for fitting stacking #
setwd("set your local directory where four files are located")
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')

# Specifying screened variables index #
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


### CV-AUC, Pearson correlation, Coefficients ###
stacking.set <- c('RF:tcell+glm:antibody','RF:tcell+glm:all','RF:all+glm:antibody','RF:all+glm:all')
pred.mat <- matrix(NA, nrow=100, ncol=4 ); colnames(pred.mat) <- stacking.set
corr.mat <- matrix(NA, nrow=100, ncol=4 ); colnames(corr.mat) <- stacking.set
coef.mat <- matrix(NA, nrow=100, ncol=2); colnames(coef.mat) <- c('GLM','RF')
for(j in 1:length(stacking.set)){
  print(j)
  if(j==1){
    # RF:tcell + glm:antibody #
    var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.tcell))
  } else if(j==2){
    # RF:tcell + glm:all #
    var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.tcell))
  } else if(j==3){
    # RF:all + glm:antibody #
    var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.all))
  } else if(j==4){
    # RF:all + glm:all #
    var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.all))
  }
  
  # One-hundred replicates #
  for( i in 1:100 ){
    print(i)
    seed <- i
    # Stacking (CV-AUC) #
    res.temp <- get.st.cvauc(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, var.index=var.index.set, method='method.NNloglik', seed=seed)
    pred.mat[i, stacking.set[j]] <- mean(do.call(rbind,res.temp)[,'est.cvauc']) # cv-auc
    corr.mat[i, stacking.set[j]] <- mean(do.call(rbind,res.temp)[,'est.corr']) # correlation
    if(j==1) {coef.mat[i,] <- apply(do.call(rbind,res.temp)[,c('coef1','coef2')], 2, mean)} # coefficients for meta-learner
  }
}
pred.mat; corr.mat; coef.mat


### Plotting prediction scores ###
var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.tcell)) # RF:tcell + glm:antibody
res.temp <- get.st.cvauc(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, var.index=var.index.set, method='method.NNloglik', seed=1, detail=TRUE) # only seed=1
res <- do.call(rbind,res.temp); res <- as.data.frame(res)

# Box plot #
myfigure(mfrow=c(1,3))
myboxplot(ranger~case, res, main="RF", names=c("controls","cases"), col=1+res$case); title(sub="CV-AUC = "%.%formatDouble(fast.auc(res[,"ranger"], res[,"case"]), 3), line=3)
myboxplot(glm~case, res, main="GLM", names=c("controls","cases"), col=1+res$case); title(sub="CV-AUC = "%.%formatDouble(fast.auc(res[,"glm"], res[,"case"]), 3), line=3)
myboxplot(stacking~case, res, main="Stacking", names=c("controls","cases"), col=1+res$case); title(sub="CV-AUC = "%.%formatDouble(fast.auc(res[,"stacking"], res[,"case"]), 3), line=3)

# Scatter plot #
myfigure(mfrow=c(1,2))
plot(glm~ranger, res, col=res[,"case"]+1, xlim=c(0,1), ylim=c(0,1), cex=.7, pch=ifelse(res$ptid %in% c(180,183), 2, 1))
abline(v=c(.45,.58), lty=2); abline(0,1)
plot(stacking~ranger, res, col=res[,"case"]+1, xlim=c(0,1), ylim=c(0,1), cex=.7, pch=ifelse(res$ptid %in% c(180,183), 2, 1))
abline(0,1); abline(v=c(.45,.58), lty=2); abline(h=c(.4), lty=2)



##### 2-4. Supplementary Materials  #####
can.set <- 'all'
screening <- c('no','lasso')[2]  # no: no screening; lasso: lasso screening
proj <- c('auc_boot','auc_train','plot_tree','plot_splitvar')[4] # (1) auc_boot: auc on bootstrapped data; (2) auc_train: auc on training data; (3) plot_tree: tree plot; (4) plot_splitvar: bar plot for split variables

# Generating data #
res.screen <- screen.dat.varidx(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set=can.set, screen.method=screening)
dat <- res.screen$dat

# Experiments #
if(startsWith(proj,'auc')){
  res.mat <- c()
  for(i in 1:1){
    print(i)
    seed <- i
    Y=dat[,'Y']; X=select(dat,-"Y"); obsWeights=weights_vaccine; method='RF'; ipw=FALSE; seed=seed # standard RF
    
    splits <- get.kfold.splits(Y=Y, X=X, k=5, seed=seed) # 5-fold CV
    res.fold <- c()
    for(k in 1:5){
      split <- splits[[k]]
      # Standard random forest (RF) #
      dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
      dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
      weights.train <- obsWeights[c(split$training$case,split$training$control)]; weights.test <- obsWeights[c(split$test$case,split$test$control)]
      
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, probability = TRUE, min.node.size = 1, keep.inbag=TRUE ) # no weights
      
      if(proj=='auc_boot'){
        # 1. bootstrapped AUC for individual trees #
        res <- c()
        for(l in 1:fit.rf$num.trees){
          inbag.info <- fit.rf$inbag.counts[[l]]
          inbag.idx <- c()
          for(a in 1:nrow(dat.train)){
            if(inbag.info[a]!=0) inbag.idx <- c(inbag.idx, rep(a,inbag.info[a]))
          }
          dat.boot <- dat.train[inbag.idx,] # bootstrapped data
          temp <- predict(fit.rf, data=dat.boot, predict.all=TRUE) 
          res[l] <- WeightedAUC(WeightedROC(guess=temp$predictions[,,l][,'1'], label=dat.boot$Y)) # auc on bootstrapped data
        }
        res.fold[k] <- mean(res)
        
      } else if(proj=='auc_train'){
        # 2. training AUC for individual trees #
        temp <- predict(fit.rf, data=dat.train, predict.all=TRUE)
        res <- c()
        for(l in 1:fit.rf$num.trees){
          res[l] <- WeightedAUC(WeightedROC(guess=temp$predictions[,,l][,'1'], label=dat.train$Y)) # auc on training data
        }
        res.fold[k] <- mean(res)
      }
    }
    res.mat[i] <- mean(res.fold)
  }
  res.mat
  
} else if(startsWith(proj,'plot')){
  seed <- 1
  Y=dat[,'Y']; X=select(dat,-"Y"); obsWeights=weights_vaccine; method='RF'; ipw=FALSE; seed=seed # standard RF
  splits <- get.kfold.splits(Y=Y, X=X, k=5, seed=seed); split <- splits[[1]] # only for the first stage of the 5-fold CV
  
  # Standard random forest (RF) #
  dat.train <- cbind(Y=Y[c(split$training$case,split$training$control)], X[c(split$training$case,split$training$control),,drop=F])
  dat.test <- cbind(Y=Y[c(split$test$case,split$test$control)], X[c(split$test$case,split$test$control),,drop=F])
  weights.train <- obsWeights[c(split$training$case,split$training$control)]; weights.test <- obsWeights[c(split$test$case,split$test$control)]
  
  set.seed( 123 )
  fit.rf <- ranger( factor(Y)~., data = dat.train, probability = TRUE, min.node.size = 1, keep.inbag=TRUE ) # no weights
  
  if(proj=='plot_tree'){
    # 3. a single tree #
    print(treeInfo(fit.rf, tree=250)) # the tree result was used to plot Figure C.2
    
  } else if(proj=='plot_splitvar'){
    # 4. split variables #
    res <- c()
    for(l in 1:fit.rf$num.trees){
      temp <- treeInfo(fit.rf, l)[,'splitvarName']
      res <- c(res, temp[!is.na(temp)])
    }
    
    if(screening=='lasso'){
      num.var <- table(res)
      names(num.var) <- c('z1','z2','z3','x51','x53','x55','x175','x397')
      barplot(num.var, main="RF with screening", angle=90, col='blue', border=NA)
    } else if(screening=='no'){
      num.var <- table(res)
      num.var <- num.var[sort(names(num.var))]
      scr.var <- c("age","bhvrisk","BMI","CD8_ANYVRCENV_IFNg_logpctpos","CD8_ANYVRCENV_IL2_logpctpos","CD8_ANYVRCENV_TNFa_logpctpos","IgG3_VRC_A_avi","IgG3w28_VRC_B_gp140_bin")
      num.var <- num.var[c(scr.var, setdiff(names(num.var), scr.var))]; names(num.var) <- NULL
      barplot(num.var, main="RF without screening", angle=90, col=c(rep('blue',8), rep('grey',398)), border=NA)
    }
  }
}
