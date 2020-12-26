##### Import helper functions #####
setwd("set your local directory where helper_rf_hvtn.R file is located")
source('helper_rf_hvtn.R')



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

# Total data (vaccine arm + placebo arm = 189 observations) #
X_markers <- dat.505 %>% select( var.super$varname, paste0(var.super$varname, '_bin') )
X_covariates <- dat.505 %>% select( age, BMI, bhvrisk )
X <- data.frame( trt = dat.505$trt, X_covariates, X_markers )
X_full_none <- data.frame( trt = dat.505$trt, X_covariates )
weights <- dat.505$wt
Y <- dat.505$case

# Vaccine arm (150 observations) #
vaccinees <- cbind.data.frame( Y, weights, X ) %>% filter( trt == 1 ) %>% select( -trt )
X_none <- cbind.data.frame( Y, weights, X_full_none ) %>% filter( trt == 1 ) %>% select( -trt, -Y, -weights )
Y_vaccine <- vaccinees$Y
weights_vaccine <- vaccinees$weights
X_vaccine <- vaccinees %>% select( -Y, -weights )



##### 2. Experiments ######
##### 2-1. Random forest (RF) #####
# Input values #
can.set <- c('no','antibody','tcell','all') # candidate set: (1) no markers, (2) antibody markers, (3) t cell markers, (4) all markers
screening <- 'no'; screening <- 'lasso' # no: no screening, lasso: lasso screening
covariates <- 'no'; covariates <- 'cvr' # no: without covariates, cvr: with coveriate

# Empty matrix for containing results #
pred.mat <- matrix(NA, nrow=10, ncol=4 ); colnames(pred.mat) <- can.set

for(j in 1:length(can.set)){
  print(j)
  # Generating data #
  res.screen <- screen.dat.index(Y=Y_vaccine, X=X_vaccine, X_markers=X_markers, obsWeights=weights_vaccine,
                                 fit.set=can.set[j], screen.dat.method=screening, screen.index.method=screening)
  if(covariates == 'no'){
    # Without coveriate #
    dat.X <- res.screen$dat
    dat.X <- list(case=select(dat.X$case, select=-c('age','BMI','bhvrisk')), control=select(dat.X$control, select=-c('age','BMI','bhvrisk'))) # no covariates 
  } else if(covariates == 'cvr'){
    # With coveriate #
    dat.X <- res.screen$dat
  }
  
  # Ten replicates #
  for( i in 1:10 ){
    print(i)
    seed <- i
    # RF #
    res.temp <- get.rf.cvauc(dat=dat.X, obsWeights=weights_vaccine, method='RF', ipw=TRUE, seed=seed)
    pred.mat[i, can.set[j]] <- mean( res.temp )
  }
}
pred.mat



##### 2-2. Generalized linear models (GLM) #####
# Input values #
can.set <- c('no','antibody','tcell','all') # candidate set: (1) no markers, (2) antibody markers, (3) t cell markers, (4) all markers
screening <- 'no'; screening <- 'lasso' # no: no screening, lasso: lasso screening
covariates <- 'no'; covariates <- 'cvr' # no: without covariates, cvr: with coveriate

# Empty matrix for containing results #
pred.mat <- matrix(NA, nrow=10, ncol=4 ); colnames(pred.mat) <- can.set

for(j in 1:length(can.set)){
  print(j)
  # Generating data #
  res.screen <- screen.dat.index(Y=Y_vaccine, X=X_vaccine, X_markers=X_markers, obsWeights=weights_vaccine,
                                 fit.set=can.set[j], screen.dat.method=screening, screen.index.method=screening)
  if(covariates == 'no'){
    # Without coveriate #
    dat.X <- res.screen$dat
    dat.X <- list(case=select(dat.X$case, select=-c('age','BMI','bhvrisk')), control=select(dat.X$control, select=-c('age','BMI','bhvrisk'))) # no covariates 
  } else if(covariates == 'cvr'){
    # With coveriate #
    dat.X <- res.screen$dat
  }
  
  # Ten replicates #
  for( i in 1:10 ){
    print(i)
    seed <- i
    # GLM #
    res.temp <- get.glm.cvauc(dat=dat.X, obsWeights=weights_vaccine, ipw=TRUE, seed=seed)
    pred.mat[i, can.set[j]] <- mean( res.temp )
  }
}
pred.mat



##### 2-3. Stacking (RF + GLM) #####
# The following four files are required for fitting stacking #
#setwd("set your local directory where five files are located")
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')

# Generating data where stacking is fitted #
dat.X <- screen.dat.index(Y=Y_vaccine, X=X_vaccine, X_markers=X_markers, obsWeights=weights_vaccine,
                          fit.set='all', screen.dat.method='no', screen.index.method='no')$dat

# Specifying screened variables index #
# For all markers #
screen.var.all <- screen.dat.index(Y=Y_vaccine, X=X_vaccine, X_markers=X_markers, obsWeights=weights_vaccine,
                                   fit.set='all', screen.dat.method='lasso', screen.index.method='lasso')$screen.index.var

# For t cell markers #
screen.var.tcell.temp <- screen.dat.index(Y=Y_vaccine, X=X_vaccine, X_markers=X_markers, obsWeights=weights_vaccine,
                                          fit.set='tcell', screen.dat.method='lasso', screen.index.method='lasso')$screen.index.var
screen.var.tcell <- rep(FALSE,ncol(X_vaccine))
for(i in 1:length(which(screen.var.tcell.temp))){
  temp <- which(names(X_vaccine) == names(which(screen.var.tcell.temp)[i]))
  screen.var.tcell[temp] <- TRUE
}
screen.var.tcell

# For antibody markers #
screen.var.antibody.temp <- screen.dat.index(Y=Y_vaccine, X=X_vaccine, X_markers=X_markers, obsWeights=weights_vaccine,
                                             fit.set='antibody', screen.dat.method='lasso', screen.index.method='lasso')$screen.index.var
screen.var.antibody <- rep(FALSE,ncol(X_vaccine))
for(i in 1:length(which(screen.var.antibody.temp))){
  temp <- which(names(X_vaccine) == names(which(screen.var.antibody.temp)[i]))
  screen.var.antibody[temp] <- TRUE
}
screen.var.antibody

# Empty matrix for containing results #
stacking.set <- c('RF:tcell+glm:antibody','RF:tcell+glm:all','RF:all+glm:antibody','RF:all+glm:all')
pred.mat <- matrix(NA, nrow=10, ncol=4 ); colnames(pred.mat) <- stacking.set
corr.vec <- rep(NA,4); names(corr.vec) <- stacking.set
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
  
  # Ten replicates #
  for( i in 1:10 ){
    print(i)
    seed <- i
    # Stacking (CV-AUC) #
    res.temp <- get.st.cvauc(dat=dat.X, obsWeights=weights_vaccine, var.index=var.index.set, method='method.NNloglik', seed=seed)
    pred.mat[i, stacking.set[j]] <- mean(do.call(rbind,res.temp)[,'est.cvauc'])
    
    # Stacking (Pearson correlation coeffieict; for specific random seed 1) #
    if(i==1){ corr.vec[j] <- mean(do.call(rbind,res.temp)[,'est.corr']) }
  }
}
pred.mat
corr.vec
