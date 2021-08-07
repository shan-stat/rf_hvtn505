source("myhelper.R")
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')

#data( 'dat.505', package = 'HVTN505' )
#suppressWarnings( data( 'var.super', package = 'HVTN505' ) )
dat.505 <- read.table( 'dat_505.txt', header = T, sep = ',' )

dat.505=subset(dat.505, trt==1)

var.super <- read.table( 'var_super.txt', header = T, sep = ',' )

# Scaling each variables: mean 0 and sd 1 #
dat.505[,c('age', 'BMI', 'bhvrisk', var.super$varname)]=scale(dat.505[,c('age', 'BMI', 'bhvrisk', var.super$varname)])

Y_vaccine <- dat.505$case
X_vaccine <- dat.505 %>% select(age, BMI, bhvrisk, var.super$varname)
weights_vaccine <- vaccinees$weights
strata_vaccine <- vaccinees$strata


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
# For all markers #
screen.var.all <- screen.dat.index(Y=Y_vaccine, X=X_vaccine, obsWeights=weights_vaccine, fit.set='all', screen.dat.method='lasso')$screen.index.var


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

### 2. Experiments ###
res=sapply(seeds, simplify="array", function (seed) {
  
  myprint(seed)
  
  # Bootstrap #
  set.seed(seed)
  idx.stratum1.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[18.4, 25)"), replace=TRUE)
  idx.stratum1.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[18.4, 25)"), replace=TRUE)
  idx.stratum2.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[25, 29.8)"), replace=TRUE)
  idx.stratum2.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[25, 29.8)"), replace=TRUE)
  #idx.stratum3.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/Blk_Hisp/[29.8, 40)"), replace=TRUE)
  idx.stratum3.case <- 134
  idx.stratum3.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/Blk_Hisp/[29.8, 40)"), replace=TRUE)
  idx.stratum4.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[18.4, 25)"), replace=TRUE)
  idx.stratum4.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[18.4, 25)"), replace=TRUE)
  idx.stratum5.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[25, 29.8)"), replace=TRUE)
  idx.stratum5.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[25, 29.8)"), replace=TRUE)
  idx.stratum6.case <- sample(which(Y_vaccine==1 & strata_vaccine=="Vacc/White/[29.8, 40)"), replace=TRUE)
  idx.stratum6.control <- sample(which(Y_vaccine==0 & strata_vaccine=="Vacc/White/[29.8, 40)"), replace=TRUE)
  
  idx.case <- c(idx.stratum1.case,idx.stratum2.case,idx.stratum3.case,idx.stratum4.case,idx.stratum5.case,idx.stratum6.case)
  idx.control <- c(idx.stratum1.control,idx.stratum2.control,idx.stratum3.control,idx.stratum4.control,idx.stratum5.control,idx.stratum6.control)
  Y_vaccine_boot <- Y_vaccine[c(idx.case, idx.control)]
  X_vaccine_boot <- X_vaccine[c(idx.case, idx.control),]
  weights_vaccine_boot <- weights_vaccine[c(idx.case, idx.control)]
  
  # Generating data where stacking is fitted #
  dat.X <- screen.dat.index(Y=Y_vaccine_boot, X=X_vaccine_boot, X_markers=X_markers, obsWeights=weights_vaccine_boot,
                            fit.set='all', screen.dat.method='no', screen.index.method='no')$dat
  
  if(can.set=='RF:tcell+glm:antibody'){
    # RF:tcell + glm:antibody #
    var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.tcell))
  } else if(can.set=='RF:tcell+glm:all'){
    # RF:tcell + glm:all #
    var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.tcell))
  } else if(can.set=='RF:all+glm:antibody'){
    # RF:all + glm:antibody #
    var.index.set <- list(glm=c(TRUE,screen.var.antibody), rf=c(TRUE,screen.var.all))
  } else if(can.set=='RF:all+glm:all'){
    # RF:all + glm:all #
    var.index.set <- list(glm=c(TRUE,screen.var.all), rf=c(TRUE,screen.var.all))
  }
  
  # Stacking (CV-AUC) #
  pred.vec <- c()
  for( i in 1:500 ){
    print(i)
    seed <- i
    # Stacking (CV-AUC) #
    res.temp <- get.st.cvauc(dat=dat.X, obsWeights=weights_vaccine, strata=strata_vaccine, var.index=var.index.set, method='method.NNloglik', seed=seed)
    pred.vec[i] <- mean(do.call(rbind,res.temp)[,'est.cvauc'])
  }
  mean(pred.vec)
})
res



# save res
foldername="res_"%.%method%.%"_"%.%can.set%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
