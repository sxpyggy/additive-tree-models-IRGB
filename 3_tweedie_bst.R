library(CASdatasets)
library(gamlss)
library(fitdistrplus)
library(nnet)
library(mboost)
library(xgboost)
library(mvtnorm)
library(tweedie)
library(statmod)
library(Matrix)

data (freMTPL2freq)
dat <- freMTPL2freq
dat$VehGas <- factor (dat$VehGas)
dat$ClaimNb <- as.numeric(dat$ClaimNb)
data (freMTPL2sev)
sev <- freMTPL2sev
sev$ClaimNb0 <- 1
dat0 <-
  aggregate (sev , by = list (IDpol = sev$IDpol), FUN = sum)[,c(1, 3:4)]
names (dat0)[2] <- "ClaimTotal"
dat <- merge (x = dat,
              y = dat0,
              by = "IDpol",
              all.x = TRUE)
dat [is.na(dat)] <- 0
cor(dat$ClaimNb, dat$ClaimNb0)
dat <- dat [which (dat$ClaimNb0 <= 5) , ]
dat <- dat[dat$Exposure>0.5, ]
dat <- dat[, -14]
dat$Exposure <- pmin (dat$Exposure , 1)
sev <- sev [which (sev$IDpol %in% dat$IDpol), c (1 , 2)]
dat$VehBrand <- factor (
  dat$VehBrand,
  levels = c (
    "B1" ,
    "B2" ,
    "B3" ,
    "B4" ,
    "B5" ,
    "B6",
    "B10" ,
    "B11" ,
    "B12" ,
    "B13" ,
    "B14"
  )
)
####
dat$AreaGLM <- as.integer (dat$Area)
dat$VehPowerGLM <- as.factor (pmin (dat$VehPower , 9))
dat$VehAgeGLM <- as.factor (cut (
  dat$VehAge ,
  c (0 , 5 , 12 , 101) ,
  labels = c ("0 -5" , "6 -12" , "12+") ,
  include.lowest = TRUE
))
dat$DrivAgeGLM <-
  as.factor (cut (
    dat$DrivAge ,
    c (18 , 20 , 25 , 30 , 40 , 50 , 70 , 101) ,
    labels = c ("18 -20" , "21 -25" , "26 -30" , "31 -40" , "41 -50" , "51 -70" , "71+") ,
    include.lowest = TRUE
  ))
dat$BonusMalusGLM <- pmin (dat$BonusMalus , 150)
dat$DensityGLM <- log (dat$Density)
dat$VehGasGLM <- as.integer(dat$VehGas)

####
dat$ClaimNb[dat$ClaimTotal==0]<-0
dat$Sev<-dat$ClaimTotal/dat$ClaimNb
dat$Y<-dat$ClaimTotal/dat$Exposure
dat<-dat[order(dat$ClaimTotal),]
dat$ind<-rep(1:5,round(nrow(dat)/5,0)+1)[1:nrow(dat)]
# dput(names(dat))
dat<-dat[ ,c("IDpol", "ind",
            "Exposure", "ClaimNb",
            "ClaimTotal", "Sev", "Y",
            "Region", "Area", "Density", "AreaGLM", "DensityGLM",
            "DrivAge", "DrivAgeGLM",
            "BonusMalus", "BonusMalusGLM",
            "VehPower", "VehAge", "VehBrand", "VehGas",
            "VehPowerGLM", "VehAgeGLM", "VehGasGLM")]
# str(dat)

# functions used

vnames_glm<-c("Region", "AreaGLM", "DrivAgeGLM", "BonusMalusGLM",
              "VehPowerGLM", "VehAgeGLM", "VehBrand", "VehGas")

condi_n<-function(Y,N,e,lambda,alpha,tau){
  # this is the conditional distribution of N given Y
    p<-(alpha+2)/(alpha+1)
    mu<-lambda*tau
    phi<-(lambda^(1-p))*(tau^(2-p))/(2-p)
    cond<-dpois(N,lambda*e)*dgamma(Y,shape = N*alpha, scale=tau/(alpha*e))/dtweedie(Y, power=p,mu=mu,phi=phi/e)
  cond
}

unit_d<-function(Y,mu,p,e){
  # this is the unit deviance in the tweedie CP model
  2*e*(Y^(2-p)/(1-p)/(2-p)-Y*mu^(1-p)/(1-p)+mu^(2-p)/(2-p))
}

DAT_AUG<-function(dat){
  dat1<-dat[dat$Y>0,]
  dat1$ClaimNb<-1
  dat2<-dat[dat$Y>0,]
  dat2$ClaimNb<-2
  dat3<-dat[dat$Y>0,]
  dat3$ClaimNb<-3
  dat4<-dat[dat$Y>0,]
  dat4$ClaimNb<-4
  dat5<-dat[dat$Y>0,]
  dat5$ClaimNb<-5
  dat_aug<-rbind(dat[dat$Y==0,],dat1,dat2,dat3,dat4,dat5)
}

EM_tweedie <- function(dat_learn, M0){
  # initial values
  dat_learn$ClaimNb<-0
  dat_learn$O<-as.numeric(dat_learn$Y>0)
  em_poi<-glm(O ~ 1 + offset(log(Exposure)), family =poisson(link="log"),data=dat_learn)
  em_gam<-glm(Y ~ 1, family = Gamma(link="log"), data=dat_learn[dat_learn$Y>0,])
  (em_alpha<-1/summary(em_gam)$dispersion)
  dat_aug<-DAT_AUG(dat_learn)
  dat_aug$em_lambda<-mean(dat_learn$O)
  dat_aug$em_tau<-mean(dat_learn$Y[dat_learn$Y>0])
  dat_aug$cond<-NULL
  dat_aug$cond[dat_aug$Y==0]<-1
  for (m in 1:M0){
    # weight updates
    for (i in which(dat_aug$Y>0)){
      dat_aug$cond[i]<-condi_n(dat_aug$Y[i],dat_aug$ClaimNb[i],dat_aug$Exposure[i],dat_aug$em_lambda[i],em_alpha,dat_aug$em_tau[i])
    }
    sum(dat_aug$cond);nrow(dat_aug)
    dat_aug$em_y_weight<-dat_aug$cond*dat_aug$ClaimNb
    em_poi<-glm(ClaimNb ~ Region + AreaGLM + DrivAgeGLM + BonusMalusGLM +
                  VehPowerGLM  + VehAgeGLM  +  VehBrand + VehGas +
                  offset(log(Exposure)), weights = cond,
                family =poisson(link="log"), data=dat_aug)
    em_gam<-glm(Y ~ Region + AreaGLM + DrivAgeGLM + BonusMalusGLM +
                  VehPowerGLM  + VehAgeGLM  +  VehBrand + VehGas +
                  offset(log(ClaimNb/Exposure)), weights=em_y_weight,
                family = Gamma(link="log"), data=dat_aug[dat_aug$Y>0,])
    dat_aug$em_lambda<-
      predict(em_poi,newdata =data.frame(dat_aug[,vnames_glm],Exposure=1), type="response")
    dat_aug$em_tau<-
      predict(em_gam, newdata = data.frame(dat_aug[,vnames_glm],Exposure=1,
                                           ClaimNb=1), type="response")
    yind<-which(dat_aug$Y>0)
    logL_alpha<-function(alpha){
      sum(dat_aug$cond[yind]*(alpha*dat_aug$ClaimNb[yind]*(log(dat_aug$Y[yind])-log(dat_aug$em_tau[yind])+log(alpha*dat_aug$Exposure[yind]))-alpha*dat_aug$Y[yind]*dat_aug$Exposure[yind]/dat_aug$em_tau[yind]-lgamma(alpha*dat_aug$ClaimNb[yind])))
    }
    em_alpha<- optimise(logL_alpha, maximum = T,interval = c(0, 10))$maximum
    em_p<-(em_alpha+2)/(em_alpha+1)
    dat_aug$em_mu<-dat_aug$em_lambda*dat_aug$em_tau
    dat_aug$em_phi <- dat_aug$em_lambda^(1-em_p)*dat_aug$em_tau^(2-em_p)/(2-em_p)
    dat_aug$logL_em<- -log(dtweedie(dat_aug$Y,xi=em_p,mu=dat_aug$em_mu,
                                    phi=dat_aug$em_phi/dat_aug$Exposure))
    loss_EM<-mean(dat_aug$logL_em[1:nrow(dat_learn)])
    print(paste("iter,p,loss:",m, round(em_p,2), round(loss_EM,8)))
  }
  list(em_alpha=em_alpha,
       em_p=em_p,
       em_poi=em_poi,
       em_gam=em_gam,
       dat_aug=dat_aug)
}

EB_tweedie <- function(dat_learn, valid_rows, vnames_bst,
                       lambda0, tau0, alpha0, M0,
                       n_tree_lambda, maxdepth_lambda, eta_lambda,
                       n_tree_tau, maxdepth_tau, eta_tau,
                       structure, trace, patience) {
  # dat_learn are the learning data;
  # the valid_rows are used to early stop the algorithm;
  # vnames_bst are the variables used in the regressions;
  # lambda0, tau0, alpha0 are the initial parameters;
  # M0 is the number of EB iterations;
  # n_tree_lambda is the number of trees in lambda boosting in each B step;
  # maxdepth_lambda is the maximum depth of tree for lambda;
  # eta_lambda is the learning rate in the boosting algorithm;
  # similar for n_tree_tau, maxdepth_tau and eta_tau
  # structure can be "both", "lambda", "tau";
  # patience is the number of non-improvements in validation loss before early stop
  dat_learn$ClaimNb<-0
  dat_learn0<-data.frame(dat_learn, lambda_hat=lambda0, tau_hat=tau0, lambda_margin=log(lambda0), tau_margin=log(tau0))
  dat_train_aug<-DAT_AUG(dat_learn0[-valid_rows,])
  dat_valid<-dat_learn0[valid_rows,]
  train_loss <- NULL; valid_loss <- NULL; learn_loss <- NULL
  train_impT <- NULL; train_impT[1]<-NULL
  valid_impT <- NULL; valid_impT[1]<-NULL
  valid_impT2<- NULL; valid_impT2[1:patience]<-NULL
  lambda_models<-list()
  tau_models<-list()
  dat_train_aug$cond<-NULL
  dat_train_aug$cond[dat_train_aug$Y==0]<-1
  alpha_hat<-alpha0
  train_mat<-sparse.model.matrix(~., data=dat_train_aug[,vnames_bst])
  train_mat_0<-as.matrix(rep(1,nrow(dat_train_aug)))
  valid_mat<-sparse.model.matrix(~., data=dat_valid[,vnames_bst])
  valid_mat_0<-as.matrix(rep(1,nrow(dat_valid)))
  learn_mat<-sparse.model.matrix(~., data=dat_learn0[,vnames_bst])
  learn_mat_0<-as.matrix(rep(1,nrow(dat_learn0)))

  for (m in 1:M0){
    # weight updates
    for (i in which(dat_train_aug$Y>0)){
      dat_train_aug$cond[i]<-condi_n(dat_train_aug$Y[i],dat_train_aug$ClaimNb[i],dat_train_aug$Exposure[i],dat_train_aug$lambda_hat[i],alpha_hat,dat_train_aug$tau_hat[i])
    }
    sum(dat_train_aug$cond);nrow(dat_train)
    dat_train_aug$y_weight<-dat_train_aug$cond*dat_train_aug$ClaimNb

    # boosting for lambda
    param<-list(max_depth=maxdepth_lambda, eta =eta_lambda, objective="count:poisson")
    if (structure == "lambda"| structure == "both") {
      dtrain<-xgb.DMatrix(data=train_mat, label = dat_train_aug$ClaimNb, weight=dat_train_aug$cond, base_margin= dat_train_aug$lambda_margin)
      dvalid<-xgb.DMatrix(data=valid_mat, base_margin= dat_valid$lambda_margin)
      dlearn<-xgb.DMatrix(data=learn_mat, base_margin= dat_learn0$lambda_margin)
    }
    if (structure == "tau"){
      dtrain<-xgb.DMatrix(data=train_mat_0, label = dat_train_aug$ClaimNb, weight=dat_train_aug$cond, base_margin=dat_train_aug$lambda_margin)
      dvalid<-xgb.DMatrix(data=valid_mat_0, base_margin= dat_valid$lambda_margin)
      dlearn<-xgb.DMatrix(data=learn_mat_0, base_margin= dat_learn0$lambda_margin)
    }
    watchlist=list(train=dtrain)
    lambda_bst <- xgb.train(param, dtrain, nrounds=n_tree_lambda, verbose = 0, watchlist, early_stopping_rounds = NULL)
    dat_train_aug$lambda_hat<-predict(lambda_bst, newdata = dtrain)
    dat_train_aug$lambda_margin<-predict(lambda_bst, newdata = dtrain, outputmargin = T)
    dat_valid$lambda_hat<-predict(lambda_bst, newdata = dvalid)
    dat_valid$lambda_margin<-predict(lambda_bst, newdata = dvalid, outputmargin = T)
    dat_learn0$lambda_hat<-predict(lambda_bst, newdata = dlearn)
    dat_learn0$lambda_margin<-predict(lambda_bst, newdata = dlearn, outputmargin = T)
    lambda_models[[m]]<-lambda_bst

    # boosting for tau
    param<-list(max_depth=maxdepth_tau, eta =eta_tau, objective="reg:gamma")
    if (structure == "tau"| structure == "both") {
      dtrain<-xgb.DMatrix(data=train_mat, label = dat_train_aug$Y, weight=dat_train_aug$y_weight, base_margin= dat_train_aug$tau_margin)
      dvalid<-xgb.DMatrix(data=valid_mat, base_margin=dat_valid$tau_margin)
      dlearn<-xgb.DMatrix(data=learn_mat, base_margin=dat_learn0$tau_margin)
    }
    if (structure == "lambda"){
      dtrain<-xgb.DMatrix(data=train_mat_0, label = dat_train_aug$Y, weight=dat_train_aug$y_weight, base_margin=dat_train_aug$tau_margin)
      dvalid<-xgb.DMatrix(data=valid_mat_0, base_margin=dat_valid$tau_margin)
      dlearn<-xgb.DMatrix(data=learn_mat_0, base_margin=dat_learn0$tau_margin)
    }
    yind<-which(dat_train_aug$Y>0)
    dtrain_y<-slice(dtrain,yind)
    watchlist=list(train=dtrain_y)
    tau_bst <- xgb.train(param, dtrain_y, nrounds=n_tree_tau, verbose = 0, watchlist, early_stopping_rounds = NULL)
    dat_train_aug$tau_hat<-predict(tau_bst, newdata = dtrain)
    dat_train_aug$tau_margin<-predict(tau_bst, newdata = dtrain, outputmargin = T)
    dat_valid$tau_hat<-predict(tau_bst, newdata = dvalid)
    dat_valid$tau_margin<-predict(tau_bst, newdata = dvalid, outputmargin = T)
    dat_learn0$tau_hat<-predict(tau_bst, newdata = dlearn)
    dat_learn0$tau_margin<-predict(tau_bst, newdata = dlearn, outputmargin = T)
    tau_models[[m]]<-tau_bst

    logL_alpha<-function(alpha){
      sum(dat_train_aug$cond[yind]*(
        alpha*dat_train_aug$ClaimNb[yind]*(
          log(dat_train_aug$Y[yind])-log(dat_train_aug$tau_hat[yind])+
            log(alpha*dat_train_aug$Exposure[yind]))-
          alpha*dat_train_aug$Y[yind]*dat_train_aug$Exposure[yind]/
          dat_train_aug$tau_hat[yind]-lgamma(alpha*dat_train_aug$ClaimNb[yind])
        ))
    }
    alpha_hat<- optimise(logL_alpha, maximum = T,interval = c(0, 10))$maximum
    p_hat<-(alpha_hat+2)/(alpha_hat+1)

    # outer train loss and validation loss
    dat_learn0$mu_hat<-dat_learn0$lambda_hat*dat_learn0$tau_hat
    dat_learn0$phi_hat <- dat_learn0$lambda_hat^(1-p_hat)*
      dat_learn0$tau_hat^(2-p_hat)/(2-p_hat)
    dat_learn0$logL_hat<- -log(dtweedie(dat_learn0$Y,xi=p_hat,mu=dat_learn0$mu_hat, phi=dat_learn0$phi_hat/dat_learn0$Exposure))
    valid_loss[m]<-mean(dat_learn0$logL_hat[valid_rows])
    train_loss[m]<-mean(dat_learn0$logL_hat[-valid_rows])
    learn_loss[m]<-mean(dat_learn0$logL_hat)

    if (m>1){
      train_impT[m]<- (train_loss[m-1]-train_loss[m]>0)
      valid_impT[m]<- (valid_loss[m-1]-valid_loss[m]>0)
    }
    if(m>patience) {
      valid_impT2[m]<- (valid_loss[m-patience]-valid_loss[m]>0)
    }
    if (trace == T) {
      # iter-m: train_ls (whether improved compared to the previous iteration),
      # valid_ls (whether improved compared to the previous iteration/
      # whether improved compared to the iteration patience ahead)
      print(paste("iter-", m, ":",
                  "learn_ls", round(learn_loss[m], 4), ",",
                  "train_ls", round(train_loss[m], 4), "(" ,train_impT[m], "),",
                  "valid_ls", round(valid_loss[m], 4), "(" ,valid_impT[m],"/",valid_impT2[m], ")", sep=""))
    }
    if (m>patience) if (valid_impT2[m]==0) {break}
  }
  list(
    lambda0 = lambda0,
    tau0 = tau0,
    alpha0 = alpha0,
    lambda_hat = dat_learn0$lambda_hat,
    tau_hat = dat_learn0$tau_hat,
    alpha_hat = alpha_hat,
    p_hat = p_hat,
    structure = structure,
    train_loss = train_loss,
    valid_loss = valid_loss,
    learn_loss = learn_loss,
    lambda_models = lambda_models,
    tau_models = tau_models
  )
}


