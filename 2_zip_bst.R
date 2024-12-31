library("xgboost")
## negative likelihood
neg_ll<-function(y,pi0,lambda,w=1){
  # the log-likelihood of zip model
  -ifelse(y==0,w*log(pi0+(1-pi0)*exp(-lambda)),
          w*(log(1-pi0)+y*log(lambda)-lambda))
}

p_linear<-function(x){log(x/(1-x))} # the log odds of success prob x

zip_sim<-function(n,seed){
  set.seed(seed)
  # covariates
  x1<-rnorm(n,0,0.5)
  x2<-runif(n,0,1)
  x3<-rgamma(n,2,0.5)
  x4<-rbinom(n,1,0.5)
  x5<-rbinom(n,1,0.2)
  # linear predictor
  # eta<--x1^2 + 0.2*log(x3) - 0.2*x4 +0.2*x4*x1 + log(0.5)
  eta<-tanh(x1^2) -tanh(x2*x4) + tanh(0.2*log(x3))
  lambda <- exp(eta)
  hist(lambda)
  hist(exp(-lambda))
  # pi_f<--2*x2^2 +x2 + 0.2*x5+0.3
  pi_f<-0
  pi0<-exp(pi_f)/(1+exp(pi_f))
  # generate data from zip
  z<-rbinom(n,1,pi0)
  y<-(1-z)*rpois(n,lambda = lambda)
  dat<-data.frame(y,x1,x2,x3,x4,x5,lambda,eta,pi0,pi_f,z)
  dat<-dat[order(dat$y),]
  dat$ind<-rep(1:5,round(n/5,0)+1)[1:n]
  dat
}

dat_aug<-function(dat){
  dat_aug<-dat
  dat_aug$ZZ<-0
  dat0<-dat
  dat0$ZZ<-1
  dat_aug<-rbind(dat_aug,dat0[dat0$y==0,])
  dat_aug
}

weight_update<-function(y, ZZ, lambda_hat, pi_hat){
  weight_vec<-rep(NA,length(y))
  for (i in 1:length(y)){
    yi<-y[i]; zi<-ZZ[i]; lambda_i<-lambda_hat[i]; pi_i<-pi_hat[i]
    if (yi>0 & zi==0) weight_vec[i]<-1
    if (yi>0 & zi==1) weight_vec[i]<-0
    if (yi==0 & zi==0) {
      weight_vec[i]<-(1-pi_i)*exp(-lambda_i)/(pi_i+(1-pi_i)*exp(-lambda_i))
    }
    if (yi==0 & zi==1) {
      weight_vec[i]<-pi_i/(pi_i+(1-pi_i)*exp(-lambda_i))
    }
  }
  weight_vec
}

mle_zip<-function(y,iter){
  # iter is the number of iterations
  # EM for null model
  n<-length(y)
  n0<-sum(y==0)
  y_bar<-mean(y)
  pi_hat<-n0/n
  lambda_hat<-y_bar
  loss<-NULL
  for (i in 1:iter){
    lambda_hat<-y_bar*(1-exp(-lambda_hat))/(1-n0/n)
    pi_hat<-1-y_bar/lambda_hat
    loss[i]<-mean(neg_ll(y,rep(pi_hat,n),rep(lambda_hat,n)))
  }
  list(pi_hat=pi_hat,lambda_hat=lambda_hat,loss=loss)
}

EM_zip <- function(dat_learn, var_names, M0, structure){
  # structure = "both", "lambda", "pi", "null"
  # dat_learn are learning data
  dat_learn_aug<-dat_aug(dat_learn[,c("y",var_names)])
  dat_learn_aug$pi_hat<-0.5
  dat_learn_aug$lambda_hat<-mean(dat_learn_aug$y)
  learn_loss<-rep(NA,M0)

  for (m in 1:M0){
  dat_learn_aug$w<-weight_update(dat_learn_aug$y, dat_learn_aug$ZZ, dat_learn_aug$lambda_hat, dat_learn_aug$pi_hat)
  if (structure=="both"){
    glm_poi<-glm(dat_learn_aug$y[dat_learn_aug$ZZ==0] ~., data=dat_learn_aug[dat_learn_aug$ZZ==0,var_names], family = poisson(link="log"), weights=dat_learn_aug$w[dat_learn_aug$ZZ==0])
    glm_ber<-glm(dat_learn_aug$ZZ ~ ., data = dat_learn_aug[,var_names], family = quasibinomial(link="logit"), weights = dat_learn_aug$w)
  }
  if (structure=="lambda"){
    glm_poi<-glm(dat_learn_aug$y[dat_learn_aug$ZZ==0] ~ ., data= dat_learn_aug[dat_learn_aug$ZZ==0,var_names], family = poisson(link="log"), weights=dat_learn_aug$w[dat_learn_aug$ZZ==0])
    glm_ber<-glm(dat_learn_aug$ZZ ~ 1, family = quasibinomial(link="logit"), weights=dat_learn_aug$w)
  }
  if (structure=="pi"){
    glm_poi<-glm(dat_learn_aug$y[dat_learn_aug$ZZ==0] ~ 1, family = poisson(link="log"), weights=dat_learn_aug$w[dat_learn_aug$ZZ==0])
    glm_ber<-glm(dat_learn_aug$ZZ ~ ., data=dat_learn_aug[,var_names], family = quasibinomial(link="logit"), weights = dat_learn_aug$w)
  }
  if (structure=="null"){
    glm_poi<-glm(dat_learn_aug$y[dat_learn_aug$ZZ==0] ~ 1, family = poisson(link="log"), weights=dat_learn_aug$w[dat_learn_aug$ZZ==0])
    glm_ber<-glm(dat_learn_aug$ZZ ~ 1, family = quasibinomial(link="logit"), weights=dat_learn_aug$w)
  }
  dat_learn_aug$lambda_hat<-predict(glm_poi, newdata = dat_learn_aug, type="response")
  dat_learn_aug$pi_hat<-predict(glm_ber, newdata = dat_learn_aug, type="response")
  dat_learn$lambda_hat<-predict(glm_poi, newdata = dat_learn, type="response")
  dat_learn$pi_hat<-predict(glm_ber, newdata = dat_learn, type="response")
  learn_loss[m]<-mean(neg_ll(dat_learn$y, dat_learn$pi_hat, dat_learn$lambda_hat))
  learn_loss_aug<-sum(neg_ll(dat_learn_aug$y, dat_learn_aug$pi_hat, dat_learn_aug$lambda_hat, dat_learn_aug$w))/sum(dat_learn_aug$w)
  print(c(m,round(c(learn_loss[m],learn_loss_aug),6)))
  }
  list(
    learn_loss=learn_loss,
    glm_poi=glm_poi,
    glm_ber=glm_ber
  )
}

EB_zip <- function(dat_learn, valid_rows, var_names, lambda0, pi0, M0,
                   n_tree_lambda, maxdepth_lambda, eta_lambda,
                   n_tree_pi, maxdepth_pi, eta_pi,
                   structure, trace, patience) {
  # dat_learn are the learning data;
  # the valid_rows are used to early stop the algorithm;
  # lambda0, pi0 are the initial parameters;
  # M0 is the number of EB iterations;
  # n_tree_lambda is the number of trees in lambda boosting in each B step;
  # maxdepth_lambda is the maximum depth of tree for lambda;
  # eta_lambda is the learning rate in the boosting algorithm;
  # similar for n_tree_pi, maxdepth_pi and eta_pi
  # structure can be "both", "p", "mu";
  # patience is the number of non-improvements in validation loss before early stop

  dat_learn0<-data.frame(dat_learn, pi_hat=pi0, lambda_hat=lambda0, pi_margin=p_linear(pi0), lambda_margin=log(lambda0))
  dat_train_aug<-dat_aug(dat_learn0[-valid_rows,])
  n_train<-nrow(dat_learn0[-valid_rows,])
  dat_valid<-dat_learn0[valid_rows,]
  train_loss <- NULL; valid_loss <- NULL
  train_impT <- NULL; train_impT[1]<-NULL
  valid_impT <- NULL; valid_impT[1]<-NULL
  valid_impT2<- NULL; valid_impT2[1:patience]<-NULL
  pi_models<-list()
  lambda_models<-list()

  for (m in 1:M0) {
    # update the weights
    dat_train_aug$w<-weight_update(dat_train_aug$y, dat_train_aug$ZZ, dat_train_aug$lambda_hat, dat_train_aug$pi_hat)
    dat_train_aug$lambda_w<-dat_train_aug$w*(1-dat_train_aug$ZZ)

    # boosting for lambda
    param<-list(max_depth=maxdepth_lambda, eta =eta_lambda, objective="count:poisson")
    if (structure == "lambda"| structure == "both") {
      dtrain<-xgb.DMatrix(data=sparse.model.matrix(~., data=dat_train_aug[,var_names]), label = dat_train_aug$y, weight=dat_train_aug$lambda_w, base_margin= dat_train_aug$lambda_margin)
      dvalid<-xgb.DMatrix(data=sparse.model.matrix(~., data=dat_valid[,var_names]), base_margin= dat_valid$lambda_margin)
      dlearn<-xgb.DMatrix(data=sparse.model.matrix(~., data=dat_learn0[,var_names]), base_margin= dat_learn0$lambda_margin)
    }
    if (structure == "pi"){
      dtrain<-xgb.DMatrix(data=as.matrix(rep(1,nrow(dat_train_aug))), label = dat_train_aug$y, weight=dat_train_aug$lambda_w, base_margin=dat_train_aug$lambda_margin)
      dvalid<-xgb.DMatrix(data=as.matrix(rep(1,nrow(dat_valid))), base_margin= dat_valid$lambda_margin)
      dlearn<-xgb.DMatrix(data=as.matrix(rep(1,nrow(dat_learn0))), base_margin= dat_learn0$lambda_margin)
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

    # boosting for pi
    param<-list(max_depth=maxdepth_pi, eta =eta_pi, objective="binary:logistic")
    if (structure == "pi"| structure == "both") {
      dtrain<-xgb.DMatrix(data=sparse.model.matrix(~., data=dat_train_aug[,var_names]), label = dat_train_aug$ZZ, weight=dat_train_aug$w, base_margin= dat_train_aug$pi_margin)
      dvalid<-xgb.DMatrix(data=sparse.model.matrix(~., data=dat_valid[,var_names]), base_margin=dat_valid$pi_margin)
      dlearn<-xgb.DMatrix(data=sparse.model.matrix(~., data=dat_learn0[,var_names]), base_margin= dat_learn0$pi_margin)
    }
    if (structure == "lambda"){
      dtrain<-xgb.DMatrix(data=as.matrix(rep(1,nrow(dat_train_aug))), label = dat_train_aug$ZZ, weight=dat_train_aug$w, base_margin=dat_train_aug$pi_margin)
      dvalid<-xgb.DMatrix(data=as.matrix(rep(1,nrow(dat_valid))), base_margin=dat_valid$pi_margin)
      dlearn<-xgb.DMatrix(data=as.matrix(rep(1,nrow(dat_learn0))), base_margin= dat_learn0$pi_margin)
    }
    watchlist=list(train=dtrain)
    pi_bst <- xgb.train(param, dtrain, nrounds=n_tree_pi, verbose = 0, watchlist, early_stopping_rounds = NULL)
    dat_train_aug$pi_hat<-predict(pi_bst, newdata = dtrain)
    dat_train_aug$pi_margin<-predict(pi_bst, newdata = dtrain, outputmargin = T)
    dat_valid$pi_hat<-predict(pi_bst, newdata = dvalid)
    dat_valid$pi_margin<-predict(pi_bst, newdata = dvalid, outputmargin = T)
    dat_learn0$pi_hat<-predict(pi_bst, newdata = dlearn)
    dat_learn0$pi_margin<-predict(pi_bst, newdata = dlearn, outputmargin = T)
    pi_models[[m]]<-pi_bst

    # outer train loss and validation loss
    train_loss[m] <- sum(neg_ll(dat_train_aug$y, dat_train_aug$pi_hat, dat_train_aug$lambda_hat, dat_train_aug$w))/sum(dat_train_aug$w)
    valid_loss[m] <- mean(neg_ll(dat_valid$y, dat_valid$pi_hat, dat_valid$lambda_hat))
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
      print(paste("iter-", m, ":", "train_ls", round(train_loss[m], 4), "(" ,train_impT[m], "),", "valid_ls", round(valid_loss[m], 4), "(" ,valid_impT[m],"/",valid_impT2[m], ")", sep=""))
    }
    if (m>patience) if (valid_impT2[m]==0) {break}
  }
  list(
    lambda0 = lambda0,
    pi0 = pi0,
    lambda_hat = dat_learn0$lambda_hat,
    pi_hat = dat_learn0$pi_hat,
    structure = structure,
    train_loss = train_loss,
    valid_loss = valid_loss,
    lambda_models = lambda_models,
    pi_models = pi_models
  )
}

