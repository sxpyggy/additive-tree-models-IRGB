library("gbm")
library("splines")
# variable names in loop
vname<-function(vname,k){paste(vname, k, sep="")}
vname2<-function(vname, K){
  temp<-paste(vname, 1, sep="")
  if (K>1) for (k in 2:K) temp<-c(temp,paste(vname, k, sep=""))
  temp
}

neg_ll <- function(y, mu, sigma, p, w=NULL) {
  # calculate the average negative log-likelihood
  # y has dimension of n; mu, sigma, p have dimension of K*n
  # w is the weights
  K<-ncol(mu)
  LL<-rep(0,nrow(mu))
  if (is.null(w)==T) w=rep(1,nrow(mu))
  for (k in 1:K){
    LL<-LL+p[, k] * dnorm(y, mean = mu[, k], sd = sigma[, k], log=F)
  }
  -sum(w*log(LL))/sum(w)
}

sim_gaussian <- function(n, seed) {
  # mixing probs are related to covariates, means are not related to covariates
  set.seed(seed)
  X1 = rnorm(n)
  X2 = rnorm(n)
  X3 = rnorm(n)
  X4 = rbinom(n,1,0.5)
  X5 = X1+rnorm(n)
  X6 = rnorm(n)

  F1 = tanh(X1^2+X2)
  F2 = tanh(X1+X2^2)
  F3 = tanh(2*X4+X3+X3^2+2*X3*X4)

  # F1 = rep(1,n)
  # F2 = rep(1,n)
  # F3 = rep(1,n)
  # pairs(cbind(F1,F2,F3,X1,X2,X3,X4))

  P <- FtoP(FF = cbind(F1, F2, F3))

  # pairs(cbind(P,X1,X2,X3,X4))
  # boxplot(P)

  MU1 <- rep(-5,n)
  MU2 <- rep(0,n)
  MU3 <- rep(5,n)
  # MU1 <- 10 + 2 * X1 + exp(0.3 * X2)
  # MU2 <- 2- 0.1 * X1 ^ 2 + 2* X3 - X1 * X2
  # MU3 <- -10 - 2 * X1 +  3*X4 + X6^2
  # # boxplot(cbind(MU1,MU2,MU3))
  dat <- data.frame(X1, X2, X3, X4, X5, X6, F1, F2, F3, P1 = P[, 1], P2 = P[, 2], P3 = P[, 3], MU1, MU2, MU3, SG1 = 0.5, SG2 = 1, SG3 = 2, Z1 = NA, Z2 = NA, Z3 = NA, Y1 = NA, Y2 = NA, Y3 = NA, Y = NA)
  for (i in 1:nrow(dat)) {
    dat[i, c("Z1", "Z2", "Z3")] <-
      t(rmultinom(1, size = 1, dat[i, c("P1", "P2", "P3")]))
    dat[i, c("Y1", "Y2", "Y3")] <-
      t(rnorm(3, mean = as.numeric(dat[i, c("MU1", "MU2", "MU3")]), sd = as.numeric(dat[i, c("SG1", "SG2", "SG3")])))
    dat$Y[i] <-
      sum(dat[i, c("Z1", "Z2", "Z3")] * dat[i, c("Y1", "Y2", "Y3")])
  }
  dat[,c("F1","F2","F3")]<-PtoF(dat[,c("P1","P2","P3")])
  dat
}

dat_augmentation<-function(dat,K){
  n<-nrow(dat)
  dat_aug<-dat
  dat_aug$ZZ<-1
  dat_aug$R0<-1:n
  for (k in 2:K){
    dat_augk<-dat
    dat_augk$ZZ<-k
    dat_augk$R0<-1:n
    dat_aug<-rbind(dat_aug,dat_augk)
  }
  ZZZ<-data.frame(model.matrix(~ 0 + as.factor(dat_aug$ZZ)))
  names(ZZZ)<-vname2("ZZ",K)
  dat_aug<-data.frame(dat_aug,ZZZ)
  dat_aug
}

dat_Norm<-function(K,n){
  # put into a suitable format for mixture of Gaussian examples
  # z is the expected indicator of component w.r.t. the conditional pmf z|y; f is pdf of components; mu is the mean of component
  # p is the mixing probabilities; plinear is  the linear predictor in mixing probs
  # sigma is the variance; weight is the weight for mu regression (same as z for mixture model)
  dat_bst<-data.frame(matrix(NA,ncol=7*K,nrow=n))
  bst_names<-NULL
  for (k in 1:K){
    vnames<-c(paste("z",k,sep=""), paste("f",k,sep=""), paste("mu",k,sep=""), paste("p",k,sep=""), paste("plinear",k,sep=""), paste("sigma",k,sep=""), paste("weight",k,sep=""))
    bst_names<-c(bst_names,vnames)
  }
  names(dat_bst)<-bst_names
  dat_bst
}

EM_gaussian<-function(X,Y,K,M0,structure,trace,patience) {
  # structure = "both", "mu", "p", "null"
  # X,Y are learning data
  par_mat <- dat_Norm(K = K, n = nrow(X))

  # initialization
  mu_models<-NULL
  p0<-rep(1/K,K); mu0<-quantile(Y,(1:K)/(K+1))
  for (k in 1:K){
    par_mat[,names(par_mat)==vname("p",k)]=p0[k]
    par_mat[,names(par_mat)==vname("mu",k)]=mu0[k]
    par_mat[,names(par_mat)==vname("sigma", k)]=sd(Y)
    mu_models[[k]]<-list()
  }
  p_models<-list()
  learn_loss <- NULL
  learn_impT <- NULL; learn_impT[1]<-T
  learn_impT2<- NULL; learn_impT2[1]<-T

  for (m in 1:M0) {
    # expectation of latent variable
    # also the weights for mean and response for mixing probabilities
    sum_fz<-rep(0,nrow(par_mat))
    for (k in 1:K){
      par_mat[,names(par_mat)==vname("f",k)]<-dnorm(Y, get(vname("mu",k), par_mat), get(vname("sigma",k),par_mat))
      sum_fz<-sum_fz+get(vname("p",k), par_mat)*get(vname("f",k), par_mat)
    }
    par_mat[,vname2("z",K)] <- par_mat[,vname2("p",K)]*par_mat[,vname2("f",K)]/sum_fz
    par_mat[,vname2("weight",K)] <- par_mat[,vname2("z",K)]

    #  mean regression
    for (k in 1:K){
      if (structure == "mu" | structure == "both") {
        mu_glm <- lm(Y ~., data = X, weights=get(vname("weight",k),par_mat))
        par_mat[,vname("mu",k)] <- predict(mu_glm, newdata = X)
      }
      if (structure == "null" | structure == "p") {
        mu_glm <- lm(Y ~ 1, weights=get(vname("weight",k),par_mat))
        par_mat[,vname("mu",k)] <- predict(mu_glm, newdata = data.frame(rep(1,nrow(par_mat))))
      }
      mu_models[[k]][[m]]<-mu_glm
    }

    # mixing probabilities regression
    if (structure == "p" | structure == "both") {
      p_glm <-nnet::multinom(as.matrix(par_mat[, vname2("z",K)]) ~ ., data = X,trace = F)
      par_mat[, vname2("p",K)] <- predict(p_glm, newdata = X, type = "probs")
    }
    if (structure=="null" | structure == "mu"){
      p_glm <- nnet::multinom(as.matrix(par_mat[,vname2("z",K)]) ~ 1,  trace=F)
      par_mat[, vname2("p",K)] <- predict(p_glm, newdata = data.frame(rep(1,nrow(par_mat))), type = "probs")
    }
    p_models[[m]]<-p_glm

    # sigma
    for (k in 1:K){
      sigma_est<-sqrt(sum(get(vname("z",k), par_mat)*(Y - get(vname("mu",k), par_mat)) ^ 2) / sum(get(vname("z",k), par_mat)))
      par_mat[,vname("sigma",k)]<-sigma_est
    }

    learn_loss[m] <- neg_ll(Y, par_mat[, vname2("mu",K)], par_mat[, vname2("sigma",K)], par_mat[, vname2("p",K)])
    if (m>1){
      learn_impT[m] <- (learn_loss[m-1]-learn_loss[m])>0
      learn_impT2[m] <- (learn_loss[m-1]-learn_loss[m])>10^-4
    }
    if (trace == T) {
      # iter-m: learning loss (whether improved compared to the previous/
      # whether improved more than 10^4 compared to the previous)
      print(paste("iter-",m, ":", "learn_ls", round(learn_loss[m], 4),"(" ,learn_impT[m],"/",learn_impT2[m],")",sep=""))
    }
    # if patience number of consecutive iteration less than 10^4, then stop
    if (m>patience) if (sum(learn_impT2[(m-patience+1):m])==0) break
  }
  par_mat[, vname2("plinear",K)] <- PtoF(par_mat[, vname2("p",K)])

  list(
    iter = m,
    learn_loss = learn_loss,
    par_mat = par_mat,
    mu_models = mu_models,
    p_models = p_models,
    sigma = apply(par_mat[,vname2("sigma",K)], 2, unique)
  )
}

EM_gaussian_aug <-
  function(X,Y,K,M0,structure,trace,patience) {
  # structure = "both", "mu", "p", "null"
  # X,Y are learning data
  dat_aug<-dat_augmentation(data.frame(X,Y),K)
  X_aug<-dat_aug[,vname2("X",6)]
  Y_aug<-dat_aug$Y
  ZZ_aug<-dat_aug[,vname2("ZZ",K)]
  par_mat <- dat_Norm(K = K, n = nrow(X_aug))

  # initialization
  mu_models<-NULL
  p0<-rep(1/K,K); mu0<-quantile(Y,(1:K)/(K+1))
  for (k in 1:K){
    par_mat[,vname("p",k)]=p0[k]
    par_mat[,vname("mu",k)]=mu0[k]
    par_mat[,vname("sigma", k)]=sd(Y)
    mu_models[[k]]<-list()
  }
  p_models<-list()
  learn_loss <- NULL
  learn_impT <- NULL; learn_impT[1]<-T
  learn_impT2<- NULL; learn_impT2[1]<-T

  for (m in 1:M0) {
    # update the weights
    sum_fz<-rep(0,nrow(par_mat))
    for (k in 1:K){
      par_mat[,vname("f",k)]<-dnorm(Y_aug, par_mat[,vname("mu",k)], par_mat[,vname("sigma",k)])
      sum_fz<-sum_fz + par_mat[,vname("p",k)]*par_mat[,vname("f",k)]
    }
    par_mat[,vname2("weight",K)] <- (par_mat[,vname2("p",K)]*par_mat[,vname2("f",K)])/sum_fz
    par_mat[,vname2("weight",K)] <- par_mat[,vname2("weight",K)]*ZZ_aug
    par_mat$weight<-apply(par_mat[,vname2("weight",K)],1,sum)

    #  mean regression
    for (k in 1:K){
      if (structure == "mu" | structure == "both") {
        mu_glm <- lm(Y_aug ~., data = X_aug, weights=par_mat[,vname("weight",k)])
        par_mat[,vname("mu",k)] <- predict(mu_glm, newdata = X_aug)
      }
      if (structure == "null" | structure == "p") {
        mu_glm <- lm(Y_aug ~ 1, weights=par_mat[,vname("weight",k)])
        par_mat[,vname("mu",k)] <- predict(mu_glm, newdata = data.frame(rep(1,nrow(par_mat))))
      }
      mu_models[[k]][[m]]<-mu_glm
    }

    # mixing probabilities regression
    if (structure == "p" | structure == "both") {
      p_glm <-nnet::multinom(as.matrix(ZZ_aug) ~ ., data = X_aug, trace = F, weights=par_mat$weight)
      par_mat[, vname2("p",K)] <- predict(p_glm, newdata = X_aug, type = "probs")
    }
    if (structure=="null" | structure == "mu"){
      p_glm <- nnet::multinom(as.matrix(ZZ_aug) ~ 1,  trace=F, weights=par_mat$weight)
      par_mat[, vname2("p",K)] <- predict(p_glm, newdata = data.frame(rep(1,nrow(par_mat))), type = "probs")
    }
    p_models[[m]]<-p_glm

    # sigma
    for (k in 1:K){
      sigma_est<-sqrt(sum(par_mat[,vname("weight",k)]*(Y_aug - par_mat[vname("mu",k)]) ^ 2) / sum(par_mat[,vname("weight",k)]))
      par_mat[,vname("sigma",k)]<-sigma_est
    }

    learn_loss[m] <- neg_ll(Y, par_mat[, vname2("mu",K)], par_mat[, vname2("sigma",K)], par_mat[, vname2("p",K)], w=par_mat$weight)
    if (m>1){
      learn_impT[m] <- (learn_loss[m-1]-learn_loss[m])>0
      learn_impT2[m] <- (learn_loss[m-1]-learn_loss[m])>10^-4
    }
    if (trace == T) {
      # iter-m: learning loss (whether improved compared to the previous/
      # whether improved more than 10^4 compared to the previous)
      print(paste("iter-",m, ":", "learn_ls", round(learn_loss[m], 4),"(" ,learn_impT[m],"/",learn_impT2[m],")",sep=""))
    }
    # if patience number of consecutive iteration less than 10^4, then stop
    if (m>patience) if (sum(learn_impT2[(m-patience+1):m])==0) break
  }
  par_mat[, vname2("plinear",K)] <- PtoF(par_mat[, vname2("p",K)])

  list(
    iter = m,
    learn_loss = learn_loss,
    par_mat = par_mat,
    mu_models = mu_models,
    p_models = p_models,
    sigma = apply(par_mat[,vname2("sigma",K)], 2, unique)
  )
}

EB_gaussian <- function(X, Y, valid_rows, mu0, p0, sigma0, K, M0,
                        n_tree_mu, maxdepth_mu, eta_mu,
                        n_tree_p, cp_p, maxdepth_p, lr_p,
                        structure, trace, patience) {
  # X, Y are the learning data;
  # the valid_rows are used to early stop the algorithm;
  # mu0, p0, sigma0 are the initial parameters;
  # K is the number of components;
  # M0 is the number of EB iterations;
  # n_tree_mu is the number of trees in mu boosting in each B step;
  # maxdepth_mu, eta_mu are the parameters in trees for mu
  # n_tree_p is the number of trees in mixing probabilities boosting in each B step.
  # cp_p, maxdepth_p, lr_p are the parameters in trees for mixing probs;
  # structure can be "both", "p", "mu";
  # patience is the number of non-improvements in validation loss before early stop

  # initialization
  train_rows<-(1:nrow(X))[-valid_rows]
  Xtrain<-X[train_rows,]; Ytrain<-Y[train_rows]
  Xvalid<-X[valid_rows,]; Yvalid<-Y[valid_rows]
  par_mat <- dat_Norm(K = K, n = nrow(X))
  mu_models<-NULL
  par_mat[,vname2("mu",K)]=mu0
  par_mat[,vname2("p",K)]=p0
  par_mat[,vname2("sigma",K)]=sigma0
  for (k in 1:K) mu_models[[k]]<-list()
  train_loss <- NULL; valid_loss <- NULL
  train_impT <- NULL; train_impT[1]<-NULL
  valid_impT <- NULL; valid_impT[1]<-NULL
  valid_impT2<- NULL; valid_impT2[1:patience]<-NULL
  p_models<-list()

  for (m in 1:M0) {
    # expectation of latent variable
    # also the weights for mean and response for mixing probabilities
    sum_fz<-rep(0,nrow(par_mat))
    for (k in 1:K){
      par_mat[,names(par_mat)==vname("f",k)]<-dnorm(Y, get(vname("mu",k), par_mat), get(vname("sigma",k),par_mat))
      sum_fz<-sum_fz+get(vname("p",k), par_mat)*get(vname("f",k), par_mat)
    }
    par_mat[,vname2("z",K)] <- par_mat[,vname2("p",K)]*par_mat[,vname2("f",K)]/sum_fz
    par_mat[,vname2("weight",K)] <- par_mat[,vname2("z",K)]

    # boosting for mu
    param<-list(max_depth=maxdepth_mu, eta =eta_mu, objective="reg:squarederror")
    for (k in 1:K){
      if (structure == "mu"| structure == "both") {
        dat_xgb<-xgb.DMatrix(data=as.matrix(X), label = Y, weight=par_mat[,vname("weight",k)], base_margin=par_mat[,vname("mu",k)])
      }
      if (structure == "p"){
        dat_xgb<-xgb.DMatrix(data=as.matrix(rep(1,length(Y))), label = Y, weight=par_mat[,vname("weight",k)],base_margin=par_mat[,vname("mu",k)])
      }
      dtrain<-xgboost::slice(dat_xgb,train_rows)
      dvalid<-xgboost::slice(dat_xgb,valid_rows)
      watchlist=list(train=dtrain, eval= dvalid)
      mu_bst <- xgb.train(param, dtrain, nrounds=n_tree_mu, verbose = 0, watchlist, early_stopping_rounds = NULL)
      par_mat[, vname("mu",k)]<-predict(mu_bst, newdata = dat_xgb)
      mu_models[[k]][[m]]<-mu_bst
    }

    # boosting for p
    if (structure == "p" |structure == "both") {Xtrain_p<-Xtrain; Xvalid_p<-Xvalid}
    if (structure == "mu"){Xtrain_p<-data.frame(rep(1,nrow(Xtrain))); Xvalid_p<-data.frame(rep(1,nrow(Xvalid)))}
    p_bst <- BST(Xtrain = Xtrain_p, Ytrain=par_mat[train_rows,vname2("z",K)], Xval=Xvalid_p, Yval=par_mat[valid_rows,vname2("z",K)], train_init=par_mat[train_rows, vname2("p",K)], valid_init = par_mat[valid_rows, vname2("p",K)], M=n_tree_p, cp=cp_p, maxdepth = maxdepth_p, lr=lr_p, trace=F, patience = 1)
    par_mat[train_rows,vname2("p",K)] <- p_bst$train_bst[,vname2("p",K)]
    par_mat[valid_rows,vname2("p",K)] <- p_bst$valid_bst[,vname2("p",K)]
    p_models[[m]]<-p_bst

    # sigma
    for (k in 1:K){
      sigma_est<-sqrt(sum(get(vname("z",k), par_mat[train_rows,]) * (Ytrain - get(vname("mu",k), par_mat[train_rows,])) ^ 2) / sum(get(vname("z",k), par_mat[train_rows,])))
      par_mat[,vname("sigma",k)]<-sigma_est
    }

    # outer train loss and validation loss
    train_loss[m] <- neg_ll(Ytrain, par_mat[train_rows, vname2("mu",K)], par_mat[train_rows, vname2("sigma",K)], par_mat[train_rows, vname2("p",K)])
    valid_loss[m] <- neg_ll(Yvalid, par_mat[valid_rows, vname2("mu",K)], par_mat[valid_rows, vname2("sigma",K)], par_mat[valid_rows, vname2("p",K)])
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
  par_mat[, vname2("plinear",K)] <- PtoF(par_mat[, vname2("p",K)])
  list(
    mu0 = mu0,
    p0 = p0,
    structure = structure,
    train_loss = train_loss,
    valid_loss = valid_loss,
    par_mat = par_mat,
    mu_models = mu_models,
    p_models = p_models,
    sigma = apply(par_mat[,vname2("sigma",K)], 2, unique)
  )
}

EB_gaussian_aug<- function(X, Y, valid_rows, mu0, p0, sigma0, K, M0,
                        n_tree_mu, maxdepth_mu, eta_mu,
                        n_tree_p, cp_p, maxdepth_p, lr_p,
                        structure, trace, patience) {
  # X, Y are the learning data;
  # the valid_rows are used to early stop the algorithm;
  # mu0, p0, sigma0 are the initial parameters;
  # K is the number of components;
  # M0 is the number of EB iterations;
  # n_tree_mu is the number of trees in mu boosting in each B step;
  # maxdepth_mu, eta_mu are the parameters in trees for mu
  # n_tree_p is the number of trees in mixing probabilities boosting in each B step.
  # cp_p, maxdepth_p, lr_p are the parameters in trees for mixing probs;
  # structure can be "both", "p", "mu";
  # patience is the number of non-improvements in validation loss before early stop
  train_rows<-(1:nrow(X))[-valid_rows]
  colnames(mu0)<-NULL; colnames(p0)<-NULL; colnames(sigma0)<-NULL
  dat_aug<-dat_augmentation(data.frame(X,Y,mu=mu0,p=p0,sigma=sigma0),K)
  par_mat <- dat_Norm(K = K, n = nrow(dat_aug))
  # initialization
  train_rows1<-which(dat_aug$R0%in%train_rows)
  valid_rows1<-which(dat_aug$R0%in%valid_rows)
  X<-dat_aug[,vname2("X",6)];Y<-dat_aug$Y;ZZ<-dat_aug[,vname2("ZZ",K)]
  Xtrain<-X[train_rows1,]; Ytrain<-Y[train_rows1]; ZZtrain<-ZZ[train_rows1,]
  Xvalid<-X[valid_rows1,]; Yvalid<-Y[valid_rows1]; ZZvalid<-ZZ[valid_rows1,]
  mu_models<-NULL
  par_mat[,vname2("mu",K)]=dat_aug[,vname2("mu.",K)]
  par_mat[,vname2("p",K)]=dat_aug[,vname2("p.",K)]
  par_mat[,vname2("sigma",K)]=dat_aug[,vname2("sigma.",K)]
  rm(dat_aug)
  for (k in 1:K) mu_models[[k]]<-list()
  train_loss <- NULL; valid_loss <- NULL
  train_impT <- NULL; train_impT[1]<-NULL
  valid_impT <- NULL; valid_impT[1]<-NULL
  valid_impT2<- NULL; valid_impT2[1:patience]<-NULL
  p_models<-list()

  for (m in 1:M0) {
    # update the weights
    sum_fz<-rep(0,nrow(par_mat))
    for (k in 1:K){
      par_mat[,vname("f",k)]<-dnorm(Y, par_mat[,vname("mu",k)], par_mat[,vname("sigma",k)])
      sum_fz<-sum_fz + par_mat[,vname("p",k)]*par_mat[,vname("f",k)]
    }
    par_mat[,vname2("weight",K)] <- (par_mat[,vname2("p",K)]*par_mat[,vname2("f",K)])/sum_fz
    par_mat[,vname2("weight",K)] <- par_mat[,vname2("weight",K)]*ZZ
    par_mat$weight<-apply(par_mat[,vname2("weight",K)],1,sum)
    par_mat[valid_rows1,vname2("weight",K)]<-1
    par_mat$weight[valid_rows1]<-1

    # boosting for mu
    param<-list(max_depth=maxdepth_mu, eta =eta_mu, objective="reg:squarederror")
    for (k in 1:K){
      if (structure == "mu"| structure == "both") {
        dat_xgb<-xgb.DMatrix(data=as.matrix(X), label = Y, weight=par_mat[,vname("weight",k)], base_margin=par_mat[,vname("mu",k)])
      }
      if (structure == "p"){
        dat_xgb<-xgb.DMatrix(data=as.matrix(rep(1,length(Y))), label = Y, weight=par_mat[,vname("weight",k)], base_margin=par_mat[,vname("mu",k)])
      }
      dtrain<-xgboost::slice(dat_xgb,train_rows1)
      dvalid<-xgboost::slice(dat_xgb,valid_rows1)
      watchlist=list(train=dtrain, eval= dvalid)
      mu_bst <- xgb.train(param, dtrain, nrounds=n_tree_mu, verbose = 0, watchlist, early_stopping_rounds = NULL)
      par_mat[,vname("mu",k)]<-predict(mu_bst, newdata = dat_xgb)
      mu_models[[k]][[m]]<-mu_bst
    }

    # boosting for p
    if (structure == "p" |structure == "both") {Xtrain_p<-Xtrain; Xvalid_p<-Xvalid}
    if (structure == "mu"){Xtrain_p<-data.frame(rep(1,nrow(Xtrain))); Xvalid_p<-data.frame(rep(1,nrow(Xvalid)))}
    p_bst <- BST(Xtrain = Xtrain_p, ZZtrain, Xval=Xvalid_p, ZZvalid, train_init=par_mat[train_rows1, vname2("p",K)], valid_init = par_mat[valid_rows1, vname2("p",K)], train_weight=par_mat$weight[train_rows1], M=n_tree_p, cp=cp_p, maxdepth = maxdepth_p, lr=lr_p, trace=0, patience = 1)
    par_mat[train_rows1,vname2("p",K)] <- p_bst$train_bst[,vname2("p",K)]
    par_mat[valid_rows1,vname2("p",K)] <- p_bst$valid_bst[,vname2("p",K)]
    p_models[[m]]<-p_bst

    # sigma
    for (k in 1:K){
      sigma_est<-sqrt(sum(par_mat[train_rows1,vname("weight",k)]*(Ytrain - par_mat[train_rows1,vname("mu",k)]) ^ 2) / sum(par_mat[train_rows1,vname("weight",k)]))
      par_mat[,vname("sigma",k)]<-sigma_est
    }

    # outer train loss and validation loss
    train_loss[m] <- neg_ll(Ytrain, par_mat[train_rows1, vname2("mu",K)], par_mat[train_rows1, vname2("sigma",K)], par_mat[train_rows1, vname2("p",K)], w=par_mat$weight[train_rows1])
    valid_loss[m] <- neg_ll(Yvalid, par_mat[valid_rows1, vname2("mu",K)], par_mat[valid_rows1, vname2("sigma",K)], par_mat[valid_rows1, vname2("p",K)])
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
par_mat[, vname2("plinear",K)] <- PtoF(par_mat[, vname2("p",K)])
par_mat <- par_mat[1:n,]

  list(
    mu0 = mu0,
    p0 = p0,
    structure = structure,
    train_loss = train_loss,
    valid_loss = valid_loss,
    par_mat = par_mat,
    mu_models = mu_models,
    p_models = p_models,
    sigma = apply(par_mat[,vname2("sigma",K)], 2, unique)
  )
}

# Data genaration
n = 10000
ntest = 2000
dat<-sim_gaussian(n,seed=1)
dat_test<-sim_gaussian(ntest,seed=7)
dat_test$true_mu<-apply((dat_test[,vname2("P",3)]*dat_test[,vname2("MU",3)]),1,sum)



