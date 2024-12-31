# load the package used in the parallel computing
library("parallel")
cores <- detectCores(logical = FALSE)
c1 <- makeCluster(cores)
clusterEvalQ(c1,library(partykit))
clusterEvalQ(c1,library(rpart))
negLL<-function(y,p,w=1){
  # the negative log-likelihood for the multinomial distribution.
  # y and p are of dimension of n*K.
  if (is.vector(y)) (w*sum(y*log(p)))
  else -mean(w*diag(as.matrix(y)%*%t(log(as.matrix(p)))))
}

FtoP<-function(FF){
  # FF has dimension of n*K.
  # From linear predictor (FF) to probs (P).
  # Different FF may lead to same P.
  # One FF lead to a P.
  EF<-exp(FF)
  SEF<-apply(EF,1,sum)
  P<-EF/SEF
  P
}

PtoF<-function(P){
  # PP has dimension of n*K.
  # From probs (P) to linear predictors (FF).
  # A P may lead to different FF, depending on the constant added to log P.
  # The constant here is the neg average of log P such that sum of FF is zero.
  logP<-log(P)
  logPa<-apply(logP,1,mean)
  FF<-logP-logPa
  FF
}

FtoP2<-function(FF){
  exp(FF)/(1+exp(FF))
}

PtoF2<-function(P){
  log(P/(1-P))
}

vname<-function(vname,k){
  # create "vnamek".
  paste(vname, k, sep="")
  }
vname2<-function(vname, K){
  # create "vname1",..."vnameK".
  temp<-paste(vname, 1, sep="")
  for (k in 2:K) temp<-c(temp,paste(vname, k, sep=""))
  temp
}

node_pre<-function(yw,K){
  # optimal updates of node in multinomial boosting.
  # see equation 32 in Friedman (2001) with modification of weights.
  # yw of size n*K, the first column contains gradients (y),
  # the second column contains weights (w).
  y<-yw[,1];w<-yw[,2]
  (K-1)/K * sum(w*y)/sum(w*abs(y)*(1-abs(y)))
}

dat_BST<-function(K,n){
  # create a dataset with suitable format ready for boosting. The columns are:
  # First:  gradients (observed - prediction, y);
  # Second: linear predictor (weaker learner, approximated gradient, f);
  # Third:  prediction of probabilities (p)
  dat_bst<-data.frame(matrix(NA,ncol=3*K,nrow=n))
  names(dat_bst)<-c(vname2("y",K),vname2("f",K),vname2("p",K))
  dat_bst
}

fit_tree <- function(Xtrain, train_bst, Xval, val_bst, train_weight,
                     maxdepth, cp, lr, K, k){
  # this function is used to perform one gradient boosting step,
  # i.e., training one weaker learner of tree for one class k.
  # Xtrain, Xval contains all the covariates;
  # train_bst, val_bst contains gradients(y), linear predictor(f), probs (p);
  # train_bst, val_bst are the objects from dat_BST function;
  # maxdepth is the depth of the weaker learners;
  # cp is the complexity parameter;
  # lr is the learning rate in the boosting;
  # K is the total number of components.
  set.seed(1)
  # weaker learner of tree for gradients (y)
  tree_rpart<-rpart(train_bst[,vname("y",k)]~., method = "anova", data = Xtrain,
               weights=train_weight,
               control = rpart.control(minbucket = 50, cp = cp, maxcompete = 4,
                                       maxsurrogate = 0, usesurrogate = 2,
                                       xval = 0, surrogatestyle = 0,
                                       maxdepth = maxdepth))
  tree_party<-as.party(tree_rpart) # node index in as.party is different from rpart
  train_node<-fitted(tree_party) # node index and prediction for gradients (y)
  names(train_node) = c("node", "response", "weights")
  # equation 32 in Friedman（2001）
  f_pred <-tapply(X=train_node[,c("response","weights")],
                  INDEX= list(train_node$node), FUN = node_pre, K = K)
  f_pred <-data.frame(f = as.vector(f_pred), node = as.numeric(names(f_pred)))
  valid_node <-predict(tree_party, newdata = Xval, type = "node")
  # update the linear predictor as f_{m-1}+ lr * f_{m}
  train_bst[, vname("f",k)] <- train_bst[, vname("f",k)] +
    lr * f_pred$f[match(train_node$node, f_pred$node)]
  val_bst[, vname("f",k)] <- val_bst[, vname("f",k)] +
    lr * f_pred$f[match(valid_node, f_pred$node)]
  return(list(train_bst=train_bst,val_bst=val_bst,
              tree_party=tree_party,tree_rpart=tree_rpart, f_pred=f_pred))
  # tree_party is the party tree used for prediction.
  # tree_rpart is the rpart tree used for relative importance.
  # tree_party and tree_rpart have different node index.
}

BST <- function(Xtrain, Ytrain, Xval, Yval, train_init, valid_init,
                train_weight=NULL, M, cp, maxdepth, lr, trace, patience) {
  # Xtrain, Xval: covariates. Ytrain, Yval: multinomial response.
  # train_init, valid_init: initial values for probs.
  # M: number of iterations.
  # cp: complex parameter.
  # maxdepth: depth of weak learners of tree.
  # lr: learning rate in the boosting.
  # trace: show loss or not.
  # patience: early stop if validation loss is not decreased.
  set.seed(1)
  K <- ncol(Ytrain)
  n <- nrow(Ytrain)
  n_val <- nrow(Yval)
  train_bst <- dat_BST(K, n)
  val_bst <- dat_BST(K, n_val)
  if (is.null(train_weight)) train_weight<-rep(1,n)
  # initialization
  train_bst[, vname2("p",K)]<-train_init
  train_bst[, vname2("f",K)]<-PtoF(train_bst[, vname2("p",K)])
  val_bst[, vname2("p",K)]<-valid_init
  val_bst[, vname2("f",K)]<-PtoF(val_bst[, vname2("p",K)])
  Train_loss <- NULL; Valid_loss <- NULL
  Train_impT<-NULL; Valid_impT<-NULL; Valid_impT2<-NULL
  Train_impT[1]<-NULL;Valid_impT[1]<-NULL;Valid_impT2[1:patience]<-NULL
  Tree_party <- list(); Tree_rpart<-list(); F_pred<-list()
  clusterExport(c1, varlist = list("node_pre","vname","vname2"))
  # boosting
  for (m in 1:M) {
    train_bst[, vname2("y",K)] <- Ytrain - train_bst[, vname2("p",K)]
    tree_party<-list()
    tree_rpart<-list()
    f_pred<-list()
    res<- clusterApplyLB(c1, x=1:K, fun=fit_tree, Xtrain=Xtrain, train_bst=train_bst, Xval=Xval, val_bst=val_bst, train_weight=train_weight, maxdepth=maxdepth, cp=cp, lr=lr, K=K)
    for (k in 1:K){
      tree_party[[k]]=res[[k]][["tree_party"]]
      tree_rpart[[k]]=res[[k]][["tree_rpart"]]
      f_pred[[k]]=res[[k]][["f_pred"]]
      train_bst[,vname("f",k)]<- res[[k]][["train_bst"]][,vname("f",k)]
      val_bst[,vname("f",k)]<- res[[k]][["val_bst"]][,vname("f",k)]
    }
    train_bst[, vname2("p",K)] <- FtoP(train_bst[, vname2("f",K)])
    val_bst[, vname2("p",K)] <- FtoP(val_bst[, vname2("f",k)])
    Train_loss[m] <- negLL(Ytrain, train_bst[, vname2("p",K)], w=train_weight)
    Valid_loss[m] <- negLL(Yval, val_bst[, vname2("p",K)])
    Tree_party[[m]]<- tree_party # K trees
    Tree_rpart[[m]]<- tree_rpart
    F_pred[[m]]<- f_pred
    # Tree_party is the party version of the rpart tree
    # Tree_party is used for predicition
    # Tree_rpart is the rpart tree with different node index from party tree
    # Tree_rpart is used for variable importance
    # F_pred is the weak learner
    if (m>1){
      Train_impT[m]<- (Train_loss[m-1]-Train_loss[m]>0)
      Valid_impT[m]<- (Valid_loss[m-1]-Valid_loss[m]>0)
    }
    if(m>patience) {
      Valid_impT2[m]<- (Valid_loss[m-patience]-Valid_loss[m]>0)
    }
    if(trace==T){
      # iter-m: train_ls (whether improved compared to the previous iteration),
      # valid_ls (whether improved compared to the previous iteration/
      # whether improved compared to the iteration patience ahead)
      print(paste("iter-",m,":","train_ls",round(Train_loss[m],4),"(",Train_impT[m],"),","valid_ls",round(Valid_loss[m],4),"(",Valid_impT[m],"/",Valid_impT2[m],")",sep=""))
    }
    if(m>patience) if(Valid_impT2[m]==0) break
  }
  list(Train_loss=Train_loss,
       Valid_loss=Valid_loss,
       train_bst=train_bst,
       valid_bst=val_bst,
       Tree_party=Tree_party,
       Tree_rpart=Tree_rpart,
       F_pred=F_pred,
       lr=lr)
}

predict_BST <- function(X, BST_fit, init, type) {
  # X: the covariates matrix;
  # BST_fit: the output from the function BST;
  # init: the initial probabilities;
  # type: response ("response") or the linear predictor ("link")
  Tree_party<-BST_fit$Tree_party
  F_pred<-BST_fit$F_pred
  lr <- BST_fit$lr
  M_best<-which.min(BST_fit$Valid_loss) # number of weak learner
  K <- unique(sapply(Tree_party, length)) # number of components
  n <- nrow(X) # number of samples
  dat_bst <- dat_BST(K, n)
  dat_bst[,vname2("p",K)]<-init
  dat_bst[,vname2("f",K)]<-PtoF(dat_bst[,vname2("p",K)])
  for (m in 1:M_best) {
    for (k in 1:K) {
      fit_party <- Tree_party[[m]][[k]]
      f_pred <- F_pred[[m]][[k]]
      dat_node <- predict(fit_party, newdata = X, type = "node")
      dat_bst[, vname("f",k)] <- dat_bst[, vname("f",k)] +
        lr * f_pred$f[match(dat_node,f_pred$node)]
    }
  }
  dat_bst[, vname2("p",K)] <- FtoP(dat_bst[, vname2("f",K)])
  if (type=="response") return(data.frame(dat_bst[,vname2("p",K)]))
  if (type=="link") return(data.frame(dat_bst[,vname2("f",K)]))
}

