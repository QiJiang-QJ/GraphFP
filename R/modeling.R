library('deSolve')
library('dplyr')
library('igraph')

addMin <- 1e-14

log_one <- function(a){
  label_0_a <- which(a==0)
  for(i in label_0_a){
    a[i] <- a[i]+addMin
  }
  new_loga <- log(a) %>% as.numeric()
  return(new_loga)
}

log_fraction <- function(a,b){
  label_0_a <- which(a==0)
  label_0_b <- which(b==0)
  for(i in label_0_a){
    a[i] <- a[i]+addMin
  }
  for(i in label_0_b){
    b[i] <- b[i]+addMin
  }
  c <- log(a/b) %>% as.numeric()
  return(c)
}

KL_div <- function(a,b){
  a <- as.numeric(a)
  b <- as.numeric(b)
  label_0_a <- which(a==0)
  label_0_b <- which(b==0)
  for(i in label_0_a){
    a[i] <- a[i]+addMin
  }
  for(i in label_0_b){
    b[i] <- b[i]+addMin
  }
  div <- 0
  for (i in 1:length(a)) {
    div <- div + a[i]*log(a[i]/b[i])
  }
  return(div)
}

KL_divs <- function(Pred,Actual){
  T <- length(Pred)
  KLs <- c(0)
  for(i in 2:f){
    KLs <- c(KLs,KL_div(Pred[i,],Actual[i,]))
  }
  return(KLs)
}

Time_array_array <- function(a,b){
  NN_row <- length(a)
  NN_col <- length(b)
  ab <- matrix(0,ncol = NN_col,nrow = NN_row)
  for (i in 1:NN_row) {
    for (j in 1:NN_col) {
      ab[i,j]<- a[i]*b[j]
    }
  }
  return(ab)
}

F <- function(Phi,W,p,beta){
  bigf <- (Phi + as.numeric(W %*% p) + beta*(log_one(p))) %>% as.numeric()
  return(bigf)
}

du_dt <- function(t,state,pars){
  u <- state %>% as.numeric()
  df_dp <- pars$df_dp
  dL_dp <- pars$dL_dp %>% as.numeric()
  du <- as.numeric(u %*% df_dp) + dL_dp
  return(list(du))
}

derive_a_dKL_dp <- function(Phi,W,pred_P,probs){
  a_dKL_dp <- list()
  array_1 <- rep(1,length(Phi))
  a_dKL_dp[[1]] <- rep(0,length(Phi))
  for(i in 2:dim(probs)[1]){
    actual_p <- probs[i,] %>% as.numeric()
    pred_prob <- pred_P[i,]
    a_dKL_dp[[i]] <- log_fraction(pred_prob,actual_p) + array_1
  }
  return(a_dKL_dp)
}

derive_dv <- function(prob,Phi,W,beta){
  Fprob <- F(Phi,W,prob,beta) %>% as.numeric()
  aaF <- (B %*% Fprob) %>% as.numeric() %>% abs()
  return(aaF)
}

derive_dp <- function(prob,Phi,W,beta){
  Fprob <- F(Phi,W,prob,beta) %>% as.numeric()
  aaF <- ((B %*% Fprob))%>% as.numeric()
  gij <- c()
  for (l in 1:num_E) {
    vertex_two <- Edgelist[l,]
    if(aaF[l]>0){
      index_output <- which(CelltypeNames==vertex_two[2])
      gij[l] <- prob[index_output]
    }else{
      index_output <- which(CelltypeNames==vertex_two[1])
      gij[l] <- prob[index_output]
    }
  }
  dp <- ((-t(B)) %*% (gij*aaF)) %>% as.numeric()
  return(dp)
}

derive_u <- function(Phi,W,pred_P,beta){
  u <- list()
  u[[length(Timepoints)]] <- rep(0,length(Phi))
  for(i in 1:(length(Timepoints)-1)){
    j <- length(Timepoints)-i
    pp <- pred_P[j+1,]
    Ap <- (A %*% (pp/2)) %>% as.numeric()
    Bp <- (B %*% (pp/2)) %>% as.numeric()
    Fp <- F(Phi,W,pp,beta)
    BF <- (B %*% Fp) %>% as.numeric()
    inverse_p <- 1/pp
    Inf_label <- which(inverse_p==Inf)
    inverse_p[Inf_label] <- 0
    BW <- B %*% W + B * as.numeric(beta*inverse_p)
    U <- (Ap*BW + BF*(A/2)) + (sign(BF)*(Bp*BW + BF*(B/2)))
    u_df_dp <- (t(B)) %*% U
    u_dL_dp <- (as.numeric(BF) %*% U)+(( (Ap*BF)  + (Bp*(abs(BF)))  ) %*% (BW))
    u_pars <- list(df_dp=u_df_dp,dL_dp=u_dL_dp)
    u_times <- c(Timepoints[j+1],Timepoints[j])
    out_u <- ode(y=u[[j+1]],func = du_dt,times = u_times,parms = u_pars)
    uu <- out_u[2,] %>% as.numeric()
    uu <- uu[-1]
    u[[j]] <- uu %>% as.numeric()
  }
  return(u)
}

derive_dPhi_dW <- function(p,Phi,W,u1,a_dKL_dp2,lambda_i,beta){
  list_d <- list()
  Fp <- F(Phi,W,p,beta)
  D1 <- (A1 %*% (p/2)) %>% as.numeric()
  D2 <- (A2 %*% (p/2)) %>% as.numeric()
  A2F <- (A2 %*% Fp) %>% as.numeric()
  D3 <- sign(A2F)
  ADA <- (t(A2)) %*% (D1 * A2)
  ADDA <- (t(A2)) %*% (as.numeric(D2 * D3) * A2)
  A2Fp1 <- Time_array_array(A2F,p)
  A2Fp2 <- Time_array_array(abs(A2F),p)
  ula <- (u1 -lambda_i*a_dKL_dp2) %>% as.numeric()
  dPhi <- as.numeric(as.numeric(A2F * D1) %*% A2)+as.numeric(as.numeric(abs(A2F) * D2) %*% A2)+as.numeric(ula %*% (ADA + ADDA))
  
  ADDAFp1 <- t(A2) %*% (D1 * A2Fp1)
  ADDAFp2 <- t(A2) %*% (D2 * A2Fp2)
  dW1 <- (ADDAFp1 + ADDAFp2)
  dW2 <- matrix(0,ncol = length(Phi),nrow = length(Phi))
  for (k in 1:length(Phi)) {
    ADAp <- Time_array_array(as.numeric(ADA[k,]),p)
    ADDAp <- Time_array_array(as.numeric(ADDA[k,]),p)
    dW2 <- dW2 + ula[k]*((ADAp + ADDAp))
  }
  dW <- (dW1+dW2)
  list_d <- list(dPhi=dPhi,dW=dW)
  return(list_d)
}


sgn_Phi <- function(X){
  n <- length(X)
  Y <- rep(0,n)
  for (i in 1:n) {
      if(X[i]>0){
        Y[i]=1
      }else if(X[i]<0){
        Y[i]=-1
      }else{
        Y[i]=0
      }
  }
  return(Y)
}

sgn_W <- function(X){
  n <- dim(X)[1]
  Y <- matrix(0,ncol=n,nrow=n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(X[i,j]>0){
        Y[i,j]=1
      }else if(X[i,j]<0){
        Y[i,j]=-1
      }else{
        Y[i,j]=0
      }
    }
  }
  return(Y)
}


Estimate_THETA <- function(lambda,beta,initial_Phi,initial_W,Graph,probs,Manual_Critical=2,Integral_step=10,Max_iter=2000,Par_step=1e-02,tol=1e-09){
  T1=Sys.time()
  
  Phi <- initial_Phi
  W <- initial_W
  Edgelist <- get.edgelist(Graph)
  critical <- (1/(max(degree(Graph))-1))/Manual_Critical
  delta_ti <- 1/Integral_step
  Par_step <- Par_step/sum(lambda)*(f-1)
  tol=tol*sum(lambda)/(f-1)
  
  
  derive_loss <- function(Phi,W,Timepoints){
    derive_dloss <- function(prob){
      Fp <- F(Phi,W,prob,beta)
      BF <- ((B %*% Fp)) %>% as.numeric()
      gij <- c()
      for (l in 1:(dim(A)[1])) {
        vertex_two <- sort(Edgelist[l,])
        if(BF[l]>0){
          label_output <- vertex_two[2]
          index_output <- which(CelltypeNames==label_output)
          gij[l] <- prob[index_output]
        }else{
          label_output <- vertex_two[1]
          index_output <- which(CelltypeNames==label_output)
          gij[l] <- prob[index_output]
        }
      }
      dloss <- (as.numeric(BF) %*% as.numeric(gij*BF)) %>% as.numeric()
      return(dloss)
    }
    Loss <- 0
    p <- probs[1,]
    Pred_P <- probs
    for(i in 1:(length(Timepoints)-1)){
      MM <- Integral_step*(Timepoints[i+1] - Timepoints[i])
      delta_t <- 1/Integral_step
      for (j in 1:MM) {
        count_time <- 0
        while((count_time)< delta_t){
          rest_time <- delta_t - count_time
          dv <- derive_dv(p,Phi,W,beta)
          new_dv <- (rest_time)*dv
          if(length(which(new_dv > critical))==0){
            count_time <- count_time + rest_time
            dloss <- derive_dloss(p)
            Loss <- Loss + rest_time*dloss
            p <- p + (rest_time*derive_dp(p,Phi,W,beta))
          }else{
            new_step <- (rest_time/(max(dv)))*critical
            count_time <- count_time + new_step
            dloss <- derive_dloss(p)
            Loss <- Loss + new_step*dloss
            p <- p + (new_step*derive_dp(p,Phi,W,beta))
          }
        }
      }
      actual_p2 <- probs[i+1,] %>% as.numeric()
      Loss <- Loss + lambda[i]*KL_div(p,actual_p2)
      Pred_P[i+1,] <- p
    }
    LLoss <- list(Loss=Loss,Pred_P=Pred_P)
    return(LLoss)
  }
  
  #### initial pred_P & U1 & A_dKL_dp2 & Loss
  LLoss1 <- derive_loss(Phi,W,Timepoints)
  Loss1 <- LLoss1$Loss
  pred_P <- LLoss1$Pred_P
  U1 <- derive_u(Phi,W,pred_P,beta)
  A_dKL_dp2 <- derive_a_dKL_dp(Phi,W,pred_P,probs)
  KLs <- KL_divs(pred_P,probs)
  Loss_int <- Loss1-((lambda %*% KLs[2:f])/(f-1))
  
  Loss_KLs <- list()
  Loss <- c(Loss1)
  Loss_KLs[[1]] <- KLs
  Loss_intergral <- c(Loss_int)
  
  #### gradient descent to estimate Phi & W
  delta_loss <- 100
  num_iter <- 1
  while(delta_loss > tol & num_iter < Max_iter){
    num_iter <- num_iter+1
    delta_Phi <- rep(0,N)
    delta_W <- matrix(0,ncol=N,nrow=N)
    pars_step <- Par_step
    p <- probs[1,]
    #### compute delta_Phi & delta_W
    for(t in 1:(f-1)){
      M <- (Timepoints[t+1] - Timepoints[t]) %/% delta_ti
      last_ti <- (Timepoints[t+1] - Timepoints[t]) %% delta_ti
      for (j in 1:M) {
        ccount_time <- 0
        while((ccount_time)<delta_ti){
          rrest_time <- delta_ti - ccount_time
          ddv <- derive_dv(p,Phi,W,beta)
          new_ddv <- ((rrest_time)*ddv)
          if(length(which(new_ddv > critical))==0){
            ccount_time <- ccount_time + rrest_time
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[t],beta)
            delta_Phi <- delta_Phi + (rrest_time*(ddd$dPhi))
            delta_W <- delta_W +(rrest_time*(ddd$dW))
            p <- p + (rrest_time*derive_dp(p,Phi,W,beta))
          }else{
            new_step <- (rrest_time/(max(ddv)))*critical
            ccount_time <- ccount_time + new_step
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[t],beta)
            delta_Phi <- delta_Phi + (new_step*(ddd$dPhi))
            delta_W <- delta_W +(new_step*(ddd$dW))
            p <- p + (new_step*derive_dp(p,Phi,W,beta))
          }
        }
      }
      if(last_ti >0){
        ccount_time <- 0
        while((ccount_time)<last_ti){
          rrest_time <- last_ti - ccount_time
          ddv <- derive_dv(p,Phi,W,beta)
          new_ddv <- ((rrest_time)*ddv)
          if(length(which(new_ddv > critical))==0){
            ccount_time <- ccount_time + rrest_time
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[f-1],beta)
            delta_Phi <- delta_Phi + (rrest_time*(ddd$dPhi))
            delta_W <- delta_W +(rrest_time*(ddd$dW))
            p <- p + (rrest_time*derive_dp(p,Phi,W,beta))
          }else{
            new_step <- (rrest_time/(max(ddv)))*critical
            ccount_time <- ccount_time + new_step
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[f-1],beta)
            delta_Phi <- delta_Phi + (new_step*(ddd$dPhi))
            delta_W <- delta_W +(new_step*(ddd$dW))
            p <- p + (new_step*derive_dp(p,Phi,W,beta))
          }
        }
      }
    }
    new_Phi <- Phi - delta_Phi*pars_step
    new_W <- W - delta_W*pars_step
    new_LLoss <- derive_loss(new_Phi,new_W,Timepoints)
    new_Loss <- new_LLoss$Loss
    if(new_Loss >= Loss1){
      repeat{
        pars_step <- pars_step/2
        print(pars_step)
        new_Phi <- Phi - delta_Phi*pars_step
        new_W <- W - delta_W*pars_step
        new_LLoss <- derive_loss(new_Phi,new_W,Timepoints)
        new_Loss <- new_LLoss$Loss
        if(new_Loss < Loss1 | (pars_step*max(abs(delta_Phi),abs(delta_W)))<1e-04){
          if(new_Loss < Loss1){
            Phi <- new_Phi
            W <- new_W
            #### new pred_P & U1 & A_dKL_dp2 & Loss1
            pred_P <- new_LLoss$Pred_P
            U1 <- derive_u(Phi,W,pred_P,beta)
            A_dKL_dp2 <- derive_a_dKL_dp(Phi,W,pred_P,probs)
            KLs <- KL_divs(pred_P,probs)
            Loss_int <- Loss1-((lambda %*% KLs[2:f])/(f-1))
            delta_loss <- Loss1 - new_Loss
            Loss1 <- new_Loss
            
            Loss <- c(Loss,new_Loss)
            Loss_KLs[[num_iter]] <- KLs
            Loss_intergral <- c(Loss_intergral,Loss_int)
          }else{
            delta_loss <- 0
          }
          break
        }
      }
    }else{
      Phi <- new_Phi
      W <- new_W
      #### new pred_P & U1 & A_dKL_dp2 & Loss
      pred_P <- new_LLoss$Pred_P
      U1 <- derive_u(Phi,W,pred_P,beta)
      A_dKL_dp2 <- derive_a_dKL_dp(Phi,W,pred_P,probs)
      KLs <- KL_divs(pred_P,probs)
      Loss_int <- Loss1-((lambda %*% KLs[2:f])/(f-1))
      delta_loss <- Loss1 - new_Loss
      Loss1 <- new_Loss
      
      Loss <- c(Loss,new_Loss)
      Loss_KLs[[num_iter]] <- KLs
      Loss_intergral <- c(Loss_intergral,Loss_int)
    }
  }
  parameter <- list(Phi=Phi,W=W)
  Intermediate_Process <- list(Loss=Loss,Loss_intergral=Loss_intergral,Loss_KLs=Loss_KLs)
  T2=Sys.time()
  list_result_ALL <- list(parameter=parameter,Intermediate_Process=Intermediate_Process,Pred_probs=pred_P,Actual_probs=probs,Time=(T2-T1))
  return(list_result_ALL)
}


Estimate_THETA_regularization <- function(lambda,beta,initial_Phi,initial_W,Graph,probs,Manual_Critical=2,Integral_step=10,Max_iter=2000,Par_step=1e-02,tol=1e-09,alpha=5){
  T1=Sys.time()
  
  Phi <- initial_Phi
  W <- initial_W
  Edgelist <- get.edgelist(Graph)
  critical <- (1/(max(degree(Graph))-1))/Manual_Critical
  delta_ti <- 1/Integral_step
  Par_step <- Par_step/sum(lambda)*(f-1)
  tol=tol*sum(lambda)/(f-1)
  
  
  derive_loss <- function(Phi,W,Timepoints){
    derive_dloss <- function(prob){
      Fp <- F(Phi,W,prob,beta)
      BF <- ((B %*% Fp)) %>% as.numeric()
      gij <- c()
      for (l in 1:(dim(A)[1])) {
        vertex_two <- sort(Edgelist[l,])
        if(BF[l]>0){
          label_output <- vertex_two[2]
          index_output <- which(CelltypeNames==label_output)
          gij[l] <- prob[index_output]
        }else{
          label_output <- vertex_two[1]
          index_output <- which(CelltypeNames==label_output)
          gij[l] <- prob[index_output]
        }
      }
      dloss <- (as.numeric(BF) %*% as.numeric(gij*BF)) %>% as.numeric()
      return(dloss)
    }
    Loss <- 0
    p <- probs[1,]
    Pred_P <- probs
    for(i in 1:(length(Timepoints)-1)){
      MM <- Integral_step*(Timepoints[i+1] - Timepoints[i])
      delta_t <- 1/Integral_step
      for (j in 1:MM) {
        count_time <- 0
        while((count_time)< delta_t){
          rest_time <- delta_t - count_time
          dv <- derive_dv(p,Phi,W,beta)
          new_dv <- (rest_time)*dv
          if(length(which(new_dv > critical))==0){
            count_time <- count_time + rest_time
            dloss <- derive_dloss(p)
            Loss <- Loss + rest_time*dloss
            p <- p + (rest_time*derive_dp(p,Phi,W,beta))
          }else{
            new_step <- (rest_time/(max(dv)))*critical
            count_time <- count_time + new_step
            dloss <- derive_dloss(p)
            Loss <- Loss + new_step*dloss
            p <- p + (new_step*derive_dp(p,Phi,W,beta))
          }
        }
      }
      actual_p2 <- probs[i+1,] %>% as.numeric()
      Loss <- Loss + lambda[i]*KL_div(p,actual_p2)
      Pred_P[i+1,] <- p
    }
    LLoss <- list(Loss=Loss,Pred_P=Pred_P)
    return(LLoss)
  }
  
  #### initial pred_P & U1 & A_dKL_dp2 & Loss
  LLoss1 <- derive_loss(Phi,W,Timepoints)
  Loss1 <- LLoss1$Loss
  pred_P <- LLoss1$Pred_P
  U1 <- derive_u(Phi,W,pred_P,beta)
  A_dKL_dp2 <- derive_a_dKL_dp(Phi,W,pred_P,probs)
  KLs <- KL_divs(pred_P,probs)
  Loss_int <- Loss1-((lambda %*% KLs[2:f])/(f-1))
  
  Loss_KLs <- list()
  Loss <- c(Loss1)
  Loss_KLs[[1]] <- KLs
  Loss_intergral <- c(Loss_int)
  
  #### gradient descent to estimate Phi & W
  delta_loss <- 100
  num_iter <- 1
  while(delta_loss > tol & num_iter < Max_iter){
    num_iter <- num_iter+1
    delta_Phi <- rep(0,N)
    delta_W <- matrix(0,ncol=N,nrow=N)
    pars_step <- Par_step
    p <- probs[1,]
    #### compute delta_Phi & delta_W
    for(t in 1:(f-1)){
      M <- (Timepoints[t+1] - Timepoints[t]) %/% delta_ti
      last_ti <- (Timepoints[t+1] - Timepoints[t]) %% delta_ti
      for (j in 1:M) {
        ccount_time <- 0
        while((ccount_time)<delta_ti){
          rrest_time <- delta_ti - ccount_time
          ddv <- derive_dv(p,Phi,W,beta)
          new_ddv <- ((rrest_time)*ddv)
          if(length(which(new_ddv > critical))==0){
            ccount_time <- ccount_time + rrest_time
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[t],beta)
            delta_Phi <- delta_Phi + (rrest_time*(ddd$dPhi))
            delta_W <- delta_W +(rrest_time*(ddd$dW))
            p <- p + (rrest_time*derive_dp(p,Phi,W,beta))
          }else{
            new_step <- (rrest_time/(max(ddv)))*critical
            ccount_time <- ccount_time + new_step
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[t],beta)
            delta_Phi <- delta_Phi + (new_step*(ddd$dPhi))
            delta_W <- delta_W +(new_step*(ddd$dW))
            p <- p + (new_step*derive_dp(p,Phi,W,beta))
          }
        }
      }
      if(last_ti >0){
        ccount_time <- 0
        while((ccount_time)<last_ti){
          rrest_time <- last_ti - ccount_time
          ddv <- derive_dv(p,Phi,W,beta)
          new_ddv <- ((rrest_time)*ddv)
          if(length(which(new_ddv > critical))==0){
            ccount_time <- ccount_time + rrest_time
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[f-1],beta)
            delta_Phi <- delta_Phi + (rrest_time*(ddd$dPhi))
            delta_W <- delta_W +(rrest_time*(ddd$dW))
            p <- p + (rrest_time*derive_dp(p,Phi,W,beta))
          }else{
            new_step <- (rrest_time/(max(ddv)))*critical
            ccount_time <- ccount_time + new_step
            ddd <- derive_dPhi_dW(p,Phi,W,U1[[t]],A_dKL_dp2[[t+1]],lambda[f-1],beta)
            delta_Phi <- delta_Phi + (new_step*(ddd$dPhi))
            delta_W <- delta_W +(new_step*(ddd$dW))
            p <- p + (new_step*derive_dp(p,Phi,W,beta))
          }
        }
      }
    }
    delta_Phi <- delta_Phi+alpha*sgn_Phi(Phi)
    delta_W <- delta_W+alpha*sgn_W(W)
    new_Phi <- Phi - delta_Phi*pars_step
    new_W <- W - delta_W*pars_step
    new_LLoss <- derive_loss(new_Phi,new_W,Timepoints)
    new_Loss <- new_LLoss$Loss
    if(new_Loss >= Loss1){
      repeat{
        pars_step <- pars_step/2
        print(pars_step)
        new_Phi <- Phi - delta_Phi*pars_step
        new_W <- W - delta_W*pars_step
        new_LLoss <- derive_loss(new_Phi,new_W,Timepoints)
        new_Loss <- new_LLoss$Loss
        if(new_Loss < Loss1 | (pars_step*max(abs(delta_Phi),abs(delta_W)))<1e-04){
          if(new_Loss < Loss1){
            Phi <- new_Phi
            W <- new_W
            #### new pred_P & U1 & A_dKL_dp2 & Loss1
            pred_P <- new_LLoss$Pred_P
            U1 <- derive_u(Phi,W,pred_P,beta)
            A_dKL_dp2 <- derive_a_dKL_dp(Phi,W,pred_P,probs)
            KLs <- KL_divs(pred_P,probs)
            Loss_int <- Loss1-((lambda %*% KLs[2:f])/(f-1))
            delta_loss <- Loss1 - new_Loss
            Loss1 <- new_Loss
            
            Loss <- c(Loss,new_Loss)
            Loss_KLs[[num_iter]] <- KLs
            Loss_intergral <- c(Loss_intergral,Loss_int)
          }else{
            delta_loss <- 0
          }
          break
        }
      }
    }else{
      Phi <- new_Phi
      W <- new_W
      #### new pred_P & U1 & A_dKL_dp2 & Loss
      pred_P <- new_LLoss$Pred_P
      U1 <- derive_u(Phi,W,pred_P,beta)
      A_dKL_dp2 <- derive_a_dKL_dp(Phi,W,pred_P,probs)
      KLs <- KL_divs(pred_P,probs)
      Loss_int <- Loss1-((lambda %*% KLs[2:f])/(f-1))
      delta_loss <- Loss1 - new_Loss
      Loss1 <- new_Loss
      
      Loss <- c(Loss,new_Loss)
      Loss_KLs[[num_iter]] <- KLs
      Loss_intergral <- c(Loss_intergral,Loss_int)
    }
  }
  parameter <- list(Phi=Phi,W=W)
  Intermediate_Process <- list(Loss=Loss,Loss_intergral=Loss_intergral,Loss_KLs=Loss_KLs)
  T2=Sys.time()
  list_result_ALL <- list(parameter=parameter,Intermediate_Process=Intermediate_Process,Pred_probs=pred_P,Actual_probs=probs,Time=(T2-T1))
  return(list_result_ALL)
}
