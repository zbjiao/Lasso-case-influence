library(lars)
library(tidyverse)
library(MASS)
library(latex2exp)
library(cowplot)
library(foreach)
library(doParallel)
library(tmvtnorm)

# case influence all samples, all fractions
CD_all_concise <- function(X,y,finesse=0.02,plot = T,threshold = F){
  centralize <- function(x, r = F){
    # centralize x with mean 0 and sd 1/sqrt(n-1)
    n = dim(x)[1]
    meanx <- drop(rep(1,n) %*% x)/n
    x = scale(x, meanx, FALSE)
    normx = sqrt(drop(rep(1,n) %*% (x^2)))
    if (r){
      return(list(m = meanx, d = scale(x,FALSE,normx), v = normx)) 
    }
    scale(x, FALSE, normx) 
  }
  get_xi <- function(w, h, n){
    # if (h>=(1-1/n-0.0000001)&(w==0)){
    #   10
    # }
    # else{
      # given omega, n and hkk, return xi
    ans = n*(1-w)/(n-1+w-n*(1-w)*h)
    if (ans<0){
      return(-ans)
    }
    else{
      return(ans)
    }
  }
  get_w <- function(xi, h, n){
    # reverse function of get_xi, 
    # given xi, hkk and n, return omega
    (n-xi*(n-1-n*h))/(n+(1+n*h)*xi)}
  
  #---------initialization and find lambda-fraction correspondence---------#
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]
  if (n<=p){
    denom = 1
    cd0 = rep(NA,n)
  }
  else{
    denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
    cd0 = cooks.distance(lm(y~X))
  }
  
  step = log(1+max(abs(t(X)%*%y))*1.001)/10
  # lambda values we take to draw the graph
  l_list = exp(0:10*step) - 1
  # record the Cook's distance of all data points given lambda
  
  l1norm = c()
  # record all gram inverses we come across during the calculation 
  XX_recorder = list()
  XX_recorder[[paste(rep(0,p),collapse = '')]] = matrix(0,nrow=p,ncol=p)
  if (n<p & p>500){
    obj = lars(X,y,type='lasso',use.Gram = F)
  }
  else{
    obj = lars(X,y,type='lasso')
  }
  ybar = mean(y)
  
  l1norm_list = c()
  beta_hat_table = matrix(nrow = p, ncol = 11)
  p_list = rep(0,11)
  for (i in 1:11){
    lambda = l_list[i]
    beta_hat = as.vector(predict.lars(obj,mode='lambda', s=lambda, type = 'coefficients')$coefficients)
    beta_hat_table[,i] = beta_hat
    p_list[i] = sum(beta_hat!=0)
    l1norm = sum(abs(beta_hat))
    l1norm_list = c(l1norm_list, l1norm)
  }
  fraction = as.numeric(l1norm_list)/l1norm_list[1]
  fraction_dif = fraction[1:10] - fraction[2:11]
  while(max(fraction_dif)>finesse){
    ind = which(fraction_dif == max(fraction_dif))[1]
    new_l = (l_list[ind] + l_list[ind+1])/2
    l_list = c(l_list[1:ind], new_l, l_list[(ind+1):length(l_list)])
    beta_hat = as.vector(predict.lars(obj,mode='lambda', s=new_l, type = 'coefficients')$coefficients)
    beta_hat_table = cbind(beta_hat_table[,1:ind],beta_hat,beta_hat_table[,(ind+1):(length(l_list)-1)])
    p_list = c(p_list[1:ind], sum(beta_hat!=0), p_list[(ind+1):length(p_list)])
    l = length(l_list)
    fraction = c(fraction[1:ind], sum(abs(beta_hat))/l1norm_list[1], fraction[(ind+1):length(fraction)])
    fraction_dif = fraction[1:(l-1)] - fraction[2:l]
  }
  # CD_Mat = matrix(nrow = n, ncol = length(fraction))
  
  # threshold_list = c()
  CD_Mat = c()
  # return(list('f' = fraction,'l'=l_list))
  for (i in 2:length(fraction)){
    # print(c(i,length(fraction)))
    lambda = l_list[i]
    
    beta_hat_backup = beta_hat_table[,i]
    s_backup = sign(beta_hat_backup)
    A_backup = s_backup!=0
    # derivative of f in terms of covariates in nonactive set
    d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
    
    A_id = paste(A_backup*1,collapse = '')
    if (A_id %in% names(XX_recorder)){
      XX_backup = XX_recorder[[A_id]]
    }
    else{
      XX_backup = solve(t(X[,A_backup])%*%X[,A_backup])
      XX_recorder[[A_id]] = XX_backup
    }
    
    if (sum(A_backup)==0){
      h_backup = matrix(0,ncol=n,nrow=n)
    }
    else{
      h_backup = X[,A_backup]%*%XX_backup%*%t(X[,A_backup])
    }
    
    # threshold_list = c(threshold_list, tmvnorm_to_f_mc(beta_hat_backup, s_backup, 0.05))
    CD_list = c()
    
    for (k in 1:n){
      beta_hat = beta_hat_backup
      s = s_backup
      A = A_backup
      
      d_hat = d_hat_backup
      XX = XX_backup
      hk = h_backup[,k]
      
      # current predictor of yk 
      yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
      
      # beta path records beta hat's value at each breakpoint
      beta_path = c(beta_hat)
      # and sign change
      s_path = c(s)
      # and omega value at each breakpoint
      hkk_path = c()
      w = 1
      while (T){
        hkk_path = c(hkk_path, hk[k])
        xi = get_xi(w, hk[k], n)
        bias = yk_hat - y[k]
        if (sum(A) == 0){
          xi_cand1 = c()
        }
        else{
          slope_beta = XX%*%X[k,A]*bias
          xi_cand1 = -beta_hat[A]/slope_beta
        }
        slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
        # xi candidates
        # xi_cand1 = -beta_hat[A]/slope_beta
        xi_cand2p = (-lambda-d_hat)/slope_d
        xi_cand2m = (lambda-d_hat)/slope_d
        xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
        xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
        # print(c(xi_cand0, xi))
        ind = which(xi_cand == xi_cand0)
        
        # update beta
        if (sum(A) > 0){
          beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
        }
        beta_path = rbind(beta_path, beta_hat)
        beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
        
        # if the xi is off the bound, stop the algorithm
        # print(get_xi(0,hk[k],n))
        # print(beta_hat[A])
        if (xi_cand0 > get_xi(0,hk[k],n)){
          hkk_path = c(hkk_path, hk[k])
          s_path = rbind(s_path, s)
          break
        }
        # if not, locate the covariate (go to or leave the active set)
        else if(ind<=sum(A)){
          coef = which(cumsum(A)==ind)[1]
          # print(c('-',coef))
          A[coef] = F
          s[coef] = 0
          # check if active set is empty 
          if (sum(A) == 0){
            hkk_path = c(hkk_path, 0, 0)
            s_path = rbind(s_path, s, s)
            beta_path = rbind(beta_path, beta_hat)
            break
          }
        }
        else if(ind>p){
          coef = which(cumsum(1-A)==(ind-p))[1]
          # print(c('+',coef))
          A[coef] = T
          s[coef] = 1
        }
        else{
          coef = which(cumsum(1-A)==(ind-sum(A)))[1]
          # print(c('+',coef))
          A[coef] = T
          s[coef] = -1
        }
        
        s_path = rbind(s_path, s)
        # update omega, XX, beta_hat, d_hat
        w = get_w(xi_cand0, hk[k], n)
        A_id = paste(A*1,collapse = '')
        if (A_id %in% names(XX_recorder)){
          XX = XX_recorder[[A_id]]
        }
        else{
          XX = solve(crossprod(X[,A]))
          # XX = solve(t(X[,A])%*%X[,A])
          XX_recorder[[A_id]] = XX
        }
        beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
        beta_hat[!A] = 0
        
        # d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
        d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X[,A]%*%beta_hat[A]
        hk = X[,A]%*%XX%*%X[k,A]
        yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
        # print(c(yk_hat,w))
      }
      
      #--------------------
      
      K = length(hkk_path)
      # print(c(k,n,K))
      A = s_path[1,]!=0
      if (sum(A)==1){
        y_tilde = X[,A]*beta_path[1,A] + ybar
      } 
      else{
        y_tilde = X[,A]%*%beta_path[1,A] + ybar
      }
      
      A = s_path[K,]!=0
      A_id = paste(A*1,collapse = '')
      if (sum(A) == 0){
        y_last = ybar + get_xi(0, hkk_path[K], n)*1/n*(ybar - y[k])
      }
      else{
        xxx =  X[,A]%*%XX_recorder[[A_id]]
        y_hat =  xxx%*%(t(X[,A])%*%y-lambda*s_path[K,A]) + ybar
        y_last = y_hat + get_xi(0, hkk_path[K], n)*(xxx%*%X[k,A]+1/n)*(y_hat[k] - y[k])
      }
      
      CD_list = c(CD_list, sum((y_last - y_tilde)**2))
    }
    # CD_Mat[,i] = CD_list/denom
    CD_Mat = cbind(CD_Mat, CD_list/denom)
  }
  CD_Mat = cbind(cd0,CD_Mat)
  # threshold_list = c(qf(0.1,p+1,n-p-1)*(p+1)/p,threshold_list)
  # threshold_list = c(tmvnorm_to_f_mc())
  threshold_table = cbind(fraction,sqrt(apply(CD_Mat,2,var)/2)*qchisq(0.95,1))
  
  # threshold_table = cbind(fraction, apply(CD_Mat,2,sd)*3)
  # threshold_table = cbind(fraction, qf(1-0.05**(1/n),p_list+1,(n-p-1))*(p_list+1)/p)
  # threshold_table = cbind(fraction, threshold_list)
  #--------------- plot case influence graph 0.42s----------#
  if (plot){
    plot_table = t(rbind(fraction, CD_Mat))
    # mean_table = data.frame('fraction' = fraction, 'meancd' = apply(CD_Mat,2,mean)*10)
    # numfrac = dim(mean_table)[1]
    # mean_table[numfrac,2] = mean_table[numfrac,2]/2
    # print(mean_table)
    colnames(plot_table) = c('fraction',1:n)
    df = as_tibble(plot_table)%>%gather('obs','val',-fraction)
    # cutting_fractions = c()
    # L = length(p_list)
    # for (index in 2:(L-1)){
    #   if (p_list[index]!=p_list[index-1]){
    #     cutting_fractions = c(cutting_fractions, (fraction[index]+fraction[index-1])/2)
    #   }
    # }
    # cutting_fractions = c(cutting_fractions,0)
    # cutting_fractions = data.frame(c = cutting_fractions)
    p = ggplot(data=df,aes(x=fraction,y=val,group=obs,color=obs)) +
      geom_line() + xlim(c(0,NA))+ labs(x = "|coef|/max|coef|")+
      theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +labs(y="Cook's distance under Lasso")+
      # geom_line(data=mean_table, aes(x=fraction, y=meancd),color='blue',linetype='dashed',alpha=1,inherit.aes = FALSE)+
      theme(panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
            panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
    # +
    #   geom_vline(data=cutting_fractions, aes(xintercept = c),color='blue',linetype='dashed',alpha=0.5)
    if (threshold){
      colnames(threshold_table)=c('fraction','val')
      df2 = as_tibble(threshold_table)%>%mutate(obs='threshold')
      p = p + geom_line(data=df2,aes(x=fraction, y=val),linetype = "dotted", color = 'red1')
    }
    return(list(CD_Mat = CD_Mat, Lambda_list = l_list, fraction = fraction, 
                threshold = threshold_table, beta_table = beta_hat_table, p=p))
  }
  else{
    return(list(CD_Mat = CD_Mat, Lambda_list = l_list, fraction = fraction, 
                threshold = threshold_table, beta_table = beta_hat_table))
  }
  
}
# case influence all samples, one fraction
CD_one_fraction <- function(X,y,lambda=1){
  centralize <- function(x, r = F){
    # centralize x with mean 0 and sd 1/sqrt(n-1)
    n = dim(x)[1]
    meanx <- drop(rep(1,n) %*% x)/n
    x = scale(x, meanx, FALSE)
    normx = sqrt(drop(rep(1,n) %*% (x^2)))
    if (r){
      return(list(m = meanx, d = scale(x,FALSE,normx), v = normx)) 
    }
    scale(x, FALSE, normx) 
  }
  get_xi <- function(w, h, n){
    # if (h>=(1-1/n-0.0000001)&(w==0)){
    #   10
    # }
    # else{
    # given omega, n and hkk, return xi
    ans = n*(1-w)/(n-1+w-n*(1-w)*h)
    if (ans<0){
      return(-ans)
    }
    else{
      return(ans)
    }
  }
  get_w <- function(xi, h, n){
    # reverse function of get_xi, 
    # given xi, hkk and n, return omega
    (n-xi*(n-1-n*h))/(n+(1+n*h)*xi)}
  
  #---------initialization and find lambda-fraction correspondence---------#
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]
  if (n<=p){
    denom = 1
  }
  else{
    denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
  }
  
  l1norm = c()
  # record all gram inverses we come across during the calculation 
  XX_recorder = list()
  XX_recorder[[paste(rep(0,p),collapse = '')]] = matrix(0,nrow=p,ncol=p)
  if (n<p & p>500){
    obj = lars(X,y,type='lasso',use.Gram = F)
  }
  else{
    obj = lars(X,y,type='lasso')
  }
  ybar = mean(y)
  # obj=lars(x,y,type='lasso',normalize = F)
  # coef(obj,s=1,mode='lambda')
  fit = glmnet(X,y,family="gaussian",lambda=lambda/n,standardize=FALSE,thresh=1e-16,intercept=T)
  beta_hat = as.vector(fit$beta)
  
  beta_hat_backup = beta_hat
  s_backup = sign(beta_hat_backup)
  A_backup = s_backup!=0
  # derivative of f in terms of covariates in nonactive set
  if(sum(!A_backup)==0){}
  d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
  
  A_id = paste(A_backup*1,collapse = '')
  if (A_id %in% names(XX_recorder)){
    XX_backup = XX_recorder[[A_id]]
  }
  else{
    XX_backup = solve(t(X[,A_backup])%*%X[,A_backup])
    XX_recorder[[A_id]] = XX_backup
  }
  
  if (sum(A_backup)==0){
    h_backup = matrix(0,ncol=n,nrow=n)
  }
  else{
    h_backup = X[,A_backup]%*%XX_backup%*%t(X[,A_backup])
  }
  
  # threshold_list = c(threshold_list, tmvnorm_to_f_mc(beta_hat_backup, s_backup, 0.05))
  CD_list = c()
  
  for (k in 1:n){
    beta_hat = beta_hat_backup
    s = s_backup
    A = A_backup
    
    d_hat = d_hat_backup
    XX = XX_backup
    hk = h_backup[,k]
    
    # current predictor of yk 
    yk_hat = drop(ybar + X[k,] %*% beta_hat)
    
    # beta path records beta hat's value at each breakpoint
    beta_path = c(beta_hat)
    # and sign change
    s_path = c(s)
    # and omega value at each breakpoint
    hkk_path = c()
    w = 1
    while (T){
      hkk_path = c(hkk_path, hk[k])
      xi = get_xi(w, hk[k], n)
      bias = yk_hat - y[k]
      if (sum(A) == 0){
        xi_cand1 = c()
      }
      else{
        slope_beta = XX%*%X[k,A]*bias
        xi_cand1 = -beta_hat[A]/slope_beta
      }
      slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
      # xi candidates
      # xi_cand1 = -beta_hat[A]/slope_beta
      xi_cand2p = (-lambda-d_hat)/slope_d
      xi_cand2m = (lambda-d_hat)/slope_d
      xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
      xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
      # print(c(xi_cand0, xi))
      ind = which(xi_cand == xi_cand0)
      
      # update beta
      if (sum(A) > 0){
        beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
      }
      beta_path = rbind(beta_path, beta_hat)
      beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
      
      # if the xi is off the bound, stop the algorithm
      # print(get_xi(0,hk[k],n))
      # print(beta_hat[A])
      if (xi_cand0 > get_xi(0,hk[k],n)){
        hkk_path = c(hkk_path, hk[k])
        s_path = rbind(s_path, s)
        break
      }
      # if not, locate the covariate (go to or leave the active set)
      else if(ind<=sum(A)){
        coef = which(cumsum(A)==ind)[1]
        # print(c('-',coef))
        A[coef] = F
        s[coef] = 0
        # check if active set is empty 
        if (sum(A) == 0){
          hkk_path = c(hkk_path, 0, 0)
          s_path = rbind(s_path, s, s)
          beta_path = rbind(beta_path, beta_hat)
          break
        }
      }
      else if(ind>p){
        coef = which(cumsum(1-A)==(ind-p))[1]
        # print(c('+',coef))
        A[coef] = T
        s[coef] = 1
      }
      else{
        coef = which(cumsum(1-A)==(ind-sum(A)))[1]
        # print(c('+',coef))
        A[coef] = T
        s[coef] = -1
      }
      
      s_path = rbind(s_path, s)
      # update omega, XX, beta_hat, d_hat
      w = get_w(xi_cand0, hk[k], n)
      A_id = paste(A*1,collapse = '')
      if (A_id %in% names(XX_recorder)){
        XX = XX_recorder[[A_id]]
      }
      else{
        XX = solve(t(X[,A])%*%X[,A])
        XX_recorder[[A_id]] = XX
      }
      beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
      beta_hat[!A] = 0
      
      # d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
      d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X[,A]%*%beta_hat[A]
      hk = X[,A]%*%XX%*%X[k,A]
      yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
      # print(c(yk_hat,w))
    }
    
    #--------------------
    
    K = length(hkk_path)
    # print(c(k,n,K))
    A = s_path[1,]!=0
    if (sum(A)==1){
      y_tilde = X[,A]*beta_path[1,A] + ybar
    } 
    else{
      y_tilde = X[,A]%*%beta_path[1,A] + ybar
    }
    
    A = s_path[K,]!=0
    A_id = paste(A*1,collapse = '')
    if (sum(A) == 0){
      y_last = ybar + get_xi(0, hkk_path[K], n)*1/n*(ybar - y[k])
    }
    else{
      xxx =  X[,A]%*%XX_recorder[[A_id]]
      y_hat =  xxx%*%(t(X[,A])%*%y-lambda*s_path[K,A]) + ybar
      y_last = y_hat + get_xi(0, hkk_path[K], n)*(xxx%*%X[k,A]+1/n)*(y_hat[k] - y[k])
    }
    
    CD_list = c(CD_list, sum((y_last - y_tilde)**2))
  }
  # CD_Mat[,i] = CD_list/denom
  CD_vec = CD_list/denom
  # threshold_list = c(qf(0.1,p+1,n-p-1)*(p+1)/p,threshold_list)
  # threshold_list = c(tmvnorm_to_f_mc())
  threshold = sqrt(var(CD_vec)/2)*qchisq(0.95,1)

  return(list(CD_vec = CD_vec, threshold = threshold, beta = beta_hat))
}
# updated version faster than the previous one (this one assume X is already normalized)
CD_one_fraction <- function(X,y,lambda=1, alpha=0.05, external = F){
  get_xi <- function(w, h, n){
    # if (h>=(1-1/n-0.0000001)&(w==0)){
    #   10
    # }
    # else{
    # given omega, n and hkk, return xi
    ans = n*(1-w)/(n-1+w-n*(1-w)*h)
    if (ans<0){
      return(-ans)
    }
    else{
      return(ans)
    }
  }
  get_w <- function(xi, h, n){
    # reverse function of get_xi, 
    # given xi, hkk and n, return omega
    (n-xi*(n-1-n*h))/(n+(1+n*h)*xi)}
  
  #---------initialization and find lambda-fraction correspondence---------#
  n = dim(X)[1]
  p = dim(X)[2]
  if (n<=p){
    denom = 1
  }
  else{
    denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
  }
  
  l1norm = c()
  # record all gram inverses we come across during the calculation 
  XX_recorder = list()
  XX_recorder[[paste(rep(0,p),collapse = '')]] = matrix(0,nrow=p,ncol=p)
  if (n<p & p>500){
    obj = lars(X,y,type='lasso',use.Gram = F)
  }
  else{
    obj = lars(X,y,type='lasso')
  }
  ybar = mean(y)
  # obj=lars(x,y,type='lasso',normalize = F)
  # coef(obj,s=1,mode='lambda')
  fit = glmnet(X,y,family="gaussian",lambda=lambda/n,standardize=FALSE,thresh=1e-16,intercept=T)
  beta_hat = as.vector(fit$beta)
  
  beta_hat_backup = beta_hat
  s_backup = sign(beta_hat_backup)
  A_backup = s_backup!=0
  # derivative of f in terms of covariates in nonactive set
  if(sum(!A_backup)==0){}
  d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
  
  A_id = paste(A_backup*1,collapse = '')
  if (A_id %in% names(XX_recorder)){
    XX_backup = XX_recorder[[A_id]]
  }
  else{
    XX_backup = solve(crossprod(X[,A_backup]))
    # XX_backup = solve(t(X[,A_backup])%*%X[,A_backup])
    XX_recorder[[A_id]] = XX_backup
  }
  
  if (sum(A_backup)==0){
    h_backup = matrix(0,ncol=n,nrow=n)
  }
  else{
    h_backup = X[,A_backup]%*%XX_backup%*%t(X[,A_backup])
  }
  yhat = ybar + X[,A_backup] %*% beta_hat[A_backup]
  # threshold_list = c(threshold_list, tmvnorm_to_f_mc(beta_hat_backup, s_backup, 0.05))
  CD_list = c()
  
  for (k in 1:n){
    beta_hat = beta_hat_backup
    s = s_backup
    A = A_backup
    
    d_hat = d_hat_backup
    XX = XX_backup
    hk = h_backup[,k]
    
    # current predictor of yk 
    yk_hat = drop(yhat[k])
    # yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
    # yk_hat = drop(ybar + X[k,] %*% beta_hat)
    
    # beta path records beta hat's value at each breakpoint
    beta_path = c(beta_hat)
    # and sign change
    s_path = c(s)
    # and omega value at each breakpoint
    hkk_path = c()
    w = 1
    while (T){
      hkk_path = c(hkk_path, hk[k])
      xi = get_xi(w, hk[k], n)
      bias = yk_hat - y[k]
      if (sum(A) == 0){
        xi_cand1 = c()
      }
      else{
        slope_beta = XX%*%X[k,A]*bias
        xi_cand1 = -beta_hat[A]/slope_beta
      }
      slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
      # xi candidates
      # xi_cand1 = -beta_hat[A]/slope_beta
      xi_cand2p = (-lambda-d_hat)/slope_d
      xi_cand2m = (lambda-d_hat)/slope_d
      xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
      xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
      # print(c(xi_cand0, xi))
      ind = which(xi_cand == xi_cand0)
      
      # update beta
      if (sum(A) > 0){
        beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
      }
      beta_path = rbind(beta_path, beta_hat)
      beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
      
      # if the xi is off the bound, stop the algorithm
      if (xi_cand0 > get_xi(0,hk[k],n)){
        hkk_path = c(hkk_path, hk[k])
        s_path = rbind(s_path, s)
        break
      }
      # if not, locate the covariate (go to or leave the active set)
      else if(ind<=sum(A)){
        coef = which(cumsum(A)==ind)[1]
        # print(c('-',coef))
        A[coef] = F
        s[coef] = 0
        # check if active set is empty 
        if (sum(A) == 0){
          hkk_path = c(hkk_path, 0, 0)
          s_path = rbind(s_path, s, s)
          beta_path = rbind(beta_path, beta_hat)
          break
        }
      }
      else if(ind>p){
        coef = which(cumsum(1-A)==(ind-p))[1]
        # print(c('+',coef))
        A[coef] = T
        s[coef] = 1
      }
      else{
        coef = which(cumsum(1-A)==(ind-sum(A)))[1]
        # print(c('+',coef))
        A[coef] = T
        s[coef] = -1
      }
      
      s_path = rbind(s_path, s)
      # update omega, XX, beta_hat, d_hat
      w = get_w(xi_cand0, hk[k], n)
      A_id = paste(A*1,collapse = '')
      if (A_id %in% names(XX_recorder)){
        XX = XX_recorder[[A_id]]
      }
      else{
        XX = solve(crossprod(X[,A]))
        # XX = solve(t(X[,A])%*%X[,A])
        XX_recorder[[A_id]] = XX
      }
      beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
      beta_hat[!A] = 0
      # d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X[,A_backup]%*%beta_hat_backup[A_backup]
      d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X[,A]%*%beta_hat[A]
      hk = X[,A]%*%XX%*%X[k,A]
      yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
      # print(c(yk_hat,w))
    }
    
    #--------------------
    
    K = length(hkk_path)
    # print(c(k,n,K))
    A = s_path[1,]!=0
    if (sum(A)==1){
      y_tilde = X[,A]*beta_path[1,A] + ybar
    } 
    else{
      y_tilde = X[,A]%*%beta_path[1,A] + ybar
    }
    
    A = s_path[K,]!=0
    A_id = paste(A*1,collapse = '')
    if (sum(A) == 0){
      y_last = ybar + get_xi(0, hkk_path[K], n)*1/n*(ybar - y[k])
    }
    else{
      xxx =  X[,A]%*%XX_recorder[[A_id]]
      y_hat =  xxx%*%(t(X[,A])%*%y-lambda*s_path[K,A]) + ybar
      y_last = y_hat + get_xi(0, hkk_path[K], n)*(xxx%*%X[k,A]+1/n)*(y_hat[k] - y[k])
    }
    
    CD_list = c(CD_list, sum((y_last - y_tilde)**2))
  }
  # CD_Mat[,i] = CD_list/denom
  CD_vec = CD_list/denom
  # threshold_list = c(qf(0.1,p+1,n-p-1)*(p+1)/p,threshold_list)
  # threshold_list = c(tmvnorm_to_f_mc())
  if (external){
    threshold = c()
    for (i in 1:n){
      threshold = c(threshold, sqrt(var(CD_vec[-i])/2)*qchisq(1-alpha,1))
    }
  }
  else{
    threshold = sqrt(var(CD_vec)/2)*qchisq(1-alpha,1)
  }
  
  
  return(list(CD_vec = CD_vec, threshold = threshold, beta = beta_hat))
}
# case influence one sample, all fraction
CD_one_concise <- function(X,y,k,finesse=0.02){
  centralize <- function(x, r = F){
    # centralize x with mean 0 and sd 1/sqrt(n-1)
    n = dim(x)[1]
    meanx <- drop(rep(1,n) %*% x)/n
    x = scale(x, meanx, FALSE)
    normx = sqrt(drop(rep(1,n) %*% (x^2)))
    if (r){
      return(list(m = meanx, d = scale(x,FALSE,normx), v = normx)) 
    }
    scale(x, FALSE, normx) 
  }
  get_xi <- function(w, h, n){
    # given omega, n and hkk, return xi
    n*(1-w)/(n-1+w-n*(1-w)*h)
  }
  get_w <- function(xi, h, n){
    # reverse function of get_xi, 
    # given xi, hkk and n, return omega
    (n-xi*(n-1-n*h))/(n+(1+n*h)*xi)}
  
  #---------initialization and find lambda-fraction correspondance 0.1s---------#
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]
  denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
  # denom = 1
  cd0 = cooks.distance(lm(y~X))
  # cd0 = rep(0,n)
  step = log(1+max(abs(t(X)%*%y))*1.001)/10
  # lambda values we take to draw the graph
  l_list = exp(0:10*step) - 1
  # record the Cook's distance of all data points given lambda
  
  l1norm = c()
  # record all gram inverses we come across during the calculation 
  XX_recorder = list()
  XX_recorder[[paste(rep(0,p),collapse = '')]] = matrix(0,nrow=p,ncol=p)
  obj = lars(X,y,type='lasso')
  ybar = mean(y)
  
  l1norm_list = c()
  beta_hat_table = matrix(nrow = p, ncol = 11)
  p_list = rep(0,11)
  for (i in 1:11){
    lambda = l_list[i]
    beta_hat = as.vector(predict.lars(obj,mode='lambda', s=lambda, type = 'coefficients')$coefficients)
    beta_hat_table[,i] = beta_hat
    p_list[i] = sum(beta_hat!=0)
    l1norm = sum(abs(beta_hat))
    l1norm_list = c(l1norm_list, l1norm)
  }
  fraction = as.numeric(l1norm_list)/l1norm_list[1]
  fraction_dif = fraction[1:10] - fraction[2:11]
  while(max(fraction_dif)>finesse){
    ind = which(fraction_dif == max(fraction_dif))[1]
    new_l = (l_list[ind] + l_list[ind+1])/2
    l_list = c(l_list[1:ind], new_l, l_list[(ind+1):length(l_list)])
    beta_hat = as.vector(predict.lars(obj,mode='lambda', s=new_l, type = 'coefficients')$coefficients)
    beta_hat_table = cbind(beta_hat_table[,1:ind],beta_hat,beta_hat_table[,(ind+1):(length(l_list)-1)])
    p_list = c(p_list[1:ind], sum(beta_hat!=0), p_list[(ind+1):length(p_list)])
    l = length(l_list)
    fraction = c(fraction[1:ind], sum(abs(beta_hat))/l1norm_list[1], fraction[(ind+1):length(fraction)])
    fraction_dif = fraction[1:(l-1)] - fraction[2:l]
  }
  # CD_Mat = matrix(nrow = n, ncol = length(fraction))
  
  # threshold_list = c()
  CD_Mat = c()
  
  for (i in 2:length(fraction)){
    lambda = l_list[i]
    
    beta_hat_backup = beta_hat_table[,i]
    s_backup = sign(beta_hat_backup)
    A_backup = s_backup!=0
    # derivative of f in terms of covariates in nonactive set
    d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X%*%beta_hat_backup
    
    A_id = paste(A_backup*1,collapse = '')
    if (A_id %in% names(XX_recorder)){
      XX_backup = XX_recorder[[A_id]]
    }
    else{
      XX_backup = solve(t(X[,A_backup])%*%X[,A_backup])
      XX_recorder[[A_id]] = XX_backup
    }
    
    if (sum(A_backup)==0){
      h_backup = matrix(0,ncol=n,nrow=n)
    }
    else{
      h_backup = X[,A_backup]%*%XX_backup%*%t(X[,A_backup])
    }
    
    # threshold_list = c(threshold_list, tmvnorm_to_f_mc(beta_hat_backup, s_backup, 0.05))
    CD_list = c()
    
  
    beta_hat = beta_hat_backup
    s = s_backup
    A = A_backup
    
    d_hat = d_hat_backup
    XX = XX_backup
    hk = h_backup[,k]
    
    # current predictor of yk 
    yk_hat = drop(ybar + X[k,] %*% beta_hat)
    
    # beta path records beta hat's value at each breakpoint
    beta_path = c(beta_hat)
    # and sign change
    s_path = c(s)
    # and omega value at each breakpoint
    hkk_path = c()
    w = 1
    while (T){
      hkk_path = c(hkk_path, hk[k])
      xi = get_xi(w, hk[k], n)
      bias = yk_hat - y[k]
      if (sum(A) == 0){
        xi_cand1 = c()
      }
      else{
        slope_beta = XX%*%X[k,A]*bias
        xi_cand1 = -beta_hat[A]/slope_beta
      }
      slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
      # xi candidates
      # xi_cand1 = -beta_hat[A]/slope_beta
      xi_cand2p = (-lambda-d_hat)/slope_d
      xi_cand2m = (lambda-d_hat)/slope_d
      xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
      xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
      ind = which(xi_cand == xi_cand0)
      
      # update beta
      if (sum(A) > 0){
        beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
      }
      beta_path = rbind(beta_path, beta_hat)
      beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
      
      # if the xi is off the bound, stop the algorithm
      if (xi_cand0 > get_xi(0,hk[k],n)){
        hkk_path = c(hkk_path, hk[k])
        s_path = rbind(s_path, s)
        break
      }
      # if not, locate the covariate (go to or leave the active set)
      else if(ind<=sum(A)){
        coef = which(cumsum(A)==ind)[1]
        A[coef] = F
        s[coef] = 0
        # check if active set is empty 
        if (sum(A) == 0){
          hkk_path = c(hkk_path, 0, 0)
          s_path = rbind(s_path, s, s)
          beta_path = rbind(beta_path, beta_hat)
          break
        }
      }
      else if(ind>p){
        coef = which(cumsum(1-A)==(ind-p))[1]
        A[coef] = T
        s[coef] = 1
      }
      else{
        coef = which(cumsum(1-A)==(ind-sum(A)))[1]
        A[coef] = T
        s[coef] = -1
      }
      
      s_path = rbind(s_path, s)
      # update omega, XX, beta_hat, d_hat
      w = get_w(xi_cand0, hk[k], n)
      A_id = paste(A*1,collapse = '')
      if (A_id %in% names(XX_recorder)){
        XX = XX_recorder[[A_id]]
      }
      else{
        XX = solve(t(X[,A])%*%X[,A])
        XX_recorder[[A_id]] = XX
      }
      beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
      beta_hat[!A] = 0
      d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X%*%beta_hat
      hk = X[,A]%*%XX%*%X[k,A]
      yk_hat = drop(ybar + X[k,] %*% beta_hat)
    }
    
    #--------------------
    
    K = length(hkk_path)
    A = s_path[1,]!=0
    if (sum(A)==1){
      y_tilde = X[,A]*beta_path[1,A] + ybar
    } 
    else{
      y_tilde = X[,A]%*%beta_path[1,A] + ybar
    }
    
    A = s_path[K,]!=0
    A_id = paste(A*1,collapse = '')
    if (sum(A) == 0){
      y_last = ybar + get_xi(0, hkk_path[K], n)*1/n*(ybar - y[k])
    }
    else{
      xxx =  X[,A]%*%XX_recorder[[A_id]]
      y_hat =  xxx%*%(t(X[,A])%*%y-lambda*s_path[K,A]) + ybar
      y_last = y_hat + get_xi(0, hkk_path[K], n)*(xxx%*%X[k,A]+1/n)*(y_hat[k] - y[k])
    }
    
    CD_list = c(CD_list, sum((y_last - y_tilde)**2))
    
    
    # CD_Mat[,i] = CD_list/denom
    CD_Mat = cbind(CD_Mat, CD_list/denom)
  }
  CD_Mat = cbind(cd0[k],CD_Mat)
  # threshold_list = c(qf(0.1,p+1,n-p-1)*(p+1)/p,threshold_list)
  # threshold_list = c(tmvnorm_to_f_mc())
  
  plot_table = t(rbind(fraction, CD_Mat))
  # threshold_table = cbind(fraction, apply(CD_Mat,2,mean)*3)
  # threshold_table = cbind(fraction, qf(1-0.05**(1/n),p_list+1,(n-p-1))*(p_list+1)/p)
  # threshold_table = cbind(fraction, threshold_list)
  colnames(plot_table) = c('fraction',k)
  # colnames(threshold_table)=c('fraction','val')
  df = as_tibble(plot_table)%>%gather('obs','val',-fraction)
  # df2 = as_tibble(threshold_table)%>%mutate(obs='threshold')
  p = ggplot(data=df,aes(x=fraction,y=val)) +
    geom_line() + xlim(c(0,NA))+ xlab(TeX("$|coef|/max|coef|$"))+
    # geom_line(data=df2,aes(x=fraction, y=val),linetype = "dotted", color = 'red1')+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +labs(y="Cook's distance")
  # print(p)
  # return(list(CD_Mat = CD_Mat, Lambda_list = l_list, fraction = fraction))
}
# case influence all samples, all fractions, parallel computing 
CD_all_concise_parallel <- function(X,y,finesse=0.02){
  #---------initialization and find lambda-fraction correspondence 0.1s---------#
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]
  denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
  # denom = 1
  cd0 = cooks.distance(lm(y~X))
  # cd0 = rep(0,n)
  step = log(1+max(abs(t(X)%*%y))*1.001)/10
  # lambda values we take to draw the graph
  l_list = exp(0:10*step) - 1
  # record the Cook's distance of all data points given lambda
  
  l1norm = c()
  # record all gram inverses we come across during the calculation 
  XX_recorder = list()
  XX_recorder[[paste(rep(0,p),collapse = '')]] = matrix(0,nrow=p,ncol=p)
  obj = lars(X,y,type='lasso')
  ybar = mean(y)
  
  l1norm_list = c()
  beta_hat_table = matrix(nrow = p, ncol = 11)
  for (i in 1:11){
    lambda = l_list[i]
    beta_hat = as.vector(predict.lars(obj,mode='lambda', s=lambda, type = 'coefficients')$coefficients)
    beta_hat_table[,i] = beta_hat
    l1norm = sum(abs(beta_hat))
    l1norm_list = c(l1norm_list, l1norm)
  }
  fraction = as.numeric(l1norm_list)/l1norm_list[1]
  fraction_dif = fraction[1:10] - fraction[2:11]
  while(max(fraction_dif)>finesse){
    ind = which(fraction_dif == max(fraction_dif))[1]
    new_l = (l_list[ind] + l_list[ind+1])/2
    l_list = c(l_list[1:ind], new_l, l_list[(ind+1):length(l_list)])
    beta_hat = as.vector(predict.lars(obj,mode='lambda', s=new_l, type = 'coefficients')$coefficients)
    beta_hat_table = cbind(beta_hat_table[,1:ind],beta_hat,beta_hat_table[,(ind+1):(length(l_list)-1)])
    l = length(l_list)
    fraction = c(fraction[1:ind], sum(abs(beta_hat))/l1norm_list[1], fraction[(ind+1):length(fraction)])
    fraction_dif = fraction[1:(l-1)] - fraction[2:l]
  }
  # CD_Mat = matrix(nrow = n, ncol = length(fraction))
  
  
  
  cfun <- function(a, b) cbind(a,b)
  lof = length(fraction)
  #---------- Cook's distance for all fractions selected 0.82s ----------#
  group = list(2:(lof%/%3), ((lof%/%3)+1):((lof%/%3)*2), (lof%/%3*2+1):lof)
  A = foreach (j=1:3, .combine = cfun) %dopar% {
    CD_Mat = c()
    get_xi <- function(w, h, n){
      # given omega, n and hkk, return xi
      n*(1-w)/(n-1+w-n*(1-w)*h)
    }
    
    get_w <- function(xi, h, n){
      # reverse function of get_xi, 
      # given xi, hkk and n, return omega
      (n-xi*(n-1-n*h))/(n+(1+n*h)*xi)}
    
    for (i in group[[j]]){
      print(i)
      # for (i in 1:length(fraction)){
      lambda = l_list[i]
      
      beta_hat_backup = beta_hat_table[,i]
      s_backup = sign(beta_hat_backup)
      A_backup = s_backup!=0
      # derivative of f in terms of covariates in nonactive set
      d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%X%*%beta_hat_backup
      
      A_id = paste(A_backup*1,collapse = '')
      if (A_id %in% names(XX_recorder)){
        XX_backup = XX_recorder[[A_id]]
      }
      else{
        XX_backup = solve(t(X[,A_backup])%*%X[,A_backup])
        XX_recorder[[A_id]] = XX_backup
      }
      
      if (sum(A_backup)==0){
        h_backup = matrix(0,ncol=n,nrow=n)
      }
      else{
        h_backup = X[,A_backup]%*%XX_backup%*%t(X[,A_backup])
      }
      
      CD_list = c()
      for (k in 1:n){
        beta_hat = beta_hat_backup
        s = s_backup
        A = A_backup
        
        d_hat = d_hat_backup
        XX = XX_backup
        hk = h_backup[,k]
        
        # current predictor of yk 
        yk_hat = drop(ybar + X[k,] %*% beta_hat)
        
        # beta path records beta hat's value at each breakpoint
        beta_path = c(beta_hat)
        # and sign change
        s_path = c(s)
        # and omega value at each breakpoint
        hkk_path = c()
        w = 1
        while (T){
          hkk_path = c(hkk_path, hk[k])
          xi = get_xi(w, hk[k], n)
          bias = yk_hat - y[k]
          if (sum(A) == 0){
            xi_cand1 = c()
          }
          else{
            slope_beta = XX%*%X[k,A]*bias
            xi_cand1 = -beta_hat[A]/slope_beta
          }
          slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
          # xi candidates
          # xi_cand1 = -beta_hat[A]/slope_beta
          xi_cand2p = (-lambda-d_hat)/slope_d
          xi_cand2m = (lambda-d_hat)/slope_d
          xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
          xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
          ind = which(xi_cand == xi_cand0)
          
          # update beta
          if (sum(A) > 0){
            beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
          }
          beta_path = rbind(beta_path, beta_hat)
          beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
          
          # if the xi is off the bound, stop the algorithm
          if (xi_cand0 > get_xi(0,hk[k],n)){
            hkk_path = c(hkk_path, hk[k])
            s_path = rbind(s_path, s)
            break
          }
          # if not, locate the covariate (go to or leave the active set)
          else if(ind<=sum(A)){
            coef = which(cumsum(A)==ind)[1]
            A[coef] = F
            s[coef] = 0
            # check if active set is empty 
            if (sum(A) == 0){
              hkk_path = c(hkk_path, 0, 0)
              s_path = rbind(s_path, s, s)
              beta_path = rbind(beta_path, beta_hat)
              break
            }
          }
          else if(ind>p){
            coef = which(cumsum(1-A)==(ind-p))[1]
            A[coef] = T
            s[coef] = 1
          }
          else{
            coef = which(cumsum(1-A)==(ind-sum(A)))[1]
            A[coef] = T
            s[coef] = -1
          }
          
          s_path = rbind(s_path, s)
          # update omega, XX, beta_hat, d_hat
          w = get_w(xi_cand0, hk[k], n)
          A_id = paste(A*1,collapse = '')
          if (A_id %in% names(XX_recorder)){
            XX = XX_recorder[[A_id]]
          }
          else{
            XX = solve(t(X[,A])%*%X[,A])
            XX_recorder[[A_id]] = XX
          }
          beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
          beta_hat[!A] = 0
          d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X%*%beta_hat
          hk = X[,A]%*%XX%*%X[k,A]
          yk_hat = drop(ybar + X[k,] %*% beta_hat)
        }
        
        #--------------------
        
        K = length(hkk_path)
        A = s_path[1,]!=0
        if (sum(A)==1){
          y_tilde = X[,A]*beta_path[1,A] + ybar
        } 
        else{
          y_tilde = X[,A]%*%beta_path[1,A] + ybar
        }
        
        A = s_path[K,]!=0
        A_id = paste(A*1,collapse = '')
        if (sum(A) == 0){
          y_last = ybar + get_xi(0, hkk_path[K], n)*1/n*(ybar - y[k])
        }
        else{
          xxx =  X[,A]%*%XX_recorder[[A_id]]
          y_hat =  xxx%*%(t(X[,A])%*%y-lambda*s_path[K,A]) + ybar
          y_last = y_hat + get_xi(0, hkk_path[K], n)*(xxx%*%X[k,A]+1/n)*(y_hat[k] - y[k])
        }
        
        CD_list = c(CD_list, sum((y_last - y_tilde)**2))
      }
      # CD_Mat[,i] = CD_list/denom
      CD_Mat = cbind(CD_Mat, CD_list/denom)
    }
    print(CD_Mat)
  }
  #--------------- plot case influence graph 0.42s----------#
  A = cbind(cd0,A)
  plot_table = t(rbind(fraction, A))
  threshold_table = cbind(fraction, apply(A,2,mean)*3)
  colnames(plot_table) = c('fraction',1:n)
  colnames(threshold_table)=c('fraction','val')
  df = as_tibble(plot_table)%>%gather('obs','val',-fraction)
  df2 = as_tibble(threshold_table)%>%mutate(obs='threshold')
  p = ggplot(data=df,aes(x=fraction,y=val,group=obs,color=obs)) +
    geom_line() + ggtitle("Cook's distance") + xlim(c(0,NA))+ labs(x = "|coef|/max|coef|")+
    geom_line(data=df2,aes(x=fraction, y=val),linetype = "dotted", color = 'red1')+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +labs(y="cook's distance")
  print(p)
  return(A)
}
