source("caseweightlasso2023.R")

data(diabetes)
attach(diabetes)

get_mse <- function(est,x_train,y_train,x_test,y_test){
  beta0 = est$a0
  beta_lasso = as.vector(est$beta)
  active_indices = which(beta_lasso!=0)
  x_train_active = x_train[,active_indices]
  beta_ols = solve(crossprod(x_train_active))%*%t(x_train_active)%*%y_train
  x_test_active = x_test[,active_indices]
  mse_ols = sum((x_test_active %*% beta_ols + beta0 - y_test)**2)
  mse_lasso = sum((x_test %*% beta_lasso + beta0 - y_test)**2)
  # print(sqrt(sum((x_train %*% beta_lasso + beta0 - y_train)**2)/400))
  return(c(mse_ols, mse_lasso))
}

set.seed(1)

test_mse = list()
for (alpha in c(0.01,0.025,0.05)){
  ans_sum=c()
  for (iter in 1:10000){
    ## total number of samples 442
    test_indices = sample(442,440*0.2)
    # train_indices = setdiff(1:442, test_indices)
    # test_indices = 391:442
    # train_indices = 1:390
    
    x_train = x[-test_indices,]
    y_train = y[-test_indices]
    x_test = x[test_indices,]
    y_test = y[test_indices]
    
    ############# regular procedure #################
    central_pack = centralize(x_train,r=T)
    x_train_c = central_pack$d
    x_test_c =  scale(x_test, central_pack$m, central_pack$v)
    
    fit = cv.glmnet(x_train_c,y_train,standardize = F)
    est = glmnet(x_train_c,y_train,lambda=fit$lambda.min,standardize=FALSE,thresh=1e-16)
    
    ans1 = get_mse(est, x_train_c,y_train,x_test_c,y_test)
    
    ############# caseinfluence procedure #################
    
    result = CD_one_fraction(x_train_c,y_train,fit$lambda.min*nrow(x_train),alpha = alpha)
    noninfluential_indices = which(result$CD_vec<=result$threshold)
    x_train_trimmed = x_train[noninfluential_indices,]
    y_train_trimmed = y_train[noninfluential_indices]
    
    central_pack = centralize(x_train_trimmed,r=T)
    x_train_trimmed_c = central_pack$d
    x_test_c =  scale(x_test, central_pack$m, central_pack$v)
    
    
    fit2 = cv.glmnet(x_train_trimmed_c,y_train_trimmed,standardize = F)
    est2 = glmnet(x_train_trimmed_c,y_train_trimmed,lambda=fit2$lambda.min,standardize=FALSE,thresh=1e-16)
    
    ans2 = get_mse(est2, x_train_trimmed_c,y_train_trimmed,x_test_c,y_test)
    
    ans_sum = rbind(ans_sum, c(ans1,ans2))
  }
  test_mse[[toString(alpha)]] = ans_sum
}

mean(test_mse$"0.05"[,2]/test_mse$"0.05"[,4])
mean(test_mse$"0.025"[,2]/test_mse$"0.025"[,4])
mean(test_mse$"0.01"[,2]/test_mse$"0.01"[,4])

set.seed(1)

result_before = c()
result_after = c()
lambda_before = c()
lambda_after = c()
x = centralize(x)

for (i in 1:1000){
  fit = cv.glmnet(x,y,standardize = F)
  est = glmnet(x,y,lambda=fit$lambda.min,standardize=FALSE,thresh=1e-16)
  
  result = CD_one_fraction(x,y,fit$lambda.min*nrow(x))
  noninfluential_indices = which(result$CD_vec<=result$threshold)
  x_trimmed = x[noninfluential_indices,]
  y_trimmed = y[noninfluential_indices]
  
  fit2 = cv.glmnet(x_trimmed,y_trimmed,standardize = F)
  est2 = glmnet(x_trimmed,y_trimmed,lambda=fit2$lambda.min,standardize=FALSE,thresh=1e-16)
  
  result_before = cbind(result_before,as.vector(coef(est)))
  result_after = cbind(result_after,as.vector(coef(est2)))
  lambda_before = c(lambda_before,fit$lambda.min*nrow(x))
  lambda_after = c(lambda_after,fit2$lambda.min*nrow(x_trimmed))
}

apply(result_before,1,sd)/apply(result_after,1,sd)
apply(result_before,1,sd)/abs(apply(result_before,1,mean))
apply(result_after,1,sd)/abs(apply(result_before,1,mean))
apply(result_after,1,sd)
apply(result_after,1,mean)


apply(selection_before,1,sd)/apply(selection_after,1,sd)


rda_ans = data.frame(mean1 = apply(result_before,1,mean),
                     mean2 = apply(result_after,1,mean),
                     sd1 = apply(result_before,1,sd),
                     sd2 = apply(result_after,1,sd),
                     nonzero_rate1 = apply(result_before, 1,function(x) mean(x!=0)),
                     nonzero_rate2 = apply(result_after, 1,function(x) mean(x!=0)),
                     pos_rate1 = apply(result_before, 1,function(x) mean(x>0)),
                     pos_rate2 = apply(result_after, 1,function(x) mean(x>0)))

rda_ans = rda_ans[2:11,]
rda_ans = rda_ans %>% mutate(sd_mean1 = sd1/abs(mean1), sd_mean2 = sd2/abs(mean1))
rownames(rda_ans) <- colnames(x)

write.csv(round(rda_ans,2), "real_data_analysis.csv")


mean(apply(result_before,2,function(x) sum(x!=0)))
mean(apply(result_after, 2,function(x) sum(x!=0)))


