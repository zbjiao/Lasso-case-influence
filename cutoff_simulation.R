n = 50

a_v = c(0,2,5,10)
b_v = c(0,2,5,10)
p_v = c(10, 20, 200,1000)
r_v = c(0,0.5)
iter = 500

results = c()

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

abpr_simu <- function(a,b,p,r){
  beta = c(5,4,3,2,1,rep(0,p-5))
  if (r!=0){
    cormat<-matrix(0,p,p)
    for(i in 1:p) {
      for(j in 1:p) {
        cormat[i,j] = r^(abs(i-j))
      }
    }
  }
  else{
    cormat = diag(p)
  }
  x = rmnorm(n=n,mean=rep(0,p),varcov=cormat)
  y = x %*% beta + rnorm(n)
  x[1,p] = x[1,p] + a
  y[1] = y[1] + b
  x = centralize(x)
  fit = cv.glmnet(x,y,standardize = F)
  result = CD_one_fraction(x,y,fit$lambda.min*n)
  inclusion = result$beta[p]!=0
  tf0 = result$CD_vec[1]>result$threshold
  tf1 = sum(result$CD_vec[2:n]>result$threshold)/(n-1)
  return(c(inclusion, tf0, tf1))
}

for (a in a_v){
  for (b in b_v){
    for (p in p_v){
      for (r in r_v){
        for(iter in 1:iter){
          result <- abpr_simu(a, b, p, r) 
          results = rbind(results,c(a,b,p,r,iter,result))
        }
        print(c(a,b,p,r))
      }
    }
  }
}

result_df = data.frame(data=results)
colnames(result_df) = c('a','b','p','r','iter','inclusion','tf0','tf1')
# start at 12:34pm
write.csv(result_df, "cv_cd_cutoff_simu.csv")

summary_df = 
  result_df %>% mutate(tf_include = ((inclusion+tf0)==2)|((inclusion+tf0)==0)) %>%
  group_by(a,b,p,r) %>% summarise(outlier = mean(tf_include), others = mean(tf1))
# write.csv(summary_df, "cv_cd_cutoff_simu_summary.csv")

summary_df2 = 
  result_df %>% group_by(a,b,p,r) %>% summarise(outlier = mean(tf0), others = mean(tf1)) %>%
  ungroup() %>% filter(r==0.5) %>% dplyr::select(-r)
write.csv(summary_df, "cv_cd_cutoff_simu_summary.csv")

# tf = c()
# inclusion = c()
# for (i in 1:100){
#   print(i)
#   x = rmnorm(n=n,mean=rep(0,p),varcov=cormat)
#   y = x %*% beta + rnorm(n)
#   
#   x[1,p] = x[1,p] + a
#   y[1] = y[1] + b
#   x = centralize(x)
#   fit = cv.glmnet(x,y,standardize = F)
#   result = CD_one_fraction(x,y,fit$lambda.min*n)
#   inclusion = c(inclusion,result$beta[p]!=0)
#   tf = c(tf,result$CD_vec[1]>result$threshold)
#   # result = CD_all_concise(x,y,finesse=0.05,plot = F)
#   # cv_result = cv.lars(x,y,plot.it = F,use.Gram=F)
#   # fraction = cv_result$index[which.min(cv_result$cv)]
#   # fraction_index = findInterval(-fraction, -result$fraction)
#   # dropped_tf = (result$beta_table[p,fraction_index]==0)
#   # print(c(result$threshold[fraction_index,2], result$CD_Mat[1,fraction_index],dropped_tf,i))
#   # tf = c(tf,result$threshold[fraction_index,2] < result$CD_Mat[1,fraction_index])
# }

# --------------------- add based on mean --------------------- #
abpr_simu <- function(a,b,p,r){
  beta = c(5,4,3,2,1,rep(0,p-5))
  if (r!=0){
    cormat<-matrix(0,p,p)
    for(i in 1:p) {
      for(j in 1:p) {
        cormat[i,j] = r^(abs(i-j))
      }
    }
  }
  else{
    cormat = diag(p)
  }
  x = rmnorm(n=n,mean=rep(0,p),varcov=cormat)
  ep = rnorm(n)
  ep[1] = 0
  y = x %*% beta + ep
  x[1,p] = a
  y[1] = y[1] + b
  x = centralize(x)
  fit = cv.glmnet(x,y,standardize = F)
  result = CD_one_fraction(x,y,fit$lambda.min*n)
  inclusion = result$beta[p]!=0
  tf0 = result$CD_vec[1]>result$threshold
  tf1 = sum(result$CD_vec[2:n]>result$threshold)/(n-1)
  return(c(inclusion, tf0, tf1))
}

n = 50

a_v = c(0,2,5,10)
b_v = c(0,2,5,10)
p_v = c(10, 20, 200,1000)
r_v = c(0.5)
iter = 500

results = c()

for (a in a_v){
  for (b in b_v){
    if ((a==0) & (b==0)){
      next 
    }
    for (p in p_v){
      for (r in r_v){
        for(iter in 1:iter){
          result <- abpr_simu(a, b, p, r) 
          results = rbind(results,c(a,b,p,r,iter,result))
        }
        print(c(a,b,p,r))
      }
    }
  }
}
result_df2 = data.frame(data=results)
colnames(result_df2) = c('a','b','p','r','iter','inclusion','tf0','tf1')
write.csv(result_df2, "cv_cd_cutoff_simu2.csv")
view(result_df2 %>% group_by(a,b,p,r) %>% summarise(outlier = mean(tf0), others = mean(tf1)) %>%
  ungroup() %>%  dplyr::select(-r))
write.csv(summary_df, "cv_cd_cutoff_simu_summary2.csv")


abpr_simu(0,0,1000, r) 

# --------------------- details about inclusion or not --------------------- #

temp = result_df2 %>% group_by(a,b,p,r) %>% 
  summarise(outlier = mean(tf0), others = mean(tf1)) %>%
  ungroup() %>%  dplyr::select(-r)

view(result_df2 %>% filter(inclusion == 1)%>% group_by(a,b,p,r) %>% 
       summarise(outlier = mean(tf0), others = mean(tf1)) %>%
       ungroup() %>%  dplyr::select(-r) %>% left_join(temp,
       by = c('a'='a','b'='b','p'='p')))

a = 2
b = 2
p = 20
r = 0.5
beta = c(5,4,3,2,1,rep(0,p-5))
cormat<-matrix(0,p,p)
for(i in 1:p) {
  for(j in 1:p) {
    cormat[i,j] = r^(abs(i-j))
  }
}

x = rmnorm(n=n,mean=rep(0,p),varcov=cormat)
ep = rnorm(n)
ep[1] = 0
y = x %*% beta + ep
x[1,p] = a
y[1] = y[1] + b
x = centralize(x)
fit = cv.glmnet(x,y,standardize = F)
result = CD_one_fraction(x,y,fit$lambda.min*n)
temp = CD_all_concise(x,y,0.02)
inclusion = result$beta[p]!=0
tf0 = result$CD_vec[1]>result$threshold
tf1 = sum(result$CD_vec[2:n]>result$threshold)/(n-1)

# --------------------- add a on x5 --------------------- #
 
result_df3 = data.frame(data=results)
colnames(result_df3) = c('a','b','n','p','q','iter','inclusion','tf0','tf1')
write.csv(result_df2, "cv_cd_cutoff_simu3.csv")
summary_df = result_df3 %>% group_by(a,b,n,p,q) %>% summarise(outlier = mean(tf0), others = mean(tf1)) %>%
       ungroup()
view(result_df3 %>%group_by(a,b,n,p,q) %>% summarise(outlier = mean(tf0), others = mean(tf1)) %>%
  ungroup())
write.csv(summary_df, "cv_cd_cutoff_simu_summary3.csv")




########################## jcgs paper replicate ##########################
lambda.int<-function(cvob,mult) {
  #Interval of "plausible" lambda vaues
  #cvob:  fitted cv.glmnet object
  #mult:  width (SE's) of interval
  
  ind<-sort.list(cvob$cvm)[1]
  upper<-min(cvob$cvm)+mult*cvob$cvsd[ind]
  temp<-cvob$lambda[which(cvob$cvm<=upper)]
  
  return(c(min(temp),max(temp)))
}
df.cvpath<-function(fit1,fit2) {
  #Change in cv-path between fit1 and fit2
  #fit1: fitted cv.glmnet
  #fit2: fitted cv.glmnet
  
  if (length(fit1$cvm) == length(fit2$cvm)) {
    return(trapz(rev(fit1$lambda),rev(abs(fit1$cvm-fit2$cvm))))
  }
  else {
    comlam = intersect(fit1$lambda, fit2$lambda)
    idx1 = fit1$lambda%in%comlam
    idx2 = fit2$lambda%in%comlam        
    x = rev(fit1$lambda[idx1])
    y = rev(abs(fit1$cvm[idx1]-abs(fit2$cvm[idx2])))
    return(trapz(x,y))
  }
}
###############
lasso.influence<-function(x,y,mult,assign) {
  #Applies measures to data
  #x:       matrix of covariates
  #y:       response vector
  #mult:    length of "plausible" lambda interval
  #assign:  assignment of obs to folds for CV
  
  n<-length(y)
  resmat<-matrix(0,n,1)
  t1 = Sys.time()
  fitfull<-cv.glmnet(x,y,standardize=TRUE,foldid=assign)
  t2 = Sys.time()
  #df.regpath and df.cvpath based on "plausible" range of lambda
  lamint<-lambda.int(cvob=fitfull,mult=mult)
  #df.regpath and df.cvpath based on "plausible" range of lambda
  lamint<-exp(seq(log(lamint[1]),log(lamint[2]),length.out=100))
  fitfull<-cv.glmnet(x,y,standardize=TRUE,foldid=assign,lambda=lamint)
  t3 = Sys.time()
  for(i in 1:n){
    fitred<-cv.glmnet(x[-i,],y[-i],standardize=TRUE,foldid=assign[-i],lambda=lamint)
    resmat[i,1]<-df.cvpath(fitfull,fitred)
  }
  t4 = Sys.time()
  print(c(t4-t3,t3-t2,t2-t1))
  return(scale(resmat))
}
###############
onesim<-function(a,b,bad,nsims,assign,cormat,cutoff=1.96,xreal,multiple) {
  
  xreal[1,bad]<-a
  ymean<-xreal%*%c(1:5,rep(0,p-5))
  
  for(i in 1:nsims) {
    
    yreal<-rnorm(n=n,mean=ymean,sd=1)
    yreal[1]<-ymean[1]+b
    res<-abs(lasso.influence(x=xreal,y=yreal,mult=multiple,assign=assign))
    
    # write.table(t(ifelse(res[1,]>cutoff,1,0)),file="bad.txt",append=T,col.names=F,row.names=F)
    # write.table(t(apply(res[-1,]>cutoff,2,sum)),file="good.txt",append=T,col.names=F,row.names=F)
    # res<-apply(res,2,sort.list)[48:50,]
    # write.table(t(apply(res==1,2,sum)),file="absolute.txt",append=T,col.names=F,row.names=F)
  }
  
}

#########################

nsims<-10
n<-50
p<-1000
rho<-0.5
multiple<-2
bad=100
set.seed(1)
assign<-whichfolds(n,10)
cormat<-matrix(0,p,p)
for(i in 1:p) {
  for(j in 1:p) {
    cormat[i,j]<-rho^(abs(i-j))
  }
}
set.seed(2)
xfixed<-rmnorm(n=n,mean=rep(0,p),varcov=cormat)

t1 = Sys.time()
onesim(a=10,b=10,bad=bad,nsims=nsims,assign=assign,cormat=cormat,cutoff=1.96,xreal=xfixed,multiple=multiple)
t2= Sys.time()
t2-t1
# 22.6s

abpr_simu <- function(a,b,p,r){
  t1 = Sys.time()
  x[1,p] = a
  y[1] = y[1] + b
  x = centralize(x)
  t2 = Sys.time()
  fit = cv.glmnet(x,y,standardize = F)
  t3 = Sys.time()
  result = CD_one_fraction(x,y,fit$lambda.min*n)
  print(result$beta[1])
  t4 = Sys.time()
  inclusion = result$beta[p]!=0
  tf0 = result$CD_vec[1]>result$threshold
  tf1 = sum(result$CD_vec[2:n]>result$threshold)/(n-1)
  t5 = Sys.time()
  print(c(t5-t4,t4-t3,t3-t2,t2-t1))
  return(c(inclusion, tf0, tf1))
}

beta = c(5,4,3,2,1,rep(0,p-5))
cormat<-matrix(0,p,p)
for(i in 1:p) {
  for(j in 1:p) {
    cormat[i,j] = r^(abs(i-j))
  }
}
x = rmnorm(n=n,mean=rep(0,p),varcov=cormat)
ep = rnorm(n)
ep[1] = 0
y = x %*% beta + ep

n = 50
a = 10
b = 10
p = 1000
r = 0.5
iter = 10
t3 = Sys.time()
for(i in 1:iter){ abpr_simu(a, b, p, r) }
t4 = Sys.time()
t4-t3
# 3.76s 

################# discussion of matrix inverse #####################
f <- function (QR) {
  ## thin Q-factor
  Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
  ## QQ'
  tcrossprod(Q)
}
set.seed(0)
X <- matrix(rnorm(1200), 40, 30)
t1 = Sys.time()
qr_linpack <- qr.default(X)
H1 <- f(qr_linpack)
t2 = Sys.time()
qr_lapack <- qr.default(X, LAPACK = TRUE)
H2 <- f(qr_lapack)
t3 = Sys.time()

SVD <- svd(X)
H3 <- tcrossprod(SVD$u)
t4 = Sys.time()

L <- t(suppressWarnings(chol(crossprod(X), pivot = TRUE)))

## compute `Q'`
r <- attr(L, "rank")
piv <- attr(L, "pivot")
Qt <- forwardsolve(L, t(X[, piv]), r)

## P = QQ'
H4 <- crossprod(Qt)
t5= Sys.time()
H5 = X%*%solve(t(X)%*%X)%*%t(X)
t6 = Sys.time()
H6 = X%*%solve(crossprod(X))%*%t(X)
t7 = Sys.time()
step1 = qr(X)
q = qr.Q(step1)
H7 = q%*%t(q)
t8 = Sys.time()
print(c(t2-t1,t3-t2,t4-t3,t5-t4,t6-t5,t7-t6,t8-t7))



t1 = Sys.time()
ans1 = X %*% solve(crossprod(X), t(X))
t2 = Sys.time()
ans2 = X%*%solve(crossprod(X))%*%t(X)
t3 = Sys.time()
an3 = 
t4 = Sys.time()
print(t3-t2)
print(t2-t1)


