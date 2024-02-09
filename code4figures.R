library(tidyverse)
require(ggpubr)
require(latex2exp)
library(mnormt)
#### figure 1 relationship between xi and omega ####
plot_xi <- function(w, h){
  # given omega and hkk, return xi
  (1-w)/(1-(1-w)*h)
}
w = 0:100/100
lev = c(0.01,0.2,0.5)
df = data.frame(cbind(c(plot_xi(w,0.01), plot_xi(w,0.2), 
                        plot_xi(w,0.5)), rep(w,length(lev)), rep(lev,each = 101)))
colnames(df) = c('xi','w','leverage')
df$leverage = factor(df$leverage)
ggplot(data=df,aes(x=w,y=xi,linetype=leverage))+geom_line()+
  xlab(TeX("Case weight $\\omega$"))+ylab(TeX("$\\xi(\\omega)$")) + 
  guides(linetype=guide_legend(title=TeX("$h_{kk}$")))+theme(panel.background = element_rect(fill = "white"),
                                                           # panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                                                           # panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                                                           panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

#### figure 2 solution path showcase using diabetes data ####
data(diabetes)
attach(diabetes)
# x = centralize(x)
beta_path(x,y,30,lambda=100,plot=2)


## figure 2 needs redraw using synthetic dataset
set.seed(100)
x=matrix(rnorm(20*5),nrow=20)
y = x%*%c(10,0,0,0,0)
x[1,] = c(30,10,10,2,1)/5
y[1] = x[1,]%*%c(-10,10,5,-8,0)-2
x = centralize(x)
beta_path(x,y,1,lambda=1.5,plot=2)


### figure for presentation 
x_values <- seq(0,1, length.out = 101)
l_list = CD_all_concise(x,y,0.01)
result_matrix = NULL
for (l in l_list){
  result = beta_path(x,y,1,lambda=l,plot=0)
  breakpoints = result$w_path
  ypiece = apply(result$s_path,1,function(x) sum(abs(x)))
  y_values_between_breakpoints = as.vector(ypiece[1:(length(ypiece)-1)])

  
  # Initialize y-values vector
  y_values <- numeric(length(x_values))
  
  # Assign y-values based on piecewise constant function
  for (i in 1:(length(breakpoints) - 1)) {
    indices <- which(x_values <= breakpoints[i] & x_values > breakpoints[i + 1])
    y_values[indices] <- y_values_between_breakpoints[i]
  }
  
  # Handle the last segment if needed
  indices <- which(x_values <= breakpoints[length(breakpoints) - 1])
  y_values[indices] <- y_values_between_breakpoints[length(y_values_between_breakpoints)]
  
  # Create a matrix with x and y values
  result_matrix <- cbind(result_matrix, y_values)
}

library(reshape2) # For melt function


# Convert to long format
colnames(result_matrix) = frac
rownames(result_matrix) = seq(0,1, length.out = 101)
mat_melt <- melt(result_matrix)
colnames(mat_melt) = c('w','fraction','value')

# Plot
ggplot(mat_melt, aes(x = fraction, y = w, fill = value)) +
  geom_tile(height = 0.01,width=0.01) +
  scale_fill_gradient(low = "blue", high = "red")+ # Gradient from white to blue
  theme_minimal() +
  theme(axis.ticks = element_blank())+labs(fill = "|A|",x = "|coef|/max|coef|",y=TeX('$\\omega$')) 
  # ggtitle(TeX('$A(\\lambda,\\omega)$'))
  ylab(TeX("$\\xi(\\omega)$")) 




#### figure 3 case influence graph mechanism ####
set.seed(1)
n = 10
root_x1 = rnorm(n,0,1)
root_x2 = rnorm(n,0,1)
root_error =  rnorm(n, 0, 1)

# scenario 1
x1 = root_x1
x1[10] = mean(root_x1[-10])
x2 = root_x2
x2[10] = 4
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = mean(y[-10])
result = CD_all_concise(X,y,0.001)
df2 = data.frame(fraction = result$fraction, threshold = sqrt(apply((result$CD_Mat)[1:9,],2,var)/2)*qchisq(0.95,1))

p1 = CD_one_concise(X,y,10,0.001)+ggtitle("Scenario II")+
  geom_line(data=df2,aes(x=fraction, y=threshold),linetype = "dotted", color = 'red1')+
  theme(panel.background = element_rect(fill = "white"),
 # panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
 # panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
 panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# +annotate("text", x = 0, y = 0, label = "Scenario 2",size=4,hjust = 0, vjust = 0)

# scenario 2
x1 = root_x1
x1[10] = mean(root_x1[-10])
x2 = root_x2
x2[10] = 4
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = lm(y[-10]~X[-10,])$coefficients%*%c(1,X[10,])
p2 = CD_one_concise(X,y,10,0.001)+ggtitle("Scenario III")+
  theme(panel.background = element_rect(fill = "white"),
  # panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
  # panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_line(data=df2,aes(x=fraction, y=threshold),linetype = "dotted", color = 'red1')


# scenario 3
x1 = root_x1
x1[10] = mean(root_x1[-10])
x2 = root_x2
x2[10] = mean(root_x2[-10])
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = 4
p3 = CD_one_concise(X,y,10,0.001)+ggtitle("Scenario I")+
  theme(panel.background = element_rect(fill = "white"),
  # panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
  # panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_line(data=df2,aes(x=fraction, y=threshold),linetype = "dotted", color = 'red1')

# scenario 4
x1 = root_x1
x1[10] = mean(root_x1[-10])+4
x2 = root_x2
x2[10] = mean(root_x2[-10])
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = -4
p4 = CD_one_concise(X,y,10,0.001)+ggtitle("Scenario IV")+
  theme(panel.background = element_rect(fill = "white"),
 # panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
 # panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
 panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
geom_line(data=df2,aes(x=fraction, y=threshold),linetype = "dotted", color = 'red1')

ggarrange(p3+rremove('xlab'), 
          p1+rremove('ylab')+rremove('xlab'), 
          p2, 
          p4+rremove('ylab'), 
          ncol = 2, nrow = 2)

#### figure 4 case influence graph for highly correlated features 3 group 9 features ###
# i=6
# i=16
# i=63
# i=73
set.seed(73)
n <- 50                                     
mu <- rep(0,6)                                   
v <- diag(6)
rho = 0.5
v[1,2] = rho
v[2,1] = rho
v[3,4] = rho
v[4,3] = rho
v[5,6] = rho
v[6,5] = rho
betac = 10
X <- centralize(x = mvrnorm(n, mu, Sigma = v))
beta = c(betac,0,betac,0,betac,0)
# n = 20
# # Sigma = matrix(c(1,0.99,0,0,0.99,1,0,0,0,0,1,0.6,0,0,0.6,1), nrow=4)
# # X = centralize(mvrnorm(n = n, rep(0,4), Sigma))
# x1 = rnorm(n)
# x2 = x1+rnorm(n,0,0.1)
# x3 = x1+rnorm(n,0,0.1)
# x4 = rnorm(n)
# x5 = x4+rnorm(n,0,0.1)
# x6 = x4+rnorm(n,0,0.1)
# x7 = rnorm(n)
# x8 = x7+rnorm(n,0,0.1)
# x9 = x7+rnorm(n,0,0.1)
# X = centralize(cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9))
# beta = c(10,0,0,10,0,0,10,0,0)

y = X %*% beta + rnorm(n,0,2)
par(mfrow=c(2,2))
# par(mfrow=c(1,1))
# mai = 0.65
par(mar = c(4, 4, 2, 2) + 0.1)
# par(mai = c(1, 0.8, 0.8, 0.4) + 0.02)
result = CD_all_concise(X,y,0.01,F)
##fig4-1
plot.lars <-
  function(x, xvar=c("norm","df","arc.length","step"), breaks = TRUE, plottype = c("coefficients", "Cp"), 
           omit.zeros = TRUE, eps = 1e-10, ...)
  {
    object <- x
    plottype <- match.arg(plottype)
    xvar <- match.arg(xvar)
    coef1 <- object$beta	### Get rid of many zero coefficients
    if(x$type!="LASSO"&&xvar=="norm")# problem with discontinuity in norm
      coef1=betabreaker(x)
    stepid=trunc(as.numeric(dimnames(coef1)[[1]]))
    coef1 <- scale(coef1, FALSE, 1/object$normx)
    if(omit.zeros) {
      c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))
      nonzeros <- c1 > eps
      cnums <- seq(nonzeros)[nonzeros]
      coef1 <- coef1[, nonzeros,drop=FALSE]
    }
    else cnums <- seq(ncol(coef1))
    s1<-switch(xvar,
               norm={
                 s1 <- apply(abs(coef1), 1, sum)
                 s1/max(s1)
               },
               df=object$df,
               arc.length=cumsum(c(0,object$arc.length)),
               step=seq(nrow(coef1))-1
    )
    xname<-switch(xvar,
                  norm="|coef|/max|coef|",
                  df="Df",
                  arc.length="Arc Length",
                  step="Step"
    )
    
    if(plottype == "Cp") {
      Cp <- object$Cp
      plot(s1, Cp, type = "b", xlab="|coef|/max|coef|", ...)
    }
    else {
      matplot(s1, coef1, xlab = '', ..., type = "b", lwd = 2,lty=1,
              pch = "*", ylab = "Coefficients")
      abline(h = 0, lty = 3)
      axis(4, at = coef1[nrow(coef1),  ], labels = paste(cnums
      ), cex.axis = 1, adj = 0)
      if(breaks) {
        axis(3, at = s1, labels = paste(stepid),cex.axis=1)
        abline(v = s1)
      }
      
    }
    invisible()
  }
plot(lars(X,y))


##fig4-2
cd = apply(result$CD_Mat,2,mean)
cd.error = apply(result$CD_Mat,2,sd)
plot(result$fraction, cd,type='l',
     xlab = '',ylab = "Mean Cook's distance",lwd=1)
# lines(result$fraction, apply(result$CD_Mat,2,mean)/result$gdf,type='l',col='blue')
# lines(result$fraction, apply(result$CD_Mat,2,mean)/apply(abs(result$beta_hat),2,sum),type='l',col='red')
abline(v = result$fraction[which.min(cd[0:100])],lty=2,col='gray',lwd=3)



#####
# error.bars(result$fraction, cd + cd, cd - cd.error, 
#            width = 1/length(result$fraction))
# error.bars <-
#   function(x, upper, lower, width = 0.02, ...)
#   {
#     xlim <- range(x)
#     barw <- diff(xlim) * width
#     segments(x, upper, x, lower, ...)
#     segments(x - barw, upper, x + barw, upper, ...)
#     segments(x - barw, lower, x + barw, lower, ...)
#     range(upper, lower)
#   }
#####

##fig4-3
matplot(result$fraction, t(result$CD_Mat),xlab = "|coef|/max|coef|",ylab = "Cook's distance under Lasso", lty = 1,type = "l")

##fig4-4
cverror = cv.lars(X,y,K=n, se=F)
abline(v = which.min(cverror$cv)/100 ,lty=2,col='gray',lwd=3)

plotCVLars <-
  function(cv.lars.object,se=TRUE){
    mode=cv.lars.object$mode
    xlab=switch(mode,
                fraction="|coef|/max|coef|",
                step="Number of steps"
    )
    index=cv.lars.object$index
    cv=cv.lars.object$cv
    cv.error=cv.lars.object$cv.error
    plot(index, cv, type = "l", ylim = range(cv, cv + cv.error, 
                                             cv - cv.error),xlab='|coef|/max|coef|',ylab="Cross-Validated MSE")
    if(se)
      error.bars(index, cv + cv.error, cv - cv.error, 
                 width = 1/length(index))
    invisible()
  }
cv.lars <-
  function(x, y, K = 10, index, 
           trace = FALSE, plot.it = TRUE, se = TRUE,type = c("lasso", "lar", "forward.stagewise", "stepwise"),
           mode=c("fraction", "step"),...)
  {
    type=match.arg(type)
    
    if(missing(mode)){
      mode=switch(type,
                  lasso="fraction",
                  lar="step",
                  forward.stagewise="fraction",
                  stepwise="step"
      )
    }
    else  mode=match.arg(mode)
    all.folds <- cv.folds(length(y), K)
    if(missing(index)){
      index=seq(from = 0, to = 1, length = 100)
      if(mode=="step"){
        fit=lars(x,y,type=type,...)
        nsteps=nrow(fit$beta)
        maxfold=max(sapply(all.folds,length))
        nsteps=min(nsteps,length(y)-maxfold)
        index=seq(nsteps)
      }
    }
    residmat <- matrix(0, length(index), K)
    for(i in seq(K)) {
      omit <- all.folds[[i]]
      fit <- lars(x[ - omit,,drop=FALSE  ], y[ - omit], trace = trace, type=type,...)
      fit <- predict(fit, x[omit,  ,drop=FALSE], mode = mode, s = index
      )$fit
      if(length(omit)==1)fit<-matrix(fit,nrow=1)
      residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
      if(trace)
        cat("\n CV Fold", i, "\n\n")
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object<-list(index = index, cv = cv, cv.error = cv.error,mode=mode)
    if(plot.it) plotCVLars(object,se=se)
    invisible(object)
  }
# abline(v = result$fraction[which.min(apply(result$CD_Mat,2,mean)/apply(abs(result$beta_hat),2,sum))],lty=2,col='blue')
# abline(v = result$fraction[which.min(apply(result$CD_Mat,2,mean)/result$gdf)],lty=2,col='red')



#### figure 5 case influence graph for 4 different measures.

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

case_influence_plotter = function(fraction, CD_Mat,title){
  plot_table = t(rbind(fraction, CD_Mat))
  # threshold_table = cbind(fraction, apply(CD_Mat,2,mean)*3)
  colnames(plot_table) = c('fraction',1:n)
  # colnames(threshold_table)=c('fraction','val')
  df = as_tibble(plot_table)%>%gather('obs','val',-fraction)
  # df2 = as_tibble(threshold_table)%>%mutate(obs='threshold')
  p = ggplot(data=df,aes(x=fraction,y=val,group=obs,color=obs)) +
    geom_line() + xlim(c(0,NA))+ labs(x = "|coef|/max|coef|")+ggtitle(title)+
    theme(plot.title = element_text(vjust = 0,size = 10,hjust = 0.5),legend.position = "none")+
    # geom_line(data=df2,aes(x=fraction, y=val),linetype = "dotted", color = 'red1')+
    labs(y="influence measure")+theme(panel.background = element_rect(fill = "white"),
                                      panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                                      panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                                      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  
}

df<-read.table("prostate.txt")
X = centralize(as.matrix(df[,1:8]))
y = as.matrix(df[,9])

p = dim(X)[2]
n = dim(X)[1]

obj = lars(X,y,type='lasso')
denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)

# CD exact
cd_exact = CD_all_concise(X,y,finesse = 0.002)

# CD approx JKSS
lambda_list = cd_exact$Lambda_list
cd_approx = c()
cd_approx2 = c()
cd_local = c()
for (l in lambda_list){
  # regular lasso estimate with lambda given
  beta_hat = predict.lars(obj,mode='lambda', s=l, type = 'coefficients')$coefficients
  index = which(beta_hat!=0)
  # W^-
  beta_hat_inv = 1/abs(beta_hat)
  beta_hat_inv[which(is.infinite(beta_hat_inv))] = 0
  
  resi = y - mean(y) - X%*%beta_hat 
  H = X%*%solve(t(X)%*%X + l*diag(beta_hat_inv))%*%t(X)
  H2 = cbind(1,X[,index])%*%solve(t(cbind(1,X[,index]))%*%cbind(1,X[,index]))%*%t(cbind(1,X[,index]))
  H_vec = diag(H)
  H_vec2 = diag(H2)
  cd_approx = cbind(cd_approx,H_vec/(1-H_vec)**2*resi**2/denom)
  cd_approx2 = cbind(cd_approx2, H_vec2/(1-H_vec2)**2*resi**2/denom)
  cd_local = cbind(cd_local, H_vec2*resi**2/denom)
}

p1 = case_influence_plotter(cd_exact$fraction,cd_exact$CD_Mat,'Exact')
p4 = case_influence_plotter(cd_exact$fraction,cd_approx,'Ridge')
p3 = case_influence_plotter(cd_exact$fraction,cd_approx2,'Approx')
p2 = case_influence_plotter(cd_exact$fraction,cd_local,'Local')

ggarrange(p1+rremove('xlab'), 
          p2+rremove('ylab')+rremove('xlab'), 
          p3, 
          p4+rremove('ylab'), 
          ncol = 2, nrow = 2)


### figure 6 needle plot of 4 measures
new_slice_pot = function(a){
  df = cbind(1:dim(cd_approx)[1], cd_exact$CD_Mat[,a], cd_approx[,a], cd_approx2[,a], cd_local[,a])
  colnames(df) = c('index','Exact','JKSS','Approx','Local')
  beta_hat = predict.lars(obj,mode='lambda', s=lambda_list[a], type = 'coefficients')$coefficients
  threshold = qf(1-0.05**(1/97),sum(beta_hat!=0)+1,97-8-1)*(sum(beta_hat!=0)+1)/8
  # df = data.frame(df) %>% gather('method', 'cd',-index)
  df = data.frame(df)
  q1 = ggplot(data=df,aes(x=index, ymax=Exact,ymin=0)) + geom_linerange()+
    ggtitle('Exact')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)+
    geom_segment(aes(x = 3 , y = 0, xend = 3, yend = df[3,2]),color='red',linewidth=1.5)+
    # geom_segment(aes(x = 58 , y = 0, xend = 58, yend = df[58,2]),color='red',linewidth=1.5)+
    geom_segment(aes(x = 69 , y = 0, xend = 69, yend = df[69,2]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 95 , y = 0, xend = 95, yend = df[95,2]),color='red',linewidth=1.5)+
    geom_segment(aes(x = 96 , y = 0, xend = 96, yend = df[96,2]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 39 , y = 0, xend = 39, yend = df[39,2]),color='blue',linewidth=1.5)+
    # annotate("text", x = 25, y = max(df$Exact)*0.9, label = "Exact")+
    geom_hline(yintercept = threshold,color='red',linetype='dashed',alpha=0.5)
  q2 = ggplot(data=df,aes(x=index, ymax=JKSS,ymin=0)) + geom_linerange()+
    ggtitle('Ridge')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)+
    geom_segment(aes(x = 69 , y = 0, xend = 69, yend = df[69,3]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 96 , y = 0, xend = 96, yend = df[96,3]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 39 , y = 0, xend = 39, yend = df[39,3]),color='blue',linewidth=1.5)+
    # annotate("text", x = 25, y = max(df$JKSS)*0.9, label = "Ridge")+
    geom_hline(yintercept = threshold,color='red',linetype='dashed',alpha=0.5)
  q3 = ggplot(data=df,aes(x=index, ymax=Approx,ymin=0)) + geom_linerange()+
    ggtitle('Approx')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)+
    geom_segment(aes(x = 3 , y = 0, xend = 3, yend = df[3,4]),color='red',linewidth=1.5)+
    geom_segment(aes(x = 95 , y = 0, xend = 95, yend = df[95,4]),color='red',linewidth=1.5)+
    # annotate("text", x = 25, y = max(df$Approx)*0.9, label = "Approx")+
    geom_hline(yintercept = threshold,color='red',linetype='dashed',alpha=0.5)
  q4 = ggplot(data=df,aes(x=index, ymax=Local,ymin=0)) + geom_linerange()+
    ggtitle('Local')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)+
    # geom_segment(aes(x = 3 , y = 0, xend = 3, yend = df[3,5]),color='red')+
    # geom_segment(aes(x = 58 , y = 0, xend = 58, yend = df[58,5]),color='red')+
    # geom_segment(aes(x = 94 , y = 0, xend = 94, yend = df[94,5]),color='red')+
    # annotate("text", x = 25, y = max(df$Local)*0.9, label = "Local")+
    geom_hline(yintercept = threshold,color='red',linetype='dashed',alpha=0.5)
  # +xlab('Observation index')
  qq = ggarrange(q1+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+theme(panel.background = element_rect(fill = "white"),
                                                                                              panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                                                                                              panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                                                                                              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 q3+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+theme(panel.background = element_rect(fill = "white"),
                                                                                              panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                                                                                              panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                                                                                              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 q4+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+theme(panel.background = element_rect(fill = "white"),
                                                                                              panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                                                                                              panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                                                                                              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 q2+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+theme(panel.background = element_rect(fill = "white"),
                                                                                              panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                                                                                              panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                                                                                              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 ncol=1,nrow=4)
                 # q4+rremove("x.grid")+rremove('ylab'), ncol = 1, nrow = 4)
  # qq = annotate_figure(qq,top = text_grob(paste0('|coef|/max|coef|=',round(cd_exact$fraction[a],2)), size = 12))
}

print(new_slice_pot(315))

### figure 8 case influence graph of the diabetes data. 
data(diabetes)
attach(diabetes)
result = CD_all_concise(x,y,finesse=0.01,plot = T,threshold = T)
x = centralize(x)
fit = cv.glmnet(x,y,standardize = F)
result = CD_one_fraction(x,y,fit$lambda.min*nrow(x))

### simulation study 
n = 50
p_v = c(10,50,100,200)
a_v = c(0,2,5,10)
b_v = c(0,2,5,10)
p_v = c(10, 20, 100, 200)


p = 1000
a = 2
b = 2
beta = c(5,4,3,2,1,rep(0,p-5))

rho<-0.5
cormat<-matrix(0,p,p)
for(i in 1:p) {
  for(j in 1:p) {
    cormat[i,j] = rho^(abs(i-j))
  }
}
tf = c()
inclusion = c()
for (i in 1:100){
  print(i)
  x = rmnorm(n=n,mean=rep(0,p),varcov=cormat)
  y = x %*% beta + rnorm(n)
  
  x[1,p] = x[1,p] + a
  y[1] = y[1] + b
  x = centralize(x)
  fit = cv.glmnet(x,y,standardize = F)
  result = CD_one_fraction(x,y,fit$lambda.min*n)
  inclusion = c(inclusion,result$beta[p]!=0)
  tf = c(tf,result$CD_vec[1]>result$threshold)
  # result = CD_all_concise(x,y,finesse=0.05,plot = F)
  # cv_result = cv.lars(x,y,plot.it = F,use.Gram=F)
  # fraction = cv_result$index[which.min(cv_result$cv)]
  # fraction_index = findInterval(-fraction, -result$fraction)
  # dropped_tf = (result$beta_table[p,fraction_index]==0)
  # print(c(result$threshold[fraction_index,2], result$CD_Mat[1,fraction_index],dropped_tf,i))
  # tf = c(tf,result$threshold[fraction_index,2] < result$CD_Mat[1,fraction_index])
}
  

###### ### ### ###   justify the chisq1 threshold ###### ### ### ###
par(mfrow=c(2,4))

abpr_simu <- function(n,p,q,r,pic){
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
  y = x %*% beta + ep 
  x = centralize(x)
  fit = cv.glmnet(x,y,standardize = F)
  result = CD_one_fraction(x,y,fit$lambda.min*n)
  hist(result$CD_vec,main=paste(TeX('Example'),pic),xlab="")
}


set.seed(100)
nppair = t(matrix(c(100,10,200,10,100,500,200,200),nrow=2))
r_v = c(0,0.5)
counter = 0
pic = 0
for (np in 1:nrow(nppair)){
  for (r in r_v){
      n = nppair[np,1]
      p = nppair[np,2]
      pic = pic + 1
      abpr_simu(n, p, q, r, pic) 
      if (pic==7){
        break
      }
  }
}
hist(rchisq(1000,1),xlab='',main=TeX('Sampling Distribution of $\\chi^2_1$'))
# x <- seq(0,150)/10
# pdf<- dchisq(x,1)
# plot(x,pdf,type = 'l')

############# real data analysis ##############
data(diabetes)
attach(diabetes)
# result = CD_all_concise(x,y,finesse=0.01,plot = T,threshold = T)


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
ans_sum=c()

for (iter in 1:1000){
  ## total number of samples 442
  test_indices = sample(442,88)
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
  
  result = CD_one_fraction(x_train_c,y_train,fit$lambda.min*nrow(x_train))
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

set.seed(1)

result_before = c()
# selection_before = c()
result_after = c()
# selection_after = c()
lambda_before = c()
lambda_after = c()
x = centralize(x)

for (i in 1:10000){
  fit = cv.glmnet(x,y,standardize = F)
  est = glmnet(x,y,lambda=fit$lambda.min,standardize=FALSE,thresh=1e-16)
  
  result = CD_one_fraction(x,y,fit$lambda.min*nrow(x))
  noninfluential_indices = which(result$CD_vec<=result$threshold)
  x_trimmed = x[noninfluential_indices,]
  y_trimmed = y[noninfluential_indices]
  
  fit2 = cv.glmnet(x_trimmed,y_trimmed,standardize = F)
  est2 = glmnet(x_trimmed,y_trimmed,lambda=fit2$lambda.min,standardize=FALSE,thresh=1e-16)
  
  result_before = cbind(result_before,as.vector(coef(est)))
  # selection_before = cbind(selection_before,as.vector(coef(est))!=0)
  result_after = cbind(result_after,as.vector(coef(est2)))
  # selection_after = cbind(selection_after,as.vector(coef(est2))!=0)
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




