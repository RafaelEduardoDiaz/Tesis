pois.HMM.confint <- function(mod, n, B = 25, alpha = 0.05){
  m <- mod$m
  lambda <- matrix(data=NA,nrow=B,ncol=m,dimnames=list(c(1:B),paste0("lambda",1:m)))
  Gamma <- matrix(data=NA,nrow=B,ncol=m*m,dimnames=list(c(1:B),
                  paste0("gamma",as.vector(outer(1:m, 1:m, paste, sep = "")))))
  for(i in 1:B){
    y <- pois.HMM.generate_sample(ns = n, mod = mod)
    model <- pois.HMM.mle(y,m, mod$lambda , mod$gamma , stationary = TRUE)
    lambda[i,] <- model$lambda
    Gamma[i,] <- c(model$gamma)}
  CI.lambda <- apply(X = lambda,MARGIN = 2, FUN = function(x)quantile(x,probs = c(alpha/2,1-alpha/2)))
  CI.Gamma <- apply(X = Gamma,MARGIN = 2, FUN = function(x)quantile(x,probs = c(alpha/2,1-alpha/2)))
  output <- rbind(cbind(mean = c(mod$gamma),t(CI.Gamma)),cbind(mean = mod$lambda,t(CI.lambda))) 
  output <- as.data.frame(output)
  output$Ancho <- output[,3] - output[,2]
  return(output)
}

bayes.stadist <- function(model, m){
  df <- as.data.frame(summary(model)$summary)
  A <- matrix(df$mean[1:(m*m)],m,m, byrow = TRUE)
  delta <- rep(1,m)%*%solve(diag(m) - A + matrix(1,m,m))
  return(c(delta))
}

bayes.viterbi <- function(x, model, m){
  n <- length(x)
  xi <- matrix(0, n, m)
  df <- as.data.frame(summary(model)$summary)
  A <- matrix(df$mean[1:(m*m)],m,m, byrow = TRUE)
  mu <- df$mean[(m*m+1):(m*m+m)]
  delta <- bayes.stadist(model, m)
  foo <- delta * dpois(x[1], mu)
  xi[1 ,] <- foo / sum(foo)
  for(i in 2:n){
    foo <- apply(xi[i - 1, ] * A, 2, max) * dpois(x[i], mu)
    xi[i, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for(i in (n - 1):1) iv[i] <- which.max(A[ ,iv[i + 1]] * xi[i, ])
  return (iv)
}

bayes.lforward <- function(x, model, m){
  n <- length(x)
  lalpha <- matrix(NA, n, m)
  df <- as.data.frame(summary(model)$summary)
  A <- matrix(df$mean[1:(m*m)],m,m, byrow = TRUE)
  mu <- df$mean[(m*m+1):(m*m+m)]
  delta <- bayes.stadist(model, m)
  foo <- delta * dpois(x[1], mu)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[1, ] <- lscale + log(foo)
  for(i in 2:n){
    foo <- foo %*% A * dpois(x[i], mu)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[i, ] <- log(foo) + lscale
  }
  return(lalpha)
}

bayes.state_prediction <- function(h = 1, x, model, m){
  n <- length(x)
  df <- as.data.frame(summary(model)$summary)
  A <- matrix(df$mean[1:(m*m)],m,m, byrow = TRUE)
  la <- t(bayes.lforward(x, model, m))
  c <- max(la[, n])
  llk <- c + log(sum(exp(la[, n] - c)))
  statepreds <- matrix(NA, ncol = h, nrow = m)
  foo <- exp(la[, n] - llk)
  for(i in 1:h){
    foo <- foo %*% A
    statepreds[, i] <- foo
  }
  return(t(statepreds))
}
