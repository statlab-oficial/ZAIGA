# Departamento de Estatística e Matemática Aplicada
# Universidade Federal do Ceará
# Orientador : Prof . Dr . Manoel Santos-Neto
# Autores : Jaiany Nunes, Luan Sousa, Marcus Vinicius

####
#Codes  obtained from Bourguignon and Gallardo (2020)
dIGAMMA2 <- function (x, mu = 1, sigma = 0.5, log = FALSE) 
{
    if (any(mu < 0)) 
        stop(paste("mu must be greater than 0", "\n", ""))
    if (any(sigma <= 0)) 
        stop(paste("sigma must be greater than 0", "\n", ""))
    if (any(x < 0)) 
        stop(paste("x must be greater than 0", "\n", ""))
    lfy <- (sigma+2) * (log(mu) + log1p(sigma)) - lgamma(sigma+2) - (sigma + 3) * log(x) - mu*(sigma+1)/x
    if (log == FALSE) 
        fy <- exp(lfy)
    else fy <- lfy
    fy
}


pIGAMMA2<-function (q, mu = 1, sigma = 0.5, lower.tail = TRUE, log.p = FALSE) 
{
    if (any(mu <= 0)) 
        stop(paste("mu must be greater than 0", "\n", ""))
    if (any(sigma <= 0)) 
        stop(paste("sigma must be greater than 0", "\n", ""))
    if (any(q < 0)) 
        stop(paste("q must be greater than 0", "\n", ""))
    lcdf <- pgamma(((mu * (1 + sigma))/q), shape = (sigma+2), lower.tail = FALSE, 
        log.p = TRUE)
    if (log.p == FALSE) 
        cdf <- exp(lcdf)
    else cdf <- lcdf
    if (lower.tail == TRUE) 
        cdf <- cdf
    else cdf <- 1 - cdf
    cdf
}

qIGAMMA2 <-function (p, mu = 1, sigma = 0.5, lower.tail = TRUE, log.p = FALSE) 
{
    if (any(mu < 0)) 
        stop(paste("mu must be positive", "\n", ""))
    if (any(sigma < 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    if (log.p == TRUE) 
        p <- exp(p)
    else p <- p
    if (lower.tail == TRUE) 
        p <- p
    else p <- 1 - p
    if (any(p < 0) | any(p > 1)) 
        stop(paste("p must be between 0 and 1", "\n", ""))
	y <- qgamma(1-p, shape= (sigma+2), rate=mu*(sigma+1))

	1/y
}

###########################################################

ZAIGA <- function(mu.link = "log", sigma.link = "log", nu.link = "logit"){
  mstats <- checklink("mu.link", "ZAIGA", substitute(mu.link), 
      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "ZAIGA", substitute(sigma.link), 
      c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "ZAIGA", substitute(nu.link), 
      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  structure(list(family = c("ZAIGA", "Zero Adjusted Inverse GA"), parameters = list(mu = TRUE, 
      sigma = TRUE, nu = TRUE), nopar = 3, type = "Mixed", 
      mu.link = as.character(substitute(mu.link)), 
      sigma.link = as.character(substitute(sigma.link)), 
      nu.link = as.character(substitute(nu.link)), 
      mu.linkfun = mstats$linkfun, 
      sigma.linkfun = dstats$linkfun, 
      nu.linkfun = vstats$linkfun, 
      mu.linkinv = mstats$linkinv, 
      sigma.linkinv = dstats$linkinv, 
      nu.linkinv = vstats$linkinv, 
      mu.dr = mstats$mu.eta, 
      sigma.dr = dstats$mu.eta, 
      nu.dr = vstats$mu.eta, 
      dldm = function(y, mu, sigma) ifelse(y == 0, 0, (sigma+2)/mu - (sigma+1)/y ), 
      d2ldm2 = function(y, mu, sigma) ifelse(y == 0, 0, -(sigma+2)/(mu^2) ), 
      dldd = function(y, mu, sigma) ifelse(y == 0, 0, 1 + log(mu) + log1p(sigma) + 1/(sigma +1)
      - digamma(sigma +2) - log(y) - mu/y), 
      d2ldd2 = function(y, sigma) ifelse(y == 0, 0, sigma/(sigma+1)^2 - psigamma(sigma+2 ,1) ), 
      dldv = function(y, nu) ifelse(y == 0, 1/nu, -1/(1 - nu)), 
      d2ldv2 = function(nu) -1/(nu*(1 - nu)), 
      d2ldmdd = function (y , mu) ifelse(y == 0, 0, 1/mu - 1/y), 
      d2ldmdv = function(y) rep(0, length(y)), 
      d2ldddv = function(y) rep(0, length(y)), 
      G.dev.incr = function(y, mu, sigma, nu, ...) -2 * dZAIGA(y, mu, sigma, nu, log = TRUE), 
      rqres = expression(rqres(pfun = "pZAIGA", type = "Mixed", mass.p = 0, prob.mp = nu, y = y, mu = mu, sigma = sigma, nu = nu)), 
      mu.initial = expression(mu <- (y + mean(y))/2), 
      sigma.initial = expression(sigma <- rep(1, length(y))), 
      nu.initial = expression(nu <- rep(0.5, length(y))), 
      mu.valid = function(mu) TRUE, 
      sigma.valid = function(sigma) all(sigma > 0), 
      nu.valid = function(nu) all(nu > 0) && all(nu < 1), 
      y.valid = function(y) all(y >= 0), 
      mean = function(mu, sigma, nu) (1 - nu) * mu, 
      variance = function(mu, sigma, nu) (1 - nu) * (mu^2)*(1/sigma + nu)), 
      class = c("gamlss.family", "family"))
}

dZAIGA <- function(x, mu = 1, sigma = 1, nu = 0.1, log = FALSE){
  if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0) | any(nu > 1)) 
      stop(paste("nu must be between 0 and 1", "\n", ""))
  log.lik <- ifelse(x == 0, log(nu), log(1-nu) +  dIGAMMA2(x, mu, sigma, log = TRUE) )
  if (log == FALSE) 
      fy <- exp(log.lik)
  else fy <- log.lik
  fy <- ifelse(x < 0, 0, fy)
  fy
}

pZAIGA <- function(q, mu = 1, sigma = 1, nu = 0.1, lower.tail = TRUE, log.p = FALSE){
  if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0) | any(nu > 1)) 
      stop(paste("nu must be between 0 and 1", "\n", ""))
  cdf <- pIGAMMA2(q = q, mu = mu, sigma = sigma)
  cdf <- ifelse((q == 0), nu, nu + (1 - nu)*cdf)
  if (lower.tail == TRUE) 
      cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
      cdf <- cdf
  else cdf <- log(cdf)
  cdf <- ifelse(q < 0, 0, cdf)
  cdf
}

qZAIGA <- function (p, mu = 1, sigma = 1, nu = 0.1, lower.tail = TRUE, 
  log.p = FALSE){
  if (any(mu <= 0)) 
      stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
      stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0) | any(nu > 1)) 
      stop(paste("nu must be between 0 and 1", "\n", ""))
  ly <- max(length(p), length(mu), length(sigma), length(nu))
  p <- rep(p, length = ly)
  if (log.p == TRUE) 
      p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
      p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p >= 1)) 
      stop(paste("p must be between 0 and 1", "\n", ""))
  if (!(length(nu) %in% c(1, length(p)))) 
      stop(paste("nu is of length", length(nu), "\n", "Must be of lenght 1 or length(p) =", 
          length(p)))
  which_zero <- which(p <= nu)
  if (length(nu) == 1) 
      p[which_zero] <- nu
  else p[which_zero] <- nu[which_zero]
  return(qIGAMMA2((p - nu)/(1 - nu), mu = mu, sigma = sigma))
}


rZAIGA <- function (n, mu = 1, sigma = 1, nu = 0.1, ...){
  if (any(mu <= 0)) 
      stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
      stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0) | any(nu > 1)) 
      stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(n <= 0)) 
      stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qZAIGA(p, mu = mu, sigma = sigma, nu = nu, ...)
  r
}

plotZAIGA <- function (n_sample, mu = 1, sigma = 1, nu = 0.1, from = 0, n = 101, 
  main = NULL){
  x <-  rZAIGA(n_sample, mu, sigma, nu)  
  xp <- x[x>0]  
  h <- hist(xp, plot = FALSE)
  nu0 <- length(x[x==0])/length(x)
  h$counts <- (1-nu0)*(h$density) 
  po <- c(0)
  y = seq(from = 0.001, to = max(x), length.out = n)
  pdf <- dZAIGA(y, mu = mu, sigma = sigma, nu = nu)
  pr0 <- c(dZAIGA(0, mu = mu, sigma = sigma, nu = nu))

  if (is.null(main)) 
    main = "Zero Adj. Inverse Gamma"  
  plot(h, main = main, ylim = c(0, max(pdf, pr0)), xlab = "Sample", ylab = "Density")
  box()
  points(po, nu0, type = "h")
  points(po, nu0, type = "p", col = "red")
    
  lines(pdf ~ y, col = "blue")
  points(po, pr0, type = "h", col = "blue")
  points(po, pr0, type = "p", col = "blue", pch = 16)
}

meanZAIGA <- function (obj){
  if (obj$family[1] != "ZAIGA") 
      stop("the object do not have a ZAIGA distribution")
  meanofY <- (1 - fitted(obj, "nu")) * fitted(obj, "mu")
  meanofY
}




