# Departamento de Estatística e Matemática Aplicada
# Universidade Federal de Campina Grande
# Orientador : Prof . Dr . Manoel Santos-Neto
# Autores : Jaiany Nunes, Luan Sousa, Marcus

# Pacotes
library(gamlss)
library(MonteCarlo)
library(tidyverse)
library(simhelpers)
library(boot)
library(LaplacesDemon)

## Carregando a funcoes de ZAIGA.R
devtools::source_url("https://github.com/statlab-oficial/ZAIGA/blob/main/ZAIGA.R?raw=TRUE")

# Funcoes auxiliares
my.gamlss <- function (...) tryCatch (expr = gamlss(...), error = function(e) NA)

est <- function(n, coef_mu, coef_sigma, coef_nu){
  n <- 100
  coef_mu <- c(0.5, 1.0, 2.5)
  coef_sigma <- c(1.1, 2.0)
  coef_nu <- 0.2
  x1 <- runif(n)
  x2 <- runif(n)
  X <- model.matrix(~x1+x2)
  Z <- model.matrix(~x2)
  eta <- as.vector(X%*%coef_mu)
  eta1 <- as.vector(Z%*%coef_sigma)
  mu   <- as.vector(exp(eta))
  sigma   <- as.vector(exp(eta1))
  nu <- invlogit(coef_nu)
  

  repeat{
    y <- mapply(rZAIGA, n = 1,  mu, sigma, nu)
    data_sim <- data.frame(y = y, x1 = x1, x2 = x2)
    conh0     <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
    fit <- gamlss(y ~ x1 + x2, 
      sigma.fo = ~ x2, 
      nu.fo = ~1, 
      family = ZAIGA(),
      data = data_sim,
      control= conh0)
  }

}