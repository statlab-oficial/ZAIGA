# Departamento de Estatística e Matemática Aplicada
# Universidade Federal do Ceará
# Orientador : Prof . Dr . Manoel Santos-Neto
# Autores : Jaiany Nunes, Luan Sousa, Marcus

# Pacotes
library(gamlss)
library(MonteCarlo)
library(tidyverse)
library(boot)
library(simhelpers)
library(LaplacesDemon) 

## Carregando a funcoes de ZAIGA.R
devtools::source_url("https://github.com/statlab-oficial/ZAIGA/blob/main/ZAIGA.R?raw=TRUE")

# Funcoes auxiliares
my.gamlss <- function (...) tryCatch (expr = gamlss(...), error = function(e) NA)

est <- function(n, nu0){
  coef_mu <- c(0.5, 1.0, 2.5)
  coef_sigma <- c(1.1, 2.0)
  coef_nu <- nu0
  
  set.seed(10)
  x1 <- runif(n)
  set.seed(200)
  x2 <- runif(n)
  set.seed(100)
  z1 <- runif(n)
  
  X <- model.matrix(~x1+x2)
  Z <- model.matrix(~z1)
  eta <- as.vector(X%*%coef_mu)
  eta1 <- as.vector(Z%*%coef_sigma)
  vector_mu   <- as.vector(exp(eta))
  vector_sigma   <- as.vector(exp(eta1))
  vector_nu <- rep(invlogit(coef_nu), n)
  
  repeat{
   y <- mapply(rZAIGA, n = 1,  vector_mu, vector_sigma, vector_nu)
   data_sim <- data.frame(yi = y, x1i = x1, x2i = x2, z1i = z1)
    conh0 <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
    
    fit <- my.gamlss(yi ~ x1i + x2i, 
      sigma.fo = ~ z1i, 
      nu.fo = ~1, 
      family = ZAIGA(),
      data = data_sim,
      control= conh0)

    #var test
    gen_test <- (all(is.na(fit)) == FALSE)
    conv_test <- ifelse(all(is.na(fit)) == FALSE, fit$converged, FALSE) 
    nu_test <- ifelse(all(is.na(fit)) == FALSE, between(fit$nu.coefficients, 0, 1), FALSE) 
    mu_test <- ifelse(all(is.na(fit)) == FALSE, all(sign(fit$mu.coefficients) == sign(coef_mu)), FALSE) 
    sigma_test <- ifelse(all(is.na(fit)) == FALSE, all(sign(fit$sigma.coefficients) == sign(coef_sigma)), FALSE) 
    
    if(all(gen_test, conv_test, nu_test, mu_test, sigma_test) == TRUE) break
  }

  mle_mu <- unname(fit$mu.coefficients)
  mle_sigma <- unname(fit$sigma.coefficients)
  mle_nu <- unname(fit$nu.coefficients)

  out <- list("mu1" = mle_mu[1],
              "mu2" = mle_mu[2],
              "mu3" = mle_mu[3],
              "sigma1" = mle_sigma[1],
              "sigma2" = mle_sigma[2],
              "nu" = mle_nu)

  out
}

n_grid <- c(50, 100, 200)
nu_grid <- c(0.20, 0.50, 0.70)

param_list <- list("n" = n_grid, "nu0" = nu_grid)

MC_result <- MonteCarlo(est, 
nrep = 100,
ncpus = 1,
param_list = param_list)

df <- MakeFrame(MC_result) |> 
  mutate(truemu1 = 0.5, 
  truemu2 = 1.0,
  truemu3 = 2.5,
    truesigma1 = 1.1,
      truesigma2 = 2.0)

result_mu1 <- df |> 
  select(n, mu1, nu0, truemu1) |>
  group_by(n, nu0, truemu1) |>
  do(calc_absolute(., estimates = mu1, true_param = truemu1, criteria = c("bias", "rmse")))
print(result_mu1)

result_mu2 <- df |> 
  select(n, mu2, nu0, truemu2) |>
  group_by(n, nu0, truemu2) |>
  do(calc_absolute(., estimates = mu2, true_param = truemu2, criteria = c("bias", "rmse")))
print(result_mu2)

result_mu3 <- df |> 
  select(n, mu3, nu0, truemu3) |>
  group_by(n, nu0, truemu3) |>
  do(calc_absolute(., estimates = mu3, true_param = truemu3, criteria = c("bias", "rmse")))
print(result_mu3)

result_sigma1 <- df |> 
  select(n, sigma1, nu0, truesigma1) |>
  group_by(n, nu0, truesigma1) |>
  do(calc_absolute(., estimates = sigma1, true_param = truesigma1, criteria = c("bias", "rmse")))
print(result_sigma1)

result_sigma2 <- df |> 
  select(n, sigma2, nu0, truesigma2) |>
  group_by(n, nu0, truesigma2) |>
  do(calc_absolute(., estimates = sigma2, true_param = truesigma2, criteria = c("bias", "rmse")))
print(result_sigma2)

result_nu <- df |> 
  select(n, nu, nu0) |>
  group_by(n, nu0) |>
  do(calc_absolute(., estimates = nu, true_param = nu0, criteria = c("bias", "rmse")))
print(result_nu)

