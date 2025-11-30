library(stats4)
library(pracma)
library(minpack.lm)

# Utility / AIC functions

get_AIC <- function(m2logL, K, N) {
  m2logL + 2 * K * N / (N - K - 1)
}

AIC_nls <- function(n, RSS, p) {
  n * log(2 * pi) + n * log(RSS / n) + n + 2 * (p + 1)
}

# DEGREE DISTRIBUTION

# Prepare output tables
param_table <- data.frame()
delta_table <- data.frame()

# Load degree sequence and compute statistics
degs <- read.table("degrees_tmax.txt", header = FALSE)$V1
N <- length(degs)
ks <- 1:N
M <- sum(degs)
M_prime <- sum(log(degs))
C <- sum(sapply(degs, function(k) if (k > 1) sum(log(2:k)) else 0))

# Initial values for parameters
lambda0 <- M / N
q0 <- N / M
gamma0 <- 2
kmax0 <- max(degs)
delta0 <- 1

# Fit models
mle_pois <- mle(function(lambda) {
    -(M * log(lambda) - N * (lambda + log(1 - exp(-lambda))) - C)
}, start = list(lambda=lambda0), method="L-BFGS-B", lower=c(1e-6))

mle_geom <- mle(function(q) {
    -((M - N) * log(1 - q) + N * log(q))
}, start = list(q=q0), method="L-BFGS-B", lower=c(1e-6), upper=c(0.9999))

mle_zeta <- mle(function(gamma) {
    -( -gamma * M_prime - N * log(zeta(gamma)) )
}, start = list(gamma=gamma0), method="L-BFGS-B", lower=c(1.000001))

mle_trunc_zeta <- mle(function(gamma, kmax) {
    H <- sum((1:kmax)^(-gamma))
    -( -gamma * M_prime - N * log(H) )
}, start = list(gamma=gamma0, kmax=kmax0), method="L-BFGS-B", lower=c(1.000001, max(degs)), upper=c(10, N))

mle_altm <- mle(function(gamma, delta) {
    Z <- sum(ks^(-gamma) * exp(-delta * ks))
    -( -N * log(Z) - gamma * M_prime - delta * M)
}, start = list(gamma=gamma0, delta=delta0), method="L-BFGS-B", lower=c(0.000001, 0.000001), upper=c(10, 10))

# m2logL values
m2logL_vals <- c(
  attributes(summary(mle_pois))$m2logL,
  attributes(summary(mle_geom))$m2logL,
  -2 * (-3 * M_prime - N * log(zeta(3))), # fixed gamma=3
  attributes(summary(mle_zeta))$m2logL,
  attributes(summary(mle_trunc_zeta))$m2logL,
  attributes(summary(mle_altm))$m2logL
)

# Number of parameters
K_vals <- c(1, 1, 0, 1, 2, 2)

# Model selection through AICc
AIC_vals <- mapply(get_AIC, m2logL_vals, K_vals, N)
AIC_best <- min(AIC_vals)
Delta_vals <- AIC_vals - AIC_best

# Extract parameters
lambda <- attributes(summary(mle_pois))$coef[1]
q <- attributes(summary(mle_geom))$coef[1]
gamma1 <- attributes(summary(mle_zeta))$coef[1]
gamma2 <- attributes(summary(mle_trunc_zeta))$coef[1]
kmax <- attributes(summary(mle_trunc_zeta))$coef[2]
gamma3 <- attributes(summary(mle_altm))$coef[1]
delta <- attributes(summary(mle_altm))$coef[2]

# Add to tables
param_table <- rbind(param_table, data.frame(
  lambda = lambda,
  q = q,
  gamma1 = gamma1,
  gamma2 = gamma2,
  kmax = kmax,
  gamma3 = gamma3,
  delta = delta
))

delta_table <- rbind(delta_table, data.frame(
  Model1 = Delta_vals[1],
  Model2 = Delta_vals[2],
  Model3 = Delta_vals[3],
  Model4 = Delta_vals[4],
  Model5 = Delta_vals[5],
  Model6 = Delta_vals[6]
))

# Print results
rownames(param_table) <- NULL
cat("Parameter estimates\n")
print(param_table)

cat("AIC differences (Δ)\n")
print(delta_table)

# -----------------------------------------------------------------------------

# SCALING OF VERTEX DEGREE OVER TIME

f0 <- function(t, a) a*t
f1 <- function(t, a) a*sqrt(t)
f2 <- function(t, a, b) a*t^b
f3 <- function(t, a, c) a*exp(c*t)
f4 <- function(t, a, d1) a*log(t + d1)
f0p <- function(t, a, d) a*t + d
f1p <- function(t, a, d) a*sqrt(t) + d
f2p <- function(t, a, b, d) a*t^b + d
f3p <- function(t, a, c, d) a*exp(c*t) + d
f4p <- function(t, a, d1, d2) a*log(t + d1) + d2

# Prepare tables
row_names <- c("t1", "t10", "t100", "t1000")
model_names <- c("M0","M1","M2","M3","M4","M0+","M1+","M2+","M3+","M4+")
param_names <- c("M0_a","M1_a","M2_a","M2_b","M3_a","M3_c","M4_a","M4_d1","M0+_a","M0+_d","M1+_a","M1+_d","M2+_a","M2+_b","M2+_d","M3+_a","M3+_c","M3+_d","M4+_a","M4+_d1","M4+_d2")

table_s <- matrix(NA, nrow = length(row_names), ncol = length(model_names), dimnames = list(row_names, model_names))
table_AIC <- matrix(NA, nrow = length(row_names), ncol = length(model_names), dimnames = list(row_names, model_names))
table_delta <- matrix(NA, nrow = length(row_names), ncol = length(model_names), dimnames = list(row_names, model_names))
table_params <- matrix(NA, nrow = length(row_names), ncol = length(param_names), dimnames = list(row_names, param_names))

# Iterate over the time series
files <- list.files(pattern = "^timeseries")
for (fname in files) {
  vertex <- sub("timeseries_(t[0-9]+)\\.txt", "\\1", fname)
  ts <- read.table(fname, header=FALSE, col.names = c("t","degree"))
  time <- ts$t
  k_obs <- ts$degree
  n <- length(time)

  # Model 0
  linear_model <- lm(degree ~ t, data = ts)
  initial_a0 <- coef(linear_model)[1]
  m0 <- nlsLM(k_obs ~ f0(time, a), start = list(a = initial_a0), trace = TRUE)
  rss <- deviance(m0)
  p <- length(coef(m0))
  table_AIC[vertex, "M0"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M0"] <- sqrt(rss/df.residual(m0))
  table_params[vertex, "M0_a"] <- coef(m0)[1]

  # Model 1
  linear_model <- lm(degree ~ sqrt(t), data = ts)
  initial_a1 <- coef(linear_model)[1]
  m1 <- nlsLM(k_obs ~ f1(time, a), start=list(a = initial_a1), trace = TRUE)
  rss <- deviance(m1)
  p <- length(coef(m1))
  table_AIC[vertex, "M1"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M1"] <- sqrt(rss/df.residual(m1))
  table_params[vertex, "M1_a"] <- coef(m1)[1]

  # Model 2
  linear_model <- lm(log(degree) ~ log(t), data = ts)
  initial_a2 <- exp(coef(linear_model)[1])
  initial_b2 <- coef(linear_model)[2]
  m2 <- nlsLM(k_obs ~ f2(time, a, b), start = list(a = initial_a2, b = initial_b2), trace = TRUE)
  rss <- deviance(m2)
  p <- length(coef(m2))
  table_AIC[vertex, "M2"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M2"] <- sqrt(rss/df.residual(m2))
  table_params[vertex, "M2_a"] <- coef(m2)[1]
  table_params[vertex, "M2_b"] <- coef(m2)[2]

  # Model 3
  linear_model <- lm(log(degree) ~ t, data = ts)
  initial_a3 <- exp(coef(linear_model)[1])
  initial_c3 <- coef(linear_model)[2]
  m3 <- nlsLM(k_obs ~ f3(time, a, c), start = list(a = initial_a3, c = initial_c3), trace = TRUE)
  rss <- deviance(m3)
  p <- length(coef(m3))
  table_AIC[vertex, "M3"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M3"] <- sqrt(rss/df.residual(m3))
  table_params[vertex, "M3_a"] <- coef(m3)[1]
  table_params[vertex, "M3_c"] <- coef(m3)[2]

  # Model 4
  initial_d4 <- 0
  linear_model <- lm(degree ~ log(t), data = ts)
  initial_a4 <- coef(linear_model)[1]
  m4 <- nlsLM(k_obs ~ f4(time, a, d1), start = list(a = initial_a4, d1 = initial_d4), trace = TRUE)
  rss <- deviance(m4)
  p <- length(coef(m4))
  table_AIC[vertex, "M4"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M4"] <- sqrt(rss/df.residual(m4))
  table_params[vertex, "M4_a"] <- coef(m4)[1]
  table_params[vertex, "M4_d1"] <- coef(m4)[1]

  # Model 0+
  m0p <- nlsLM(k_obs ~ f0p(time, a, d), start = list(a = initial_a0, d = 0), trace = TRUE)
  rss <- deviance(m0p)
  p <- length(coef(m0p))
  table_AIC[vertex, "M0+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M0+"] <- sqrt(rss/df.residual(m0p))
  table_params[vertex, "M0+_a"] <- coef(m0p)[1]
  table_params[vertex, "M0+_d"] <- coef(m0p)[2]

  # Model 1+
  m1p <- nlsLM(k_obs ~ f1p(time, a, d), start = list(a = initial_a1, d = 0), trace = TRUE)
  rss <- deviance(m1p)
  p <- length(coef(m1p))
  table_AIC[vertex, "M1+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M1+"] <- sqrt(rss/df.residual(m1p))
  table_params[vertex, "M1+_a"] <- coef(m1p)[1]
  table_params[vertex, "M1+_d"] <- coef(m1p)[2]

  # Model 2+
  m2p <- nlsLM(k_obs ~ f2p(time, a, b, d), start = list(a = initial_a2, b = initial_b2, d = 0), control = nls.lm.control(maxiter=500), trace = TRUE)
  rss <- deviance(m2p)
  p <- length(coef(m2p))
  table_AIC[vertex, "M2+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M2+"] <- sqrt(rss/df.residual(m2p))
  table_params[vertex, "M2+_a"] <- coef(m2p)[1]
  table_params[vertex, "M2+_b"] <- coef(m2p)[2]
  table_params[vertex, "M2+_d"] <- coef(m2p)[3]

  # Model 3+
  m3p <- nlsLM(k_obs ~ f3p(time, a, c, d), start = list(a = initial_a3, c = initial_c3, d = 0), control = nls.lm.control(maxiter=500), trace = TRUE)
  rss <- deviance(m3p)
  p <- length(coef(m3p))
  table_AIC[vertex, "M3+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M3+"] <- sqrt(rss/df.residual(m3p))
  table_params[vertex, "M3+_a"] <- coef(m3p)[1]
  table_params[vertex, "M3+_c"] <- coef(m3p)[2]
  table_params[vertex, "M3+_d"] <- coef(m3p)[3]

  # Model 4+
  m4p <- nlsLM(k_obs ~ f4p(time, a, d1, d2), start = list(a = initial_a4, d1 = initial_d4, d2 = 0), trace = TRUE)
  rss <- deviance(m4p)
  p <- length(coef(m4p))
  table_AIC[vertex, "M4+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M4+"] <- sqrt(rss/df.residual(m4p))
  table_params[vertex, "M4+_a"] <- coef(m4p)[1]
  table_params[vertex, "M4+_d1"] <- coef(m4p)[2]
  table_params[vertex, "M4+_d2"] <- coef(m4p)[3]
  
  # Compute Δ (delta) = AIC - min AIC in that vertex
  rowA <- table_AIC[vertex, ]
  minA <- min(rowA, na.rm = TRUE)
  table_delta[vertex, ] <- rowA - minA
}

df_homo <- as.data.frame(table_homo)
df_s <- as.data.frame(table_s)
df_AIC <- as.data.frame(table_AIC)
df_delta <- as.data.frame(table_delta)
df_params <- as.data.frame(table_params)

cat("\n=== s ===\n")
print(df_s)
#digits <- c(0,2,2,2,2,2,2,2,2,2)
#latex_table <- xtable(df_s, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_s.tex", include.rownames = TRUE)

cat("\n=== AICs ===\n")
print(df_AIC)
#latex_table <- xtable(df_AIC, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_AIC.tex", include.rownames = TRUE)

cat("\n=== Δ AIC ===\n")
print(df_delta)
#latex_table <- xtable(df_delta, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_delta.tex", include.rownames = TRUE)

cat("\n=== Parameter estimates ===\n")
print(df_params)
#digits <- c(0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
#latex_table <- xtable(df_params, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_params.tex", include.rownames = TRUE)
========
library(stats4)
library(pracma)
library(minpack.lm)
library(dplyr)
library(ggplot2)

# SCALING OF VERTEX DEGREE OVER TIME

AIC_nls <- function(n, RSS, p) {
  n * log(2 * pi) + n * log(RSS / n) + n + 2 * (p + 1)
}

f0 <- function(t, a) a*t
f1 <- function(t, a) a*sqrt(t)
f2 <- function(t, a, b) a*t^b
f3 <- function(t, a, c) a*exp(c*t)
f4 <- function(t, a, d1) a*log(t + d1)
f0p <- function(t, a, d) a*t + d
f1p <- function(t, a, d) a*sqrt(t) + d
f2p <- function(t, a, b, d) a*t^b + d
f3p <- function(t, a, c, d) a*exp(c*t) + d
f4p <- function(t, a, d1, d2) a*log(t + d1) + d2

# Prepare tables
row_names <- c("1", "10", "100", "1000")
model_names <- c("M0","M1","M2","M3","M4","M0+","M1+","M2+","M3+","M4+")
param_names <- c("M0_a","M1_a","M2_a","M2_b","M3_a","M3_c","M4_a","M4_d1","M0+_a","M0+_d","M1+_a","M1+_d","M2+_a","M2+_b","M2+_d","M3+_a","M3+_c","M3+_d","M4+_a","M4+_d1","M4+_d2")

table_s <- matrix(NA, nrow = length(row_names), ncol = length(model_names), dimnames = list(row_names, model_names))
table_AIC <- matrix(NA, nrow = length(row_names), ncol = length(model_names), dimnames = list(row_names, model_names))
table_delta <- matrix(NA, nrow = length(row_names), ncol = length(model_names), dimnames = list(row_names, model_names))
table_params <- matrix(NA, nrow = length(row_names), ncol = length(param_names), dimnames = list(row_names, param_names))
kp_list <- list()

# Iterate over the time series
files <- list.files(pattern = "^timeseries")
msub0 <- 3 # AQUESTS S'HAURIEN D'OBTENIR DEL NOM DEL FITXER O DE LA CARPETA
nsub0 <- 10
for (fname in files) {
  vertex <- sub("timeseries_t([0-9]+)\\.txt", "\\1", fname)
  ts <- read.table(fname, header=FALSE, col.names = c("t","degree"))
  time <- ts$t
  k_obs <- ts$degree
  n <- length(time)
  
  ti <- as.numeric(vertex)
  kp_list[[vertex]] <- data.frame(
    t = time,
    k1 = sqrt(ti)*k_obs,
    k2 = k_obs + msub0*log(nsub0 + ti - 1) - msub0,
    k3 = k_obs,
    vertex = vertex
  )
  
  # Model 0
  linear_model <- lm(degree ~ t, data = ts)
  initial_a0 <- coef(linear_model)[1]
  m0 <- nlsLM(k_obs ~ f0(time, a), start = list(a = initial_a0), trace = TRUE)
  rss <- deviance(m0)
  p <- length(coef(m0))
  table_AIC[vertex, "M0"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M0"] <- sqrt(rss/df.residual(m0))
  table_params[vertex, "M0_a"] <- coef(m0)[1]
  
  # Model 1
  linear_model <- lm(degree ~ sqrt(t), data = ts)
  initial_a1 <- coef(linear_model)[1]
  m1 <- nlsLM(k_obs ~ f1(time, a), start=list(a = initial_a1), trace = TRUE)
  rss <- deviance(m1)
  p <- length(coef(m1))
  table_AIC[vertex, "M1"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M1"] <- sqrt(rss/df.residual(m1))
  table_params[vertex, "M1_a"] <- coef(m1)[1]
  
  # Model 2
  linear_model <- lm(log(degree) ~ log(t), data = ts)
  initial_a2 <- exp(coef(linear_model)[1])
  initial_b2 <- coef(linear_model)[2]
  m2 <- nlsLM(k_obs ~ f2(time, a, b), start = list(a = initial_a2, b = initial_b2), trace = TRUE)
  rss <- deviance(m2)
  p <- length(coef(m2))
  table_AIC[vertex, "M2"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M2"] <- sqrt(rss/df.residual(m2))
  table_params[vertex, "M2_a"] <- coef(m2)[1]
  table_params[vertex, "M2_b"] <- coef(m2)[2]
  
  # Model 3
  linear_model <- lm(log(degree) ~ t, data = ts)
  initial_a3 <- exp(coef(linear_model)[1])
  initial_c3 <- coef(linear_model)[2]
  m3 <- nlsLM(k_obs ~ f3(time, a, c), start = list(a = initial_a3, c = initial_c3), trace = TRUE)
  rss <- deviance(m3)
  p <- length(coef(m3))
  table_AIC[vertex, "M3"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M3"] <- sqrt(rss/df.residual(m3))
  table_params[vertex, "M3_a"] <- coef(m3)[1]
  table_params[vertex, "M3_c"] <- coef(m3)[2]
  
  # Model 4
  initial_d4 <- 0
  linear_model <- lm(degree ~ log(t), data = ts)
  initial_a4 <- coef(linear_model)[1]
  m4 <- nlsLM(k_obs ~ f4(time, a, d1), start = list(a = initial_a4, d1 = initial_d4), trace = TRUE)
  rss <- deviance(m4)
  p <- length(coef(m4))
  table_AIC[vertex, "M4"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M4"] <- sqrt(rss/df.residual(m4))
  table_params[vertex, "M4_a"] <- coef(m4)[1]
  table_params[vertex, "M4_d1"] <- coef(m4)[1]
  
  # Model 0+
  m0p <- nlsLM(k_obs ~ f0p(time, a, d), start = list(a = initial_a0, d = 0), trace = TRUE)
  rss <- deviance(m0p)
  p <- length(coef(m0p))
  table_AIC[vertex, "M0+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M0+"] <- sqrt(rss/df.residual(m0p))
  table_params[vertex, "M0+_a"] <- coef(m0p)[1]
  table_params[vertex, "M0+_d"] <- coef(m0p)[2]
  
  # Model 1+
  m1p <- nlsLM(k_obs ~ f1p(time, a, d), start = list(a = initial_a1, d = 0), trace = TRUE)
  rss <- deviance(m1p)
  p <- length(coef(m1p))
  table_AIC[vertex, "M1+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M1+"] <- sqrt(rss/df.residual(m1p))
  table_params[vertex, "M1+_a"] <- coef(m1p)[1]
  table_params[vertex, "M1+_d"] <- coef(m1p)[2]
  
  # Model 2+
  m2p <- nlsLM(k_obs ~ f2p(time, a, b, d), start = list(a = initial_a2, b = initial_b2, d = 0), control = nls.lm.control(maxiter=500), trace = TRUE)
  rss <- deviance(m2p)
  p <- length(coef(m2p))
  table_AIC[vertex, "M2+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M2+"] <- sqrt(rss/df.residual(m2p))
  table_params[vertex, "M2+_a"] <- coef(m2p)[1]
  table_params[vertex, "M2+_b"] <- coef(m2p)[2]
  table_params[vertex, "M2+_d"] <- coef(m2p)[3]
  
  # Model 3+
  m3p <- nlsLM(k_obs ~ f3p(time, a, c, d), start = list(a = initial_a3, c = initial_c3, d = 0), control = nls.lm.control(maxiter=500), trace = TRUE)
  rss <- deviance(m3p)
  p <- length(coef(m3p))
  table_AIC[vertex, "M3+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M3+"] <- sqrt(rss/df.residual(m3p))
  table_params[vertex, "M3+_a"] <- coef(m3p)[1]
  table_params[vertex, "M3+_c"] <- coef(m3p)[2]
  table_params[vertex, "M3+_d"] <- coef(m3p)[3]
  
  # Model 4+
  m4p <- nlsLM(k_obs ~ f4p(time, a, d1, d2), start = list(a = initial_a4, d1 = initial_d4, d2 = 0), trace = TRUE)
  rss <- deviance(m4p)
  p <- length(coef(m4p))
  table_AIC[vertex, "M4+"] <- AIC_nls(rss, n, p)
  table_s[vertex, "M4+"] <- sqrt(rss/df.residual(m4p))
  table_params[vertex, "M4+_a"] <- coef(m4p)[1]
  table_params[vertex, "M4+_d1"] <- coef(m4p)[2]
  table_params[vertex, "M4+_d2"] <- coef(m4p)[3]
  
  # Compute Δ (delta) = AIC - min AIC in that vertex
  rowA <- table_AIC[vertex, ]
  minA <- min(rowA, na.rm = TRUE)
  table_delta[vertex, ] <- rowA - minA
}

df_s <- as.data.frame(table_s)
df_AIC <- as.data.frame(table_AIC)
df_delta <- as.data.frame(table_delta)
df_params <- as.data.frame(table_params)

cat("\n=== s ===\n")
print(df_s)
#digits <- c(0,2,2,2,2,2,2,2,2,2)
#latex_table <- xtable(df_s, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_s.tex", include.rownames = TRUE)

cat("\n=== AICs ===\n")
print(df_AIC)
#latex_table <- xtable(df_AIC, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_AIC.tex", include.rownames = TRUE)

cat("\n=== Δ AIC ===\n")
print(df_delta)
#latex_table <- xtable(df_delta, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_delta.tex", include.rownames = TRUE)

cat("\n=== Parameter estimates ===\n")
print(df_params)
#digits <- c(0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
#latex_table <- xtable(df_params, type = "latex", row.names = TRUE, digits = digits)
#print(latex_table, file = "table_params.tex", include.rownames = TRUE)

#Plots
kp_all <- do.call(rbind, kp_list)
t_max <- max(kp_all$t)
t_theo <- seq(10, t_max, length.out = 1000)

theory_df <- data.frame(
  t = t_theo,
  k1 = msub0 * sqrt(t_theo),
  k2 = msub0*log(msub0 + t_theo - 1),
  k3 = 2*msub0*t_theo/nsub0
)
################ IF BARABÁSI-ALBERT ################ 
ggplot() +
  geom_line(data = subset(kp_all, t >= 10), aes(x = t, y = k1, color = vertex), size = 1) +
  geom_line(data = subset(theory_df, t >= 10), aes(x = t, y = k1), color = "black", size = 1.2, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Growth of vertex degree k'_i(t) with theoretical curve",
    x = "t",
    y = "k'_i(t)"
  )

################ IF RANDOM ATTACHMENT ################
ggplot() +
  geom_line(data = subset(kp_all, t >= 10), aes(x = t, y = k2, color = vertex), size = 1) +
  geom_line(data = subset(theory_df, t >= 10), aes(x = t, y = k2), color = "black", size = 1.2, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Growth of vertex degree k''_i(t) with theoretical curve",
    x = "t",
    y = "k''_i(t)"
  )

################ IF NO GROWTH ################
ggplot() +
  geom_line(data = subset(kp_all, t >= 10), aes(x = t, y = k3, color = vertex), size = 1) +
  geom_line(data = subset(theory_df, t >= 10), aes(x = t, y = k3), color = "black", size = 1.2, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Growth of vertex degree k_i(t) with theoretical curve",
    x = "t",
    y = "k_i(t)"
  )
