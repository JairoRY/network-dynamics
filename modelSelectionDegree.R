library(stats4)
library(pracma)
library(minpack.lm)
library(dplyr)
library(ggplot2)

# DEGREE DISTRIBUTION

get_AIC <- function(m2logL, K, N) {
  m2logL + 2 * K * N / (N - K - 1)
}

# Prepare output tables
param_table <- data.frame()
delta_table <- data.frame()

# Load degree sequence and compute statistics
degs <- read.table("degree_sequence.txt", header = FALSE)$V1
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

# Plots
distribution_df <- data.frame(degree = degs) %>%
  count(degree, name = "frequency") %>%
  mutate(p_k = frequency / sum(frequency))

best_model <- which.min(delta_table[1, ])
model_names <- c(
  "Displaced Poisson",
  "Displaced Geometric",
  "Zeta (γ=3)",
  "Zeta (fitted γ)",
  "Right-truncated Zeta",
  "Altmann"
)
best_model_name <- model_names[best_model]

params <- param_table[1, ]
lambda <- params$lambda
q      <- params$q
gamma1 <- params$gamma1
gamma2 <- params$gamma2
kmax   <- params$kmax
gamma3 <- params$gamma3
delta  <- params$delta

p_pois <- function(k) dpois(k, lambda) / (1 - exp(-lambda))
p_geom <- function(k) (1 - q)^(k - 1) * q
p_zeta3 <- function(k) k^(-3) / zeta(3)
p_zeta <- function(k) k^(-gamma1) / zeta(gamma1)
p_trunc <- function(k) k^(-gamma2) / sum((1:kmax)^(-gamma2))
p_altmann <- function(k) {
  Z <- sum((1:N)^(-gamma3) * exp(-delta * (1:N)))
  k^(-gamma3) * exp(-delta * k) / Z
}

model_list <- list(
  p_pois,
  p_geom,
  p_zeta3,
  p_zeta,
  p_trunc,
  p_altmann
)

fit_fun <- model_list[[best_model]]

distribution_df <- distribution_df %>%
  mutate(p_fit = fit_fun(degree))

# POTSER CALDRIA AFEGIR LA DISTRIBUCIÓ ESPERADA SI NO COINCIDEIX AMB L'OBTINGUDA
ggplot(distribution_df, aes(x = degree)) +
  geom_point(aes(y = p_k), color = "black", size = 2) +
  geom_line(aes(y = p_fit), color = "red", linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = paste("Empirical vs Best-Fit Model:", best_model_name),
    x = "Degree k",
    y = "P(k)"
  ) +
  theme_minimal(base_size = 14)
