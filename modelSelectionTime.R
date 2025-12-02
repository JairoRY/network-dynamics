library(stats4)
library(pracma)
library(minpack.lm)
library(dplyr)
library(ggplot2)
library(xtable)

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

dirpath <- "results/barabasi_albert_n10_m3"

dirname <- basename(dirpath)
nsub0 <- as.numeric(sub(".*_n([0-9]+)_m[0-9]+.*", "\\1", dirname))
msub0 <- as.numeric(sub(".*_n[0-9]+_m([0-9]+).*", "\\1", dirname))
files <- list.files(dirpath, pattern = "^timeseries_t[0-9]+\\.txt$", full.names = TRUE)
for (fname in files) {
  vertex <- sub(".*timeseries_t([0-9]+)\\.txt", "\\1", fname)
  ts <- read.table(fname, header = FALSE, col.names = c("t", "degree"))
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
digits <- c(0,2,2,2,2,2,2,2,2,2,2)
latex_table <- xtable(df_s, type = "latex", row.names = TRUE, digits = digits)
print(latex_table, file = "table_s.tex", include.rownames = TRUE)

cat("\n=== AICs ===\n")
print(df_AIC)
latex_table <- xtable(df_AIC, type = "latex", row.names = TRUE, digits = digits)
print(latex_table, file = "table_AIC.tex", include.rownames = TRUE)

cat("\n=== Δ AIC ===\n")
print(df_delta)
latex_table <- xtable(df_delta, type = "latex", row.names = TRUE, digits = digits)
print(latex_table, file = "table_delta2.tex", include.rownames = TRUE)

cat("\n=== Parameter estimates ===\n")
print(df_params)
digits <- c(0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
latex_table <- xtable(df_params, type = "latex", row.names = TRUE, digits = digits)
print(latex_table, file = "table_params2.tex", include.rownames = TRUE)

#Plots
kp_all <- do.call(rbind, kp_list)
t_max <- max(kp_all$t)
t_theo <- seq(10, t_max, length.out = 1000)

theory_df <- data.frame(
  t = t_theo,
  k1 = msub0 * sqrt(t_theo),
  k2 = msub0*log(msub0 + t_theo - 1)
)
################ IF BARABÁSI-ALBERT ################ 
ggplot() +
  geom_line(data = subset(kp_all, t >= 10), aes(x = t, y = k1, color = vertex), size = 1) +
  geom_line(data = subset(theory_df, t >= 10), aes(x = t, y = k1), color = "black", size = 1.2, linetype = "dashed") +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_log10(labels = scales::label_number()) +
  theme_bw() +
  labs(
    x = "t",
    y = "k'_i(t)"
  )

################ IF RANDOM ATTACHMENT ################
ggplot() +
  geom_line(data = subset(kp_all, t >= 10), aes(x = t, y = k2, color = vertex), size = 1) +
  geom_line(data = subset(theory_df, t >= 10), aes(x = t, y = k2), color = "black", size = 1.2, linetype = "dashed") +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_log10(labels = scales::label_number()) +
  theme_bw() +
  labs(
    x = "t",
    y = "k''_i(t)"
  )

################ IF NO GROWTH ################
t_theo <- seq(1, t_max, length.out = 1000)
theory_df$t <- t_theo
theory_df$k3 <- (2*msub0/nsub0)*t_theo

ggplot() +
  geom_line(data = kp_all, aes(x = t, y = k3, color = vertex), size = 1) +
  geom_line(data = theory_df, aes(x = t, y = k3), color = "black", size = 1.2, linetype = "dashed") +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  theme_bw() +
  labs(
    x = "t",
    y = "k_i(t)"
  )