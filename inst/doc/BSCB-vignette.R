## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BSCB)

## -----------------------------------------------------------------------------
# Simulate data from a quadratic model using a D-optimal covariate design
theta_true <- c(-6, -3, 0.25)
alpha <- 0.05
a <- -0.5
b <- 0.5
L <- 5000 
p <- 2
n <- 20
e_sd <- 0.2

sim_data <- generate_simulation_data(
   p           = p,
   n           = n,
   e_sd        = e_sd,
   theta_true  = theta_true,
   a           = a,
   b           = b,
   replication = 1,
   design_index = 2,
   center_index = 1
)

X <- sim_data$X
x <- sim_data$X[,2]
Y <- as.numeric(sim_data$Y.list[[1]])

## -----------------------------------------------------------------------------
# You could also generate data in the simple way
# set.seed(123)
# n <- 20
# x <- seq(-0.5, 0.5, length.out = n)
# X <- cbind(1, x, x^2)
# theta_true <- c(-6, -3, 0.25)
# Y <- as.numeric(X %*% theta_true + rnorm(n, sd = 0.2))

## -----------------------------------------------------------------------------
# --- BSCB-C: Bayesian simultaneous credible bands under the Normal-Gamma conjugate prior ---
fit_c <- compute_bscb_conjugate(
  X              = X,
  Y              = Y,
  alpha          = alpha,
  a              = a,
  b              = b,
  L              = L,
  theta_true     = theta_true,
  hyperparameter = "g_prior",   # "empirical", "unit_info", or "g_prior"
  optimize_type  = "P"          # "P" = polyroot (recommended)
)

cat("Critical constant (BSCB-C):", fit_c$lambda, "\n")
cat("Posterior mean of theta:\n")
print(round(fit_c$mu_star, 4))

## -----------------------------------------------------------------------------
# --- BSCB-I-J: Bayesian simultaneous credible bands under the Independent Jeffreys prior ---
fit_j <- compute_bscb_ind_jeffreys(
  X     = X,
  Y     = Y,
  alpha = alpha,
  a     = a,
  b     = b,
  theta_true = theta_true,
  L     = L
)
cat("Critical constant (BSCB-J):", fit_j$lambda, "\n")

## -----------------------------------------------------------------------------
# mod <- instantiate::stan_package_model(
#   name    = "HMC_model",
#   package = "BSCB",
#   compile = TRUE
# )

## -----------------------------------------------------------------------------
# --- BSCB-H-C: BSCB under the normal-half-Cauchy prior(0,2) implemented via HMC ---

# fit_h <- compute_bscb_hmc(
#   X     = X,
#   Y     = Y,
#   V     = diag(n),
#   alpha = alpha,
#   a     = a,
#   b     = b,
#   theta_true = theta_true,
#   prior_type = "normal_half_cauchy",
#   L     = L,
#   draw_num = 10000
# )
# 
# cat("Critical constant (BSCB-H-C):", fit_h$lambda, "\n")

## -----------------------------------------------------------------------------
# --- BPCB-I-J: Bayesian pointwise credible bands under the Independent Jeffreys prior ---
fit_p <- compute_bpcb_ind_jeffreys(
  X     = X,
  Y     = Y,
  alpha = alpha,
  a     = a,
  b     = b,
  theta_true = theta_true
)


## ----fig.width=8, fig.height=6, out.width="100%"------------------------------
library(ggplot2)
x_seq  <- seq(-0.5, 0.5, length.out = 500)
y_true <- as.numeric(cbind(1, x_seq, x_seq^2) %*% theta_true)
df_obs <- data.frame(x = x, Y = as.numeric(Y))

# Collect all band boundaries into a single data frame
df_bands <- data.frame(
  x       = rep(x_seq, 3),
  lower   = c(as.numeric(fit_c$lower_bound(x_seq)),
              as.numeric(fit_j$lower_bound(x_seq)),
              as.numeric(fit_p$lower_bound(x_seq))),
  upper   = c(as.numeric(fit_c$upper_bound(x_seq)),
              as.numeric(fit_j$upper_bound(x_seq)),
              as.numeric(fit_p$upper_bound(x_seq))),
  method  = rep(c("BSCB-C-G",
                  "BSCB-I-J",
                  "BPCB-I-J"),
                each = length(x_seq))
)

df_true <- data.frame(x = x_seq, y = y_true)

band_colours <- c(
  "BSCB-C-G" = "#4DAF4A",
  "BSCB-I-J" = "#E41A1C",
  "BPCB-I-J" = "#377EB8"
)

band_linetypes <- c(
  "BSCB-C-G" = "F1",
  "BSCB-I-J" = "dotdash",
  "BPCB-I-J" = "solid"
)

ggplot() +
  # Shaded credible regions
  geom_ribbon(
    data    = df_bands,
    mapping = aes(x = x, ymin = lower, ymax = upper,
                  fill = method),
    alpha   = 0.10
  ) +
  # Band boundaries
  geom_line(
    data    = df_bands,
    mapping = aes(x = x, y = lower,
                  colour   = method,
                  linetype = method),
    linewidth = 0.8
  ) +
  geom_line(
    data    = df_bands,
    mapping = aes(x = x, y = upper,
                  colour   = method,
                  linetype = method),
    linewidth = 0.8
  ) +
  # True regression curve
  geom_line(
    data      = df_true,
    mapping   = aes(x = x, y = y),
    colour    = "navyblue",
    linewidth = 0.7,
    linetype  = "solid"
  ) +
  # Observed data
  geom_point(
    data    = df_obs,
    mapping = aes(x = x, y = Y),
    colour  = "gray50",
    size    = 1.5
  ) +
  scale_colour_manual(
  values = band_colours,
  name   = "Method",
  breaks = c("BSCB-C-G", "BSCB-I-J", "BPCB-I-J")) +
  scale_fill_manual(
  values = band_colours,
  name   = "Method",
  breaks = c("BSCB-C-G", "BSCB-I-J", "BPCB-I-J")) +
  scale_linetype_manual(
  values = band_linetypes,
  name   = "Method",
  breaks = c("BSCB-C-G", "BSCB-I-J", "BPCB-I-J")) +
  labs(
    title = "95% BSCB-C-G, BSCB-I-J and BPCB-I-J for Quadratic Regression",
    x     = "x",
    y     = "y"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title      = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------

escr_j <- coverage_ESCR(fit_j, optimize_type = "P", verbose = TRUE)
escr_c <- coverage_ESCR(fit_c, optimize_type = "P", verbose = TRUE)
escr_p <- coverage_ESCR(fit_p, optimize_type = "P", verbose = TRUE)

cat("ESCR (BSCB-I-J):", escr_j, "\n")
cat("ESCR (BSCB-C-G):", escr_c, "\n")
cat("ESCR (BPCB-I-J):", escr_p, "\n")

## -----------------------------------------------------------------------------

pscp_j <- coverage_PSCP(fit_j, draw_num = 10000,
                         optimize_type = "P", verbose = TRUE)
pscp_c <- coverage_PSCP(fit_c, draw_num = 10000,
                         optimize_type = "P", verbose = TRUE)
pscp_p <- coverage_PSCP(fit_p, draw_num = 10000,
                         optimize_type = "P", verbose = TRUE)

cat("PSCP (BSCB-I-J):", round(pscp_j, 4), "\n")
cat("PSCP (BSCB-C-G):", round(pscp_c, 4), "\n")
cat("PSCP (BPCB-I-J):", round(pscp_p, 4), "\n")

## -----------------------------------------------------------------------------
summary_tab <- data.frame(
  Method   = c("BSCB-C-G", "BSCB-I-J"),
  Lambda   = round(c(fit_c$lambda, fit_j$lambda), 4),
  ESCR     = c(escr_c, escr_j),
  PSCP     = round(c(pscp_c, pscp_j), 4)
)
knitr::kable(summary_tab, caption = "Coverage summary for one simulated dataset")

