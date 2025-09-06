# =============================================================================
# Copyright 2025. Somjit Roy and Pritam Dey.
# This program implements runtimes of HierBOSSS and iBART for the simulated
# example as developed in:
#   Roy, S., Dey, P., Pati, D., and Mallick, B.K.
#   'Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions'.
# Authors:
#   Somjit Roy (sroy_123@tamu.edu) and Pritam Dey (pritam.dey@tamu.edu)
# =============================================================================

## storing runtimes for HierBOSSS when applied to the simulated example
## with K=2 symbolic trees over 25 repetitions

source("MCMC.R")

data_gen = function(n=1000, p=3, seed = 10, sigma_sq=1.5) {
  set.seed(seed)
  X = matrix(c(rnorm(n, 4, 1),
               rnorm(n, 6, 1),
               rnorm(n, 8, 1)), ncol = p, nrow = n, byrow = F)
  y = 5 * (X[, 1] + X[, 2]) * X[, 3] + rnorm(n, 0, sqrt(sigma_sq))
  data.list = list(y = y, X = X)
  return(data.list)
}

runtime_HierBOSSS = rep(0, 25)   # to record elapsed time

for(d in 1:25) {
  sigma_sq = 1.5
  data = data_gen(1000, 3, seed=d, sigma_sq)
  X = data$X
  y = data$y
  K1 = 2 # no. of symbolic trees
  n = nrow(X)
  p = ncol(X)
  niter = 1000
  Ops = c("+", "*")
  Op_type = c(2, 2)
  Op_weights = rep(1/length(Ops), length(Ops))
  Ft_weights = rep(1/p, p)
  delta_0 = 1.2
  max_depth = 3
  
  conc_op = 1
  conc_ft = 1
  nu = 1
  lambda = 1
  
  w_op_prop = rep(1/length(Ops), length(Ops))
  w_ft_prop = rep(1/p, p)
  Ops_prop = c("+", "*")
  Ops_type_prop = c(2, 2)
  
  nfeature = p
  
  alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
  alpha_ft = matrix(0, ncol = K1, nrow = p)
  for(j in 1:K1) {
    alpha_op[, j] = Op_weights * conc_op
    alpha_ft[, j] = Ft_weights * conc_ft
  }
  
  mu_beta = rep(1.0, K1)
  Sigma_beta = diag(1, K1)
  
  MCMCenv = listenv::listenv()
  
  ## measure run time for each repetition
  t0 = proc.time()
  
  results = HierBOSSS(y, X, K1, nfeature, Ops, Op_type,
                      beta_0 = delta_0, max_depth,
                      alpha_op, alpha_ft,
                      mu_beta, Sigma_beta, nu, lambda,
                      w_op_prop, w_ft_prop,
                      Ops_prop, Ops_type_prop, p_grow = 0.5,
                      niter, sigma_sq = sigma_sq, MCMCenv)
  
  t1 = proc.time()
  runtime_HierBOSSS[d] = (t1 - t0)[["elapsed"]]
}

## summarize runtime
mean_time = mean(runtime_HierBOSSS)
print(runtime_HierBOSSS)
cat("Average runtime per repetition:", mean_time, "seconds\n")

## saving the runtime for HierBOSSS
write.csv(runtime_HierBOSSS, file="runtimes/simulated_example/runtime_HierBOSSS.csv")

# -------------------------------------------------------------------------

## storing runtimes for iBART when applied to the simulated example
## with K=2 symbolic trees over 25 repetitions

# install.packages("rJava", INSTALL_opts = "--no-multriarch")
library(rJava)
# install.packages("bartMachine", INSTALL_opts = "--no-multiarch")
library(bartMachine)
library(remotes)  # if not already installed
# remotes::install_url("https://cran.r-project.org/src/contrib/Archive/bartMachineJARs/bartMachineJARs_1.1.tar.gz")
library(bartMachineJARs)
# devtools::install_version("bartMachine", version = "1.2.6")
library(bartMachine)

set.seed(10)
options(java.parameters = "-Xmx10g") # Allocate 10GB of memory for Java
library(iBART)

runtime_iBART = rep(0, 25)

for(d in 1:25) {
  n = 250
  p = 3
  data = data_gen(n, 3, seed=d, 1.5)
  y = data$y
  X = data$X
  
  colnames(X) <- paste("x.", seq(from = 1, to = p, by = 1), sep = "")
  
  t0 = proc.time()
  
  iBART_results <- iBART(X = X, y = y,
                         head = colnames(X),
                         unit = NULL,                         # no unit information for simulation
                         opt = c("binary", "unary", "binary"), # unary operator first
                         sin_cos = FALSE,                      # add sin and cos to operator set
                         apply_pos_opt_on_neg_x = FALSE,      # e.g. do not apply log() on negative x
                         Lzero = TRUE,                        # best subset selection
                         K = 2,                               # at most 4 predictors in best subset model
                         standardize = FALSE,                 # don't standardize input matrix X
                         hold = 2)
  
  t1 = proc.time()
  
  runtime_iBART[d] = (t1-t0)[["elapsed"]]
  
  print(d)
}

## saving the runtime for iBART
write.csv(runtime_iBART, file="runtimes/simulated_example/runtime_iBART.csv")

## boxplot of the runtimes of HierBOSSS and competing methods across 25 repetitions

runtime_HierBOSSS = read.csv("runtimes/simulated_example/runtime_HierBOSSS.csv")
runtime_HierBOSSS = runtime_HierBOSSS$x
runtime_iBART = read.csv("runtimes/simulated_example/runtime_iBART.csv")
runtime_iBART = runtime_iBART$x
runtime_BSR = read.csv("runtimes/simulated_example/runtime_BSR.csv")
runtime_BSR = runtime_BSR$runtime
runtime_QLattice = read.csv("runtimes/simulated_example/runtime_QLattice.csv")
runtime_QLattice = runtime_QLattice$runtime

library(ggplot2)

# build data (each vector length = 25)
rt_df <- data.frame(
  runtime = c(as.numeric(runtime_HierBOSSS),
              as.numeric(runtime_BSR),
              as.numeric(runtime_iBART),
              as.numeric(runtime_QLattice)),
  method = factor(rep(c("HierBOSSS","BSR","iBART","QLattice"),
                      each = length(runtime_HierBOSSS)),
                  levels = c("HierBOSSS","BSR","iBART","QLattice"))
)

# pastel palette (as in your image)
pal <- c(
  "BSR"       = "#66C2A5",  # teal
  "HierBOSSS" = "#FC8D62",  # orange
  "iBART"     = "#8DA0CB",  # periwinkle
  "QLattice"  = "#E78AC3"   # pink
)

main_title <- "Runtime across methods (n = 25 each)"  # change to your exact title

p <- ggplot(rt_df, aes(x = method, y = runtime, fill = method, color = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9) +
  # black median line overlay
  stat_summary(fun = median, geom = "crossbar", width = 0.6,
               color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.55, size = 1.8) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  scale_y_continuous(trans="log",labels = scales::number_format(accuracy = 1)) +
  labs(
    #title = expression(paste(
    #  "Runtime (in log-scale) for learning: ",
    #  y == 5*(x[1] + x[2])*x[3] + epsilon, ", ",
    #  sigma^2 == 1.50
    #)),
    x = "Method",
    y = "Runtime (in log-seconds)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p)

# Save figure
ggsave("runtimes/simulated_example/runtime_boxplot.png", p, dpi = 300)
