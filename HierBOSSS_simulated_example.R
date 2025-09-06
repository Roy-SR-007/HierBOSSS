# =============================================================================
# Copyright 2025. Somjit Roy and Pritam Dey.
# This program implements HierBOSSS simulations as developed in:
#   Roy, S., Dey, P., Pati, D., and Mallick, B.K.
#   'Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions'.
# Authors:
#   Somjit Roy (sroy_123@tamu.edu) and Pritam Dey (pritam.dey@tamu.edu)
# =============================================================================

# Simulation for the simulated example in Section 5.1

## Fitting K1 = (2, 3, 4, 5) trees, across 10000 MCMC iterations
## For each of the K1 trees, we consider plotting the RMSE box plot
source("MCMC.R")

set.seed(10)

K1 = c(2, 3, 4, 5) # no. of symbolic trees
n = 1000
p = 3
niter = 10000 # no. of posterior sampling iterations

sigma_sq = 1.5 # model noise
X = matrix(c(rnorm(n, 4, 1), rnorm(n, 6, 1), rnorm(n, 8, 1)), ncol = p, nrow = n, byrow = F)
y = 5 * (X[, 1] + X[ , 2])*X[ , 3] + rnorm(n, 0, sqrt(sigma_sq)) 

Ops = c("+", "*") # set of operators
Op_type = c(2, 2) # type of operators, unary:1 and binary:2
Op_weights = rep(1/length(Ops), length(Ops)) # w_op initial values
Ft_weights = rep(1/p, p) # w_ft initial values
delta_0 = 1.2 # split-probability parameter
max_depth = 3 # maximum depth of each symbolic tree

## hyper-parameter initial specification for Dirichlet prior over weights
## alpha_op and alpha_ft chosen inside the tree-loop over k
conc_op = 1
conc_ft = 1

## hyper-parameter initial specification for NIG prior over model parameters
## mu_beta and Sigma_beta chosen inside the tree-loop over k
nu = 1
lambda = 1

## choice of proposal weights and operator set
w_op_prop = rep(1/length(Ops), length(Ops))
w_ft_prop = rep(1/p, p)
Ops_prop = c("+", "*")
Ops_type_prop = c(2, 2)

nfeature = p
RMSE_mat = matrix(0, ncol = length(K1), nrow = niter)

for(k in 1:length(K1)) {
  
  alpha_op = matrix(0, ncol = K1[k], nrow = length(Ops))
  alpha_ft = matrix(0, ncol = K1[k], nrow = p)
  for(j in 1:K1[k]) {
    alpha_op[, j] = Op_weights * conc_op
    alpha_ft[, j] = Ft_weights * conc_ft
  }
  
  mu_beta = rep(1.0, K1[k])
  Sigma_beta = diag(1, K1[k])
  
  MCMCenv = listenv::listenv()
  
  ## running HierBOSSS
  results = HierBOSSS(y, X, K1[k], nfeature, Ops, Op_type, beta_0 = delta_0, max_depth, alpha_op, alpha_ft,
                      mu_beta, Sigma_beta, nu, lambda,
                      w_op_prop, w_ft_prop, 
                      Ops_prop, Ops_type_prop, p_grow = 0.5,
                      niter, sigma_sq = sigma_sq, MCMCenv)
  
  ## determining the symbolic expressions
  expressions_sim1 = matrix("", niter, K1[k])
  for(i in 1:niter){
    for(j in 1:k)
    {
      expressions_sim1[i,j] = Express(MCMCenv$forests[[i]][[j]])
    }
  }
  
  ## computing mean-squared-error (MSE)
  MSE = rep(0, niter)
  for(i in 1:niter) {
    tree_val_y = rep(0, n)
    for(j in 1:K1[k]) {
      tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X)
    }
    
    MSE[i] = mean((y - tree_val_y) ^ 2)
  }
  
  RMSE_mat[ , k] = MSE
}

RMSE_mat = sqrt(RMSE_mat)

## box plots of RMSE values over different values of K

library(reshape2)
library(ggplot2)

colnames(RMSE_mat) = c("K=2", "K=3", "K=4", "K=5")
df_long <- melt(RMSE_mat, variable.name = "X", value.name = "y")
df_long = df_long[, -1]
colnames(df_long) = c("K", "RMSE")

# Compute the minimum RMSE values for each K
min_RMSE_values <- aggregate(RMSE ~ K, data = df_long, FUN = min)
min_RMSE_values$vjust <- ifelse(min_RMSE_values$K %in% c("K=2", "K=3"), -9, -4)

# Create box plot
box_plot <- ggplot(df_long, aes(x = K, y = RMSE, fill = K)) +
  geom_boxplot(color = "black", alpha = 0.6, outlier.shape = NA) +
  geom_hline(yintercept = sqrt(1.5), linetype = "dashed", color = "red", size = 1) +
  geom_point(data = min_RMSE_values, aes(x = K, y = RMSE), color = "blue", size = 3, shape = 18) +  # Mark min y values
  geom_text(data = min_RMSE_values, aes(x = K, y = RMSE, label = round(RMSE, 5)), vjust = min_RMSE_values$vjust, color = "blue", size = 4, fontface = "bold") +  # Label min values
  theme_minimal() +
  labs(title = "RMSE across K",
       x = "Number of symbolic trees (K)",
       y = "RMSE") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA)   # White outer background
  )
ggsave("results/simulated_example/boxplot_RMSE_simulated_example.png", plot = box_plot, width = 8, height = 6, dpi = 300)

## fitting K=2 symbolic trees and noting the results
alpha_op = matrix(0, ncol = K1[1], nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1[1], nrow = p)
for(j in 1:K1[1]) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}

mu_beta = rep(1.0, K1[1])
Sigma_beta = diag(1, K1[1])

MCMCenv = listenv::listenv()

## running HierBOSSS
results = HierBOSSS(y, X, K1[1], nfeature, Ops, Op_type, beta_0 = delta_0, max_depth, alpha_op, alpha_ft,
                    mu_beta, Sigma_beta, nu, lambda,
                    w_op_prop, w_ft_prop, 
                    Ops_prop, Ops_type_prop, p_grow = 0.5,
                    niter, sigma_sq = sigma_sq, MCMCenv)

## determining the symbolic expressions
expressions_sim1 = matrix("", niter, K1[1])
for(i in 1:niter){
  for(j in 1:K1[1])
  {
    expressions_sim1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

## computing mean-squared-error (MSE)
MSE = rep(0, niter)
for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1[1]) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X)
  }
  
  MSE[i] = mean((y - tree_val_y) ^ 2)
}

indices = order(MSE, decreasing = F)
MSE_trees = cbind(Tree1.expr = expressions_sim1[indices, 1],
                  Tree2.expr = expressions_sim1[indices, 2],
                  RMSE = sqrt(MSE[indices]), 
                  beta1.est = MCMCenv$beta[indices, 1],
                  beta2.est = MCMCenv$beta[indices, 2])
View(MSE_trees)

indices_JMP = order(MCMCenv$JMP_trees[1:niter], decreasing = T)
indices_JMP = indices_JMP[-1]
JMP_trees = cbind(Tree1.expr = expressions_sim1[indices_JMP, 1],
                  Tree2.expr = expressions_sim1[indices_JMP, 2],
                  JMP = MCMCenv$JMP_trees[indices_JMP],
                  RMSE = sqrt(MSE[indices_JMP]),
                  #Roshoman_ratio = MCMCenv$JMP_trees[indices_JMP]-MCMCenv$JMP_trees[indices_JMP[1]],
                  beta1.est = MCMCenv$beta[indices_JMP, 1],
                  beta2.est = MCMCenv$beta[indices_JMP, 2])
View(JMP_trees)

indices_marg_lik = order(MCMCenv$marg_lik, decreasing = T)
indices_marg_lik = indices_marg_lik[-1]
marg_lik_trees = cbind(Tree1.expr = expressions_sim1[indices_marg_lik, 1],
                       Tree2.expr = expressions_sim1[indices_marg_lik, 2],
                       marg_lik = MCMCenv$marg_lik[indices_marg_lik],
                       RMSE = sqrt(MSE[indices_marg_lik]),
                       #Roshoman_ratio = exp(MCMCenv$marg_lik[indices_marg_lik]-MCMCenv$marg_lik[indices_marg_lik[1]]),
                       beta1.est = MCMCenv$beta[indices_marg_lik, 1],
                       beta2.est = MCMCenv$beta[indices_marg_lik, 2])
View(marg_lik_trees)

## plotting the predicted y values against the original for K=2 trees
fitted_y = 4.9964 * (X[,3] * X[,2]) + 5.0045 * (X[,1] * X[,3])

fitted = data.frame(y = y, y_hat = fitted_y)
model <- lm(y_hat ~ y, data = fitted)  # Linear model fitting
r_squared <- summary(model)$r.squared  # Extract R-squared

fitted_original_plot = ggplot(fitted, aes(x = y, y = y_hat)) +
  geom_point(color = "blue", size = 2) +  # Plot points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add line of best fit
  labs(title = "Fitted vs Original Values, HierBOSSS with K = 2",
       x = "y",
       y = bquote(hat(y))) +
  annotate("text", x = min(fitted$y), y = max(fitted$y_hat), label = paste("R² =", round(r_squared, 4)), 
           color = "black", size = 5, hjust = 0, vjust = 1) + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("results/simulated_example/fitted_original_simulated_example.png", plot = fitted_original_plot, width = 8, height = 6, dpi = 300)

# -------------------------------------------------------------------------

## Data generation.
source("MCMC.R")

data_gen = function(n=1000, p=3, seed = 10, sigma_sq=1.5) {
  set.seed(seed)
  X = matrix(c(rnorm(n, 4, 1), rnorm(n, 6, 1), rnorm(n, 8, 1)), ncol = p, nrow = n, byrow = F)
  y = 5 * (X[, 1] + X[ , 2])*X[ , 3] + rnorm(n, 0, sqrt(sigma_sq))
  data.list = list(y = y, X = X)
  return(data.list)
}

rmse_HierBOSSS = rep(0, 25)

for(d in 1:1) {
  sigma_sq = 1.5
  data = data_gen(1000, 3, seed=d, sigma_sq)
  X = data$X
  y = data$y
  K1 = 2 # no. of symbolic trees
  n = nrow(X)
  p = ncol(X)
  niter = 1000
  Ops = c("+", "*") # set of operators
  Op_type = c(2, 2) # type of operators, unary:1 and binary:2
  Op_weights = rep(1/length(Ops), length(Ops)) # w_op initial values
  Ft_weights = rep(1/p, p) # w_ft initial values
  delta_0 = 1.2 # split-probability parameter
  max_depth = 3 # maximum depth of each symbolic tree
  
  ## hyper-parameter initial specification for Dirichlet prior over weights
  ## alpha_op and alpha_ft chosen inside the tree-loop over k
  conc_op = 1
  conc_ft = 1
  
  ## hyper-parameter initial specification for NIG prior over model parameters
  ## mu_beta and Sigma_beta chosen inside the tree-loop over k
  nu = 1
  lambda = 1
  
  ## choice of proposal weights and operator set
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
  
  ## running HierBOSSS
  results = HierBOSSS(y, X, K1, nfeature, Ops, Op_type, beta_0 = delta_0, max_depth, alpha_op, alpha_ft,
                      mu_beta, Sigma_beta, nu, lambda,
                      w_op_prop, w_ft_prop, 
                      Ops_prop, Ops_type_prop, p_grow = 0.5,
                      niter, sigma_sq = sigma_sq, MCMCenv)
  
  # ## determining the symbolic expressions
  # expressions_sim1 = matrix("", niter, K1[k])
  # for(i in 1:niter){
  #   for(j in 1:k)
  #   {
  #     expressions_sim1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  #   }
  # }
  
  ## computing mean-squared-error (MSE)
  MSE = rep(0, niter)
  complexity_iter = rep(0, niter)
  for(i in 1:niter) {
    complexity = 0
    tree_val_y = rep(0, n)
    for(j in 1:K1) {
      tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X)
      complexity = complexity + getNum(MCMCenv$forests[[i]][[j]])
    }
    complexity_iter[i] = complexity  
    MSE[i] = mean((y - tree_val_y) ^ 2)
  }
  
  indices_JMP = order(MCMCenv$JMP_trees[1:niter], decreasing = T)
  indices_JMP = indices_JMP[-1]
  rmse_HierBOSSS[d] = (sqrt(MSE[indices_JMP[1]]))
}

rmse_HierBOSSS_25 = rmse_HierBOSSS

saveRDS(rmse_HierBOSSS_25, "results/simulated_example/rmse_HierBOSSS_25.RDS")

## boxplot for comparison

rmse_BSR_25 = read.csv("results/simulated_example/rmse_BSR_25.csv")
rmse_BSR_25 = as.numeric(rmse_BSR_25$rmse)
rmse_QLattice_25 = read.csv("results/simulated_example/rmse_QLattice_25.csv")
rmse_QLattice_25 = as.numeric(rmse_QLattice_25$rmse)
rmse_HierBOSSS_25 = readRDS("results/simulated_example/rmse_HierBOSSS_25.RDS")

# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Combine into a tidy data frame
rmse_df <- data.frame(
  rmse = c(rmse_HierBOSSS_25, rmse_BSR_25, rmse_QLattice_25),
  method = factor(rep(c("BSR", "HierBOSSS", "QLattice"), each = 25))
)

# Boxplot using ggplot2
p = ggplot(rmse_df, aes(x = method, y = rmse, fill = method)) +
  geom_boxplot(width = 0.6, outlier.color = "red", outlier.shape = NA) +
  geom_hline(yintercept = sqrt(1.50), color = "blue", linetype = "dashed", size = 0.5) +
  annotate("text", x = 2.5, y = sqrt(1.50) - 0.02, 
           label = expression(sigma == sqrt(1.50)), 
           parse = TRUE, color = "blue") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = expression("RMSE in learning " ~ y == 5*(x[1] + x[2])*x[3] + epsilon * ", " ~ sigma^2 == 1.50),
    x = "Method",
    y = "RMSE"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + ylim(min(rmse_df$rmse), 1.5)

ggsave("results/simulated_example/rmse_boxplot_25_HierBOSSS_BSR_QLattice.png", plot = p, width = 8, height = 6, dpi = 300)

BSR_complexity = c(10, 10, 16, 14, 14, 20, 16, 20, 12, 28)
HierBOSSS_complexity = c(6, 6, 6, 6, 6, 6, 10, 8, 6, 6)

# Create data frame
df_complexity <- data.frame(
  Replication = rep(1:10, 2),
  Complexity = c(BSR_complexity, HierBOSSS_complexity),
  Method = rep(c("BSR", "HierBOSSS"), each = 10)
)

# Plot
p1<-ggplot(df_complexity, aes(x = Replication, y = Complexity, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Method), fill = "white", size = 3, stroke = 1) +  # Circles at each point
  scale_color_manual(values = c("BSR" = "firebrick", 
                                "HierBOSSS" = "steelblue")) +
  scale_shape_manual(values = c("BSR" = 18,  # triangle
                                "HierBOSSS" = 21)) +  # circle
  scale_x_continuous(breaks = 1:10)+
  labs(
    title = "Symbolic tree complexity across first 10 regenerations",
    x = "Regenerations",
    y = "Symbolic tree complexity",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
ggsave("results/simulated_example/Expression_Complexity_HierBOSSS_BSR_QLattice.png", plot = p1, width = 8, height = 6, dpi = 300)


#### Testing iBART on the simulated example

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

sigma_sq = 1.5 # model noise
n=250
p=3
X = matrix(c(rnorm(n, 4, 1), rnorm(n, 6, 1), rnorm(n, 8, 1)), ncol = p, nrow = n, byrow = F)
y = 5 * (X[, 1] + X[ , 2])*X[ , 3] + rnorm(n, 0, sqrt(sigma_sq)) 
colnames(X) <- paste("x.", seq(from = 1, to = p, by = 1), sep = "")

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

## boxplot values for iBART

rmse_iBART = rep(0, 25)
options(java.parameters = "-Xmx10g") # Allocate 10GB of memory for Java

for(d in 1:25) {
  n = 250
  p = 3
  data = data_gen(n, 3, seed=d, 1.5)
  y = data$y
  X = data$X
  
  colnames(X) <- paste("x.", seq(from = 1, to = p, by = 1), sep = "")
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
  rmse_iBART[d] = iBART_results$iBART_in_sample_RMSE
  print(d)
}


saveRDS(rmse_iBART, "results/simulated_example/rmse_iBART_25.RDS")

rmse_iBART_25 = readRDS("results/simulated_example/rmse_iBART_25.RDS")

# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Combine into a tidy data frame
rmse_df <- data.frame(
  rmse = c(rmse_HierBOSSS_25, rmse_BSR_25, rmse_QLattice_25, rmse_iBART_25),
  method = factor(rep(c("BSR", "HierBOSSS", "QLattice", "iBART"), each = 25))
)

# Boxplot using ggplot2
p = ggplot(rmse_df, aes(x = rmse, y = method, fill = method)) +
  geom_boxplot(width = 0.6, outlier.color = "red", outlier.shape = NA) +
  geom_vline(xintercept = sqrt(1.50), color = "blue", linetype = "dashed", size = 0.5) +
  annotate("text", y = 2.5, x = sqrt(1.50) + 0.10, 
           label = expression(sigma == sqrt(1.50)), 
           parse = TRUE, color = "blue") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = expression("RMSE in learning " ~ y == 5*(x[1] + x[2])*x[3] + epsilon * ", " ~ sigma^2 == 1.50),
    x = "RMSE",
    y = "Method"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) + xlim(min(rmse_df$rmse), max(rmse_df$rmse))

ggsave("results/simulated_example/rmse_boxplot_25_HierBOSSS_BSR_QLattice_iBART_horizontal.png", plot = p, width = 8, height = 6, dpi = 300)

# ####
# 
# library(rJava)
# options(java.parameters = "-Xmx20g") # Allocate 10GB of memory for Java
# library(bartMachine)
# library(iBART)
# 
# I_12_2_data = readRDS("Data/AI_Feynman/I_12_2_data.RDS")
# n = 250
# p = I_12_2_data$p # p
# X = as.matrix(I_12_2_data$X)[1:n, ] # initial feature matrix
# y = as.numeric(I_12_2_data$Y)[1:n] + rnorm(n, 0, sqrt(0.25))
# colnames(X) = c("q1^2/r1^2", "q2/q1")
# 
# iBART_results_1 <- iBART(X = X, y = y,
#                        head = colnames(X),
#                        unit = NULL,                         # no unit information for simulation
#                        opt = c("unary", "binary", "unary"), # unary operator first
#                        sin_cos = FALSE,                     # add sin and cos to operator set
#                        apply_pos_opt_on_neg_x = FALSE,      # e.g. do not apply log() on negative x
#                        Lzero = TRUE,                        # best subset selection
#                        K = 3,                               # at most 4 predictors in best subset model
#                        standardize = FALSE,                 # don't standardize input matrix X
#                        hold = 2)
# 
# library(rJava)
# options(java.parameters = "-Xmx20g") # Allocate 10GB of memory for Java
# library(bartMachine)
# library(iBART)
# 
# I_12_11_data = readRDS("Data/AI_Feynman/I_12_11_data.RDS")
# n = 250
# p = I_12_2_data$p # p
# X = as.matrix(I_12_11_data$X)[1:n, ] # initial feature matrix
# y = as.numeric(I_12_11_data$Y)[1:n]
# colnames(X) = c("Ef*q", "Bv/Ef", "theta")
# 
# iBART_results_1 <- iBART(X = X, y = y,
#                          head = colnames(X),
#                          unit = NULL,                         # no unit information for simulation
#                          opt = c("binary", "unary", "binary"), # unary operator first
#                          sin_cos = TRUE,                     # add sin and cos to operator set
#                          apply_pos_opt_on_neg_x = FALSE,      # e.g. do not apply log() on negative x
#                          Lzero = TRUE,                        # best subset selection
#                          K = 2,                               # at most 4 predictors in best subset model
#                          standardize = FALSE,                 # don't standardize input matrix X
#                          hold = 2)
# 




