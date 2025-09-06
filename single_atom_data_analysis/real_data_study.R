# =============================================================================
# Copyright 2025. Somjit Roy and Pritam Dey.
# This program implements HierBOSSS' application to single-atom catalysis dataset 
# as developed in:
#   Roy, S., Dey, P., Pati, D., and Mallick, B.K.
#   'Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions'.
# Authors:
#   Somjit Roy (sroy_123@tamu.edu) and Pritam Dey (pritam.dey@tamu.edu)
# =============================================================================

### Single-Atom catalysis data

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y

## First stage analysis

### K = 2

source("MCMC.R")

X_stage1 = cbind(X$Hfo, X$Oxv)
y_stage1 = y

K1 = 2
n = length(y_stage1)
niter = 10000

nfeature = ncol(X_stage1)

Ops = c("*", "+", "inv", "neg")
Op_type = c(2, 2, 1, 1)
Op_weights = rep(1/length(Ops), length(Ops))
Ft_weights = rep(1/nfeature, nfeature)
beta_0 = 1.2
max_depth = 5

conc_op = 1
conc_ft = 1

alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1, nrow = nfeature)

for(j in 1:K1) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}


mu_beta = rep(0, K1)
Sigma_beta = diag(1, K1)
nu = 1
lambda = 1

w_op_prop = Op_weights
w_ft_prop = Ft_weights
Ops_prop = Ops
Ops_type_prop = Op_type

MCMCenv = listenv::listenv()

results_catlysis = HierBOSSS(y_stage1, X_stage1, K1, nfeature, Ops, Op_type, beta_0 = beta_0, max_depth, alpha_op, alpha_ft,
                      mu_beta, Sigma_beta, nu, lambda,
                      w_op_prop, w_ft_prop, 
                      Ops_prop, Ops_type_prop, p_grow = 0.5,
                      niter, sigma_sq = 0.25, MCMCenv)

expressions_catalysis_1 = matrix("", niter, K1)
for(i in 1:niter){
  for(j in 1:K1)
  {
    expressions_catalysis_1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

pred_error = rep(0, niter)

for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X_stage1)
  }
  
  pred_error[i] = mean((y - tree_val_y) ^ 2)
}

indices = order(pred_error, decreasing = F)
RMSE_tab = data.frame(Tree1.expr = expressions_catalysis_1[indices, 1],
           Tree2.expr = expressions_catalysis_1[indices, 2],
           ##Tree3.expr = expressions_catalysis_1[indices, 3], 
           PredError = pred_error[indices],
           #JMP = MCMCenv$JMP_trees[indices],
           beta1.est = MCMCenv$beta[indices, 1],
           beta2.est = MCMCenv$beta[indices, 2])
           ##beta3.est = MCMCenv$beta[indices, 3]))
#uJMP = sort(unique(RMSE_tab$JMP), decreasing = T)[-1]
#r = which(uJMP == RMSE_tab$JMP[which.min(RMSE_tab$PredError)][1])


indices_JMP = order(MCMCenv$JMP_trees, decreasing = T)
View(cbind(Tree1.expr = expressions_catalysis_1[indices_JMP, 1],
           Tree2.expr = expressions_catalysis_1[indices_JMP, 2],
           ##Tree3.expr = expressions_catalysis_1[indices_JMP, 3], 
           JMP = MCMCenv$JMP_trees[indices_JMP],
           RMSE = pred_error[indices_JMP],
           beta1.est = MCMCenv$beta[indices_JMP, 1],
           beta2.est = MCMCenv$beta[indices_JMP, 2]))
           #beta3.est = MCMCenv$beta[indices_JMP, 3]))


### K = 3

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y

## First stage analysis

source("MCMC.R")

X_stage1 = cbind(X$Hfo, X$Oxv)
y_stage1 = y

K1 = 3
n = length(y_stage1)
niter = 20000

nfeature = ncol(X_stage1)

Ops = c("*", "+", "inv", "neg", "square")
Op_type = c(2, 2, 1, 1, 1)
Op_weights = rep(1/length(Ops), length(Ops))
Ft_weights = rep(1/nfeature, nfeature)
beta_0 = 1.2
max_depth = 5

conc_op = 1
conc_ft = 1

alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1, nrow = nfeature)

for(j in 1:K1) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}


mu_beta = rep(0, K1)
Sigma_beta = diag(1, K1)
nu = 1
lambda = 1

w_op_prop = Op_weights
w_ft_prop = Ft_weights
Ops_prop = Ops
Ops_type_prop = Op_type

MCMCenv = listenv::listenv()

results_catlysis = HierBOSSS(y_stage1, X_stage1, K1, nfeature, Ops, Op_type, beta_0 = beta_0, max_depth, alpha_op, alpha_ft,
                             mu_beta, Sigma_beta, nu, lambda,
                             w_op_prop, w_ft_prop, 
                             Ops_prop, Ops_type_prop, p_grow = 0.5,
                             niter, sigma_sq = 0.25, MCMCenv)
niter = 2000
expressions_catalysis_1 = matrix("", niter, K1)
for(i in 1:niter){
  for(j in 1:K1)
  {
    expressions_catalysis_1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

pred_error = rep(0, niter)

for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X_stage1)
  }
  
  pred_error[i] = mean((y - tree_val_y) ^ 2)
}

indices = order(pred_error, decreasing = F)
View(cbind(Tree1.expr = expressions_catalysis_1[indices, 1],
           Tree2.expr = expressions_catalysis_1[indices, 2],
           Tree3.expr = expressions_catalysis_1[indices, 3], 
           PredError = pred_error[indices], 
           beta1.est = MCMCenv$beta[indices, 1],
           beta2.est = MCMCenv$beta[indices, 2],
           beta3.est = MCMCenv$beta[indices, 3]))

indices_JMP = order(MCMCenv$JMP_trees, decreasing = T)
View(cbind(Tree1.expr = expressions_catalysis_1[indices_JMP, 1],
           Tree2.expr = expressions_catalysis_1[indices_JMP, 2],
           Tree3.expr = expressions_catalysis_1[indices_JMP, 3], 
           JMP = MCMCenv$JMP_trees[indices_JMP],
           RMSE = pred_error[indices_JMP],
           beta1.est = MCMCenv$beta[indices_JMP, 1],
           beta2.est = MCMCenv$beta[indices_JMP, 2],
           beta3.est = MCMCenv$beta[indices_JMP, 3]))

### K = 4

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y

## First stage analysis

source("MCMC.R")

X_stage1 = cbind(X$Hfo, X$Oxv)
y_stage1 = y

K1 = 4
n = length(y_stage1)
niter = 2000

nfeature = ncol(X_stage1)

Ops = c("*", "+", "inv", "neg", "square")
Op_type = c(2, 2, 1, 1, 1)
Op_weights = rep(1/length(Ops), length(Ops))
Ft_weights = rep(1/nfeature, nfeature)
beta_0 = 1.2
max_depth = 5

conc_op = 1
conc_ft = 1

alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1, nrow = nfeature)

for(j in 1:K1) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}


mu_beta = rep(0, K1)
Sigma_beta = diag(1, K1)
nu = 1
lambda = 1

w_op_prop = Op_weights
w_ft_prop = Ft_weights
Ops_prop = Ops
Ops_type_prop = Op_type

MCMCenv = listenv::listenv()

results_catlysis = HierBOSSS(y_stage1, X_stage1, K1, nfeature, Ops, Op_type, beta_0 = beta_0, max_depth, alpha_op, alpha_ft,
                             mu_beta, Sigma_beta, nu, lambda,
                             w_op_prop, w_ft_prop, 
                             Ops_prop, Ops_type_prop, p_grow = 0.5,
                             niter, sigma_sq = 0.25, MCMCenv)

expressions_catalysis_1 = matrix("", niter, K1)
for(i in 1:niter){
  for(j in 1:K1)
  {
    expressions_catalysis_1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

pred_error = rep(0, niter)

for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X_stage1)
  }
  
  pred_error[i] = mean((y - tree_val_y) ^ 2)
}

indices = order(pred_error, decreasing = F)
View(cbind(Tree1.expr = expressions_catalysis_1[indices, 1],
           Tree2.expr = expressions_catalysis_1[indices, 2],
           Tree3.expr = expressions_catalysis_1[indices, 3],
           Tree4.expr = expressions_catalysis_1[indices, 4],
           PredError = pred_error[indices], 
           beta1.est = MCMCenv$beta[indices, 1],
           beta2.est = MCMCenv$beta[indices, 2],
           beta3.est = MCMCenv$beta[indices, 3],
           beta4.est = MCMCenv$beta[indices, 4]))

indices_JMP = order(MCMCenv$JMP_trees, decreasing = T)
View(cbind(Tree1.expr = expressions_catalysis_1[indices_JMP, 1],
           Tree2.expr = expressions_catalysis_1[indices_JMP, 2],
           Tree3.expr = expressions_catalysis_1[indices_JMP, 3], 
           JMP = MCMCenv$JMP_trees[indices_JMP],
           RMSE = pred_error[indices_JMP],
           beta1.est = MCMCenv$beta[indices_JMP, 1],
           beta2.est = MCMCenv$beta[indices_JMP, 2],
           beta3.est = MCMCenv$beta[indices_JMP, 3]))

### K = 5

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y

## First stage analysis

source("MCMC.R")

X_stage1 = cbind(X$Hfo, X$Oxv)
y_stage1 = y

K1 = 5
n = length(y_stage1)
niter = 2000

nfeature = ncol(X_stage1)

Ops = c("*", "+", "inv", "neg", "square")
Op_type = c(2, 2, 1, 1, 1)
Op_weights = rep(1/length(Ops), length(Ops))
Ft_weights = rep(1/nfeature, nfeature)
beta_0 = 1.2
max_depth = 5

conc_op = 1
conc_ft = 1

alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1, nrow = nfeature)

for(j in 1:K1) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}


mu_beta = rep(0, K1)
Sigma_beta = diag(1, K1)
nu = 1
lambda = 1

w_op_prop = Op_weights
w_ft_prop = Ft_weights
Ops_prop = Ops
Ops_type_prop = Op_type

MCMCenv = listenv::listenv()

results_catlysis = HierBOSSS(y_stage1, X_stage1, K1, nfeature, Ops, Op_type, beta_0 = beta_0, max_depth, alpha_op, alpha_ft,
                             mu_beta, Sigma_beta, nu, lambda,
                             w_op_prop, w_ft_prop, 
                             Ops_prop, Ops_type_prop, p_grow = 0.5,
                             niter, sigma_sq = 0.25, MCMCenv)

expressions_catalysis_1 = matrix("", niter, K1)
for(i in 1:niter){
  for(j in 1:K1)
  {
    expressions_catalysis_1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

pred_error = rep(0, niter)

for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X_stage1)
  }
  
  pred_error[i] = mean((y - tree_val_y) ^ 2)
}

indices = order(pred_error, decreasing = F)
View(cbind(Tree1.expr = expressions_catalysis_1[indices, 1],
           Tree2.expr = expressions_catalysis_1[indices, 2],
           Tree3.expr = expressions_catalysis_1[indices, 3],
           Tree4.expr = expressions_catalysis_1[indices, 4],
           Tree5.expr = expressions_catalysis_1[indices, 5],
           PredError = pred_error[indices], 
           beta1.est = MCMCenv$beta[indices, 1],
           beta2.est = MCMCenv$beta[indices, 2],
           beta3.est = MCMCenv$beta[indices, 3],
           beta4.est = MCMCenv$beta[indices, 4],
           beta5.est = MCMCenv$beta[indices, 5]))

indices_JMP = order(MCMCenv$JMP_trees, decreasing = T)
View(cbind(Tree1.expr = expressions_catalysis_1[indices_JMP, 1],
           Tree2.expr = expressions_catalysis_1[indices_JMP, 2],
           Tree3.expr = expressions_catalysis_1[indices_JMP, 3], 
           JMP = MCMCenv$JMP_trees[indices_JMP],
           RMSE = pred_error[indices_JMP],
           beta1.est = MCMCenv$beta[indices_JMP, 1],
           beta2.est = MCMCenv$beta[indices_JMP, 2],
           beta3.est = MCMCenv$beta[indices_JMP, 3]))

####

# boxplot
source("MCMC.R")

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y
X = cbind(X$Hfo, X$Oxv)
y = y

nreps = 25
ooRMSE = matrix(0, ncol = 4, nrow = nreps)
train_RMSE = matrix(0, ncol = 4, nrow = nreps)
niter = 2000
nfeature = ncol(X)

Ops = c("*", "+", "inv", "neg", "square")
Op_type = c(2, 2, 1, 1, 1)
Op_weights = rep(1/length(Ops), length(Ops))
Ft_weights = rep(1/nfeature, nfeature)
beta_0 = 1.2
max_depth = 5

conc_op = 1
conc_ft = 1

nu = 1
lambda = 1

w_op_prop = Op_weights
w_ft_prop = Ft_weights
Ops_prop = Ops
Ops_type_prop = Op_type

for(rep in 1:nreps){
   for(K1 in 2:5) {
    s = sample(1:length(y), floor(0.90*length(y)), replace=F)
    s = sort(s)
    train_X = X[s, ]
    train_y = y[s]
    test_X = X[-s, ]
    test_y = y[-s]
    n = length(train_y)
    alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
    alpha_ft = matrix(0, ncol = K1, nrow = nfeature)
    for(j in 1:K1) {
      alpha_op[, j] = Op_weights * conc_op
      alpha_ft[, j] = Ft_weights * conc_ft
    }
    mu_beta = rep(0, K1)
    Sigma_beta = diag(1, K1)
    
    MCMCenv = listenv::listenv()
    
    results_catlysis = HierBOSSS(train_y, train_X, K1, nfeature, Ops, Op_type, beta_0 = beta_0, max_depth, alpha_op, alpha_ft,
                                 mu_beta, Sigma_beta, nu, lambda,
                                 w_op_prop, w_ft_prop, 
                                 Ops_prop, Ops_type_prop, p_grow = 0.5,
                                 niter, sigma_sq = 0.25, MCMCenv)
    
    pred_error = rep(0, niter)
    
    for(i in 1:niter) {
      tree_val_y = rep(0, n)
      for(j in 1:K1) {
        tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], train_X)
      }
      
      pred_error[i] = mean((train_y - tree_val_y) ^ 2)
    }
    
    index_min_RMSE = which.min(pred_error)
    train_RMSE[rep, K1-1] = pred_error[index_min_RMSE]
    
    tree_val_y_test = rep(0, length(test_y))
    for(j in 1:K1) {
      tree_val_y_test = tree_val_y_test + MCMCenv$beta[index_min_RMSE, j] * allcal(MCMCenv$forests[[index_min_RMSE]][[j]], test_X)
    }
    
    ooRMSE[rep, K1-1] = mean((test_y - tree_val_y_test) ^ 2)
    
   }
  save(ooRMSE, train_RMSE, file="single_atom_data_analysis/intermediate.RData")
}

### iBART

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y

X_stage1 = cbind(X$Hfo, X$Oxv)
colnames(X_stage1) = c("Hfo", "Oxv")
y_stage1 = y

options(java.parameters = "-Xmx10g")
library(iBART)

iBART_results = iBART(X = X_stage1, y = y_stage1,
                      head = colnames(X_stage1),  # colnames of X
                      opt = c("binary", "unary", "binary"), # binary operator first
                      out_sample = FALSE,
                      Lzero = TRUE,
                      K = 2, # maximum number of descriptors in l-zero model
                      standardize = FALSE,
                      seed = 888, hold = 2)

iBART_results_1 = iBART(X = X_stage1, y = y_stage1,
                      head = colnames(X_stage1),  # colnames of X
                      opt = c("binary", "unary", "binary"), # binary operator first
                      out_sample = FALSE,
                      Lzero = TRUE,
                      K = 3, # maximum number of descriptors in l-zero model
                      standardize = FALSE, hold = 2)

iBART_results_2 = iBART(X = X_stage1, y = y_stage1,
                        head = colnames(X_stage1),  # colnames of X
                        opt = c("binary", "unary", "binary"), # binary operator first
                        out_sample = FALSE,
                        Lzero = TRUE,
                        K = 4, # maximum number of descriptors in l-zero model
                        standardize = FALSE, hold = 2)

iBART_results_3 = iBART(X = X_stage1, y = y_stage1,
                        head = colnames(X_stage1),  # colnames of X
                        opt = c("binary", "unary", "binary"), # binary operator first
                        out_sample = FALSE,
                        Lzero = TRUE,
                        K = 5, # maximum number of descriptors in l-zero model
                        standardize = FALSE, hold = 2)

###

## This is run on the Arseven Statistics Computing Cluster

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y

X_stage1 = cbind(X$Hfo, X$Oxv)
colnames(X_stage1) = c("Hfo", "Oxv")
y_stage1 = y

options(java.parameters = "-Xmx10g")
library(iBART)

data_rep = 25
K_tree = c(2, 3, 4, 5)
train_RMSE_iBART = matrix(0, ncol = length(K_tree), nrow = data_rep)
ooRMSE_iBART = matrix(0, ncol = length(K_tree), nrow = data_rep)

for(d in 1:1) {
  message(paste0("Starting ", d, " replicated of data!"))
  for(k in K_tree) {
    message(paste0("Starting ", k, " tree!"))
    BART_results = iBART(X = X_stage1, y = y_stage1,
                         head = colnames(X_stage1),  # colnames of X
                         opt = c("binary", "unary", "binary"), # binary operator first
                         out_sample = TRUE,
                         Lzero = TRUE,
                         train_ratio = 0.9,
                         K = k, # maximum number of descriptors in l-zero model
                         standardize = FALSE, hold = 2)
    train_RMSE_iBART[d, k-1] = BART_results$iBART_in_sample_RMSE
    ooRMSE_iBART[d, k-1] = BART_results$iBART_out_sample_RMSE
    message(paste0("Ending ", k, " tree!"))
  }
  message(paste0("Data replication ", d, " done!"))
  save(ooRMSE_iBART, train_RMSE_iBART, file="single_atom_data_analysis/iBART_intermediate.RData")
}

#####

## Plot of fitted binding energy versus original binding energy, for HierBOSSS(K=4)

load("single_atom_data_analysis/catalysis.rda")

X = as.data.frame(catalysis$X)
y = catalysis$y
X = cbind(X$Hfo, X$Oxv)
y = y


s = sample(1:length(y), floor(0.90*length(y)), replace=F)
s = sort(s)
train_X = X[s, ]
train_y = y[s]
test_X = X[-s, ]
test_y = y[-s]

fitted_y = -1.81 * ((test_X[,1]/test_X[,2]) * (1/test_X[,2])) - 0.0009*(test_X[,1]*test_X[,2]*test_X[,2]) -0.43*(test_X[,1]/test_X[,2]) + 0.01*(-test_X[,2] - test_X[,1])

library(ggplot2)
library(ggpmisc)

# Create a data frame
df <- data.frame(
  Fitted = fitted_y,
  Actual = test_y
)

# Plot with ggplot
be = ggplot(df, aes(x = Actual, y = Fitted)) +
  geom_point(color="blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x, 
    parse = TRUE,
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  labs(
    x = "Original Biniding Energy",
    y = "Fitted Binding Energy",
    title = expression(bold("Fitted vs. Original Binding Energy values, HierBOSSS with K=4"))
  ) +
  theme_minimal()

ggsave("single_atom_data_analysis/fitted_original_BE.png", plot = be, width = 8, height = 6, dpi = 300)

