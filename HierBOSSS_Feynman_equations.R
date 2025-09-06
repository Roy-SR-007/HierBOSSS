# =============================================================================
# Copyright 2025. Somjit Roy and Pritam Dey.
# This program implements HierBOSSS application to Feynman equations as developed in:
#   Roy, S., Dey, P., Pati, D., and Mallick, B.K.
#   'Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions'.
# Authors:
#   Somjit Roy (sroy_123@tamu.edu) and Pritam Dey (pritam.dey@tamu.edu)
# =============================================================================

## F = k * q1q2/r^2: Coulomb's law
## F: electric force, k: Coulumb constant, q1,q2: charges, r: distance of separation
## k = 1/(4*pi*epsilon)

source("MCMC.R")

I_12_2_data = readRDS("Data/AI_Feynman/I_12_2_data.RDS")

K1 = 2 # no. of symbolic trees
n = 1000
p = I_12_2_data$p
niter = 1000 # no. of posterior sampling iterations

sigma_sq = 0.25 # model noise
X = I_12_2_data$X[1:n, ] # x1: q1^2/r^2; x2: q2/q1
y = I_12_2_data$Y[1:n] + rnorm(n, 0, sqrt(sigma_sq)) 

Ops = c("*", "+") # set of operators
Op_type = c(2, 2) # type of operators, unary:1 and binary:2
Op_weights = rep(1/length(Ops), length(Ops)) # w_op initial values
Ft_weights = rep(1/p, p) # w_ft initial values
delta_0 = 1.2 # split-probability parameter
max_depth = 3 # maximum depth of each symbolic tree

## hyper-parameter initial specification for Dirichlet prior over weights
conc_op = 1
conc_ft = 1
alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1, nrow = p)
for(j in 1:K1) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}

## hyper-parameter initial specification for NIG prior over model parameters
mu_beta = rep(1.0, K1)
Sigma_beta = diag(1, K1)
nu = 1
lambda = 1

## choice of proposal weights and operator set
w_op_prop = rep(1/length(Ops), length(Ops))
w_ft_prop = rep(1/p, p)
Ops_prop = c("*", "+")
Ops_type_prop = c(2, 2)

nfeature = p

MCMCenv = listenv::listenv()

results = HierBOSSS(y, X, K1, nfeature, Ops, Op_type, beta_0 = delta_0, max_depth, alpha_op, alpha_ft,
                    mu_beta, Sigma_beta, nu, lambda,
                    w_op_prop, w_ft_prop, 
                    Ops_prop, Ops_type_prop, p_grow = 0.5,
                    niter, sigma_sq = sigma_sq, MCMCenv)

expressions_sim1 = matrix("", niter, K1)
for(i in 1:niter){
  for(j in 1:K1)
  {
    expressions_sim1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

MSE = rep(0, niter)

for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X)
  }
  
  MSE[i] = mean((y - tree_val_y) ^ 2)
}

# indices = order(MSE, decreasing = F)
# MSE_trees = cbind(Tree1.expr = expressions_sim1[indices, 1],
#                   Tree2.expr = expressions_sim1[indices, 2],
#                   RMSE = sqrt(MSE[indices]), 
#                   beta1.est = MCMCenv$beta[indices, 1],
#                   beta2.est = MCMCenv$beta[indices, 2])
# View(MSE_trees)

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

# indices_marg_lik = order(MCMCenv$marg_lik, decreasing = T)
# indices_marg_lik = indices_marg_lik[-1]
# marg_lik_trees = cbind(Tree1.expr = expressions_sim1[indices_marg_lik, 1],
#                        Tree2.expr = expressions_sim1[indices_marg_lik, 2],
#                        marg_lik = MCMCenv$marg_lik[indices_marg_lik], 
#                        #Roshoman_ratio = exp(MCMCenv$marg_lik[indices_marg_lik]-MCMCenv$marg_lik[indices_marg_lik[1]]),
#                        beta1.est = MCMCenv$beta[indices_marg_lik, 1],
#                        beta2.est = MCMCenv$beta[indices_marg_lik, 2])
# View(marg_lik_trees)

## F = q(E_f + Bvsin(theta)): Lorentz force
## F: electromagentic force, E_f: electric field, q:charge, B: magnetic field, v: velocity, theta: angle

source("MCMC.R")

I_12_11_data = readRDS("Data/AI_Feynman/I_12_11_data.RDS")

K1 = 2 # no. of symbolic trees
n = 1000
p = I_12_11_data$p
niter = 1000 # no. of posterior sampling iterations

sigma_sq = 0.25 # model noise
X = I_12_11_data$X[1:n, ] # x1: Ef*q; x2: Bv/Ef; x3: theta
y = I_12_11_data$Y[1:n] + rnorm(n, 0, sqrt(sigma_sq))

Ops = c("*", "sin") # set of operators
Op_type = c(2, 1) # type of operators, unary:1 and binary:2
Op_weights = c(0.8, 0.2)#rep(1/length(Ops), length(Ops)) # w_op initial values
Ft_weights = rep(1/p, p) # w_ft initial values
delta_0 = 1.2 # split-probability parameter
max_depth = 3 # maximum depth of each symbolic tree

## hyper-parameter initial specification for Dirichlet prior over weights
conc_op = 1
conc_ft = 1
alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1, nrow = p)
for(j in 1:K1) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}

## hyper-parameter initial specification for NIG prior over model parameters
mu_beta = rep(1.0, K1)
Sigma_beta = diag(1, K1)
nu = 1
lambda = 1

## choice of proposal weights and operator set
w_op_prop = c(0.8, 0.2)# rep(1/length(Ops), length(Ops))
w_ft_prop = rep(1/p, p)
Ops_prop = c("*", "sin")
Ops_type_prop = c(2, 1)

nfeature = p

MCMCenv = listenv::listenv()

results = HierBOSSS(y, X, K1, nfeature, Ops, Op_type, beta_0 = delta_0, max_depth, alpha_op, alpha_ft,
                    mu_beta, Sigma_beta, nu, lambda,
                    w_op_prop, w_ft_prop, 
                    Ops_prop, Ops_type_prop, p_grow = 0.5,
                    niter, sigma_sq = sigma_sq, MCMCenv)

expressions_sim1 = matrix("", niter, K1)
for(i in 1:niter){
  for(j in 1:K1)
  {
    expressions_sim1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

MSE = rep(0, niter)

for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X)
  }
  
  MSE[i] = mean((y - tree_val_y) ^ 2)
}

# indices = order(MSE, decreasing = F)
# MSE_trees = cbind(Tree1.expr = expressions_sim1[indices, 1],
#                   Tree2.expr = expressions_sim1[indices, 2],
#                   RMSE = sqrt(MSE[indices]), 
#                   beta1.est = MCMCenv$beta[indices, 1],
#                   beta2.est = MCMCenv$beta[indices, 2])
# View(MSE_trees)

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

## U = m1m2(1/r2 - 1/r1): Change in Gravitational Potential Energy
## F: gravitational potential energy, m1, m2: point masses, 
## r1, r2: distance of separation

source("MCMC.R")

I_13_12_data = readRDS("Data/AI_Feynman/I_13_12_data.RDS")

K1 = 5 # no. of symbolic trees
n = 1000
p = I_13_12_data$p # p
niter = 1000 # no. of posterior sampling iterations

sigma_sq = 0.25 # model noise
X = as.matrix(I_13_12_data$X)[1:n, ] # x1: m1^2/r1; x2: m2/m1; x3: r1/r2
y = as.numeric(I_13_12_data$Y)[1:n] # response vector

Ops = c("*", "inv", "+") # set of operators
Op_type = c(2, 1, 2) # type of operators, unary:1 and binary:2
Op_weights = rep(1/length(Ops), length(Ops)) # w_op initial values
Ft_weights = rep(1/p, p) # w_ft initial values
delta_0 = 1.2 # split-probability parameter
max_depth = 2 # maximum depth of each symbolic tree

## hyper-parameter initial specification for Dirichlet prior over weights
conc_op = 1
conc_ft = 1
alpha_op = matrix(0, ncol = K1, nrow = length(Ops))
alpha_ft = matrix(0, ncol = K1, nrow = p)
for(j in 1:K1) {
  alpha_op[, j] = Op_weights * conc_op
  alpha_ft[, j] = Ft_weights * conc_ft
}

## hyper-parameter initial specification for NIG prior over model parameters
mu_beta = rep(1.0, K1)
Sigma_beta = diag(1, K1)
nu = 1
lambda = 1

## choice of proposal weights and operator set
w_op_prop = rep(1/length(Ops), length(Ops))
w_ft_prop = rep(1/p, p)
Ops_prop = c("*", "inv", "+")
Ops_type_prop = c(2, 1, 2)

nfeature = p

MCMCenv = listenv::listenv()

results = HierBOSSS(y, X, K1, nfeature, Ops, Op_type, beta_0 = delta_0, max_depth, alpha_op, alpha_ft,
                    mu_beta, Sigma_beta, nu, lambda,
                    w_op_prop, w_ft_prop, 
                    Ops_prop, Ops_type_prop, p_grow = 0.5,
                    niter, sigma_sq = sigma_sq, MCMCenv)

expressions_sim1 = matrix("", niter, K1)
for(i in 1:niter){
  for(j in 1:K1)
  {
    expressions_sim1[i,j] = Express(MCMCenv$forests[[i]][[j]])
  }
}

MSE = rep(0, niter)

for(i in 1:niter) {
  tree_val_y = rep(0, n)
  for(j in 1:K1) {
    tree_val_y = tree_val_y + MCMCenv$beta[i, j] * allcal(MCMCenv$forests[[i]][[j]], X)
  }
  
  MSE[i] = mean((y - tree_val_y) ^ 2)
}

# indices = order(MSE, decreasing = F)
# MSE_trees = cbind(Tree1.expr = expressions_sim1[indices, 1],
#                   Tree2.expr = expressions_sim1[indices, 2],
#                   RMSE = sqrt(MSE[indices]), 
#                   beta1.est = MCMCenv$beta[indices, 1],
#                   beta2.est = MCMCenv$beta[indices, 2])
# View(MSE_trees)

indices_JMP = order(MCMCenv$JMP_trees[1:niter], decreasing = T)
indices_JMP = indices_JMP[-1]
JMP_trees = cbind(Tree1.expr = expressions_sim1[indices_JMP, 1],
                  Tree2.expr = expressions_sim1[indices_JMP, 2],
                  Tree3.expr = expressions_sim1[indices_JMP, 3],
                  Tree4.expr = expressions_sim1[indices_JMP, 4],
                  Tree5.expr = expressions_sim1[indices_JMP, 5],
                  JMP = MCMCenv$JMP_trees[indices_JMP],
                  RMSE = sqrt(MSE[indices_JMP]),
                  #Roshoman_ratio = MCMCenv$JMP_trees[indices_JMP]-MCMCenv$JMP_trees[indices_JMP[1]],
                  beta1.est = MCMCenv$beta[indices_JMP, 1],
                  beta2.est = MCMCenv$beta[indices_JMP, 2],
                  beta3.est = MCMCenv$beta[indices_JMP, 3],
                  beta4.est = MCMCenv$beta[indices_JMP, 4],
                  beta5.est = MCMCenv$beta[indices_JMP, 5])
View(JMP_trees)