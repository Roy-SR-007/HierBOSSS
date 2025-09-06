# =============================================================================
# Copyright 2025. Somjit Roy and Pritam Dey.
# This program implements HierBOSSS as developed in:
#   Roy, S., Dey, P., Pati, D., and Mallick, B.K.
#   'Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions'.
# Authors:
#   Somjit Roy (sroy_123@tamu.edu) and Pritam Dey (pritam.dey@tamu.edu)
# =============================================================================

###
# This R script performs the Metropolis-within-collapsed Gibbs sampling algorithm
# for sampling from the HierBOSSS-induced posterior
###

# required imports of necessary R scripts/modules
source("tree_functions.R")
library(mvtnorm)

# inverting a matrix using Cholesky decomposition
solve_chol = function(A)
{
  return(chol2inv(chol(A)))
}

# drawing from Dirichlet distribution
rdirichlet = function(K, alpha) {
  Y = rep(0, K)
  
  for(i in 1:K) {
    
    Y[i] = rgamma(1, shape = alpha[i], scale = 1)
  }
  
  X = Y/sum(Y)
  return(X)
}

# updating the tree parameters
update_tree_parameters = function(alpha_op, alpha_ft, K1, forests, nfeature, Ops) {
  
  W_op = matrix(0, nrow = length(Ops), ncol = K1)
  W_ft = matrix(0, nrow = nfeature, ncol = K1)
  
  for(j in 1:K1) {
    
    tree = forests[[j]]
    lik = tree_lik(tree, nfeature, Ops)
    
    W_op[ , j] = rdirichlet(length(Ops), alpha_op[ , j] + lik$m_op)
      
    W_ft[ , j] = rdirichlet(nfeature, alpha_ft[ , j] + lik$n_ft)
  }
  
  res_lst = list(W_op = W_op, W_ft = W_ft)
  
  return(res_lst)
}

# updating the HierBOSSS model parameters
update_model_parameters = function(y, X, forests, K1, mu_beta, Sigma_beta, nu, lambda) {
  
  n = nrow(X)
  
  T.evaluated = matrix(0, nrow = n, ncol = K1)
  
  for(j in 1:length(forests)) {
    T.evaluated[ , j] = allcal(forests[[j]], X)
  }
  
  Sigma_beta.star = solve_chol(solve_chol(Sigma_beta) + crossprod(T.evaluated))
  mu_beta.star = Sigma_beta.star %*% (solve_chol(Sigma_beta) %*% mu_beta + crossprod(T.evaluated, y))
  lambda.star = lambda + n
  nu.star = crossprod(mu_beta, solve_chol(Sigma_beta) %*% mu_beta) + crossprod(y) + nu - crossprod(mu_beta.star, solve_chol(Sigma_beta.star) %*% mu_beta.star)
  
  beta = rmvnorm(1, mu_beta.star, Sigma_beta.star)[1 , ]
  sigma_square = 1 / rgamma(1, shape = lambda.star/2, rate = nu.star/2)
  
  res_lst = list(beta = beta, sigma_square = sigma_square)
  
  return(res_lst)
  
}

# computing the log normalizing constant
log_normalizing_constant = function(y, X, forests, K1, mu_beta, Sigma_beta, nu, lambda) {
  
  n = nrow(X)
  
  T.evaluated = matrix(0, nrow = n, ncol = K1)
  
  for(i in 1:length(forests)) {
    T.evaluated[ , i] = allcal(forests[[i]], X)
  }
  
  Sigma_beta.star = solve_chol(solve_chol(Sigma_beta) + crossprod(T.evaluated))
  mu_beta.star = Sigma_beta.star %*% (solve_chol(Sigma_beta) %*% mu_beta + crossprod(T.evaluated, y))
  lambda.star = lambda + n
  nu.star = crossprod(mu_beta, solve_chol(Sigma_beta) %*% mu_beta) + crossprod(y) + nu - crossprod(mu_beta.star, solve_chol(Sigma_beta.star) %*% mu_beta.star)
  
  NC = ((lambda.star / 2) * log(nu.star / 2)) - ((K1 / 2) * log(2 * pi)) - (0.5 * log(det(Sigma_beta.star))) - lgamma(lambda.star / 2)
  NC_return = -1 * NC
  
  return(NC_return)
}

# updates the symbolic tree structures using the Metropolis-Hasting algorithm
update_tree_MH = function(W_op, W_ft, w_op_prop, w_ft_prop, K1, X, y, trees, p_grow, nfeature, Ops_prop, Ops_type_prop, beta_0, max_depth,
                          mu_beta, Sigma_beta, nu, lambda) {
  forests_res = list()
  forests_temp = list()
  log_MH = rep(0, K1)
  
  naccept = 0
  
  for(j in 1:K1) {
    forests_res[[j]] = deepcopy(trees[[j]])
    forests_temp[[j]] = deepcopy(trees[[j]])
  }
  
  for(j in 1:K1) {
    step_prob = runif(1)
    if(step_prob <= p_grow) {
      grow_ind = 1
      terminal_nodes = get_all_terminal(forests_res[[j]], "")
      terminal_sel = sample(1:length(terminal_nodes), 1)
      change_address = terminal_nodes[[terminal_sel]]
      forests_temp[[j]] = deepcopy(forests_res[[j]])
      grow(access_terminal(change_address, forests_temp[[j]]), nfeature, Ops_prop, w_op_prop, w_ft_prop, Ops_type_prop, beta_0, 
                               max_depth)
    }else {
      grow_ind = 0
      nonterminal_nodes = get_all_nonterminal(forests_res[[j]], "")
      nonterminal_sel = sample(1:length(nonterminal_nodes), 1)
      change_address = nonterminal_nodes[[nonterminal_sel]]
      forests_temp[[j]] = deepcopy(forests_res[[j]])
      shrink(access_nonterminal(change_address, forests_temp[[j]]), w_ft_prop)
    }
    
    log_MH[j] = (log_normalizing_constant(y, X, forests_temp, K1, mu_beta, Sigma_beta, nu, lambda) +
      tree_log_likelihood_val(forests_temp[[j]], W_op[ , j], W_ft[ , j], beta_0, nfeature, Ops) - 
      log_normalizing_constant(y, X, forests_res, K1, mu_beta, Sigma_beta, nu, lambda) - 
      tree_log_likelihood_val(forests_res[[j]], W_op[ , j], W_ft[ , j], beta_0, nfeature, Ops) -
      tree_log_proposal_likelihood(forests_temp[[j]], forests_res[[j]], grow_ind, change_address, 
                                   w_op_prop, w_ft_prop, p_grow, nfeature, Ops_prop, beta_0) +
      tree_log_proposal_likelihood(forests_res[[j]], forests_temp[[j]], 1 - grow_ind, change_address,
                                   w_op_prop, w_ft_prop, p_grow, nfeature, Ops_prop, beta_0))[1, 1]
    
    if(log_MH[j] >= log(runif(1))) {
      forests_res[[j]] = deepcopy(forests_temp[[j]])
      naccept = naccept + 1
    }else {
      forests_temp[[j]] = deepcopy(forests_res[[j]])
    }
  }
  
  return(list(forests_res = forests_res, log_MH = log_MH, naccept = naccept))
}

# main MCMC module to sample from the HierBOSSS-induced posterior
HierBOSSS = function(y, X, K1, nfeature, Ops, Op_type, beta_0, max_depth,
                     alpha_op, alpha_ft, mu_beta, Sigma_beta, nu, lambda,
                     w_op_prop = NULL, w_ft_prop = NULL, Ops_prop, Ops_type_prop, p_grow,
                     niter, sigma_sq, env_, init_params = NULL) {
  
  message("Starting main_MCMC() for HierBOSSS!")
  message("Please Wait!")
  
  env_$log_MH = matrix(0, nrow = niter, ncol = K1)
  env_$beta = matrix(0, nrow = niter, ncol = K1)
  env_$sigma_square = rep(0, niter)
  env_$W_op = array(0, dim = c(length(Ops), K1, niter))
  env_$W_ft = array(0, dim = c(nfeature, K1, niter))
  
  env_$forests = list(list())
  
  env_$JMP_trees = rep(0, niter)
  
  env_$marg_lik = rep(0, niter)
  
  env_$beta[1 , ] = rmvnorm(1, mu_beta, sigma_sq * Sigma_beta)
  
  sigma_square_init = 1 / rgamma(1, shape = lambda/2, rate = nu/2)
  env_$sigma_square[1] = sigma_square_init
  
  for(j in 1:K1) {
    env_$W_op[ , j, 1] = rdirichlet(length(Ops), alpha_op)
    env_$W_ft[ , j, 1] = rdirichlet(nfeature, alpha_ft)
  }
  
  if(is.null(init_params)) {
    
    #beta_init = rmvnorm(1, mu_beta, Sigma_beta)
    #print(beta_init)
    env_$beta[1 , ] = rmvnorm(1, mu_beta, sigma_sq * Sigma_beta)
    
    sigma_square_init = 1 / rgamma(1, shape = lambda/2, rate = nu/2)
    env_$sigma_square[1] = sigma_square_init
    
    for(j in 1:K1) {
      env_$W_op[ , j, 1] = rdirichlet(length(Ops), alpha_op)
      env_$W_ft[ , j, 1] = rdirichlet(nfeature, alpha_ft)
      
      root_node = Node$new(depth = 0)
      options(expressions = 5000)
      grow(root_node, nfeature, Ops, env_$W_op[ , j, 1], env_$W_ft[ , j, 1], Op_type, beta_0,
           max_depth)
      
      env_$forests[[1]][[j]] = root_node
    }
  } else {
    #beta_init = rmvnorm(1, mu_beta, Sigma_beta)
    #print(beta_init)
    #env_$beta[1 , ] = init_params$beta
    
    #env_$sigma_square[1] = init_params$sigma_square
    
    for(j in 1:K1) {
      #env_$W_op[ , j, 1] = init_params$W_op[ , j]
      #env_$W_ft[ , j, 1] = init_params$W_ft[ , j]
      
      env_$forests[[1]][[j]] = deepcopy(init_params$forests[[j]])
    }
  }
  
  naccept = 0
  
  for(i in 2:niter) {
    
    # env_$forests[[i]] = env_$forests[[i-1]]
    # forests_tree = forests[[i]]
    #for(j in 1:K1) {
    #  forests[[i]][[j]] = deepcopy(forests[[i-1]][[j]])
    #}
    
    tree_params = update_tree_parameters(alpha_op, alpha_ft, K1, env_$forests[[i - 1]], nfeature, Ops)
    env_$W_op[ , , i] = tree_params$W_op
    env_$W_ft[ , , i] = tree_params$W_ft
    #for(j in 1:K1) {
    temp = update_tree_MH(env_$W_op[ , , i], env_$W_ft[ , , i], w_op_prop, w_ft_prop, K1, X, y,
                          env_$forests[[i - 1]], p_grow, nfeature, Ops_prop, Ops_type_prop, 
                                         beta_0, max_depth, mu_beta, Sigma_beta, nu, lambda)
    env_$forests[[i]] = temp$forests_res
    env_$log_MH[i , ] = temp$log_MH
    naccept = naccept + temp$naccept
    
    model_params = update_model_parameters(y, X, env_$forests[[i]], K1, mu_beta, Sigma_beta, nu, lambda)
    env_$beta[i, ] = model_params$beta
    env_$sigma_square[i] = model_params$sigma_square
    #}
    JMP_trees = log_normalizing_constant(y, X, env_$forests[[i]], K1, mu_beta, Sigma_beta, nu, lambda)
    env_$marg_lik[i] = JMP_trees
    for(j in 1:K1) {
      temp = tree_lik(env_$forests[[i]][[j]], nfeature, Ops)
      Dir_ft = alpha_ft + temp$n_ft
      Dir_op = alpha_op + temp$m_op
      
      s = sum(lgamma(Dir_ft)) + sum(lgamma(Dir_op)) - lgamma(sum(Dir_op)) - lgamma(sum(Dir_ft))
      
      JMP_trees = JMP_trees + s
    }
    
    env_$JMP_trees[i] = JMP_trees
    
    message(paste0("Iteration ", i, " done...................."))
  }
  
  # params_MCMC = list(beta_MCMC = beta, sigma_square_MCMC = sigma_square, W_op_MCMC = W_op, W_ft_MCMC = W_ft, forests_MCMC = forests,
                     #log_MH = log_MH)
  
  
  message(paste0("Acceptance: ", naccept/(niter*K1)))
  return(message("main_MCMC() for HierBOSSS done!"))
}


