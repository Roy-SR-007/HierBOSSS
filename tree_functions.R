# =============================================================================
# Copyright 2025. Somjit Roy and Pritam Dey.
# This program implements the coding for the symbolic tree structures in 
# HierBOSSS as developed in:
#   Roy, S., Dey, P., Pati, D., and Mallick, B.K.
#   'Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions'.
# Authors:
#   Somjit Roy (sroy_123@tamu.edu) and Pritam Dey (pritam.dey@tamu.edu)
# =============================================================================

###
# R code for symbolic tree structure development and related tree operations
###

# Operator reference class

Operator <- setRefClass(
  "Operator",
  fields = list(
    name = "character",
    func = "function",
    arity = "numeric"
  ),
  methods = list(
    initialize = function(name, func, arity) {
      .self$name <- name
      .self$func <- func
      .self$arity <- arity
    }
  )
)

# Node reference class

Node <- setRefClass(
  "Node",
  fields = list(
    type = "numeric", # -1, 0, 1, or 2 referring to the type of Node
    order = "numeric",
    left = "ANY", # can be NULL or another Node
    right = "ANY", # can be NULL or another Node
    depth = "numeric", # depth of a particular Node
    parent = "ANY", # can be NULL or another Node
    operator = "character", # a possible operator, binary or unary, based on the bag of operators
    op_ind = "ANY", # index of the operator (can be NULL or numeric)
    data = "ANY", # data associated with the Node
    feature = "ANY" # feature index from the input data (can be NULL or numeric)
    #a = "ANY",
    #b = "ANY"
  ),
  
  # initialize the fields of Node class
  methods = list(
    initialize = function(depth) {
      .self$type <- -1 # newly grown Node
      .self$order <- 0
      .self$left <- NULL
      .self$right <- NULL
      .self$depth <- depth
      .self$parent <- NULL
      .self$operator <- " "
      .self$op_ind <- NULL
      .self$data <- NULL
      .self$feature <- NULL
      #.self$a <- NULL
      #.self$b <- NULL
    },
    
    # print the information of the Node
    inform = function() {
      cat("order:", .self$order, "\n")
      cat("type:", .self$type, "\n")
      cat("depth:", .self$depth, "\n")
      cat("operator:", .self$operator, "\n")
      cat("data:", .self$data, "\n")
      cat("feature:", .self$feature, "\n")
      #if (.self$operator == "ln") {
      # cat("ln_a:", .self$a, "\n")
      # cat("ln_b:", .self$b, "\n")
      #}
    }
  )
)

# deepcopy() function

## The 'deepcopy' function takes in a tree (node) and creates a hard (deep)
## copy of the same

deepcopy = function(node) {
  node_copy = Node$new(depth = node$depth + 1)
  
  node_copy$type = node$type
  node_copy$order = node$order
  node_copy$depth = node$depth
  # node_copy$parent = node$parent
  node_copy$operator = node$operator
  node_copy$op_ind = node$op_ind
  node_copy$data = node$data
  node_copy$feature = node$feature
  
  if(node$type == 1) {
    node_copy$left = deepcopy(node$left)
  }else if(node$type == 2) {
    node_copy$left = deepcopy(node$left)
    node_copy$right = deepcopy(node$right)
  }
  
  return(node_copy)
}

# grow() function: growing a tree

## The 'grow' function takes in a node and grows on it, using the corresponding
## operator weights and feature weights. The 'max depth' argument controls the
## maximum depth of a tree to which it is grown

grow <- function(node, nfeature, Ops, Op_weights, Ft_weights, Op_type, beta, 
                 #sigma_a, sigma_b, 
                 max_depth = 10) {
  # node: particular node of the tree
  # nfeature: total number of features in the input data
  # Ops: the bag of operators considered
  # Op_weights: weights corresponding to each operator considered
  # Ft_weights: weights corresponding to each feature considered
  # Op_type: the type of operators considered, 1: unary and 2: binary
  # beta: controls the grow?split probability of the node
  # max_depth: controls the maximum depth that a node can have
  
  depth <- node$depth # depth of the node
  
  # stop recursion if max_depth is reached
  if (depth >= max_depth) {
    node$type <- 0  # set to terminal node
    node$feature <- sample(0:(nfeature - 1), 1, prob = Ft_weights) # assigning a feature
    return()
  }
  
  # deciding the number of child nodes
  #if (depth > 0) {
    prob <- 1 / ((1 + depth) ^ (beta)) # this is the probability of growing/splitting, i.e., pt = 1 - prob, refer to CGM
    
    test <- runif(1, 0, 1)
    if (test > prob) {  # terminal node
      node$feature <- sample(0:(nfeature - 1), 1, prob = Ft_weights) # assigning a feature uniformly from the set of features according as (1/nfeature, ..., 1/nfeature)
      node$type <- 0 # set to terminal node
      return()  # terminate recursion here if it's a terminal node
    } else {
      op_ind <- sample(1:length(Ops), 1, prob = Op_weights) # picking up an operator index uniformly according as the weights Op_weights
      node$operator <- Ops[op_ind]
      node$type <- Op_type[op_ind]
      node$op_ind <- op_ind
    }
 #}
  #else {  # root node, sure to split
  #  op_ind <- sample(1:length(Ops), 1, prob = Op_weights)
  #  node$operator <- Ops[op_ind]
  #  node$type <- Op_type[op_ind]
  #  node$op_ind <- op_ind
  #}
  
  # use the grow() function recursively on the nodes
  if (node$type == 0) { # terminal node
    node$feature <- sample(0:(nfeature - 1), 1, prob = Ft_weights)
    return()  # terminate recursion at terminal node
  } else if (node$type == 1) { # non-terminal node, having one child (left) node
    node$left <- Node$new(depth = node$depth + 1)
    # node$left$parent <- node
    #if (node$operator == 'ln') {  # linear parameters
    # node$a <- rnorm(1, mean = 1, sd = sqrt(sigma_a))
    # node$b <- rnorm(1, mean = 0, sd = sqrt(sigma_b))
    #}
    grow(node$left, nfeature, Ops, Op_weights, Ft_weights, Op_type, beta, #sigma_a, sigma_b,
         max_depth)
    
  } else if (node$type == 2) { # non-terminal node, having two children (left and right) nodes
    node$left <- Node$new(depth = node$depth + 1)
    # node$left$parent <- node
    node$right <- Node$new(depth = node$depth + 1)
    #node$right$parent <- node
    grow(node$left, nfeature, Ops, Op_weights, Ft_weights, Op_type, beta, #sigma_a, sigma_b, 
         max_depth)
    grow(node$right, nfeature, Ops, Op_weights, Ft_weights, Op_type, beta, #sigma_a, sigma_b, 
         max_depth)
  }
}


# genList() function: to generate list storing the nodes in the tree

## The 'genList' function stores the nodes of a particular tree as a list

genList <- function(node) {
  lst <- list()
  
  # terminal node
  if (is.null(node$left)) {
    lst <- append(lst, list(node))
  } else {
    if (is.null(node$right)) {  # Node with one child
      lst <- append(lst, list(node))
      lst <- c(lst, genList(node$left))
    } else {  # Node with two children
      lst <- append(lst, list(node))
      lst <- c(lst, genList(node$left))
      lst <- c(lst, genList(node$right))
    }
  }
  
  # assigning order to each node
  for (i in seq_along(lst)) {
    lst[[i]]$order <- i - 1  # Setting the order (0-based index)
  }
  
  return(lst)
  
}

# shrink() function: delete all the child nodes of the current node

## The 'shrink' function takes in a node and turns it into a terminal node
## by assigning a feature to it according to feature weights

shrink <- function(node, Ft_weights) {
  if (is.null(node$left)) {
    cat("Already a terminal node!\n")
  } else {
    node$left <- NULL
    node$right <- NULL
    node$type <- 0
    node$operator <- ""  # delete operator
    node$feature <- sample(0:(nfeature - 1), 1, prob = Ft_weights)
    #node$a <- NULL  # delete parameter 'a'
    #node$b <- NULL  # delete parameter 'b'
  }
  return()
}

# upgOd() function: upgrades the order attribute of nodes in the tree

upgOd <- function(Tree) {
  for (i in seq_along(Tree)) {
    Tree[[i]]$order <- i - 1  # Assign order starting from 0 (0-based index)
  }
  return()
}


# allcal() function: computes the tree output from the nodes

## The 'allcal' function considers evaluating the entire symbolic tree expression
## with respect to the input data

allcal <- function(node, indata) { # take in the node and data as the inputs
  if (node$type == 0) {  # terminal node
    if (!is.null(indata)) {
      node$data <- as.matrix(indata[, node$feature + 1])  # Adjust for 1-based indexing in R
    }
  } else if (node$type == 1) {  # one child node
    #if (node$operator == 'ln') {
     # node$data <- node$a * allcal(node$left, indata) + node$b
    #}
    
    # exponential
    if (node$operator == 'exp') {
      node$data <- allcal(node$left, indata)
      for (i in 1:nrow(node$data)) {
        if (node$data[i, 1] <= 200) {
          node$data[i, 1] <- exp(node$data[i, 1])
        } else {
          node$data[i, 1] <- 1e+10
        }
      }
    } else if (node$operator == 'inv') { # inverse
      node$data <- allcal(node$left, indata)
      for (i in 1:nrow(node$data)) {
        if (node$data[i, 1] == 0) {
          node$data[i, 1] <- 0
        } else {
          node$data[i, 1] <- 1 / node$data[i, 1]
        }
      }
    } else if (node$operator == 'neg') { # negative
      node$data <- -1 * allcal(node$left, indata)
    } else if (node$operator == 'sin') { # sin
      node$data <- sin(allcal(node$left, indata))
    } else if (node$operator == 'cos') { # cos
      node$data <- cos(allcal(node$left, indata))
    } else if (node$operator == 'square') {  # square
      node$data <- allcal(node$left, indata)^2
    } else if (node$operator == 'cubic') {  # cubic
      node$data <- allcal(node$left, indata)^3
    } else if(node$operator == 'id') {
      node$data <- allcal(node$left, indata)} 
    else {
      cat("No matching type and operator!\n")
    }
  } else if (node$type == 2) {  # two child nodes
    if (node$operator == '+') { # +, addition
      node$data <- allcal(node$left, indata) + allcal(node$right, indata)
    } else if (node$operator == '*') { # *, multiplication
      node$data <- allcal(node$left, indata) * allcal(node$right, indata)
    } else {
      cat("No matching type and operator!\n")
    }
  } else if (node$type == -1) {  # not grown
    cat("Not a grown tree!\n")
  } else {
    cat("No legal node type!\n")
  }
  
  return(node$data)
}

# display() function: to display the tree structure

## The 'display' function considers displaying or representing the entire
## symbolic tree structure, where the input tree is stored as a list using
## the 'genList' function

display <- function(Tree) {
  # Find the maximum depth of the tree
  tree_depth <- -1
  for (i in seq_along(Tree)) {
    if (Tree[[i]]$depth > tree_depth) {
      tree_depth <- Tree[[i]]$depth
    }
  }
  
  # Create lists to store nodes by depth
  dlists <- vector("list", tree_depth + 1)
  for (d in seq(0, tree_depth)) {
    dlists[[d + 1]] <- list()
  }
  
  # Assign nodes to their respective depth lists
  for (i in seq_along(Tree)) {
    depth <- Tree[[i]]$depth + 1  # Adjust for 1-based indexing
    dlists[[depth]] <- append(dlists[[depth]], list(Tree[[i]]))
  }
  
  # Display the tree structure by depth
  for (d in seq_along(dlists)) {
    st <- ""
    for (i in seq_along(dlists[[d]])) {
      node <- dlists[[d]][[i]]
      if (node$type > 0) {
        st <- paste0(st, node$operator, " ")
      } else {
        st <- paste0(st, node$feature, " ")
      }
    }
    cat(st, "\n")
  }
}

# getHeight() function: to compute the height (depth) of a node

getHeight <- function(node) {
  if (node$type == 0) {
    return(0)  # Terminal node
  } else if (node$type == 1) {  # One child
    return(getHeight(node$left) + 1)
  } else {  # Two children
    lheight <- getHeight(node$left)
    rheight <- getHeight(node$right)
    return(max(lheight, rheight) + 1)
  }
}

# getNum() function: computes the total number of nodes in a tree 

getNum <- function(node) {
  if (node$type == 0) {
    return(1)  # Terminal node
  } else if (node$type == 1) {  # One child
    return(getNum(node$left) + 1)
  } else {  # Two children
    lnum <- getNum(node$left)
    rnum <- getNum(node$right)
    return(lnum + rnum + 1)
  }
}

# upDepth() function: updates the depth of the nodes in a tree

upDepth <- function(Root) {
  if (is.null(Root$parent)) {
    Root$depth <- 0  # Root node has depth 0
  } else {
    Root$depth <- Root$parent$depth + 1  # Depth is parent depth + 1
  }
  
  # Update depth of left child, if present
  if (!is.null(Root$left)) {
    upDepth(Root$left)
    
    # Update depth of right child, if present
    if (!is.null(Root$right)) {
      upDepth(Root$right)
    }
  }
}


# Express() function: expression of the tree

## The 'Express' function gives back the symbolic expression represented by the
## symbolic tree

Express <- function(node) {
  expr <- ""
  
  if (node$type == 0) {  # Terminal node
    expr <- paste0("x", node$feature)
    return(expr)
  } else if (node$type == 1) {  # One child
    if (node$operator == 'exp') {
      expr <- paste0("exp(", Express(node$left), ")")
    } #else if (node$operator == 'ln') {
      #expr <- paste0(round(node$a, 4), "*(", Express(node$left), ") + ", round(node$b, 4))
    #} 
    else if (node$operator == 'inv') {
      expr <- paste0("1/(", Express(node$left), ")")
    } else if (node$operator == 'sin') {
      expr <- paste0("sin(", Express(node$left), ")")
    } else if (node$operator == 'cos') {
      expr <- paste0("cos(", Express(node$left), ")")
    } else if (node$operator == 'square') {
      expr <- paste0("(", Express(node$left), ")^2")
    } else if (node$operator == 'cubic') {
      expr <- paste0("(", Express(node$left), ")^3")
    } else if (node$operator == 'id'){
      expr <- paste0(Express(node$left))}
    else {  # Default: negative
      expr <- paste0("-(", Express(node$left), ")")
    }
  } else {  # Two children (node$type == 2)
    if (node$operator == '+') {
      expr <- paste0(Express(node$left), " + ", Express(node$right))
    } else {
      expr <- paste0("(", Express(node$left), ")*(", Express(node$right), ")")
    }
  }
  
  return(expr)
}

# tree_lik() function for computing the tree likelihood

## The 'tree_lik' computes the components required for the calculation of the 
## tree likelihood value such as, (a_1, a_2, m_op, n_ft), where m_op and n_ft
## is required in the computation of the log-tree likelihood value which is
## performed by the 'tree_log_likelihood_val' function

tree_lik = function(node, nfeature, Ops) {
  
  a = c(0, 0)
  m_op = rep(0, length(Ops))
  n_ft = rep(0, nfeature)
  
  if(node$type == 0) {
    
    a = c(1, 0)
    n_ft[node$feature + 1] = 1 
  
  } else if(node$type == 1) {
    
    a = c(0, 1)
    m_op[node$op_ind] = 1
    
    lchild_lst = tree_lik(node$left, nfeature, Ops)
    
    a = a + lchild_lst$a
    m_op = m_op + lchild_lst$m_op
    n_ft = n_ft + lchild_lst$n_ft
    
  } else if(node$type == 2) {
    
    a = c(0, 1)
    m_op[node$op_ind] = 1
    
    lchild_lst = tree_lik(node$left, nfeature, Ops)
    rchild_lst = tree_lik(node$right, nfeature, Ops)
    
    a = a + lchild_lst$a + rchild_lst$a
    m_op = m_op + lchild_lst$m_op + rchild_lst$m_op
    n_ft = n_ft + lchild_lst$n_ft + rchild_lst$n_ft
    
  }
  
  res_lst = list(a = a, m_op = m_op, n_ft = n_ft)
  return(res_lst)
  
}


# get_all_terminal() function: to get all the terminal nodes of a tree

## The 'get_all_terminal' function takes in a symbolic tree and returns all the 
## path to the terminal nodes of that particular symbolic tree

get_all_terminal = function(node, current = "") {
  
  terminal_nodes = list()
  
  if(node$type == 0) {
    terminal_nodes[[1]] = current
  }else if(node$type == 1) {
    lterm = get_all_terminal(node$left, paste0(current, "1"))
    terminal_nodes = c(terminal_nodes, lterm)
  }else {
    lterm = get_all_terminal(node$left, paste0(current, "1"))
    rterm = get_all_terminal(node$right, paste0(current, "2"))
    terminal_nodes = c(terminal_nodes, lterm, rterm)
  }
  
  return(terminal_nodes)
  
}

# access_terminal(): access the terminal nodes
access_terminal = function(get_all_terminal_pos, node) {
  
  x = as.numeric(strsplit(get_all_terminal_pos, "")[[1]])
  if(length(x) == 0) {
    return(node)
  }
  for(i in 1:length(x)) {
    if(x[i] == 1) {
      node = node$left
    }else {
      node = node$right
    }
  }
  return(node)
  
}

# get_all_nonterminal(): acquires all the nonterminal nodes
get_all_nonterminal = function(node, current = "") {
  
  nonterminal_nodes = list()
  
  if(node$type == 0) {
    return(nonterminal_nodes)
  }
  
  nonterminal_nodes[[1]] = current
  
  if(node$type == 1) {
    lnterm = get_all_nonterminal(node$left, paste0(current, "1"))
    nonterminal_nodes = c(nonterminal_nodes, lnterm)
  }else {
    lnterm = get_all_nonterminal(node$left, paste0(current, "1"))
    rnterm = get_all_nonterminal(node$right, paste0(current, "2"))
    nonterminal_nodes = c(nonterminal_nodes, lnterm, rnterm)
  }
  
  return(nonterminal_nodes)
}

# access_nonterminal(): access all the nonterminal nodes
access_nonterminal = function(get_all_nonterminal_pos, node) {
  
  x = as.numeric(strsplit(get_all_nonterminal_pos, "")[[1]])
  
  if(length(x) == 0) {
    return(node)
  }
  
  for(i in 1:length(x)) {
    if(x[i] == 1) {
      node = node$left
    }else {
      node = node$right
    }
  }
  return(node)
}

# tree_log_likelihood_val(): computes the log-likelihood of the tree structure
tree_log_likelihood_val = function(node, w_op, w_ft, beta_0, nfeature, Ops) {
  
  m_op = tree_lik(node, nfeature, Ops)$m_op
  n_ft = tree_lik(node, nfeature, Ops)$n_ft
  
  log_likelihood = sum(m_op * log(w_op)) + sum(n_ft * log(w_ft))
  
  node_terminal = get_all_terminal(node, "")
  if(length(node_terminal) > 0) {
    tab_terminal = table(sapply(node_terminal, nchar))
    counts_terminal = unname(tab_terminal)
    d_terminal = as.numeric(names(tab_terminal))
    terminal_log_likelihood = sum(counts_terminal * log(1 - (1 + d_terminal) ^ (-beta_0)))
  }else {
    terminal_log_likelihood = 0
  }
  
  
  node_nonterminal = get_all_nonterminal(node, "")
  if(length(node_nonterminal) > 0) {
    tab_nonterminal = table(sapply(node_nonterminal, nchar))
    counts_nonterminal = unname(tab_nonterminal)
    d_nonterminal = as.numeric(names(tab_nonterminal))
    non_terminal_log_likelihood = sum((-beta_0 * counts_nonterminal) * log(1 + d_nonterminal))
  }else {
    non_terminal_log_likelihood = 0
  }
  
  
  log_likelihood = log_likelihood + non_terminal_log_likelihood + terminal_log_likelihood
  
  return(log_likelihood)
}

# tree_log_proposal_likelihood(): computes the log proposal likelihood for the tree structure
tree_log_proposal_likelihood = function(node_prop, node_curr, grow, node_address, w_op_prop, w_ft_prop, p_grow, nfeature, Ops, beta_0) {
  
  if(grow == 0) {
    ft_sel = access_terminal(node_address, node_prop)$feature
    return(log(1 - p_grow) - log(length(get_all_nonterminal(node_curr, ""))) + log(w_ft_prop[ft_sel + 1]))
  }else {
    return(log(p_grow) - log(length(get_all_terminal(node_curr, ""))) + 
             tree_log_likelihood_val(access_nonterminal(node_address, node_prop), w_op_prop, w_ft_prop, beta_0, nfeature, Ops))
  }
  
}
