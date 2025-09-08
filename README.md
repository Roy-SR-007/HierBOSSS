# HierBOSSS: Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions

[![R](https://img.shields.io/badge/R-4.4.1-blue)](https://cran.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.13.5-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Forks](https://img.shields.io/github/forks/Roy-SR-007/HierBOSSS)](https://github.com/Roy-SR-007/HierBOSSS/network)
[![Repo Size](https://img.shields.io/github/repo-size/Roy-SR-007/HierBOSSS)](https://github.com/Roy-SR-007/HierBOSSS)
[![Last Commit](https://img.shields.io/github/last-commit/Roy-SR-007/HierBOSSS)](https://github.com/Roy-SR-007/HierBOSSS/commits/main)
[![Issues](https://img.shields.io/github/issues/Roy-SR-007/HierBOSSS)](https://github.com/Roy-SR-007/HierBOSSS/issues)
[![Pull Requests](https://img.shields.io/github/issues-pr/Roy-SR-007/HierBOSSS)](https://github.com/Roy-SR-007/HierBOSSS/pulls)

<p align="center">
  <img src="hierbosss_tree.gif" alt="HierBOSSS_logo" width="800"/>
</p>

This repository holds the source code and implementation of **HierBOSSS** for Bayesian structural learning of scientific symbolic expressions proposed in Roy, S., Dey, P., Pati, D., & Mallick, B. K. (2025), *Hierarchical Bayesian Operator-induced Symbolic Regression Trees for Structural Learning of Scientific Expressions*.

---

## Developers and Maintainers

**Somjit Roy**  
Department of Statistics  
Texas A&M University, College Station, TX, USA  

📧 Email: [sroy_123@tamu.edu](mailto:sroy_123@tamu.edu)  
🌐 Website: [https://roy-sr-007.github.io](https://roy-sr-007.github.io)

**Pritam Dey**  
Department of Statistics  
Texas A&M University, College Station, TX, USA  

📧 Email: [pritam.dey@tamu.edu](mailto:pritam.dey@tamu.edu)  
🌐 Website: [https://pritamdey.github.io](https://pritamdey.github.io)

---

## NEWS

- This is the first official release of `HierBOSSS v1.0.0` on `GitHub`.
- Explore HierBOSSS' ability in Bayesian scientific machine learning for discovering and learning various scientific symbolic expressions.
- HierBOSSS has been applied for learning a suite of canonical Feynman equations and discovering important descriptors for single-atom catalysis.
- Comparisons have been presented against state-of-the-art competing symbolic regression modules.

---

## Overview

**Scientific Machine Learning** (SciML) is an interdisciplinary field that combines scientific computing with modern machine learning techniques to solve complex scientific and engineering problems ([Zhong et al., 2020](https://www.nature.com/articles/s41586-020-2242-8); [Hart et al., 2021](https://www.nature.com/articles/s41578-021-00340-w); [Zhang et al., 2025](https://www.mdpi.com/2073-4433/16/1/82); [Yang et al., 2024](https://www.sciencedirect.com/science/article/pii/S1352231024004722)). Across multidisciplinary scientific fields such as *biology* ([Boadu et al., 2025](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/abs/10.1002/pmic.202300471)), *physics* ([Udrescu and Tegmark, 2020](https://www.science.org/doi/abs/10.1126/sciadv.aay2631); [Raissi et al., 2019](https://www.sciencedirect.com/science/article/pii/S0021999118307125)), *climate science* ([Rasp et al., 2018](https://www.pnas.org/doi/abs/10.1073/pnas.1810286115)), and *materials science* ([Schmidt et al., 2019](https://www.nature.com/articles/s41524-019-0221-0)); SciML has received tremendous success for accurate and efficient problem-solving. A core challenge in the above STEM disciplines is deriving interpretable scientific formulas, which are combinations of observable physical quantities with fundamental mathematical operations. **Symbolic regression** is a powerful SciML-driven approach that automatically uncovers scientific truths by capturing (or, learning) symbolic expressions from observed data. It has been an active research area for more than two decades. Applications of symbolic regression extends to *sparse identification for nonlinear dynamical systems* ([Brunton et al., 2016](https://www.pnas.org/doi/abs/10.1073/pnas.1517384113)), *materials discovery* ([Wang et al., 2024](https://onlinelibrary.wiley.com/doi/abs/10.1002/idm2.12180)), and *learning natural scientific laws from experimental datasets* ([Schmidt and Lipson, 2009](https://www.science.org/doi/abs/10.1126/science.1165893)).

Despite the large-scale popularity of symbolic regression in learning hidden scientific mysteries, contributions to symbolic regression majorly centers around algorithmic and heuristic approaches (*genetic programming* ([Willis et al., 1997](https://ieeexplore.ieee.org/document/681044); [Davidson et al, 1999, 2003](https://doi.org/10.2166/hydro.1999.0010)), *deep symbolic regression* (DSR) ([Cranmer et al., 2020](https://proceedings.neurips.cc/paper_files/paper/2020/file/c9f2f917078bd2db12f23c3b413d9cba-Paper.pdf); [Petersen et al., 2021](https://openreview.net/forum?id=m5Qsh0kBQG)), *AI-Feynman* ([Udrescu and Tegmark, 2020](https://www.science.org/doi/abs/10.1126/sciadv.aay2631)), *QLattice* ([Brolos et al., 2021](https://arxiv.org/abs/2104.05417)), *SISSO* ([Ouyang et al., 2018](https://link.aps.org/doi/10.1103/PhysRevMaterials.2.083802)), and *iBART* ([Shengbin et al., 2024](https://doi.org/10.1080/01621459.2023.2294527))). The current literature still lacks a fully model-based framework that enables principled model selection and uncertainty quantification. To address this gap, we develop **Hier**archical **B**ayesian **O**perator-induced **S**ymbolic regression trees for **S**tructural learning of **S**cientific expressions (HierBOSSS), which advances Bayesian SciML and performs symbolic regression through a Bayesian lens by using symbolic tree representation of scientific symbolic expressions (see Figure 1 below) within a hierarchical framework.

<br><br>

<figure align="center">
  <img src="assets/symbolic_tree_representation.png" alt="symbolic_tree_representation" width="800"/>
  <figcaption><em>Figure 1: Symbolic tree representation of scientific expressions.</em></figcaption>
</figure>

<br><br>


**HierBOSSS** models symbolic expressions through an operator-induced sum-of-symbolic trees. Conjugate priors are assigned to model regression parameters, while a regularizing prior is designed for the individual symbolic tree structures. To perform inference from the HierBOSSS-induced posterior distribution, we develop an efficient Metropolis-within-partially-collapsed Gibbs Markov chain Monte Carlo (MCMC) sampling algorithm. This GitHub repository showcases the success of HierBOSSS in discovering interpretable scientific laws. Specifically, we demonstrate HierBOSSS' ability to recover and learn well-known physics-based Feynman equations and identify meaningful descriptors in single-atom catalysis. Moreover, HierBOSSS consistently outperforms state-of-the-art symbolic regression methods (which often suffer from training data size, noise, and overly complicated symbolic expressions), offering advantages both in symbolic expression discovery and computational efficiency.

<br><br>

<figure align="center">
  <img src="assets/HierBOSSS_comparison.png" alt="HierBOSSS_comparison" width="800"/>
  <figcaption><em>Figure 2: Comparing HierBOSSS with competing symbolic regression methods in learning the Coulomb's law.</em></figcaption>
</figure>

---

## Installation and Dependencies

To get started with HierBOSSS, where the initial version (`HierBOSSS v1.0.0`) is programmed using the `R` software, clone the current GitHub repository and install the required dependencies:

* `listenv`
* `dplyr`
* `ggplot2`
* `reshape2`
* `tidyr`
* `mvtnorm`
* `scales`
* `tidyverse`
* `ggpmisc`

```bash
# using SSH on bash
git clone git@github.com:Roy-SR-007/HierBOSSS.git

# or, using HTTPS on bash
git clone https://github.com/Roy-SR-007/HierBOSSS.git
```

For running the competing methods:

1. `iBART`: is moduled on `R`, which can be installed using `install.packages('iBART')` along with the following dependencies:
- `bartMachine`
- `bartMachineJARs`
- `rJava`
- `remotes`

Check [https://github.com/mattsheng/iBART](https://github.com/mattsheng/iBART) for detailed setup instructions.

2. `BSR`: *Bayesian symbolic regression* (BSR) ([Jin et al., 2020](https://arxiv.org/abs/1910.0889)) has been built on `Python`. The folder [BSR](./BSR) contains the entire setup and experimental settings used across different applications considered.

3. `QLattice`: is available through the `Python` package `feyn`. See [https://docs.abzu.ai/](https://docs.abzu.ai/) and [QLattice](./QLattice) for detailed instructions on setup and explore QLattice's application to various experiments and numerical studies.

---

## Functionality and Components of the HierBOSSS Module

The backend of `HierBOSSS` is mainly moduled in the `R` scripts, [MCMC.R](./MCMC.R) and [tree_functions.R](./tree_functions.R).

* [tree_functions.R](./tree_functions.R): This `R` script implements the core symbolic tree structural representation of scientific symbolic expressions modeled by the HierBOSSS framework. It defines reference classes for operators and nodes, and provides a comprehensive set of functions to perform different tree operations as listed below.

<div align="center">

| **Component**                              | **Type**            | **Functionality**                                                                                                                                                                                    |
| ------------------------------------------ | ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Symbolic Tree Utilities ([tree_functions.R](./tree_functions.R))**    |                     |                                                                                                                                                                                                      |
| **Operator (class)**                       | Reference class     | Defines an operator with name, function, and arity (unary/binary).                                                                                                                                   |
| **Node (class)**                           | Reference class     | Represents a symbolic tree node (type, depth, operator, feature, parent/children, data).                                                                                                             |
| `deepcopy(node)`                           | Helper              | Creates a deep copy of a symbolic tree from a given node.                                                                                                                                            |
| `grow(node, …)`                            | Recursive           | Recursively grows a tree using operator/feature weights, splitting probability, and max depth.                                                                                                       |
| `genList(node)`                            | Tree traversal      | Stores tree nodes in a list and assigns order indices.                                                                                                                                               |
| `shrink(node, …)`                          | Tree edit           | Converts a non-terminal node into a terminal node.                                                                                                                                                   |
| `upgOd(Tree)`                              | Helper              | Updates the order indices of nodes in a tree.                                                                                                                                                        |
| `allcal(node, indata)`                     | Evaluation          | Computes the output of a symbolic tree for given input data.                                                                                                                                         |
| `display(Tree)`                            | Visualization       | Prints the symbolic tree structure by depth (operators/features).                                                                                                                                    |
| `getHeight(node)`                          | Property            | Returns the maximum depth (height) of a tree.                                                                                                                                                        |
| `getNum(node)`                             | Property            | Returns the total number of nodes in a tree.                                                                                                                                                         |
| `upDepth(Root)`                            | Recursive           | Updates depth attributes of all nodes starting from root.                                                                                                                                            |
| `Express(node)`                            | Expression builder  | Generates algebraic expression string of the symbolic tree.                                                                                                                                          |
| `tree_lik(node, …)`                        | Likelihood helper   | Computes counts of operator uses and feature assignments.                                                                                                                                            |
| `get_all_terminal(node)`                   | Traversal           | Returns paths to all terminal nodes.                                                                                                                                                                 |
| `access_terminal(path, node)`              | Accessor            | Accesses a specific terminal node by path string.                                                                                                                                                    |
| `get_all_nonterminal(node)`                | Traversal           | Returns paths to all non-terminal nodes.                                                                                                                                                             |
| `access_nonterminal(path, node)`           | Accessor            | Accesses a specific non-terminal node by path string.                                                                                                                                                |
| `tree_log_likelihood_val(node, …)`         | Likelihood          | Computes the log-likelihood of a tree structure under prior weights.                                                                                                                                 |
| `tree_log_proposal_likelihood(…)`          | Proposal            | Computes proposal likelihood for grow/shrink moves (used in MH steps).                                                                                                                               |

</div>

* [MCMC](./MCMC.R): This `R` script implements the Metropolis-within-partially-collapsed Gibbs sampling algorithm to sample from the HierBOSSS-induced posterior distribution. It integrates the symbolic tree utilities from [tree_functions.R](./tree_functions.R) and provides the full posterior sampling routine. The components of [MCMC.R](./MCMC.R) and their corresponding functionalities are highlighted below.

<div align="center">

| **Component**                              | **Type**            | **Functionality**                                                                                                                                                                                    |
| ------------------------------------------ | ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Posterior Sampler ([MCMC.R](./MCMC.R))** |                     |                                                                                                                                                                                                      |
| `update_tree_parameters(…)`                | Gibbs update        | Updates operator weights (`W_op`) and feature weights (`W_ft`) via Dirichlet distributions.                                                                                                          |
| `update_model_parameters(…)`               | Collapsed Gibbs     | Updates regression coefficients $\beta$ and variance $\sigma^2$ using Normal–Inverse-Gamma conjugacy.                                                                                                      |
| `log_normalizing_constant(…)`              | Marginal term       | Computes log normalizing constant after integrating out $\beta$ and $\sigma^2$.                                                                                                                                    |
| `update_tree_MH(…)`                        | Metropolis–Hastings | Proposes grow/shrink moves on symbolic trees and accepts/rejects based on Metropolis-Hastings ratio.                                                                                                           |
| `HierBOSSS(…)`                             | Main driver         | Runs the **Metropolis-within-partially-collapsed Gibbs sampler**: initializes forests, updates weights, performs MH tree updates, updates $\beta$ and $\sigma^2$, computes marginal likelihood, and stores chains in `env_`. |

</div>

---

## HierBOSSS in Action

We present three broad applications of HierBOSSS and investigate its ability in learning and discovering scientific symbolic expressions across diverse scientific domains.

1. **A Simulated Example**

* Data-generating symbolic expression: $\mathbf{y} = 5(x_1 + x_2) + \epsilon,\epsilon \sim \mathcal{N}_n(0_n, \sigma^{2} I_n)$; with $x_1 \sim \mathcal{N}_n(4\mathbf{1}_n, I_n), x_2\sim \mathcal{N}_n(6\mathbf{1}_n, I_n), x_3\sim \mathcal{N}_n(8\mathbf{1}_n, I_n)$.

* HierBOSSS accurately recovers the true data-generating expression while maintaining a low RMSE (see Figure 3).

* Outperforms competing symbolic regression methods (Bayesian symbolic regression (BSR), QLattice, and iBART; see Figures 4 and 5) by consistently identifying simpler, interpretable symbolic trees guided by the joint marginal posterior model selection (Occam's window).

* In terms of computational efficiency, HierBOSSS significantly records lower runtimes when comapared to BSR and iBART while comparing well with the machine learning module QLattice (see Figure 6).

For detailed analysis of the simulated example, see:  

<div align="center">

| **Category**        | **Files / Notebooks** |
|:-------------------:|:----------------------|
| HierBOSSS Implementation | [HierBOSSS_simulated_example.R](./HierBOSSS_simulated_example.R) |
| Results Output      | [results/simulated_example](./results/simulated_example) |
| Competing Methods   | [BSR/Simulated_example.ipynb](./BSR/Simulated_Example.ipynb), [QLattice.ipynb](./QLattice/QLattice.ipynb) |
| Runtime Comparisons | [runtimes/simulated_example](./runtimes/simulated_example) |
| mGED Computation | [GED](./GED) |

</div>

<br><br>

<figure align="center">
  <img src="results/simulated_example/rmse_boxplot_25_HierBOSSS_BSR_QLattice_iBART_horizontal.png" alt="rmse_boxplot" width="800"/>
  <figcaption><em>Figure 3: RMSE of over 25 data regenrations.</em></figcaption>
</figure>

<br><br>

<figure align="center">
  <img src="results/simulated_example/Expression_Complexity_HierBOSSS_BSR_QLattice.png" alt="expression_complexity" width="800"/>
  <figcaption><em>Figure 4: Expression complexity.</em></figcaption>
</figure>

<br><br>

<figure align="center">
  <img src="assets/GED_Boxplot.png" alt="GED_boxplot" width="800"/>
  <figcaption><em>Figure 5: Structural proximity with the true simulated example.</em></figcaption>
</figure>

<br><br>

<figure align="center">
  <img src="runtimes/simulated_example/runtime_boxplot.png" alt="runtime_comparison" width="800"/>
  <figcaption><em>Figure 6: Runtime comparisons.</em></figcaption>
</figure>

<br><br>

2. **Learning Feynman Equations**

* Applied to popular physical laws such as **Coulomb's Law** $(F=(q_1\times q_2)/r_*^2)$, **Lorentz Force** $(F = q\times (E_f + B\times v\times \sin \theta))$, and **Gravitational Potential Energy** $(U = m_1\times m_2\times (1/r_2 - 1/r_1))$.

* HierBOSSS consistently identifies (learns) the correct symbolic forms of these physics laws (structural proximity studied using mGED), even under introduction of different levels of noise in the data.

* Competing methods either overfit producing unduly complex expressions (Bayesian symbolic regression (BSR)), fail under noise (QLattice), or encounter computational issues (iBART).

* HierBOSSS achieves superior predictive accuracy while preserving interpretability.

[HierBOSSS_Feynman_equations.R](./HierBOSSS_Feynman_equations.R) comprises the full implementation of HierBOSSS and iBART under different experimental settings when applied to learn the Feynman equations. Demonstration of the competing methods are available in [BSR](./BSR/) for Bayesian symbolic regression (BSR) and [QLattice](./QLattice/) for QLattice.

<br><br>

<figure align="center">
  <img src="assets/Feynman_expressions.png" alt="Feynman_expressions" width="800"/>
  <figcaption><em>Figure 7: HierBOSSS learns the symbolic tree structures of the scientific laws defined by Feynman equations.</em></figcaption>
</figure>

<br><br>

3. **Single-atom Catalysis Data Study**

* The central task is to identify meaningful descriptors linking physical/chemical properties of metal-oxide pairs to binding energy ([O'Connor et al., 2018](https://www.nature.com/articles/s41929-018-0094-5)).

* HierBOSSS efficiently learns interpretable descriptors that align with domain knowledge while reducing RMSE compared to LASSO*, Bayesian symbolic regression (BSR), QLattice, and iBART.

Explore [single_atom_data_analysis](./single_atom_data_analysis/), [BSR/Single-Atom Catalysis for BSR.ipynb](./BSR/Single-Atom%20Catalysis%20for%20BSR%20test.ipynb), and [QLattice/Single-Atom catalysis.ipynb](./QLattice/Single-Atom%20catalysis.ipynb) for the complete data study.

<br><br>

<figure align="center">
  <img src="single_atom_data_analysis/boxplot_ooRMSE.png" alt="ooRMSE_boxplot" width="800"/>
  <figcaption><em>Figure 8:  Out-of-sample RMSE values across different symbolic trees K.</em></figcaption>
</figure>

<br><br>

<figure align="center">
    <img src="assets/single_atom_expressions.png" alt="single_atom_expressions" width="800"/>
    <figcaption><em>Figure 9: Some important descriptors learned by HierBOSSS influencing the binding energy in single-atom catalysis dataset.</em></figcaption>
</figure>

---

