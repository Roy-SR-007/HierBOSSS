###

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

for(d in 1:data_rep) {
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
}

print('ooRMSE_iBART:')
ooRMSE_iBART
print('train_RMSE_iBART:')
train_iBART