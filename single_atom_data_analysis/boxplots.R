library(tidyverse)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

load("single_atom_data_analysis/intermediate.RData")
QLattice = read.csv("single_atom_data_analysis/test_rmse_QLattice_Single_Atom.csv")
BSR = as.matrix(read.csv("single_atom_data_analysis/test_rmse_BSR_Single_Atom.csv"))
load("single_atom_data_analysis/iBART_intermediate.RData")


HierBOSSS_oormse = train_RMSE
QLattice_oormse = QLattice$oormse
BSR_oormse = BSR
iBART_oormse = ooRMSE_iBART

oormse_data = as.matrix(cbind(HierBOSSS_oormse, BSR_oormse, iBART_oormse))
colnames(oormse_data) <- c("HierBOSSS_K2", "HierBOSSS_K3", "HierBOSSS_K4", "HierBOSSS_K5",
                          "BSR_K2", "BSR_K3", "BSR_K4", "BSR_K5",
                          "iBART_K2", "iBART_K3", "iBART_K4", "iBART_K5")


# Convert to data frame
oormse_df <- as.data.frame(oormse_data)
# Add repetition ID
oormse_df$Rep <- 1:nrow(oormse_df)
# Reshape from wide to long format
tidy_oormse <- pivot_longer(oormse_df, 
                            cols = -Rep,
                            names_to = c("Method", "K"),
                            names_sep = "_K",
                            values_to = "oormse") %>%
  mutate(K = as.factor(K))


ggplot(tidy_oormse, aes(x = K, y = oormse, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8))+
  
  labs(title = "Out-of-Sample RMSE by K and Method",
       x = "Number of trees (K)",
       y = "Out-of-Sample RMSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

# Create a data frame for QLattice (single boxplot)
qlattice_df <- data.frame(
  Method = "QLattice",  # Method name for QLattice
  K = factor(rep(3.5, 25)),  # Using K=3.5 to center QLattice
  oormse = QLattice_oormse  # The 25 RMSE values for QLattice
)

# Combine QLattice data with your existing tidy_oormse data
combined_df <- bind_rows(tidy_oormse, qlattice_df)

# Plot the combined data with the boxplot for QLattice overlaid
bp <- ggplot(combined_df, aes(x = K, y = oormse, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape=NA) +
  
  # Add a boxplot for QLattice values
  geom_boxplot(data = qlattice_df, aes(x = K, y = oormse, fill = "QLattice"), 
               alpha = 0.3, color = "black", position = position_dodge(width = 0.8), outlier.shape=NA) +
  
  labs(title = "Out-of-sample RMSE by K and Method across 25 re-samplings of the dataset",
       x = "Number of descriptors (K) & RMSE of QLattice optimal model across 25 re-samplings",
       y = "Out-of-sample RMSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  #scale_fill_manual(values = c("QLattice" = "firebrick")) +
  scale_x_discrete(labels = c("2" = "K=2", "3" = "K=3", "4" = "K=4", "5" = "K=5", "3.5" = "RMSE of QLattice"))  # Set x-axis labels
ggsave("single_atom_data_analysis/boxplot_ooRMSE.png", plot = bp, width = 8, height = 6, dpi = 300)
