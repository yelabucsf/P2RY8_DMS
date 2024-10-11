library(tidyverse)
# 
# Upload DMS results. Missense_var_z.csv is same content as Sheet 5 of Table S2. 
Missense_var_z <- read_csv("~/Missense_var_z.csv", col_names = T)

# Add position column 
Missense_var_z$pos <- as.integer(sub("[^0-9]", "", Missense_var_z$full_var))

# Add domain definitions
NTD_pos <- c(1:18)
TM_pos <- c(19:51, 56:85, 92:127, 136:160, 186:222, 228:263, 271:296)
IC_pos <- c(52:55, 128:135, 223:227)
EC_pos <- c(86:91, 161:185, 264:270)
H8_pos <- c(297:311)
CTD_pos <- c(312:359)

Missense_var_z <- Missense_var_z %>%
  mutate(domain = case_when(pos %in% NTD_pos ~ "NTD", 
                            pos %in% TM_pos ~ "TM",
                            pos %in% IC_pos ~ "IC",
                            pos %in% EC_pos ~ "EC",
                            pos %in% H8_pos ~ "H8",
                            pos %in% CTD_pos ~ "CTD"))

# Upload regression scores
Regressed_AMquant2_Exp <- read_csv("~/df_amquant_ohllr2_exp.csv", col_names = T)
Regressed_AMquant2_Exp <- Regressed_AMquant2_Exp[,-1]
Regressed_AMquant4_Exp <- read_csv("~/df_amquant_ohllr4_exp.csv", col_names = T)
Regressed_AMquant4_Exp <- Regressed_AMquant4_Exp[,-1]
Regressed_AMquant8_Exp <- read_csv("~/df_amquant_ohllr8_exp.csv", col_names = T)
Regressed_AMquant8_Exp <- Regressed_AMquant8_Exp[,-1]
Regressed_AMquant12_Exp <- read_csv("~/df_amquant_ohllr12_exp.csv", col_names = T)
Regressed_AMquant12_Exp <- Regressed_AMquant12_Exp[,-1]
Regressed_AMquant16_Exp <- read_csv("~/df_amquant_ohllr16_exp.csv", col_names = T)
Regressed_AMquant16_Exp <- Regressed_AMquant16_Exp[,-1]
Regressed_AMquant2_Mig <- read_csv("~/df_amquant_ohllr2_migr.csv", col_names = T)
Regressed_AMquant2_Mig <- Regressed_AMquant2_Mig[,-1]
Regressed_AMquant4_Mig <- read_csv("~/df_amquant_ohllr4_migr.csv", col_names = T)
Regressed_AMquant4_Mig <- Regressed_AMquant4_Mig[,-1]
Regressed_AMquant8_Mig <- read_csv("~/df_amquant_ohllr8_migr.csv", col_names = T)
Regressed_AMquant8_Mig <- Regressed_AMquant8_Mig[,-1]
Regressed_AMquant12_Mig <- read_csv("~/df_amquant_ohllr12_migr.csv", col_names = T)
Regressed_AMquant12_Mig <- Regressed_AMquant12_Mig[,-1]
Regressed_AMquant16_Mig <- read_csv("~/df_amquant_ohllr16_migr.csv", col_names = T)
Regressed_AMquant16_Mig <- Regressed_AMquant16_Mig[,-1]
Regressed_AMquant2_Prol <- read_csv("~/df_amquant_ohllr2_prol.csv", col_names = T)
Regressed_AMquant2_Prol <- Regressed_AMquant_Prol[,-1]
Regressed_AMquant4_Prol <- read_csv("~/df_amquant_ohllr4_prol.csv", col_names = T)
Regressed_AMquant4_Prol <- Regressed_AMquant_Prol[,-1]
Regressed_AMquant8_Prol <- read_csv("~/df_amquant_ohllr8_prol.csv", col_names = T)
Regressed_AMquant8_Prol <- Regressed_AMquant_Prol[,-1]
Regressed_AMquant12_Prol <- read_csv("~/df_amquant_ohllr12_prol.csv", col_names = T)
Regressed_AMquant12_Prol <- Regressed_AMquant_Prol[,-1]
Regressed_AMquant16_Prol <- read_csv("~/df_amquant_ohllr16_prol.csv", col_names = T)
Regressed_AMquant16_Prol <- Regressed_AMquant_Prol[,-1]

Regressed_ESMquant2_Exp <- read_csv("~/df_esmquant_ohllr2_exp.csv", col_names = T)
Regressed_ESMquant2_Exp <- Regressed_ESMquant2_Exp[,-1]
Regressed_ESMquant4_Exp <- read_csv("~/df_esmquant_ohllr4_exp.csv", col_names = T)
Regressed_ESMquant4_Exp <- Regressed_ESMquant4_Exp[,-1]
Regressed_ESMquant8_Exp <- read_csv("~/df_esmquant_ohllr8_exp.csv", col_names = T)
Regressed_ESMquant8_Exp <- Regressed_ESMquant8_Exp[,-1]
Regressed_ESMquant12_Exp <- read_csv("~/df_esmquant_ohllr12_exp.csv", col_names = T)
Regressed_ESMquant12_Exp <- Regressed_ESMquant12_Exp[,-1]
Regressed_ESMquant16_Exp <- read_csv("~/df_esmquant_ohllr16_exp.csv", col_names = T)
Regressed_ESMquant16_Exp <- Regressed_ESMquant16_Exp[,-1]
Regressed_ESMquant2_Mig <- read_csv("~/df_esmquant_ohllr2_migr.csv", col_names = T)
Regressed_ESMquant2_Mig <- Regressed_ESMquant2_Mig[,-1]
Regressed_ESMquant4_Mig <- read_csv("~/df_esmquant_ohllr4_migr.csv", col_names = T)
Regressed_ESMquant4_Mig <- Regressed_ESMquant4_Mig[,-1]
Regressed_ESMquant8_Mig <- read_csv("~/df_esmquant_ohllr8_migr.csv", col_names = T)
Regressed_ESMquant8_Mig <- Regressed_ESMquant8_Mig[,-1]
Regressed_ESMquant12_Mig <- read_csv("~/df_esmquant_ohllr12_migr.csv", col_names = T)
Regressed_ESMquant12_Mig <- Regressed_ESMquant12_Mig[,-1]
Regressed_ESMquant16_Mig <- read_csv("~/df_esmquant_ohllr16_migr.csv", col_names = T)
Regressed_ESMquant16_Mig <- Regressed_ESMquant16_Mig[,-1]
Regressed_ESMquant2_Prol <- read_csv("~/df_esmquant_ohllr2_prol.csv", col_names = T)
Regressed_ESMquant2_Prol <- Regressed_ESMquant_Prol[,-1]
Regressed_ESMquant4_Prol <- read_csv("~/df_esmquant_ohllr4_prol.csv", col_names = T)
Regressed_ESMquant4_Prol <- Regressed_ESMquant_Prol[,-1]
Regressed_ESMquant8_Prol <- read_csv("~/df_esmquant_ohllr8_prol.csv", col_names = T)
Regressed_ESMquant8_Prol <- Regressed_ESMquant_Prol[,-1]
Regressed_ESMquant12_Prol <- read_csv("~/df_esmquant_ohllr12_prol.csv", col_names = T)
Regressed_ESMquant12_Prol <- Regressed_ESMquant_Prol[,-1]
Regressed_ESMquant16_Prol <- read_csv("~/df_esmquant_ohllr16_prol.csv", col_names = T)
Regressed_ESMquant16_Prol <- Regressed_ESMquant_Prol[,-1]

# Merge regression scores with some columns (full_var, Expr_z_score, Migr_z_score, Prolif_z_score, pos,
# and domain) from Missense_var_z
Missense_subset <- Missense_var_z[,c(3,4,5,6,9,10)]

AM_ExpQuant2_Regr <- left_join(Missense_subset, Regressed_AMquant2_Exp, by = "full_var")
AM_MigQuant2_Regr <- left_join(Missense_subset, Regressed_AMquant2_Mig, by = "full_var")
AM_ProlQuant2_Regr <- left_join(Missense_subset, Regressed_AMquant2_Prol, by = "full_var")
AM_ExpQuant4_Regr <- left_join(Missense_subset, Regressed_AMquant4_Exp, by = "full_var")
AM_MigQuant4_Regr <- left_join(Missense_subset, Regressed_AMquant4_Mig, by = "full_var")
AM_ProlQuant4_Regr <- left_join(Missense_subset, Regressed_AMquant4_Prol, by = "full_var")
AM_ExpQuant8_Regr <- left_join(Missense_subset, Regressed_AMquant8_Exp, by = "full_var")
AM_MigQuant8_Regr <- left_join(Missense_subset, Regressed_AMquant8_Mig, by = "full_var")
AM_ProlQuant8_Regr <- left_join(Missense_subset, Regressed_AMquant8_Prol, by = "full_var")
AM_ExpQuant12_Regr <- left_join(Missense_subset, Regressed_AMquant12_Exp, by = "full_var")
AM_MigQuant12_Regr <- left_join(Missense_subset, Regressed_AMquant12_Mig, by = "full_var")
AM_ProlQuant12_Regr <- left_join(Missense_subset, Regressed_AMquant12_Prol, by = "full_var")
AM_ExpQuant16_Regr <- left_join(Missense_subset, Regressed_AMquant16_Exp, by = "full_var")
AM_MigQuant16_Regr <- left_join(Missense_subset, Regressed_AMquant16_Mig, by = "full_var")
AM_ProlQuant16_Regr <- left_join(Missense_subset, Regressed_AMquant16_Prol, by = "full_var")

ESM_ExpQuant2_Regr <- left_join(Missense_subset, Regressed_AMquant2_Exp, by = "full_var")
ESM_MigQuant2_Regr <- left_join(Missense_subset, Regressed_AMquant2_Mig, by = "full_var")
ESM_ProlQuant2_Regr <- left_join(Missense_subset, Regressed_AMquant2_Prol, by = "full_var")
ESM_ExpQuant4_Regr <- left_join(Missense_subset, Regressed_AMquant4_Exp, by = "full_var")
ESM_MigQuant4_Regr <- left_join(Missense_subset, Regressed_AMquant4_Mig, by = "full_var")
ESM_ProlQuant4_Regr <- left_join(Missense_subset, Regressed_AMquant4_Prol, by = "full_var")
ESM_ExpQuant8_Regr <- left_join(Missense_subset, Regressed_AMquant8_Exp, by = "full_var")
ESM_MigQuant8_Regr <- left_join(Missense_subset, Regressed_AMquant8_Mig, by = "full_var")
ESM_ProlQuant8_Regr <- left_join(Missense_subset, Regressed_AMquant8_Prol, by = "full_var")
ESM_ExpQuant12_Regr <- left_join(Missense_subset, Regressed_AMquant12_Exp, by = "full_var")
ESM_MigQuant12_Regr <- left_join(Missense_subset, Regressed_AMquant12_Mig, by = "full_var")
ESM_ProlQuant12_Regr <- left_join(Missense_subset, Regressed_AMquant12_Prol, by = "full_var")
ESM_ExpQuant116_Regr <- left_join(Missense_subset, Regressed_AMquant16_Exp, by = "full_var")
ESM_MigQuant116_Regr <- left_join(Missense_subset, Regressed_AMquant16_Mig, by = "full_var")
ESM_ProlQuant16_Regr <- left_join(Missense_subset, Regressed_AMquant16_Prol, by = "full_var")

#Subset from AM k = 2 and k = 4 to have only first 100 columns 
AM_ExpQuant2_Regr100 <- AM_ExpQuant2_Regr[,1:106]
AM_MigQuant2_Regr100 <- AM_MigQuant2_Regr[,1:106]
AM_ProlQuant2_Regr100 <- AM_ProlQuant2_Regr[,1:106] 
AM_ExpQuant4_Regr100 <- AM_ExpQuant4_Regr[,1:106]
AM_MigQuant4_Regr100 <- AM_MigQuant4_Regr[,1:106]
AM_ProlQuant4_Regr100 <- AM_ProlQuant4_Regr[,1:106]

#Function to calculate Spearman correlation and output into new tibble 
calc_SpCorr <- function(data, pred_column, var) {
  # Calculate for the whole dataset
  all_res <- cor(data[[pred_column]], data[, var], method = "spearman", use = "complete.obs")
  
  # Calculate for each domain
  ntd <- cor(data[data$domain == "NTD", pred_column], data[data$domain == "NTD", var], method = "spearman", use = "complete.obs")
  tm <- cor(data[data$domain == "TM", pred_column], data[data$domain == "TM", var], method = "spearman", use = "complete.obs")
  ec <- cor(data[data$domain == "EC", pred_column], data[data$domain == "EC", var], method = "spearman", use = "complete.obs")
  ic <- cor(data[data$domain == "IC", pred_column], data[data$domain == "IC", var], method = "spearman", use = "complete.obs")
  h8 <- cor(data[data$domain == "H8", pred_column], data[data$domain == "H8", var], method = "spearman", use = "complete.obs")
  ctd <- cor(data[data$domain == "CTD", pred_column], data[data$domain == "CTD", var], method = "spearman", use = "complete.obs")
  
  # Return as a named vector (this will become a row in the final tibble)
  c(all_res = all_res, NTD = ntd, TM = tm, EC = ec, IC = ic, H8 = h8, CTD = ctd)
}

pred_col <- grep("pred_y_", names(AM_ExpQuant2_Regr100), value = TRUE)

# New tibbles with Spearman correlations for these iterations of regression scores 
Sp_AM_ExpQuant2_Regr100 <- map_dfr(pred_col, ~ calc_SpCorr(AM_ExpQuant2_Regr100, .x, "Expr_z_score"), .id = "column")
Sp_AM_MigQuant2_Regr100 <- map_dfr(pred_col, ~ calc_SpCorr(AM_MigQuant2_Regr100, .x, "Migr_z_score"), .id = "column")
Sp_AM_ProlQuant2_Regr100 <- map_dfr(pred_col, ~ calc_SpCorr(AM_ProlQuant2_Regr100, .x, "Proli_z_score"), .id = "column")
Sp_AM_ExpQuant4_Regr100 <- map_dfr(pred_col, ~ calc_SpCorr(AM_ExpQuant4_Regr100, .x, "Expr_z_score"), .id = "column")
Sp_AM_MigQuant4_Regr100 <- map_dfr(pred_col, ~ calc_SpCorr(AM_MigQuant4_Regr100, .x, "Migr_z_score"), .id = "column")
Sp_AM_ProlQuant4_Regr100 <- map_dfr(pred_col, ~ calc_SpCorr(AM_ProlQuant4_Regr100, .x, "Prolif_z_score"), .id = "column")
Sp_AM_ExpQuant8_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ExpQuant8_Regr, .x, "Expr_z_score"), .id = "column")
Sp_AM_MigQuant8_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_MigQuant8_Regr, .x, "Migr_z_score"), .id = "column")
Sp_AM_ProlQuant8_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ProlQuant12_Regr, .x, "Prolif_z_score"), .id = "column")
Sp_AM_ExpQuant12_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ExpQuant12_Regr, .x, "Expr_z_score"), .id = "column")
Sp_AM_MigQuant12_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_MigQuant12_Regr, .x, "Migr_z_score"), .id = "column")
Sp_AM_ProlQuant12_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ProlQuant16_Regr, .x, "Prolif_z_score"), .id = "column")
Sp_AM_ExpQuant16_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ExpQuant16_Regr, .x, "Expr_z_score"), .id = "column")
Sp_AM_MigQuant16_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_MigQuant16_Regr, .x, "Migr_z_score"), .id = "column")
Sp_AM_ProlQuant16_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ProlQuant16_Regr, .x, "Prolif_z_score"), .id = "column")

Sp_ESM_ExpQuant2_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ExpQuant2_Regr, .x, "Expr_z_score"), .id = "column")
Sp_ESM_MigQuant2_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_MigQuant2_Regr, .x, "Migr_z_score"), .id = "column")
Sp_ESM_ProlQuant2_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ProlQuant2_Regr, .x, "Prolif_z_score"), .id = "column")
Sp_ESM_ExpQuant4_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ExpQuant2_Regr, .x, "Expr_z_score"), .id = "column")
Sp_ESM_MigQuant4_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_MigQuant2_Regr, .x, "Migr_z_score"), .id = "column")
Sp_ESM_ProlQuant4_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ProlQuant2_Regr, .x, "Prolif_z_score"), .id = "column")
Sp_ESM_ExpQuant8_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ExpQuant2_Regr, .x, "Expr_z_score"), .id = "column")
Sp_ESM_MigQuant8_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_MigQuant2_Regr, .x, "Migr_z_score"), .id = "column")
Sp_ESM_ProlQuant8_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ProlQuant2_Regr, .x, "Prolif_z_score"), .id = "column")
Sp_ESM_ExpQuant12_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ExpQuant2_Regr, .x, "Expr_z_score"), .id = "column")
Sp_ESM_MigQuant12_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_MigQuant2_Regr, .x, "Migr_z_score"), .id = "column")
Sp_ESM_ProlQuant12_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ProlQuant2_Regr, .x, "Prolif_z_score"), .id = "column")
Sp_ESM_ExpQuant16_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ExpQuant2_Regr, .x, "Expr_z_score"), .id = "column")
Sp_ESM_MigQuant16_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_MigQuant2_Regr, .x, "Migr_z_score"), .id = "column")
Sp_ESM_ProlQuant16_Regr <- map_dfr(pred_col, ~ calc_SpCorr(ESM_ProlQuant2_Regr, .x, "Prolif_z_score"), .id = "column")

# Now using full 1000 iterations from k = 2 and k = 4
pred_col <- grep("pred_y_", names(AM_ExpQuant2_Regr), value = TRUE)
Sp_AM_ExpQuant2_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ExpQuant2_Regr, .x, "Expr_z_score"), .id = "column")
Sp_AM_MigQuant2_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_MigQuant2_Regr, .x, "Migr_z_score"), .id = "column")
Sp_AM_ProlQuant2_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ProlQuant2_Regr, .x, "Proli_z_score"), .id = "column")
Sp_AM_ExpQuant4_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ExpQuant4_Regr, .x, "Expr_z_score"), .id = "column")
Sp_AM_MigQuant4_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_MigQuant4_Regr, .x, "Migr_z_score"), .id = "column")
Sp_AM_ProlQuant4_Regr <- map_dfr(pred_col, ~ calc_SpCorr(AM_ProlQuant4_Regr, .x, "Prolif_z_score"), .id = "column")

#Min, 2.5%ile, median, 97.5%ile, and max for each, overall
quantile(Sp_AM_ExpQuant2_Regr100$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant2_Regr100$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant2_Regr100$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant4_Regr100$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant4_Regr100$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant4_Regr100$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant8_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant8_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant8_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant12_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant12_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant12_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant16_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant16_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant16_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))

quantile(Sp_ESM_ExpQuant2_Regr100$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_MigQuant2_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ProlQuant2_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ExpQuant4_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_MigQuant4_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ProlQuant4_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ExpQuant8_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_MigQuant8_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ProlQuant8_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ExpQuant12_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_MigQuant12_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ProlQuant12_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ExpQuant16_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_MigQuant16_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_ESM_ProlQuant16_Regr$all_res, probs = c(0,0.025,0.5,0.975,1))

#Domain specific min, 2.5%ile, median, 97.5%ile, and max for each, overall
quantile(Sp_AM_ExpQuant2_Regr$NTD, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant2_Regr$TM, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant2_Regr$EC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant2_Regr$IC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant2_Regr$H8, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant2_Regr$CTD, probs = c(0,0.025,0.5,0.975,1))

quantile(Sp_AM_MigQuant2_Regr$NTD, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant2_Regr$TM, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant2_Regr$EC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant2_Regr$IC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant2_Regr$H8, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant2_Regr$CTD, probs = c(0,0.025,0.5,0.975,1))

quantile(Sp_AM_ProlQuant2_Regr$NTD, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant2_Regr$TM, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant2_Regr$EC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant2_Regr$IC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant2_Regr$H8, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant2_Regr$CTD, probs = c(0,0.025,0.5,0.975,1))

quantile(Sp_AM_ExpQuant4_Regr$NTD, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant4_Regr$TM, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant4_Regr$EC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant4_Regr$IC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant4_Regr$H8, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ExpQuant4_Regr$CTD, probs = c(0,0.025,0.5,0.975,1))

quantile(Sp_AM_MigQuant4_Regr$NTD, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant4_Regr$TM, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant4_Regr$EC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant4_Regr$IC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant4_Regr$H8, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_MigQuant4_Regr$CTD, probs = c(0,0.025,0.5,0.975,1))

quantile(Sp_AM_ProlQuant4_Regr$NTD, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant4_Regr$TM, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant4_Regr$EC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant4_Regr$IC, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant4_Regr$H8, probs = c(0,0.025,0.5,0.975,1))
quantile(Sp_AM_ProlQuant4_Regr$CTD, probs = c(0,0.025,0.5,0.975,1))
