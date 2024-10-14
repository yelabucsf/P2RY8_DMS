library(tidyverse)
library(stringr)

# Load Enrich2 scores. These csv files are equivalent to sheets 2-4 in Table_S2, columns 1 and 11
# but with col names changed to full_var from value and to Expr_score, Migr_score, or Prolif_score 
# from score. That was done in Excel. 
Expr_scores <-read_csv(
  "~/C042_Expr_scores.csv", 
  col_names = T)

Migr_scores <-read_csv(
  "~/C042_Migr_scores.csv", 
  col_names = T)

Prolif_scores <-read_csv(
  "~/C042_Prolif_scores.csv", 
  col_names = T)

# Enrich2 scores were generated such that a positive score indicated enriched in low expression bin, 
# but to ease interpretation, altered to negative score indicates low expressing 
Expr_scores$Inv_expr_score <- (-1 * Expr_scores$Expr_score)

# From full_var (form Q2V eg) extract WT and var amino acids 
Expr_scores$wt_aa = str_sub(Expr_scores$full_var, 1, 1)
Expr_scores$var_aa = str_sub(Expr_scores$full_var, -1, -1)

# New table with just synonymous variant data
Expr_scores_synon <- Expr_scores %>%
  filter(wt_aa == var_aa)

Expr_mean <- mean(Expr_scores_synon$Inv_expr_score)
Expr_sd <- sd(Expr_scores_synon$Inv_expr_score)

Expr_mean
Expr_sd

# By design, the Enrich2 score is meant to have 0 = WT phenotype (i.e. geometric mean of synon var counts)
# Expr_mean = -0.0082 z(),SD = 0.5131, mean sufficiently close to 0 that 
# was felt using 0 as center of WT was appropriate and did not subtract mean to calc z-score
Expr_scores$Expr_z_score <- (Expr_scores$Inv_expr_score/Expr_sd)

# Same approach now for proliferation and migration 
Migr_scores$wt_aa = str_sub(Migr_scores$full_var, 1, 1)
Migr_scores$var_aa = str_sub(Migr_scores$full_var, -1, -1)

Migr_scores_synon <- Migr_scores %>%
  filter(wt_aa == var_aa)

Migr_mean <- mean(Migr_scores_synon$Migr_score)
Migr_sd <- sd(Migr_scores_synon$Migr_score)

Migr_mean
Migr_sd

#> Migr_mean
#[1] 0.0007343935
#> Migr_sd
#[1] 0.1813575

# Mean = 0.0007, SD = 0.1814
Migr_scores$Migr_z_score <- (Migr_scores$Migr_score/Migr_sd)

Prolif_scores$wt_aa = str_sub(Prolif_scores$full_var, 1, 1)
Prolif_scores$var_aa = str_sub(Prolif_scores$full_var, -1, -1)

Prolif_scores_synon <- Prolif_scores %>%
  filter(wt_aa == var_aa)

Prolif_mean <- mean(Prolif_scores_synon$Prolif_score)
Prolif_sd <- sd(Prolif_scores_synon$Prolif_score)

Prolif_mean
Prolif_sd

# Mean = -0.0002, SD = 0.1594
Prolif_scores$Prolif_z_score <- (Prolif_scores$Prolif_score/Prolif_sd)