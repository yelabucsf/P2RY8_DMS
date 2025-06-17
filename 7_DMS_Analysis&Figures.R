library(tidyverse)

################
# DMS Score Analysis 

#created three CSV file containing Column A ("full_var") and Column K (overall score) data from 
#Supplemental Table 2, tab 2 (2_Enrich2_Expression), 3 (3_Enrich2_Migration), or 4
#(4_Enrich2_Proliferation), referred to as FilePath1, FilePath2, and FilePath3

FilePath1 <- "~/Cyster-Ye/Comp Bio/P2RY8_SatMut/C042_PSM_2/Enrich2_Output/C042_Expr_scores.csv"
FilePath2 <- "~/Cyster-Ye/Comp Bio/P2RY8_SatMut/C042_PSM_2/Enrich2_Output/C042_Migr_scores.csv"
FilePath3 <- "~/Cyster-Ye/Comp Bio/P2RY8_SatMut/C042_PSM_2/Enrich2_Output/C042_Prolif_scores.csv"

Expr_scores <-read_csv(FilePath1, col_names = T)

Migr_scores <-read_csv(FilePath2, col_names = T)

Prolif_scores <-read_csv(FilePath3, col_names = T)

# Import ESM1b and AlphaMissense variant scores. FilePath5 refers to a CSV file created in 
#Microsoft Excel that contains two columns: full_var and ESM1b_score. Modified from downloadable
#file available at https://huggingface.co/spaces/ntranoslab/esm_variants, P2RY8
#Uniprot ID Q86VZ1. Similarly, FilePath5 refers to a CSV file with full_var and am_pathogenicity
#scores and am_class. These values can be obtained at https://alphamissense.hegelab.org. 
ESM1b_tb <- read_csv(FilePath4, col_names = T)
P2RY8_AlphaMissense <- read_csv(FilePath5, col_names = T)

# Multiplicative inverse of Enrich2 expression score, such that higher score now indicates
# the variant has higher surface expression 
Expr_scores$Inv_expr_score <- (-1 * Expr_scores$Expr_score)

#Create columns of WT and variant amino acids 
Expr_scores$wt_aa = str_sub(Expr_scores$full_var, 1, 1)
Expr_scores$var_aa = str_sub(Expr_scores$full_var, -1, -1)
Migr_scores$wt_aa = str_sub(Expr_scores$full_var, 1, 1)
Migr_scores$var_aa = str_sub(Expr_scores$full_var, -1, -1)
Prolif_scores$wt_aa = str_sub(Expr_scores$full_var, 1, 1)
Prolif_scores$var_aa = str_sub(Expr_scores$full_var, -1, -1)

#Generate z-scores based on synonymous variant results. The anticipated mean is 0 based
#on the methodology of Enrich2. 

Expr_scores_synon <- Expr_scores %>%
  filter(wt_aa == var_aa)
Expr_mean <- mean(Expr_scores_synon$Inv_expr_score)
Expr_sd <- sd(Expr_scores_synon$Inv_expr_score)

Expr_mean
#[1] -0.008197562
Expr_sd
#[1] 0.5131151
# Treated mean as 0 given <0.01, less than 0.02 SD
Expr_scores$Expr_z_score <- (Expr_scores$Inv_expr_score/Expr_sd)

Migr_scores_synon <- Migr_scores %>%
  filter(wt_aa == var_aa)
Migr_mean <- mean(Migr_scores_synon$Migr_score)
Migr_sd <- sd(Migr_scores_synon$Migr_score)

Migr_mean
#[1] 0.0007343935
Migr_sd
#[1] 0.1813575
# Again, treated mean as 0 given < 0.01 SD
Migr_scores$Migr_z_score <- (Migr_scores$Migr_score/Migr_sd)

Prolif_scores_synon <- Prolif_scores %>%
  filter(wt_aa == var_aa)
Prolif_mean <- mean(Prolif_scores_synon$Prolif_score)
Prolif_sd <- sd(Prolif_scores_synon$Prolif_score)

Prolif_mean
#[1] -0.0002118805
Prolif_sd
#[1] 0.1593699
# Again, treated mean as 0 given < 0.01 SD 
Prolif_scores$Prolif_z_score <- (Prolif_scores$Prolif_score/Prolif_sd)

#Compiled into a single tibble. 
All_var <- left_join(Expr_scores, select(Migr_scores, full_var, Migr_score, Migr_z_score), 
                       by = "full_var")
All_var <- left_join(All_var, select(Prolif_scores, full_var, Prolif_score, Prolif_z_score), 
                       by = "full_var")

#Calculation of expression-adjusted migration and proliferation scores
#First, best-fit line for expression vs migration for synonymous variants 
Synon_var <- All_var %>%
  filter(wt_aa == var_aa)
lm_Synon_E_M <- lm(Synon_var$Migr_score ~ Synon_var$Inv_expr_score)

summary(lm_Synon_E_M)

# Call:
#   lm(formula = Synon_var$Migr_score ~ Synon_var$Inv_expr_score)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.65768 -0.05568  0.00638  0.06751  0.26099 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.001677   0.005478  -0.306     0.76    
# Synon_var$Inv_expr_score -0.294141   0.010691 -27.513   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1007 on 336 degrees of freedom
# Multiple R-squared:  0.6926,	Adjusted R-squared:  0.6917 
# F-statistic:   757 on 1 and 336 DF,  p-value: < 2.2e-16
#Then, best-fit line for expression vs proliferation for synonymous variants

lm_Synon_E_P <- lm(Synon_var$Prolif_score ~ Synon_var$Inv_expr_score)
summary(lm_Synon_E_P)

# Call:
#   lm(formula = Synon_var$Prolif_score ~ Synon_var$Inv_expr_score)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.23841 -0.06974 -0.01311  0.04884  0.62722 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.001989   0.006217   -0.32    0.749    
# Synon_var$Inv_expr_score -0.216824   0.012132  -17.87   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1143 on 336 degrees of freedom
# Multiple R-squared:  0.4873,	Adjusted R-squared:  0.4858 
# F-statistic: 319.4 on 1 and 336 DF,  p-value: < 2.2e-16

#Then determine expression-adjusted migration and proliferation scores by distance between variant score
#and these best-fit lines. Return -distance such that points above best-fit line have
#positive values. 
Eadj_function <- function(x0, y0, slope, intercept) {
  distance <- (slope * x0 - y0 + intercept) / sqrt(slope^2 + 1)
  return(-distance)
}

All_var <- All_var %>%
  mutate(E_adj_M_score = Eadj_function(Inv_expr_score, Migr_score, 
                                       coef(lm_Synon_E_M)[2],
                                       coef(lm_Synon_E_M)[1])) %>%
  mutate(E_adj_P_score = Eadj_function(Inv_expr_score, Prolif_score, 
                                       coef(lm_Synon_E_P)[2],
                                       coef(lm_Synon_E_P)[1]))


plot(Expr_scores_synon$Inv_expr_score, Migr_scores_synon$Migr_score, 
     main = "Linear Regression", xlab = "Expression", ylab = "Migration")
abline(lm_Synon_E_M, col = "red")  

plot(Missense_var$Inv_expr_score, Missense_var$Migr_score, 
     main = "Linear Regression", xlab = "Expression", ylab = "Migration")
abline(lm_Synon_E_M, col = "red") 

plot(Expr_scores_synon$Inv_expr_score, Prolif_scores_synon$Prolif_score, 
     main = "Linear Regression", xlab = "Expression", ylab = "Proliferation")
abline(lm_Synon_E_P, col = "red")  

plot(Missense_var$Inv_expr_score, Missense_var$Prolif_score, 
     main = "Linear Regression", xlab = "Expression", ylab = "Proliferation")
abline(lm_Synon_E_P, col = "red") 

cor(All_var$Migr_z_score, All_var$Prolif_z_score, method = "pearson")
#[1] 0.7036408
cor(All_var$Migr_z_score, All_var$Prolif_z_score, method = "spearman")
#[1] 0.7245299
plot(All_var$E_adj_M_score, All_var$E_adj_P_score)
cor(All_var$E_adj_M_score, All_var$E_adj_P_score, method = "pearson")
#[1] 0.6321562
cor(All_var$E_adj_M_score, All_var$E_adj_P_score, method = "spearman")
#[1] 0.5128761

# Merge in ESM1b data and AlphaMissense data
ESM1b_tb$ESM1b_score = ESM1b_tb$score
ESM1b_tb <- ESM1b_tb %>%
  select(full_var, pos, ESM1b_score)
All_var <- left_join(ESM1b_tb, All_var, by = "full_var")
All_var <- All_var %>%
  mutate(wt_aa = str_sub(full_var, 1, 1)) %>%
  mutate(var_aa = str_sub(full_var, -1, -1))
All_var <- left_join(All_var, P2RY8_AlphaMissense, by = "full_var")

# New column specifying position 
All_var <- All_var %>%
  mutate(pos = as.double(str_sub(full_var, 2, -2)))

# Additional classifications
hydrophob_al <- c("M", "V", "L", "I", "A")
hydrophob_ar <- c("W", "Y", "F")
h_bonding <- c("Q", "N", "T", "S")
charged_pos <- c("R", "K", "H")
charged_neg <- c("E", "D")

All_var$var_aa <- factor(All_var$var_aa, 
                                  levels = c("C","P", "G", charged_neg,
                                             charged_pos, h_bonding,
                                             hydrophob_ar, hydrophob_al))
All_var <- All_var %>%
  mutate(var_type = ifelse(wt_aa == var_aa, "Synonymous", "Missense"))

#Define protein segments
NTD_pos <- c(1:18)
TM_pos <- c(19:51, 56:85, 92:127, 136:160, 186:222, 228:263, 271:296)
IC_pos <- c(52:55, 128:135, 223:227)
EC_pos <- c(86:91, 161:185, 264:270)
H8_pos <- c(297:311)
CTD_pos <-c(312:359)

#Make tibble summarizing data by position 
Missense_var <- All_var %>% filter(var_type == "Missense")
Results_by_pos <- Missense_var %>%
  group_by(pos) %>%
  summarize(
    Mean_expr_z_score = mean(Expr_z_score),
    Mean_migr_z_score = mean(Migr_z_score),
    Mean_prolif_z_score = mean(Prolif_z_score),
    Mean_EaM_score = mean(E_adj_M_score),
    Mean_EaP_score = mean(E_adj_P_score),
    Mean_ESM1b_score = mean(ESM1b_score),
    Mean_am_pathogenicity = mean(am_pathogenicity),
  )

Results_by_pos <- Results_by_pos %>%
  mutate(domain = case_when(pos %in% NTD_pos ~ "NTD", 
                            pos %in% TM_pos ~ "TM",
                            pos %in% IC_pos ~ "IC",
                            pos %in% EC_pos ~ "EC",
                            pos %in% H8_pos ~ "H8",
                            pos %in% CTD_pos ~ "CTD"))

Results_by_pos$domain <- factor(Results_by_pos$domain, 
                                levels = c("NTD", "TM", "EC", "IC", "H8", "CTD"))

All_var <- All_var %>%
  mutate(domain = case_when(pos %in% NTD_pos ~ "NTD", 
                            pos %in% TM_pos ~ "TM",
                            pos %in% IC_pos ~ "IC",
                            pos %in% EC_pos ~ "EC",
                            pos %in% H8_pos ~ "H8",
                            pos %in% CTD_pos ~ "CTD"))

All_var$domain <- factor(All_var$domain, 
                         levels = c("NTD", "TM", "EC", "IC", "H8", "CTD"))

Missense_var <- All_var %>% filter(var_type == "Missense")


# Output results for use in creating supplemental tables, analysis in 
# Prism, etc. 
write_csv(All_var, 
          "OutputPath/All_var.csv", 
          col_names = T)

write_csv(Missense_var, 
          "OutputPath/Missense_var.csv", 
          col_names = T)

write_csv(Results_by_pos, 
          "OutputPath/Results_by_pos.csv", 
          col_names = T)

########################
# Figures

# Fig. 1a created using Inkscape (https://gitlab.com/inkscape/inkscape)

#Fig. 1b - Heatmaps and average value lines; assembled further in Inkscape 
Expr_heatmap <- ggplot(All_var) +
  geom_tile(mapping = aes(x = pos, y = var_aa, fill = Expr_z_score)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red3", na.value="grey") +
  labs(x = "position", y = "substituted AA", fill = "Expression z-score") + 
  scale_x_continuous(breaks = seq(0, 359, by = 50)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) 
Expr_heatmap

Expr_z_score_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_expr_z_score)), size = 0.6) +
  theme_bw()
Expr_z_score_Line

Migr_heatmap <- ggplot(All_var) +
  geom_tile(mapping = aes(x = pos, y = var_aa, fill = Migr_z_score)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red3", na.value="grey") +
  labs(x = "position", y = "substituted AA", fill = "Migration z-score") + 
  scale_x_continuous(breaks = seq(0, 359, by = 50)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) 
Migr_heatmap

Migr_z_score_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_migr_z_score)), size = 0.6) +
  theme_bw()
Migr_z_score_Line

Migr_adj_heatmap <- ggplot(All_var) +
  geom_tile(mapping = aes(x = pos, y = var_aa, fill = E_adj_M_score)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red3", na.value="grey") +
  labs(x = "position", y = "substituted AA", fill = "Expression-adjusted migration score") + 
  scale_x_continuous(breaks = seq(0, 359, by = 50)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) 
Migr_adj_heatmap

EaM_score_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_EaM_score)), size = 0.6) +
  theme_bw()
EaM_score_Line

Prolif_heatmap <- ggplot(All_var) +
  geom_tile(mapping = aes(x = pos, y = var_aa, fill = Prolif_z_score)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red3", na.value="grey") +
  labs(x = "position", y = "substituted AA", fill = "Proliferation z-score") + 
  scale_x_continuous(breaks = seq(0, 359, by = 50)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) 
Prolif_heatmap

Prolif_z_score_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_prolif_z_score)), size = 0.6) +
  theme_bw()
Prolif_z_score_Line

Prolif_adj_heatmap <- ggplot(All_var) +
  geom_tile(mapping = aes(x = pos, y = var_aa, fill = E_adj_P_score)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red3", na.value="grey") +
  labs(x = "position", y = "substituted AA", fill = "Expression-adjusted proliferation score") + 
  scale_x_continuous(breaks = seq(0, 359, by = 50)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) 
Prolif_adj_heatmap

EaP_score_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_EaP_score)), size = 0.6) +
  theme_bw()
EaP_score_Line

#Fig. 1c: expression vs migration dotplot, with density along side, including lines at -2, 2 for x and y 
# (threshold between LoF, ~WT, GoF). Elements assembled in Inkscape. 
Expr_by_Migr_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Expr_z_score, y = Migr_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "Expression z-score", y = "Migration z-score")+
  ggtitle(label = "Expr & Migr by Variant")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))+
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, lineteype = "dashed")+
  geom_hline(yintercept = -2, linetype = "dashed")+
  geom_hline(yintercept = 2, linetype = "dashed")
Expr_by_Migr_plot

Expr_z_score_density <- ggplot(Missense_var) +
  geom_density(mapping = aes(x = Expr_z_score), fill="grey")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim((floor(min(Missense_var$Expr_z_score, na.rm = T)/0.5)*0.5), 
       (ceiling(max(Missense_var$Expr_z_score, na.rm = T)/0.5)*0.5))
Expr_z_score_density

Migr_z_score_density <- ggplot(Missense_var) +
  geom_density(mapping = aes(x = Migr_z_score), fill="grey")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim((floor(min(Missense_var$Migr_z_score, na.rm = T)/0.5)*0.5), 
       (ceiling(max(Missense_var$Migr_z_score, na.rm = T)/0.5)*0.5))
Migr_z_score_density

#Number of variants by phenotype 
nrow(Missense_var %>% filter(Expr_z_score >= 2 & Migr_z_score >=2))
nrow(Missense_var %>% filter(Expr_z_score >= 2 & Migr_z_score <2 & Migr_z_score > -2))
nrow(Missense_var %>% filter(Expr_z_score >= 2 & Migr_z_score <= -2))
nrow(Missense_var %>% filter(Expr_z_score <2 & Expr_z_score > -2 & Migr_z_score >=2))
nrow(Missense_var %>% filter(Expr_z_score <2 & Expr_z_score > -2 & Migr_z_score <2 & Migr_z_score > -2))
nrow(Missense_var %>% filter(Expr_z_score <2 & Expr_z_score > -2 & Migr_z_score <= -2))
nrow(Missense_var %>% filter(Expr_z_score <= -2 & Migr_z_score >=2))
nrow(Missense_var %>% filter(Expr_z_score <= -2 & Migr_z_score <2 & Migr_z_score > -2))
nrow(Missense_var %>% filter(Expr_z_score <= -2 & Migr_z_score <= -2))
#Results: 65,237,50,304,2942,77,2545,486,1

#Spearman correlation between expression z-score and migration z-score 
cor.test(Missense_var$Expr_z_score, Missense_var$Migr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$Expr_z_score and Missense_var$Migr_z_score
# S = 9.1443e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.8185082 

#Fig 1d: Expression by migration by position dot plot 
Expr_by_Migr_Pos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_expr_z_score, y = Mean_migr_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "Mean(expression z-score)", y = "Mean(migration z-score)")+
  ggtitle(label = "Expr & Migr by Position")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
Expr_by_Migr_Pos_plot

#Spearmant correlation between mean expression z-score and mean migration z-score 
cor.test(Results_by_pos$Mean_expr_z_score, Results_by_pos$Mean_migr_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_expr_z_score and Results_by_pos$Mean_migr_z_score
# S = 13470066, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.8373857 

#Fig. 1e created using Graphpad Prism

#Fig. 2a, ESM1b score by expression dot plot, some assembly in Inkscape
ESM1b_by_Expr_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = ESM1b_score, y = Expr_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "ESM1b score", y = "Expression z-score")+
  ggtitle(label = "ESM1b & Expr by Variant")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
ESM1b_by_Expr_plot

# Expression density as in Fig. 1c 

ESM1b_density <- ggplot(Missense_var) +
  geom_density(mapping = aes(x = ESM1b_score), fill="grey")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim((floor(min(Missense_var$ESM1b_score, na.rm = T)/0.5)*0.5), 
       (ceiling(max(Missense_var$ESM1b_score, na.rm = T)/0.5)*0.5))
ESM1b_density

# Spearman correlation between ESM1b score and expression z-score
cor.test(Missense_var$ESM1b_score, Missense_var$Expr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$ESM1b_score and Missense_var$Expr_z_score
# S = 2.0798e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5863984 

#Fig. 2b - ESM1b by Migration score dot plot, some assembly in Inkscape; 
# density curves as from Fig. 2a and Fig. 1c 
ESM1b_by_Migr_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = ESM1b_score, y = Migr_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "ESM1b score", y = "Migration z-score")+
  ggtitle(label = "ESM1b & Migr by Variant")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
ESM1b_by_Migr_plot

# Spearman correlation between ESM1b and migration z-score 
cor.test(Missense_var$ESM1b_score, Missense_var$Migr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$ESM1b_score and Missense_var$Migr_z_score
# S = 8.0088e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# -0.592705 

#Fig. 2c - ESM1b by proliferation z-score dot plot; some assembly in Inkscape;
#density plots as in Fig. 2a, Fig. S2a
ESM1b_by_Prolif_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = ESM1b_score, y = Prolif_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "ESM1b score", y = "Proliferation z-score")+
  ggtitle(label = "ESM1b & Prolif by Variant")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
ESM1b_by_Prolif_plot

# Spearman correlation ESM1b and proliferation z-score 
cor.test(Missense_var$ESM1b_score, Missense_var$Prolif_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$ESM1b_score and Missense_var$Prolif_z_score
# S = 7.5447e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.5004129 

#Fig. 2d - ESM1b by Expression by Position 
ESM1b_by_Expr_byPos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_ESM1b_score, y = Mean_expr_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "Mean(ESM1b score)", y = "Mean(Expression z-score)")+
  ggtitle(label = "ESM1b & Expr by Position")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
ESM1b_by_Expr_byPos_plot

#Spearman correlation mean ESM1b and mean expression z-score
cor.test(Results_by_pos$Mean_ESM1b_score, Results_by_pos$Mean_expr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_ESM1b_score and Results_by_pos$Mean_expr_z_score
# S = 2667386, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.636155 

#Fig. 2e - ESM1b by Migration by Position 
ESM1b_by_Migr_byPos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_ESM1b_score, y = Mean_migr_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "Mean(ESM1b score)", y = "Mean(Migration z-score)")+
  ggtitle(label = "ESM1b & Migr by Position")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
ESM1b_by_Migr_byPos_plot

# Spearman correlation mean ESM1b and mean migration z-score 
cor.test(Results_by_pos$Mean_ESM1b_score, Results_by_pos$Mean_migr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_ESM1b_score and Results_by_pos$Mean_migr_z_score
# S = 12073274, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.6468562 

#Fig. 2f - ESM1b by Proliferation by Position 
ESM1b_by_Prolif_byPos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_ESM1b_score, y = Mean_prolif_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "Mean(ESM1b score)", y = "Mean(Proliferation z-score)")+
  ggtitle(label = "ESM1b & Proliferation by Position")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
ESM1b_by_Prolif_byPos_plot

# Spearman correlation mean ESM1b and mean proliferation z-score 
cor.test(Results_by_pos$Mean_ESM1b_score, Results_by_pos$Mean_prolif_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_ESM1b_score and Results_by_pos$Mean_prolif_z_score
# S = 11710831, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.5974172 

#Fig. 2g-m, created using Prism (Graphpad)

#Fig 3a - Expression by Migration highlighting elevated expr-adjusted migration 

# Set threshold for elevated expression-adjusted migration: more than 2x
# the interquartile range above the median 
EaM_25percent = quantile(Results_by_pos$Mean_EaM_score, 0.25, na.rm = T) 
EaM_median = quantile(Results_by_pos$Mean_EaM_score, 0.5, na.rm = T) 
EaM_75percent = quantile(Results_by_pos$Mean_EaM_score, 0.75, na.rm = T) 
EaM_threshold = 2*((EaM_75percent - EaM_25percent)) + EaM_median
EaM_threshold
# 0.1530474 

Expr_by_Migr_byPos_highlighted_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_expr_z_score, y = Mean_migr_z_score,
                           color = Mean_EaM_score > EaM_threshold ))+
  theme_bw()+
  labs(x = "mean(expression z-score)", y = "mean(migration z-score)")+
  ggtitle(label = "Expr & Migr by Position")+
  scale_colour_manual(values = setNames(c('red',"grey"),c(T, F))) 
Expr_by_Migr_byPos_highlighted_plot

#Fig. 3b - Migration by Proliferation with highlighting 
# Define migration-adjusted proliferation scores 
lm_Synon_M_P <- lm(Synon_var$Prolif_score ~ Synon_var$Migr_score)
summary(lm_Synon_M_P)
# Call:
#   lm(formula = Synon_var$Prolif_score ~ Synon_var$Migr_score)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.29632 -0.07014 -0.01825  0.04332  0.86385 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -0.0005425  0.0074555  -0.073    0.942    
# Synon_var$Migr_score  0.4502332  0.0411701  10.936   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1371 on 336 degrees of freedom
# Multiple R-squared:  0.2625,	Adjusted R-squared:  0.2603 
# F-statistic: 119.6 on 1 and 336 DF,  p-value: < 2.2e-16

Missense_var <- Missense_var %>%
  mutate(M_adj_P_score = Eadj_function(Migr_score, Prolif_score, 
                                       coef(lm_Synon_M_P)[2],
                                       coef(lm_Synon_M_P)[1])) 

temp <- Missense_var %>%
  group_by(pos) %>%
  summarize(
    Mean_MaP_score = mean(M_adj_P_score)
  )

Results_by_pos <- left_join(Results_by_pos, temp, by = "pos")

# Set threshold of median + 2x interquartile range
MaP_25percent = quantile(Results_by_pos$Mean_MaP_score, 0.25, na.rm = T) 
MaP_median = quantile(Results_by_pos$Mean_MaP_score, 0.5, na.rm = T) 
MaP_75percent = quantile(Results_by_pos$Mean_MaP_score, 0.75, na.rm = T) 
MaP_threshold = 2*((MaP_75percent - MaP_25percent)) + MaP_median
MaP_threshold
#0.2099343

#Plot 
Migr_by_Prolif_byPos_highlighted_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_migr_z_score, y = Mean_prolif_z_score,
                           color = Mean_MaP_score > MaP_threshold ))+
  theme_bw()+
  labs(x = "mean(migration z-score)", y = "mean(proliferation z-score)")+
  ggtitle(label = "Migr & Prolif by Position")+
  scale_colour_manual(values = setNames(c('red',"grey"),c(T, F))) 
Migr_by_Prolif_byPos_highlighted_plot

# Fig. 3c - highlighting validated variants 
ValidatedVariants_list <- c("V29W", "A34F", "S44C", "T68K", "F103L", "Y104I",
                            "G124E", "K180I", "W181M", "L220Y", "G237N", "A240T", "T330N", "W181D", "W181S")

ValVar1_list <- c("F103L", "G124E", "L220Y", "W181D")
ValVar2_list <- c("A34F", "S44C", "G237N", "A240T")
ValVar3_list <- c("T68K", "Y104I", "K180I")
ValVar4_list <- c("V29W", "W181M", "W181S")
ValVar5_list <- c("T330N")

ValidatedVariants_EbyM <- ggplot(Missense_var) +
  geom_point(mapping = aes(x=Expr_z_score, y = Migr_z_score), color = "lightgrey", size = 0.6)+
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar1_list), mapping = aes(x = Expr_z_score, y = Migr_z_score),
             color = "black", fill = "yellow3", pch = 21, size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar2_list), mapping = aes(x = Expr_z_score, y = Migr_z_score),
             color = "turquoise4", size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar3_list), mapping = aes(x = Expr_z_score, y = Migr_z_score),
             color = "deeppink2", size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar4_list), mapping = aes(x = Expr_z_score, y = Migr_z_score),
             color = "purple3", size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar5_list), mapping = aes(x = Expr_z_score, y = Migr_z_score),
             color = "black", size=1.5) +
  theme_bw()+
  labs(x = "Expression z-score", y = "Migration z-score") 
ValidatedVariants_EbyM

ValidatedVariants_MbyP <- ggplot(Missense_var) +
  geom_point(mapping = aes(x=Migr_z_score, y = Prolif_z_score), color = "lightgrey", size = 0.6)+
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar1_list), mapping = aes(x = Migr_z_score, y = Prolif_z_score),
             color = "black", fill = "yellow3", pch = 21, size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar2_list), mapping = aes(x = Migr_z_score, y = Prolif_z_score),
             color = "turquoise4", size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar3_list), mapping = aes(x = Migr_z_score, y = Prolif_z_score),
             color = "deeppink2", size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar4_list), mapping = aes(x = Migr_z_score, y = Prolif_z_score),
             color = "purple3", size=1.5) +
  geom_point(Missense_var_z %>% filter(full_var %in% ValVar5_list), mapping = aes(x = Migr_z_score, y = Prolif_z_score),
             color = "black", size=1.5) +
  theme_bw()+
  labs(x = "Migration z-score", y = "Proliferation z-score") 
ValidatedVariants_MbyPAlt

#Fig. 3d-i made in Prism (Graphpad)

#Fig. 4a-h made using ChimeraX 

#Fig. 4i - Expr by migration by position with Ga-interface residues highlighted
# Contacts determined from cryo-EM structure 
Ga13_contacts <- c(52, 58, 59, 62, 120, 121, 124, 125, 126, 127, 128, 
                   129, 130, 132, 133, 135, 136, 138, 213, 216, 219, 
                   220, 223, 225, 227, 228, 229, 231, 232, 234, 235,
                   238, 239, 297, 298, 299, 300)
  
Expr_by_Migr_byPos_Ga13HighL_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_expr_z_score, y = Mean_migr_z_score),
             color = "grey")+
  theme_bw()+
  labs(x = "mean(expression z-score)", y = "mean(migration z-score)")+
  ggtitle(label = "Expr & Migr by Position")+
  geom_point(Results_by_pos %>% filter(pos %in% Ga13_contacts),
             mapping = aes(x = Mean_expr_z_score, y = Mean_migr_z_score),
             color = "red")
Expr_by_Migr_byPos_Ga13HighL_plot

#Fig. 5a,b made using ChimeraX 

#Fig. 5c - Expression by migration by position, ligand contacts highlighted 
# GGG contacts determined using cryo-EM structure
GGG_contacts = c(8, 9, 10, 104, 108, 155, 158, 165, 176, 177, 178, 179, 180, 
                 181, 190, 193, 194, 197, 198, 201, 202, 257, 260, 261, 264, 
                 265, 272, 275)

Expr_by_Migr_byPos_GGGHighL_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_expr_z_score, y = Mean_migr_z_score),
             color = "grey")+
  theme_bw()+
  labs(x = "mean(expression z-score)", y = "mean(migration z-score)")+
  ggtitle(label = "Expr & Migr by Position")+
  geom_point(Results_by_pos %>% filter(pos %in% GGG_contacts),
             mapping = aes(x = Mean_expr_z_score, y = Mean_migr_z_score),
             color = "red")
Expr_by_Migr_byPos_GGGHighL_plot

#Fig. 5d,e made in Prism (Graphpad)

#Fig. 6a. Figure made in Prism (Graphpad). 
# Here, non-linear regression and statistical analysis
For_Migr_ANOVA <- read_csv("~/Cyster-Ye/Comp Bio/P2RY8_SatMut/C042_PSM_2/For_Migr_ANCOVA.csv", 
         col_names = T)

logistic_model_fixed_A <- function(x, B, C) {
  A <- 1  # Fixed value for A
  A / (1 + exp(-B * (x - C)))
}

WT_mig_nls_model <- nls(Migration_index ~ logistic_model_fixed_A(log_GGG_Conc, B, C), 
                        start = list(B = -1, C = -8),
                        data = (For_Migr_ANOVA %>% filter(Genotype == "WT")))
summary(WT_mig_nls_model)
# Parameters:
#   Estimate Std. Error  t value Pr(>|t|)    
# B -4.69256    1.97381   -2.377   0.0211 *  
#   C -8.11642    0.06395 -126.913   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.28 on 52 degrees of freedom
# 
# Number of iterations to convergence: 8 
# Achieved convergence tolerance: 1.487e-06

V29W_mig_nls_model <- nls(Migration_index ~ logistic_model_fixed_A(log_GGG_Conc, B, C), 
                           start = list(B = -1, C = -8),
                           data = (For_Migr_ANOVA %>% filter(Genotype == "V29W")))
summary(V29W_mig_nls_model)
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# # B  -3.0197     1.7060   -1.77   0.0869 .  
# C  -7.8897     0.1516  -52.04   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3799 on 30 degrees of freedom
# 
# Number of iterations to convergence: 7 
# Achieved convergence tolerance: 5.112e-06

K180I_mig_nls_model <- nls(Migration_index ~ logistic_model_fixed_A(log_GGG_Conc, B, C), 
                        start = list(B = -1, C = -8),
                        data = (For_Migr_ANOVA %>% filter(Genotype == "K180I")))
summary(K180I_mig_nls_model)
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# B  -1.9635     0.5797  -3.387  0.00211 ** 
#   C  -7.5292     0.1672 -45.018  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3016 on 28 degrees of freedom
# 
# Number of iterations to convergence: 7 
# Achieved convergence tolerance: 7.716e-06

W181M_mig_nls_model <- nls(Migration_index ~ logistic_model_fixed_A(log_GGG_Conc, B, C), 
                           start = list(B = -1, C = -8),
                           data = (For_Migr_ANOVA %>% filter(Genotype == "W181M")))
summary(W181M_mig_nls_model)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# B   -1.796      1.234  -1.455    0.156    
# C   -8.222      0.296 -27.778   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.473 on 30 degrees of freedom
# 
# Number of iterations to convergence: 8 
# Achieved convergence tolerance: 5.466e-06

IC50s <- data.frame(
  genotype = factor(c("WT", "V29W", "K180I", "W181M")),
  estimate = c(summary(WT_mig_nls_model)$parameters["C", "Estimate"],
  summary(V29W_mig_nls_model)$parameters["C", "Estimate"],
summary(K180I_mig_nls_model)$parameters["C", "Estimate"],
summary(W181M_mig_nls_model)$parameters["C", "Estimate"]),
SE = c(summary(WT_mig_nls_model)$parameters["C", "Std. Error"],
summary(V29W_mig_nls_model)$parameters["C", "Std. Error"],
summary(K180I_mig_nls_model)$parameters["C", "Std. Error"],
summary(W181M_mig_nls_model)$parameters["C", "Std. Error"])
)

#pairwise z-tests against WT: 
Z_WT_V29W <- ((IC50s$estimate[IC50s$genotype == "V29W"])-(IC50s$estimate[IC50s$genotype == "WT"]))/
  (sqrt((IC50s$SE[IC50s$genotype == "V29W"])^2 + (IC50s$SE[IC50s$genotype == "WT"])^2))
#[1] 1.378202
Z_WT_K180I <- ((IC50s$estimate[IC50s$genotype == "K180I"])-(IC50s$estimate[IC50s$genotype == "WT"]))/
  (sqrt((IC50s$SE[IC50s$genotype == "K180I"])^2 + (IC50s$SE[IC50s$genotype == "WT"])^2))
# [1] 3.279412
Z_WT_W181M <- ((IC50s$estimate[IC50s$genotype == "W181M"])-(IC50s$estimate[IC50s$genotype == "WT"]))/
  (sqrt((IC50s$SE[IC50s$genotype == "W181M"])^2 + (IC50s$SE[IC50s$genotype == "WT"])^2))
# [1] -0.3479934

#Determine p-value, then adjust for multiple testing by BH method 
pvals <- c(2*pnorm(-abs(Z_WT_V29W)), 2*pnorm(-abs(Z_WT_K180I)), 2*pnorm(-abs(Z_WT_W181M)))
# [1] 0.168140878 0.001040237 0.727845123

pvals_adj <- p.adjust(pvals, method = "BH")
pvals_adj
# [1] 0.252211318 0.003120712 0.727845123

#Fig. 6b, Figure made in Prism (Graphpad)
# Here, non-linear regression and statistical analysis
For_Prolif_ANOVA <- read_csv("~/Cyster-Ye/Comp Bio/P2RY8_SatMut/C042_PSM_2/For_Prolif_ANCOVA.csv", 
                            col_names = T)

WT_prol_nls_model <- nls(Proliferation_index ~ SSlogis(log_GGG_Conc, Asym, xmid, scal),
                         data = (For_Prolif_ANOVA %>% filter(Genotype == "WT")))

summary(WT_prol_nls_model)
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# Asym  0.64383    0.02992  21.518  < 2e-16 ***
#   xmid -6.27654    0.17284 -36.315  < 2e-16 ***
#   scal -1.27101    0.25114  -5.061 6.27e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07587 on 49 degrees of freedom
# 
# Number of iterations to convergence: 5 
# Achieved convergence tolerance: 2.338e-06

V29W_prol_nls_model <- nls(Proliferation_index ~ SSlogis(log_GGG_Conc, Asym, xmid, scal),
                         data = (For_Prolif_ANOVA %>% filter(Genotype == "V29W")))

summary(V29W_prol_nls_model)
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# Asym  0.91619    0.03055  29.989  < 2e-16 ***
#   xmid -6.38893    0.10933 -58.438  < 2e-16 ***
#   scal -0.74164    0.12473  -5.946 2.84e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1204 on 49 degrees of freedom
# 
# Number of iterations to convergence: 5 
# Achieved convergence tolerance: 3.624e-06

K180I_prol_nls_model <- nls(Proliferation_index ~ SSlogis(log_GGG_Conc, Asym, xmid, scal),
                           data = (For_Prolif_ANOVA %>% filter(Genotype == "K180I")))

summary(K180I_prol_nls_model)
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# Asym  0.79776    0.03705  21.534   <2e-16 ***
#   xmid -5.78129    0.21502 -26.888   <2e-16 ***
#   scal -0.70279    0.28407  -2.474   0.0199 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1193 on 27 degrees of freedom
# 
# Number of iterations to convergence: 2 
# Achieved convergence tolerance: 9.137e-07

W181M_prol_nls_model <- nls(Proliferation_index ~ SSlogis(log_GGG_Conc, Asym, xmid, scal),
                            data = (For_Prolif_ANOVA %>% filter(Genotype == "W181M")))

summary(W181M_prol_nls_model)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# Asym  1.12267    0.06415  17.500 2.92e-16 ***
#   xmid -7.01963    0.17635 -39.805  < 2e-16 ***
#   scal -0.75400    0.16136  -4.673 7.34e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1666 on 27 degrees of freedom
# 
# Number of iterations to convergence: 0 
# Achieved convergence tolerance: 1.939e-06

IC50s_prol <- data.frame(
  genotype = factor(c("WT", "V29W", "K180I", "W181M")),
  estimate = c(summary(WT_prol_nls_model)$parameters["xmid", "Estimate"],
               summary(V29W_prol_nls_model)$parameters["xmid", "Estimate"],
               summary(K180I_prol_nls_model)$parameters["xmid", "Estimate"],
               summary(W181M_prol_nls_model)$parameters["xmid", "Estimate"]),
  SE = c(summary(WT_prol_nls_model)$parameters["xmid", "Std. Error"],
         summary(V29W_prol_nls_model)$parameters["xmid", "Std. Error"],
         summary(K180I_prol_nls_model)$parameters["xmid", "Std. Error"],
         summary(W181M_prol_nls_model)$parameters["xmid", "Std. Error"]),
  Asym_estimate = c(summary(WT_prol_nls_model)$parameters["Asym", "Estimate"],
               summary(V29W_prol_nls_model)$parameters["Asym", "Estimate"],
               summary(K180I_prol_nls_model)$parameters["Asym", "Estimate"],
               summary(W181M_prol_nls_model)$parameters["Asym", "Estimate"]),
  Asym_SE = c(summary(WT_prol_nls_model)$parameters["Asym", "Std. Error"],
         summary(V29W_prol_nls_model)$parameters["Asym", "Std. Error"],
         summary(K180I_prol_nls_model)$parameters["Asym", "Std. Error"],
         summary(W181M_prol_nls_model)$parameters["Asym", "Std. Error"])
)

#pairwise z-tests against WT, first IC50 then asymptote [Emax]: 
Z_WT_V29W_prol <- ((IC50s_prol$estimate[IC50s_prol$genotype == "V29W"])-(IC50s_prol$estimate[IC50s_prol$genotype == "WT"]))/
  (sqrt((IC50s_prol$SE[IC50s_prol$genotype == "V29W"])^2 + (IC50s_prol$SE[IC50s_prol$genotype == "WT"])^2))
#[1] -0.5495354
Z_WT_K180I_prol <- ((IC50s_prol$estimate[IC50s_prol$genotype == "K180I"])-(IC50s_prol$estimate[IC50s_prol$genotype == "WT"]))/
  (sqrt((IC50s_prol$SE[IC50s_prol$genotype == "K180I"])^2 + (IC50s_prol$SE[IC50s_prol$genotype == "WT"])^2))
# [1] 1.795213
Z_WT_W181M_prol <- ((IC50s_prol$estimate[IC50s_prol$genotype == "W181M"])-(IC50s_prol$estimate[IC50s_prol$genotype == "WT"]))/
  (sqrt((IC50s_prol$SE[IC50s_prol$genotype == "W181M"])^2 + (IC50s_prol$SE[IC50s_prol$genotype == "WT"])^2))
# [1] -3.00939
Z_WT_V29W_prolA <- ((IC50s_prol$Asym_estimate[IC50s_prol$genotype == "V29W"])-(IC50s_prol$Asym_estimate[IC50s_prol$genotype == "WT"]))/
  (sqrt((IC50s_prol$Asym_SE[IC50s_prol$genotype == "V29W"])^2 + (IC50s_prol$Asym_SE[IC50s_prol$genotype == "WT"])^2))
# [1] 6.369288
Z_WT_K180I_prolA <- ((IC50s_prol$Asym_estimate[IC50s_prol$genotype == "K180I"])-(IC50s_prol$Asym_estimate[IC50s_prol$genotype == "WT"]))/
  (sqrt((IC50s_prol$Asym_SE[IC50s_prol$genotype == "K180I"])^2 + (IC50s_prol$Asym_SE[IC50s_prol$genotype == "WT"])^2))
# [1] 3.232425
Z_WT_W181M_prolA <- ((IC50s_prol$Asym_estimate[IC50s_prol$genotype == "W181M"])-(IC50s_prol$Asym_estimate[IC50s_prol$genotype == "WT"]))/
  (sqrt((IC50s_prol$Asym_SE[IC50s_prol$genotype == "W181M"])^2 + (IC50s_prol$Asym_SE[IC50s_prol$genotype == "WT"])^2))
# [1] 6.764629

#Determine p-value, then adjust for multiple testing by BH method 
pvals_prol <- c(2*pnorm(-abs(Z_WT_V29W_prol)), 2*pnorm(-abs(Z_WT_K180I_prol)), 2*pnorm(-abs(Z_WT_W181M_prol)),
                2*pnorm(-abs(Z_WT_V29W_prolA)), 2*pnorm(-abs(Z_WT_K180I_prolA)), 2*pnorm(-abs(Z_WT_W181M_prolA)))
# [1] 5.826381e-01 7.261969e-02 2.617729e-03 1.899073e-10 1.227444e-03 1.336509e-11

pvals_adj_prol <- p.adjust(pvals_prol, method = "BH")
pvals_adj_prol
# [1] 5.826381e-01 8.714363e-02 3.926593e-03 5.697219e-10 2.454888e-03 8.019052e-11

#Fig. 7a - gnomAD variant DMS results highlighted 

#First upload gnomAD (v4.0.0) data in CSV file at FilePath8. 
gnomad_missense <- read_csv(FilePath8, col_names = T)

#Modify variant format to align with the A10G format used throughout
gnomad_missense$full_var <- str_sub(gnomad_missense$full_var, 3, -1)

aa_short <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", 
              "Q", "R", "S", "T", "V", "W", "Y")
aa_long <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
             "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")

for (i in c(1:20)) {
  gnomad_missense$full_var <- str_replace_all(gnomad_missense$full_var, 
                                              aa_long[i], aa_short[i])
}

#There are some duplicate missense variants in the original CSV
gnomad_missense_nodup <- gnomad_missense %>%
  filter(!duplicated(full_var))

nrow(gnomad_missense_nodup)
#[1] 476

#Bring allele frequencies into missense_var
Missense_var <- left_join(Missense_var, 
                          select(gnomad_missense_nodup, full_var, allele_freq),
                          by = "full_var")

#Plot all variants with those in gnomAD highlighted by domain
Expr_by_Migr_gnomAD_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Expr_z_score, y = Migr_z_score
                           ), size = 0.6, color = "grey")+
  theme_bw()+
  labs(x = "Expression z-score", y = "Migration z-score")+
  ggtitle(label = "Expression & Migr by Variant, gnomAD")+
  geom_point(Missense_var %>% filter(!is.na(Missense_var$allele_freq)),
             mapping = aes(x = Expr_z_score, y = Migr_z_score, 
                           color = domain), size = 0.6)+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
Expr_by_Migr_gnomAD_plot

Migr_by_Prolif_gnomAD_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Migr_z_score, y = Prolif_z_score
  ), size = 0.6, color = "grey")+
  theme_bw()+
  labs(x = "Migration z-score", y = "Proliferation z-score")+
  ggtitle(label = "Migration & Proliferation by Variant, gnomAD")+
  geom_point(Missense_var %>% filter(!is.na(Missense_var$allele_freq)),
             mapping = aes(x = Migr_z_score, y = Prolif_z_score, 
                           color = domain), size = 0.6)+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
Migr_by_Prolif_gnomAD_plot

#Fig. 7b, made in Prism (Graphpad)

#Fig. 7c - DLBCL and Burkitt lymphoma variants DMS phenotype plots
#Upload CSV file of 48 variant names; 32 from COSMIC and additional 16 from Muppidi et al, Nature, 2014.
All_var_D_B_nodup <- read_csv(FilePath10, col_names = T)

#Pull in expression data
All_var_D_B_nodup <- left_join(All_var_D_B_nodup, select(Missense_var, full_var, Expr_z_score, 
                                               Migr_z_score, Prolif_z_score))

Expr_by_Migr_AllDB_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Expr_z_score, y = Migr_z_score
  ), size = 0.6, color = "grey")+
  theme_bw()+
  labs(x = "Expression z-score", y = "Migration z-score")+
  geom_point(All_var_D_B_nodup, mapping = aes(x = Expr_z_score, y = Migr_z_score),
             color = "red", size = 0.9)
Expr_by_Migr_AllDB_plot

Migr_by_Prolif_AllDB_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Migr_z_score, y = Prolif_z_score
  ), size = 0.6, color = "grey")+
  theme_bw()+
  labs(x = "Migration z-score", y = "Proliferation z-score")+
  geom_point(All_var_D_B_nodup, mapping = aes(x = Migr_z_score, y = Prolif_z_score),
             color = "red", size = 0.9)
Migr_by_Prolif_AllDB_plot

#Fig. 7d - non-hematologic COSMIC missense variants DMS phenotype plot
#Upload CSV file of 159 variant names
Cosmic_nonHL <- read_csv(FilePath9, col_names = T)

#Pull in expression data
Cosmic_nonHL <- left_join(Cosmic_nonHL, select(Missense_var, full_var, Expr_z_score, 
                                              Migr_z_score, Prolif_z_score))

Expr_by_Migr_CosmicnonHL_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Expr_z_score, y = Migr_z_score
  ), size = 0.6, color = "grey")+
  theme_bw()+
  labs(x = "Expression", y = "Migration z-score")+
  geom_point(Cosmic_nonHL, mapping = aes(x = Expr_z_score, y = Migr_z_score),
                          color = "blue", size = 0.9)
Expr_by_Migr_CosmicnonHL_plot

Migr_by_Prolif_CosmicnonHL_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Migr_z_score, y = Prolif_z_score
  ), size = 0.6, color = "grey")+
  theme_bw()+
  labs(x = "Migration z-score", y = "Proliferation z-score")+
  geom_point(Cosmic_nonHL, mapping = aes(x = Migr_z_score, y = Prolif_z_score),
             color = "blue", size = 0.9)
Migr_by_Prolif_CosmicnonHL_plot

#Fig. 7e - CDFs across expression, migration, proliferation

Expr_CDF <- ggplot(Missense_var, aes(Expr_z_score))+
                     stat_ecdf(size = 0.7, geom = "step")+
       stat_ecdf(Missense_var %>% filter(!is.na(Missense_var$allele_freq)), 
                                    mapping = aes(Expr_z_score), color = "green4", size = 0.7, 
            geom = "step")+
  stat_ecdf(Cosmic_nonHL , 
            mapping = aes(Expr_z_score), color = "blue", size = 0.7,
            geom = "step")+
  stat_ecdf(All_var_D_B_nodup,
            mapping = aes(Expr_z_score), color = "red", size = 0.7,
            geom = "step")+
  labs(x = "Expression z-score", y = "Cumulative density")+
  theme_bw()
Expr_CDF

Migr_CDF <- ggplot(Missense_var, aes(Migr_z_score))+
  stat_ecdf(size = 0.7, geom = "step")+
  stat_ecdf(Missense_var %>% filter(!is.na(Missense_var$allele_freq)), 
            mapping = aes(Migr_z_score), color = "green4", size = 0.7, 
            geom = "step")+
  stat_ecdf(Cosmic_nonHL , 
            mapping = aes(Migr_z_score), color = "blue", size = 0.7,
            geom = "step")+
  stat_ecdf(All_var_D_B_nodup,
            mapping = aes(Migr_z_score), color = "red", size = 0.7,
            geom = "step")+
  labs(x = "Migration z-score", y = "Cumulative density")+
  theme_bw()
Migr_CDF

Prolif_CDF <- ggplot(Missense_var, aes(Prolif_z_score))+
  stat_ecdf(size = 0.7, geom = "step")+
  stat_ecdf(Missense_var %>% filter(!is.na(Missense_var$allele_freq)), 
            mapping = aes(Prolif_z_score), color = "green4", size = 0.7, 
            geom = "step")+
  stat_ecdf(Cosmic_nonHL , 
            mapping = aes(Prolif_z_score), color = "blue", size = 0.7,
            geom = "step")+
  stat_ecdf(All_var_D_B_nodup,
            mapping = aes(Prolif_z_score), color = "red", size = 0.7,
            geom = "step")+
  labs(x = "Proliferation z-score", y = "Cumulative density")+
  theme_bw()
Prolif_CDF

Cosmic_nonHL <- left_join(Cosmic_nonHL, select(Missense_var, full_var, E_adj_M_score))
All_var_D_B_nodup <- left_join(All_var_D_B_nodup, select(Missense_var, full_var, E_adj_M_score))

EaM_CDF <- ggplot(Missense_var, aes(E_adj_M_score))+
  stat_ecdf(size = 0.7, geom = "step")+
  stat_ecdf(Missense_var %>% filter(!is.na(Missense_var$allele_freq)), 
            mapping = aes(E_adj_M_score), color = "green4", size = 0.7, 
            geom = "step")+
  stat_ecdf(Cosmic_nonHL , 
            mapping = aes(E_adj_M_score), color = "blue", size = 0.7,
            geom = "step")+
  stat_ecdf(All_var_D_B_nodup,
            mapping = aes(E_adj_M_score), color = "red", size = 0.7,
            geom = "step")+
  labs(x = "Expr-adjusted migration score", y = "Cumulative density")+
  theme_bw()
EaM_CDF

#####################
# Supplemental Figures
#Fig. S1a, from FACSDiva 

#Figure S1b-d: Replicate comparisons. Filepath6 refers to a CSV file assembled in Microsoft Excel
#consisting of full_var and the 12 sets of individual Enrich2 scores, 4 each for expression,
#migration, and proliferation.  Further assembly in Inkscape
Replicate_Scores <- read_csv(FilePath6, col_names = T)

C36Expr_C37Expr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C36_Expr, y = C37_Expr)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 2") 
C36Expr_C37Expr_plot

C36Expr_C37Expr_CorP <- cor(Replicate_Scores$C36_Expr, Replicate_Scores$C37_Expr, 
                            use = "complete.obs", method = "pearson")

C36Expr_C39Expr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C36_Expr, y = C39_Expr)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 3") 
C36Expr_C39Expr_plot

C36Expr_C39Expr_CorP <- cor(Replicate_Scores$C36_Expr, Replicate_Scores$C39_Expr, 
                            use = "complete.obs", method = "pearson")

C36Expr_C41Expr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C36_Expr, y = C41_Expr)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 4") 
C36Expr_C41Expr_plot

C36Expr_C41Expr_CorP <- cor(Replicate_Scores$C36_Expr, Replicate_Scores$C41_Expr, 
                            use = "complete.obs", method = "pearson")

C37Expr_C39Expr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C37_Expr, y = C39_Expr)) +
  theme_bw()+
  labs(x = "Rep.2", y = "Rep. 3") 
C37Expr_C39Expr_plot

C37Expr_C39Expr_CorP <- cor(Replicate_Scores$C37_Expr, Replicate_Scores$C39_Expr, 
                            use = "complete.obs", method = "pearson")

C37Expr_C41Expr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C37_Expr, y = C41_Expr)) +
  theme_bw()+
  labs(x = "Rep. 2", y = "Rep. 4") 
C37Expr_C41Expr_plot

C37Expr_C41Expr_CorP <- cor(Replicate_Scores$C37_Expr, Replicate_Scores$C41_Expr, 
                            use = "complete.obs", method = "pearson")

C39Expr_C41Expr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C39_Expr, y = C41_Expr)) +
  theme_bw()+
  labs(x = "Rep. 3", y = "Rep. 4") 
C39Expr_C41Expr_plot

C39Expr_C41Expr_CorP <- cor(Replicate_Scores$C39_Expr, Replicate_Scores$C41_Expr, 
                            use = "complete.obs", method = "pearson")

C35Migr_C37Migr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C35_Migr, y = C37_Migr)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 2") 
C35Migr_C37Migr_plot

C35Migr_C37Migr_CorP <- cor(Replicate_Scores$C35_Migr, Replicate_Scores$C37_Migr, 
                            use = "complete.obs", method = "pearson")

C35Migr_C39Migr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C35_Migr, y = C39_Migr)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 3") 
C35Migr_C39Migr_plot

C35Migr_C39Migr_CorP <- cor(Replicate_Scores$C35_Migr, Replicate_Scores$C39_Migr, 
                            use = "complete.obs", method = "pearson")

C35Migr_C40Migr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C35_Migr, y = C40_Migr)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 4") 
C35Migr_C40Migr_plot

C35Migr_C40Migr_CorP <- cor(Replicate_Scores$C35_Migr, Replicate_Scores$C40_Migr, 
                            use = "complete.obs", method = "pearson")

C37Migr_C39Migr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C37_Migr, y = C39_Migr)) +
  theme_bw()+
  labs(x = "Rep. 2", y = "Rep. 3") 
C37Migr_C39Migr_plot

C37Migr_C39Migr_CorP <- cor(Replicate_Scores$C37_Migr, Replicate_Scores$C39_Migr, 
                            use = "complete.obs", method = "pearson")

C37Migr_C40Migr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C37_Migr, y = C40_Migr)) +
  theme_bw()+
  labs(x = "Rep. 2", y = "Rep. 4") 
C37Migr_C40Migr_plot

C37Migr_C40Migr_CorP <- cor(Replicate_Scores$C37_Migr, Replicate_Scores$C40_Migr, 
                            use = "complete.obs", method = "pearson")

C39Migr_C40Migr_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C39_Migr, y = C40_Migr)) +
  theme_bw()+
  labs(x = "Rep. 3", y = "Rep. 4") 
C39Migr_C40Migr_plot

C39Migr_C40Migr_CorP <- cor(Replicate_Scores$C39_Migr, Replicate_Scores$C40_Migr, 
                            use = "complete.obs", method = "pearson")

C34Prol_C35Prol_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C34_Prolif, y = C35_Prolif)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 2") 
C34Prol_C35Prol_plot

C34Prol_C35Prol_CorP <- cor(Replicate_Scores$C34_Prolif, Replicate_Scores$C35_Prolif, 
                            use = "complete.obs", method = "pearson")

C34Prol_C39Prol_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C34_Prolif, y = C39_Prolif)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 3") 
C34Prol_C39Prol_plot

C34Prol_C39Prol_CorP <- cor(Replicate_Scores$C34_Prolif, Replicate_Scores$C39_Prolif, 
                            use = "complete.obs", method = "pearson")

C34Prol_C40Prol_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C34_Prolif, y = C40_Prolif)) +
  theme_bw()+
  labs(x = "Rep. 1", y = "Rep. 4") 
C34Prol_C40Prol_plot

C34Prol_C40Prol_CorP <- cor(Replicate_Scores$C34_Prolif, Replicate_Scores$C40_Prolif, 
                            use = "complete.obs", method = "pearson")

C35Prol_C39Prol_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C35_Prolif, y = C39_Prolif)) +
  theme_bw()+
  labs(x = "Rep. 2", y = "Rep. 3") 
C35Prol_C39Prol_plot

C35Prol_C39Prol_CorP <- cor(Replicate_Scores$C35_Prolif, Replicate_Scores$C39_Prolif, 
                            use = "complete.obs", method = "pearson")

C35Prol_C40Prol_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C35_Prolif, y = C40_Prolif)) +
  theme_bw()+
  labs(x = "Rep. 2", y = "Rep. 4") 
C35Prol_C40Prol_plot

C35Prol_C40Prol_CorP <- cor(Replicate_Scores$C35_Prolif, Replicate_Scores$C40_Prolif, 
                            use = "complete.obs", method = "pearson")

C39Prol_C40Prol_plot <- ggplot(Replicate_Scores) +
  geom_point(mapping = aes(x = C39_Prolif, y = C40_Prolif)) +
  theme_bw()+
  labs(x = "Rep. 3", y = "Rep. 4") 
C39Prol_C40Prol_plot

C39Prol_C40Prol_CorP <- cor(Replicate_Scores$C39_Prolif, Replicate_Scores$C40_Prolif, 
                            use = "complete.obs", method = "pearson")

#Fig. S2a, expression by proliferation dot plot; as well as number of variants
# by phenotype, further assembled in Inkscape
Expr_by_Prolif_plot <- ggplot(Missense_var) +
  geom_point(mapping = aes(x = Expr_z_score, y = Prolif_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "Expression z-score", y = "Proliferation z-score")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))+
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, lineteype = "dashed")+
  geom_hline(yintercept = -2, linetype = "dashed")+
  geom_hline(yintercept = 2, linetype = "dashed")
Expr_by_Prolif_plot

#Expression score density re-used from Fig. 1c

Prolif_score_density <- ggplot(Missense_var) +
  geom_density(mapping = aes(x = Prolif_z_score), fill="grey")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim((floor(min(Missense_var$Prolif_z_score, na.rm = T)/0.5)*0.5), 
       (ceiling(max(Missense_var$Prolif_z_score, na.rm = T)/0.5)*0.5))
Prolif_score_density

#Number of variants by phenotype 
nrow(Missense_var %>% filter(Expr_z_score >= 2 & Prolif_z_score >=2))
nrow(Missense_var %>% filter(Expr_z_score >= 2 & Prolif_z_score <2 & Prolif_z_score > -2))
nrow(Missense_var %>% filter(Expr_z_score >= 2 & Prolif_z_score <= -2))
nrow(Missense_var %>% filter(Expr_z_score <2 & Expr_z_score > -2 & Prolif_z_score >=2))
nrow(Missense_var %>% filter(Expr_z_score <2 & Expr_z_score > -2 & Prolif_z_score <2 & Prolif_z_score > -2))
nrow(Missense_var %>% filter(Expr_z_score <2 & Expr_z_score > -2 & Prolif_z_score <= -2))
nrow(Missense_var %>% filter(Expr_z_score <= -2 & Prolif_z_score >=2))
nrow(Missense_var %>% filter(Expr_z_score <= -2 & Prolif_z_score <2 & Prolif_z_score > -2))
nrow(Missense_var %>% filter(Expr_z_score <= -2 & Prolif_z_score <= -2))
#Results: 113,221,18,334,2926,63,1438,1593,1

# Spearman correlation between expression z-score and proliferation z-score
cor.test(Missense_var$Expr_z_score, Missense_var$Prolif_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$Expr_z_score and Missense_var$Prolif_z_score
# S = 8.0384e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.5985789 

#Fig S2b: migration by proliferation by dot plot, further assembled in Inkscape
Migr_by_Prolif_plot <- ggplot(Missense_var) +
  geom_point(mapping = aes(x = Migr_z_score, y = Prolif_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "Migration z-score", y = "Proliferation z-score")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))+
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, lineteype = "dashed")+
  geom_hline(yintercept = -2, linetype = "dashed")+
  geom_hline(yintercept = 2, linetype = "dashed")
Migr_by_Prolif_plot

# Number of variants by phenotype 
nrow(Missense_var %>% filter(Migr_z_score >= 2 & Prolif_z_score >=2))
nrow(Missense_var %>% filter(Migr_z_score >= 2 & Prolif_z_score <2 & Prolif_z_score > -2))
nrow(Missense_var %>% filter(Migr_z_score >= 2 & Prolif_z_score <= -2))
nrow(Missense_var %>% filter(Migr_z_score <2 & Migr_z_score > -2 & Prolif_z_score >=2))
nrow(Missense_var %>% filter(Migr_z_score <2 & Migr_z_score > -2 & Prolif_z_score <2 & Prolif_z_score > -2))
nrow(Missense_var %>% filter(Migr_z_score <2 & Migr_z_score > -2 & Prolif_z_score <= -2))
nrow(Missense_var %>% filter(Migr_z_score <= -2 & Prolif_z_score >=2))
nrow(Missense_var %>% filter(Migr_z_score <= -2 & Prolif_z_score <2 & Prolif_z_score > -2))
nrow(Missense_var %>% filter(Migr_z_score <= -2 & Prolif_z_score <= -2))
# Results: 1513,1401,0,367,3223,75,5,116,7

# Spearman correlation between migration z-score and proliferation z-score 
cor.test(Missense_var$Migr_z_score, Missense_var$Prolif_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$Migr_z_score and Missense_var$Prolif_z_score
# S = 1.4302e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7155697 

#Fig. S2c: Expression z-score by expression-adjusted migration score dot plot
Expr_by_EaM_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = Expr_z_score, y = E_adj_M_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "Expression z-score", y = "Expression-adjusted migration score")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
Expr_by_EaM_plot

#Fig. S2d, expression z-score by expression-adjusted proliferation score dot plot
Expr_by_EaP_plot <- ggplot(Missense_var) +
  geom_point(mapping = aes(x = Expr_z_score, y = E_adj_P_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "Expression z-score", y = "Expression-adjusted proliferation score")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
Expr_by_EaP_plot

#Fig. S2e, expression-adjusted migration score by expression-adjusted proliferation
EaM_by_EaP_plot <- ggplot(Missense_var) +
  geom_point(mapping = aes(x = E_adj_M_score, y = E_adj_P_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "Expression-adjusted migration score", y = "Expression-adjusted proliferation score")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
EaM_by_EaP_plot

# Pearson correlation between expression-adjusted migration and expression-
# adjusted proliferation scores 
cor.test(Missense_var$E_adj_M_score, Missense_var$E_adj_P_score, 
         use = "complete.obs", method = "pearson")
# data:  Missense_var$E_adj_M_score and Missense_var$E_adj_P_score
# t = 68.259, df = 6705, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.6259657 0.6542138
# sample estimates:
#   cor 
# 0.6403062 

#Fig S2f: Expression by Proliferation by position dot plot
Expr_by_Prolif_Pos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_expr_z_score, y = Mean_prolif_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "Mean expression z-score", y = "Mean proliferation z-score")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
Expr_by_Prolif_Pos_plot

#Spearman correlation between mean expression z-score and mean proliferation z-score
cor.test(Results_by_pos$Mean_expr_z_score, Results_by_pos$Mean_prolif_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_expr_z_score and Results_by_pos$Mean_prolif_z_score
# S = 12279070, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.6749278 

#Fig S2g: Migration by Proliferation by position dot plot
Migr_by_Prolif_Pos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_migr_z_score, y = Mean_prolif_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "Mean migration z-score", y = "Mean proliferation z-score")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
Migr_by_Prolif_Pos_plot

# Spearman correlation between mean migration z-score and mean proliferation z-score
cor.test(Results_by_pos$Mean_migr_z_score, Results_by_pos$Mean_prolif_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_migr_z_score and Results_by_pos$Mean_prolif_z_score
# S = 1083598, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8521917 

#Fig. S2h,i: created in Prism (Graphpad)

#Fig. S2j: Plot of Mean +/- SD of z-scores by position, some assembly in Inkscape 

temp <- Missense_var %>%
  group_by(pos) %>%
  summarize(
    SD_expr_z_score = sd(Expr_z_score),
    SD_migr_z_score = sd(Migr_z_score),
    SD_prolif_z_score = sd(Prolif_z_score)
  )
Results_by_pos <- left_join(Results_by_pos, temp, by = "pos")

ExprSpread_byPos_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_expr_z_score + SD_expr_z_score)), 
            size = 0.6) +
  geom_line(mapping = aes(x = pos, y = (Mean_expr_z_score - SD_expr_z_score)),
            size = 0.6) +
  geom_ribbon(mapping = aes(x = pos, ymin = (Mean_expr_z_score - SD_expr_z_score),
                            ymax =  (Mean_expr_z_score + SD_expr_z_score)),
              fill = "grey")+
  theme_bw()
ExprSpread_byPos_Line

MigrSpread_byPos_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_migr_z_score + SD_migr_z_score)), 
            size = 0.6) +
  geom_line(mapping = aes(x = pos, y = (Mean_migr_z_score - SD_migr_z_score)),
            size = 0.6) +
  geom_ribbon(mapping = aes(x = pos, ymin = (Mean_migr_z_score - SD_migr_z_score),
                            ymax =  (Mean_migr_z_score + SD_migr_z_score)),
              fill = "grey")+
  theme_bw()
MigrSpread_byPos_Line

ProlifSpread_byPos_Line <- ggplot(Results_by_pos) +
  geom_line(mapping = aes(x = pos, y = (Mean_prolif_z_score + SD_prolif_z_score)), 
            size = 0.6) +
  geom_line(mapping = aes(x = pos, y = (Mean_prolif_z_score - SD_prolif_z_score)),
            size = 0.6) +
  geom_ribbon(mapping = aes(x = pos, ymin = (Mean_prolif_z_score - SD_prolif_z_score),
                            ymax =  (Mean_prolif_z_score + SD_prolif_z_score)),
              fill = "grey")+
  theme_bw()
ProlifSpread_byPos_Line

#Fig. S2k,l, made in Prism (Graphpad)
# But some analysis underlying that as follows:
Missense_var <- Missense_var %>%
  mutate(var_aa_class = case_when(var_aa %in% hydrophob_al ~ "hydrophob_al", 
                                  var_aa %in% hydrophob_ar ~ "hydrophob_ar",
                                  var_aa %in% h_bonding ~ "hydrophilic",
                                  var_aa %in% charged_pos ~ "charged_pos",
                                  var_aa %in% charged_neg ~ "charged_neg",
                                  var_aa %in% "C" ~ "C",
                                  var_aa %in% "P" ~ "P",
                                  var_aa %in% "G" ~ "G")) %>%
  mutate(wt_aa_class = case_when(wt_aa %in% hydrophob_al ~ "hydrophob_al", 
                                 wt_aa %in% hydrophob_ar ~ "hydrophob_ar",
                                 wt_aa %in% h_bonding ~ "hydrophilic",
                                 wt_aa %in% charged_pos ~ "charged_pos",
                                 wt_aa %in% charged_neg ~ "charged_neg",
                                 wt_aa %in% "C" ~ "C",
                                 wt_aa %in% "P" ~ "P",
                                 wt_aa %in% "G" ~ "G"))

Temp <- Missense_var %>%
  group_by(pos) %>%
  summarize(wt_aa_class = first(wt_aa_class)
  )

Results_by_pos <- left_join(Results_by_pos, Temp, by = "pos")
Domain_SD_means <- Results_by_pos %>%
  group_by(domain) %>%
  summarize(Mean_SD_expr = mean(SD_expr_z_score, na.rm = T))
WT_AA_SD_means <- Results_by_pos %>%
  group_by(wt_aa_class) %>%
  summarize(Mean_SD_expr = mean(SD_expr_z_score, na.rm = T))

Results_by_pos$domain <- factor(Results_by_pos$domain, 
                                levels = c("NTD", "TM", "EC", "IC", "H8", "CTD"))

Results_by_pos$wt_aa_class <- factor(Results_by_pos$wt_aa_class, 
                                     levels = c("C","P", "G", "charged_neg",
                                                "charged_pos", "hydrophilic",
                                                "hydrophob_ar", "hydrophob_al"))

Hi_SD_Expr <- Results_by_pos %>%
  filter(SD_expr_z_score >= 2)
Lo_SD_Lo_Expr <- Results_by_pos %>%
  filter(SD_expr_z_score < 2 & Mean_expr_z_score <= -2)
Lo_SD_Hi_Expr <- Results_by_pos %>%
  filter(SD_expr_z_score < 2 & Mean_expr_z_score > -2)
nrow(Hi_SD_Expr)
# [1] 154
nrow(Lo_SD_Lo_Expr)
# [1] 76
nrow(Lo_SD_Hi_Expr)
# [1] 123

cor.test(Results_by_pos$SD_expr_z_score, Results_by_pos$SD_migr_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$SD_expr_z_score and Results_by_pos$SD_migr_z_score
# S = 1772894, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7581682 
cor.test(Results_by_pos$SD_expr_z_score, Results_by_pos$SD_prolif_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$SD_expr_z_score and Results_by_pos$SD_prolif_z_score
# S = 4401682, p-value = 3.774e-15
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3995881 
cor.test(Results_by_pos$SD_migr_z_score, Results_by_pos$SD_prolif_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$SD_migr_z_score and Results_by_pos$SD_prolif_z_score
# S = 3287134, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5516181 

table(Hi_SD_Expr$domain)
# NTD  TM  EC  IC  H8 CTD 
# 2 138  12   0   2   0 
table(Lo_SD_Lo_Expr$domain)
# NTD  TM  EC  IC  H8 CTD 
# 2  55  14   2   3   0 
table(Lo_SD_Hi_Expr$domain)
# NTD  TM  EC  IC  H8 CTD 
# 13  28  10  15  10  47 

table(Hi_SD_Expr$domain, Hi_SD_Expr$wt_aa_class)
# C  P  G charged_neg charged_pos hydrophilic hydrophob_ar hydrophob_al
# NTD  0  0  0           0           0           0            0            2
# TM   6  2  2           2           7          25           19           75
# EC   0  1  0           0           3           1            3            4
# IC   0  0  0           0           0           0            0            0
# H8   0  0  0           0           0           1            0            1
# CTD  0  0  0           0           0           0            0            0
table(Lo_SD_Lo_Expr$domain, Lo_SD_Lo_Expr$wt_aa_class)
# C  P  G charged_neg charged_pos hydrophilic hydrophob_ar hydrophob_al
# NTD  0  0  0           1           0           0            0            1
# TM   3  7  1           0           3          11           10           20
# EC   1  0  2           1           1           3            2            4
# IC   0  1  0           0           1           0            0            0
# H8   0  0  0           1           0           1            1            0
# CTD  0  0  0           0           0           0            0            0
table(Lo_SD_Hi_Expr$domain, Lo_SD_Hi_Expr$wt_aa_class)
# C  P  G charged_neg charged_pos hydrophilic hydrophob_ar hydrophob_al
# NTD  0  2  1           0           1           7            0            2
# TM   1  0  3           2           7           2            2           11
# EC   0  1  0           0           3           1            2            3
# IC   0  1  1           2           4           3            1            3
# H8   1  0  1           1           4           0            1            2
# CTD  0  2  4           7           9          12            2           11

#Fig. S3a, made in Prism (Graphpad)

#Fig. S3b - representative fine-tuned ESM1b by Expression
# FilePath7 refers to csv with fie-tuned ESM1b scores for which Spearman correlation is closest to
# median spearman for 2 var/pos training: iteration 35.

FTESM1b <- read_csv(FilePath7, col_names = T)

Missense_var_FT <- left_join(Missense_var, FTESM1b, by = "full_var")

FTESM1b_by_Expr_Color <- ggplot(Missense_var_FT)+
  geom_point(mapping = aes(x = `pred_Expr_z_score:N2:iter35:FTR-LL3`, y = Expr_z_score,
                           color = domain),
             size = 0.6) +
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey", "red", "goldenrod1")) +    
  theme_bw()+
  labs(x = "Fine-tuned ESM1b score", y = "Expression z-score") + 
  ggtitle(label = "Fine-tuned ESM1b and Expression by variant")
FTESM1b_by_Expr_Color

FT_ESM1b_density <- ggplot(Missense_var_FT) +
  geom_density(mapping = aes(x = `pred_Expr_z_score:N2:iter35:FTR-LL3`), fill="grey")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim((floor(min(Missense_var_FT$`pred_Expr_z_score:N2:iter35:FTR-LL3`, na.rm = T)/0.5)*0.5), 
       (ceiling(max(Missense_var_FT$`pred_Expr_z_score:N2:iter35:FTR-LL3`, na.rm = T)/0.5)*0.5))
FT_ESM1b_density

#Spearman correlation 
cor.test(Missense_var$`pred_Expr_z_score:N2:iter35:FTR-LL3`, Missense_var$Expr_z_score, 
         use = "complete.obs", method = "spearman")

#Fig. S3c,d, made in Prism (Graphpad)

#Fig. S4a, AlphaMissense pathogenicity by expression plot; using logit(am_pathogenicity) to facilitate
# visualization of AM pathogenicity scores 

AM_by_Expr_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = logit(am_pathogenicity), y = Expr_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "logit(am_pathogenicity)", y = "Expression z-score")+
  ggtitle(label = "AM & Expr by Variant")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
AM_by_Expr_plot

AM_density <- ggplot(Missense_var) +
  geom_density(mapping = aes(x = logit(am_pathogenicity)), fill="grey")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim((floor(min(logit(Missense_var_FT$am_pathogenicity), na.rm = T)/0.5)*0.5), 
       (ceiling(max(logit(Missense_var_FT$am_pathogenicity), na.rm = T)/0.5)*0.5))
AM_density

#Spearman correlation AM pathogenicity and expression-score 
cor.test(logit(Missense_var$am_pathogenicity), Missense_var$Expr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  logit(Missense_var$am_pathogenicity) and Missense_var$Expr_z_score
# S = 8.0519e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.6012783 

# Fig. S4b - AM pathogenicity by migration
AM_by_Migr_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = logit(am_pathogenicity), y = Migr_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "logit(am_pathogenicity)", y = "Migration z-score")+
  ggtitle(label = "AM & Migr by Variant")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
AM_by_Migr_plot

# Spearman AM pathogenicity, migration 
cor.test(logit(Missense_var$am_pathogenicity), Missense_var$Migr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  logit(Missense_var$am_pathogenicity) and Missense_var$Migr_z_score
# S = 1.506e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7005085 

#Fig. S4c - AM pathogenicity by proliferation dot plot 
AM_by_Prolif_plot <- ggplot(Missense_var)+
  geom_point(mapping = aes(x = logit(am_pathogenicity), y = Prolif_z_score,
                           color = domain), size = 0.6)+
  theme_bw()+
  labs(x = "logit(am_pathogenicity)", y = "Proliferation z-score")+
  ggtitle(label = "AM & Prolif by Variant")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
AM_by_Prolif_plot

#Spearman test, AM pathogenicity and prolif z-score
cor.test(logit(Missense_var$am_pathogenicity), Missense_var$Prolif_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  logit(Missense_var$am_pathogenicity) and Missense_var$Prolif_z_score
# S = 1.9546e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6112986 

#Fig S4d - AM pathogenicity and expression by position 
AM_by_Expr_byPos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = logit(Mean_am_pathogenicity), y = Mean_expr_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "logit(am_pathogenicity)", y = "mean(Expression z-score)")+
  ggtitle(label = "AM & Expr by Position")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
AM_by_Expr_byPos_plot

#Spearman correlation mean AM pathogenicity, mean expression z-score 
cor.test(Results_by_pos$Mean_am_pathogenicity, Results_by_pos$Mean_expr_z_score,
         use = "complete.obs", method = "spearman")
# 	Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_am_pathogenicity and Results_by_pos$Mean_expr_z_score
# S = 12146330, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.6568214 

#Fig S4e - AM pathogenicity and migration by position 
AM_by_Migr_byPos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = logit(Mean_am_pathogenicity), y = Mean_migr_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "logit(am_pathogenicity)", y = "mean(Migration z-score)")+
  ggtitle(label = "AM & Migr by Position")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
AM_by_Migr_byPos_plot

# Spearman correlation, mean AM pathogenicity, mean migration z-score 
cor.test(Results_by_pos$Mean_am_pathogenicity, Results_by_pos$Mean_migr_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_am_pathogenicity and Results_by_pos$Mean_migr_z_score
# S = 1108966, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8487314 

#Fig S4f - AM pathogenicity and proliferation by position 
AM_by_Prolif_byPos_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = logit(Mean_am_pathogenicity), y = Mean_prolif_z_score,
                           color = domain))+
  theme_bw()+
  labs(x = "logit(am_pathogenicity)", y = "mean(Proliferation z-score)")+
  ggtitle(label = "AM & Prolif by Position")+
  scale_color_manual(values = c("deeppink4", "blue", "lightgreen", "grey",
                                "red", "goldenrod1"))
AM_by_Prolif_byPos_plot

#Spearman correlation, mean AM pathogenicity, mean proliferation z-score 
cor.test(Results_by_pos$Mean_am_pathogenicity, Results_by_pos$Mean_prolif_z_score,
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Results_by_pos$Mean_am_pathogenicity and Results_by_pos$Mean_prolif_z_score
# S = 1507432, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7943786 

#Fig. S4g-k, created in Prism (Graphpad)

# Fig. S5a - Expression by proliferation, highlighting high Expr-adj migr score

# Set threshold for elevated expression-adjusted proliferation: more than 2x
# the interquartile range above the median 
EaP_25percent = quantile(Results_by_pos$Mean_EaP_score, 0.25, na.rm = T) 
EaP_median = quantile(Results_by_pos$Mean_EaP_score, 0.5, na.rm = T) 
EaP_75percent = quantile(Results_by_pos$Mean_EaP_score, 0.75, na.rm = T) 
EaP_threshold = 2*((EaP_75percent - EaP_25percent)) + EaP_median
EaP_threshold
# 0.2632164

Expr_by_Prolif_byPos_highlighted_plot <- ggplot(Results_by_pos)+
  geom_point(mapping = aes(x = Mean_expr_z_score, y = Mean_prolif_z_score,
                           color = Mean_EaP_score > EaP_threshold ))+
  theme_bw()+
  labs(x = "Mean expression z-score", y = "Mean proliferation z-score")+
  scale_colour_manual(values = setNames(c('blue',"grey"),c(T, F)))+ 
  geom_point(Results_by_pos %>% filter(Mean_EaM_score > EaM_threshold),
             mapping = aes(x = Mean_expr_z_score, y = Mean_prolif_z_score),
             color = "red")
Expr_by_Prolif_byPos_highlighted_plot

#Fig. S6-9 made without R

#Fig. S10a,b made in Prism (Graphpad)

#Fig S10c - gnomAD variant allele frequency by expression plot
gnomAD_Expr_plot <- ggplot(Missense_var) +
  geom_point(mapping = aes(x = log(allele_freq, 10), y = Expr_z_score))+
  theme_bw()+
  labs(x = "log(allele frequency)", y = "Expression z-score")+
  ggtitle(label = "Allele freq & Expr by Variant")
gnomAD_Expr_plot

cor.test(Missense_var$allele_freq, Missense_var$Expr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$allele_freq and Missense_var$Expr_z_score
# S = 14582572, p-value = 0.000395
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1626161 

#Fig S10d - gnomAD variant allele frequency by migration plot
gnomAD_Migr_plot <- ggplot(Missense_var) +
  geom_point(mapping = aes(x = log(allele_freq, 10), y = Migr_z_score))+
  theme_bw()+
  labs(x = "log(allele frequency)", y = "Migration z-score")+
  ggtitle(label = "Allele freq & Migr by Variant")
gnomAD_Migr_plot

cor.test(Missense_var$allele_freq, Missense_var$Migr_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$allele_freq and Missense_var$Migr_z_score
# S = 20207819, p-value = 0.0004747
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1604059 

cor.test(Missense_var$allele_freq, Missense_var$E_adj_M_score, 
         use = "complete.obs", method = "spearman")

#Fig S10e - gnomAD variant allele frequency by proliferation plot 
gnomAD_Prolif_plot <- ggplot(Missense_var) +
  geom_point(mapping = aes(x = log(allele_freq, 10), y = Prolif_z_score))+
  theme_bw()+
  labs(x = "log(allele frequency)", y = "Proliferation z-score")+
  ggtitle(label = "Allele freq & Prolif by Variant")
gnomAD_Prolif_plot

cor.test(Missense_var$allele_freq, Missense_var$Prolif_z_score, 
         use = "complete.obs", method = "spearman")
# Spearman's rank correlation rho
# 
# data:  Missense_var$allele_freq and Missense_var$Prolif_z_score
# S = 20426823, p-value = 0.0001615
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1729819 

#Fig. S10f, Kolmogorov-Smirnov tests for all missense, gnomAD, DLBCL/BL, and
#non-hematopoietic cancer variants, for expression, migration, and proliferation 
gnomAD_temp <- Missense_var %>% filter(!is.na(allele_freq))
ks.test(Missense_var$Expr_z_score, gnomAD_temp$Expr_z_score)
ks.test(Missense_var$Expr_z_score, Cosmic_nonHL$Expr_z_score)
ks.test(Missense_var$Expr_z_score, All_var_D_B_nodup$Expr_z_score)
ks.test(gnomAD_temp$Expr_z_score, Cosmic_nonHL$Expr_z_score)
ks.test(gnomAD_temp$Expr_z_score, All_var_D_B_nodup$Expr_z_score)
ks.test(Cosmic_nonHL$Expr_z_score, All_var_D_B_nodup$Expr_z_score)

ks.test(Missense_var$Migr_z_score, gnomAD_temp$Migr_z_score)
ks.test(Missense_var$Migr_z_score, Cosmic_nonHL$Migr_z_score)
ks.test(Missense_var$Migr_z_score, All_var_D_B_nodup$Migr_z_score)
ks.test(gnomAD_temp$Migr_z_score, Cosmic_nonHL$Migr_z_score)
ks.test(gnomAD_temp$Migr_z_score, All_var_D_B_nodup$Migr_z_score)
ks.test(Cosmic_nonHL$Migr_z_score, All_var_D_B_nodupMigr_z_score)

ks.test(Missense_var$Prolif_z_score, gnomAD_temp$Prolif_z_score)
ks.test(Missense_var$Prolif_z_score, Cosmic_nonHL$Prolif_z_score)
ks.test(Missense_var$Prolif_z_score, All_var_D_B_nodup$Prolif_z_score)
ks.test(gnomAD_temp$Prolif_z_score, Cosmic_nonHL$Prolif_z_score)
ks.test(gnomAD_temp$Prolif_z_score, All_var_D_B_nodup$Prolif_z_score)
ks.test(Cosmic_nonHL$Prolif_z_score, All_var_D_B_nodup$Prolif_z_score)

###End
  
