# Pipeline Analysis General Population

library(here)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggforce)
library(patchwork)
library(factoextra)
library(psych)
library(MASS)
library(stargazer)
library(pscl)
library(sjPlot)

##################################################################################################################
# Charge and visualize dataset
##################################################################################################################
# Repository definition for R project. All files called from this repo.
here()

df_raw <- read.csv(file = here("data", "2026-controls-CDS.csv"), 
               header = TRUE, sep = ";", dec = ",", na.strings = c("NaN", "Na", "#NUL!", " "))
head(df_raw)
colnames(df_raw)
nrow(df_raw)

# Transform df to exclude rows from VESTICOR study
df_raw2 <- df_raw[1:665, ]
nrow(df_raw2)
as.data.frame(table(df_raw2$raison_exclusion))
as.data.frame(table(df_raw2$Difference_sum))

inclus <- df_raw2$INCLUSION_CDS_VALIDATION==1
df <- df_raw2[inclus,]
nrow(df)
as.data.frame(table(df$Difference_sum))

df$MIGRAINE[df$MIGRAINE == ""] <- NA
df$MIGRAINE <- factor(df$MIGRAINE)

table(df$MIGRAINE, useNA = "ifany")
levels(df$MIGRAINE)

df$SEXE <- relevel(factor(df$SEXE), ref="h")
levels(df$SEXE)


##################################################################################################################
# Preprocessing
##################################################################################################################
# NaN inspection
##################################################################################################################

# Define interest columns
cols_int <- c("SEXE", "LATERALITE", "AGE", "PROFESSION", "FAMILLE", "FUMEUR", "MIGRAINE", "ALCOOL", "aucuntroublevisuel",
              "myopie", "astigmatisme", "maladierétine", "glaucome", "trouble_vision_global", "lunettes", "lentilles", "ETUDES", "CDS1", "CDS2",
              "CDS3", "CDS4", "CDS5", "CDS6", "CDS7", "CDS8", "CDS9", "CDS10", "CDS11", "CDS12", "CDS13", "CDS14", 
              "CDS15","CDS16", "CDS17", "CDS18", "CDS19", "CDS20", "CDS21", "CDS22", "CDS23", "CDS24", "CDS25", "CDS26",
              "CDS27","CDS28", "CDS29", "OBE", "ANXIETE_recoded", "DEPRESSION_recoded", "CDStotal_matlab", "CDS_total_sum", "FREQUENCYall", "DURATIONall")

# Count number of missing/empty values for each interest column
colSums(is.na(df[cols_int]) | df[cols_int]=="" |is.null(df[cols_int]))

# Display row number of missing/empty values for each variable
missing <- sapply(df[cols_int], function(x) which(is.na(x) | x == "" | is.null(x)) + 1)
missing

##################################################################################################################
# Identify aberrant values (> total score
##################################################################################################################
# For CDS items
CDS_items <- c("CDS1", "CDS2","CDS3", "CDS4", "CDS5", "CDS6", "CDS7", "CDS8", "CDS9", "CDS10", "CDS11", "CDS12",
               "CDS13", "CDS14", "CDS15","CDS16", "CDS17", "CDS18", "CDS19", "CDS20", "CDS21", "CDS22", "CDS23", 
               "CDS24", "CDS25", "CDS26", "CDS27","CDS28", "CDS29")
for (col in CDS_items) {
  av <- sum(df[[col]]  > 10, na.rm = TRUE)
  cat (col, " : Aberrant values (>10):", av, "\n")
}


# For CDS total score
av <- sum(df$CDS_total_sum > 290, na.rm = TRUE)
cat ("CDS total : Aberrant values (>290):", av, "\n")


# For HADS A and D
anxdep <- c("ANXIETE_recoded", "DEPRESSION_recoded")
for (col in anxdep) {
  av <- sum(df[[col]]  > 21, na.rm = TRUE)
  cat (col, ": Aberrant values (>21):", av, "\n")
}


##################################################################################################################
# Visualization of variables distribution
##################################################################################################################
num_cols <- c("AGE", "CDS1", "CDS2","CDS3", "CDS4", "CDS5", "CDS6", "CDS7", "CDS8", "CDS9", "CDS10", "CDS11", "CDS12",
              "CDS13", "CDS14", "CDS15","CDS16", "CDS17", "CDS18", "CDS19", "CDS20", "CDS21", "CDS22", "CDS23", 
              "CDS24", "CDS25", "CDS26", "CDS27","CDS28", "CDS29", "ANXIETE", "DEPRESSION", "CDStotal_matlab", "CDS_total_sum", 
              "FREQUENCYall", "DURATIONall")

df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

summary(df[num_cols])

describe(df$CDS_total_sum)

print(paste("Number of participants above CDS cutoff (70) :", sum(df$CDS_total_sum > 70)))
print(paste("Proportion of participants above CDS cutoff (70) :", round((sum(df$CDS_total_sum > 70)/nrow(df))*100, 2), "%"))

####### Plot numerical variables distribution #######

#####Classic distribution plots (histogram + density) ######
distrib_plot_list <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x= .data[[col]])) +
  geom_histogram(aes(y=after_stat(density)), bins=30, fill="purple4", color="white") +
  geom_density(fill="grey", alpha = 0.5) +
  geom_vline(xintercept = mean(df[[col]], na.rm=TRUE), linetype="dashed", color="turquoise", linewidth=1) +
  geom_vline(xintercept = median(df[[col]], na.rm=TRUE), linetype="dotdash", color="red", linewidth=1) +
  theme_minimal()
})
names(distrib_plot_list) <- num_cols

# Combine all plots
wrap_plots(distrib_plot_list, ncol=6)

# Save figure
ggsave(here("figures", "numeric_distribution_plots.png"),
       wrap_plots(distrib_plot_list, ncol=6),
       width = 30, height = 20, dpi = 300)


###### Visualize (boxplots with outliers) ######
box_plot_list <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x = "", y= .data[[col]])) +
    geom_boxplot(outliers = TRUE, outlier.size = 0.7, outlier.color = "grey42", fill="skyblue") +
    labs(x=NULL, y=col) +
    theme_light()
  return(p)
})
names(box_plot_list) <- num_cols

# Combine all plots
wrap_plots(box_plot_list, ncol=6)

# Save figure
ggsave(here("figures", "numeric_box_plots.png"),
       wrap_plots(box_plot_list, ncol=6),
       width = 30, height = 20, dpi = 300)



##### Visualize (violin plots) #####
violin <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x = "", y= .data[[col]])) +
    geom_violin(stat="ydensity",
                quantile.linetype = "dashed",
                quantile.linewidth = 0.7,
                fill="dodgerblue3", 
                alpha=0.6, 
                quantile.colour = "black") +
    labs(x=NULL, y=col) +
    theme_light()
  return(p)
})
names(violin) <- num_cols

# Combine all plots
wrap_plots(violin, ncol=6)

# Save figure
ggsave(here("figures", "violin_plots.png"),
       wrap_plots(violin, ncol=6),
       width = 30, height = 20, dpi = 300)


##### Visualize (boxplots + violin plots) #####
box_violin <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x = "", y= .data[[col]])) +
    geom_boxplot(outliers = FALSE, fill="lightskyblue") +
    geom_violin(fill="midnightblue", alpha=0.4) +
    labs(x=NULL, y=col) +
    theme_light()
  return(p)
})
names(box_violin) <- num_cols

# Combine all plots
wrap_plots(box_violin, ncol=6)


###### Visualize (boxplots + violin plots_ Alternative) #####
box_violin2 <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x =" ", y= .data[[col]])) +
    geom_violin(fill="lightsteelblue2", alpha=0.7) +
    geom_boxplot(outliers=FALSE, width = 0.1, fill="cornflowerblue") +
    labs(x=NULL, y=col) +
    theme_light()
  return(p)
})
names(box_violin2) <- num_cols

# Combine all plots
wrap_plots(box_violin2, ncol=6)


# Save figure
ggsave(here("figures", "box_violin_plots2.png"),
       wrap_plots(box_violin2, ncol=6),
       width = 30, height = 20, dpi = 300)



###### Visualize (boxplots + dots) ######
box_dots <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x = "", y= .data[[col]])) +
    geom_boxplot(outliers=FALSE, fill="lightskyblue") +
    geom_sina(color= "royalblue2", alpha = 0.5, size= 1) +
    labs(x=NULL, y=col) +
    theme_light()
  return(p)
})
names(box_dots) <- num_cols

# Combine all plots
wrap_plots(box_dots, ncol=6)

# Save figure
ggsave(here("figures", "box_dots_plots.png"),
       wrap_plots(box_dots, ncol=6),
       width = 30, height = 20, dpi = 300)


##### Visualize (boxplots + violin plots + dots) #####
box_violin_dots <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x = "", y= .data[[col]])) +
    geom_boxplot(outliers=FALSE, fill="lightskyblue") +
    geom_violin(fill="royalblue4", alpha=0.6) +
    geom_sina(color= "royalblue2", alpha = 0.5, size= 1) +
    labs(x=NULL, y=col) +
    theme_light()
  return(p)
})
names(box_violin_dots) <- num_cols

# Combine all plots
wrap_plots(box_violin_dots, ncol=6)

# Save figure
ggsave(here("figures", "box_violin_dots_plots.png"),
       wrap_plots(box_violin_dots, ncol=6),
       width = 30, height = 20, dpi = 300)



########## Plot categorical variables ###############
cat_cols <- c("SEXE", "LATERALITE", "PROFESSION", "FAMILLE", "FUMEUR", "MIGRAINE", "ALCOOL", "aucuntroublevisuel",
              "myopie", "astigmatisme", "maladierétine", "glaucome", "lunettes", "lentilles", "trouble_vision_global", "ETUDES", 
              "OBE")

df[cat_cols] <- lapply(df[cat_cols], factor)

# Description categorical variables
sum_cat <- lapply(df[cat_cols], function(x) as.data.frame(table(x)))
sum_cat

# Barplot of categorical variables distribution
cat_dist_plots <- lapply(cat_cols, function(col) {
  p<- ggplot(df %>% filter(!is.na(.data[[col]]), !.data[[col]] %in% c("", "NaN", "#NUL!")), aes(x = .data[[col]])) +
    geom_bar(fill="darkblue") +
    geom_text(stat="count", aes(label=after_stat(count)), vjust=-0.5, size = 3) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.15))) +
    theme_light() +
    theme(legend.text=element_text(size=6), 
          plot.title=element_text(size=8),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8))
  return(p)
})
names(cat_dist_plots) <- cat_cols

wrap_plots(cat_dist_plots, ncol=6)

ggsave(here("figures", "catge_distribution_plots.png"),
       wrap_plots(cat_dist_plots, ncol=6),
       width = 30, height = 20, dpi = 300)



##################################################################################################################
# Outliers Identification
##################################################################################################################
# Count and identify values of outliers
for (col in num_cols) {
  out_values <- boxplot.stats(df[[col]])$out
  cat("\n=== Colonne:", col, "===\n")
  cat("Nombre d'outliers:", length(out_values), "\n")
  if (length(out_values) > 0) {
    cat("Valeurs:", paste(out_values, collapse = ", "), "\n")
  }
}




#############################################################################################################################
# Exploratory Factor Analysis #
#############################################################################################################################
items <- df %>% 
  dplyr::select(all_of(CDS_items)) %>% 
  filter(complete.cases(.))
summary(items)
sum(is.na(items))

item_names <- names(items)

###### Preliminary inspections ######
describe(items)
describe(items)[, c("skew", "kurtosis")]

# Number of observations per item category (score)
sapply(items, table)
# Percentage of zeros for each item
pct_zero <- sapply(items, function(x) mean(x == 0, na.rm = TRUE) * 100)
pct_zero_df <- data.frame(item = names(pct_zero), pct_zero = round(pct_zero, 1)) |>
  arrange(desc(pct_zero))
print(pct_zero_df)
# Percentage of zeros extremely high (>95%) for items 19, 12, 20 and 27, with really restricted variance and extremely asymetric distribution (skewness)

# Number of categories used in each item
n_categories <- sapply(items, function(x) length(unique(na.omit(x))))
print(n_categories)

min_nonzero_n <- sapply(items, function(x) {
  tab <- table(x)
  tab <- tab[names(tab) != "0"]
  if (length(tab) == 0) return(NA)
  min(tab)
})
print(sort(min_nonzero_n))


##### Assumptions ######

cor_matrix <- cor(items, use = "pairwise.complete.obs", method = "spearman")

# Multicolinearity
corr_values <- cor_matrix[lower.tri(cor_matrix)]
range(corr_values)

mean(cor_matrix[lower.tri(cor_matrix)], na.rm = TRUE)


library(EFAtools)

# EFA assumptions and factor selecton : 
# Gives KMO and Bartlett's test for assumtpions.
# Gives criteria (parallel, MAP, HULL, EKC) for factor number
EFA_assumptions <- efa_retain(items,
                              cor_method = "spearman",
                              estimator =  "ULS")
EFA_assumptions


# Kaiser criterion
eigenvalues <- eigen(cor_matrix)$values
sum(eigenvalues > 1)   # number of eigenvalues > 1
# 7 possible factors

# Scree plot
pca <- prcomp(items, scale. = TRUE)
fviz_screeplot(pca, addlabels = TRUE, ncp = 10)  # show first 10 components
ggsave(here("figures", "Screeplot.png"), width=8, height=6, dpi=300)



####### EFA ########
### Comparison of EFA models with 1 to 3 factors for Pearson and Spearman correlations
# Correlation matrices
mat_spearman <- cor(items, method = "spearman")
mat_pearson  <- cor(items, method = "pearson")
n_sujets     <- nrow(items)

compare_table <- data.frame()

# Loop to have fit indices for all models
for (type in c("Spearman", "Pearson")) {
  matrix <- if(type == "Spearman") mat_spearman else mat_pearson
  
  for (k in 1:3) {
    modele <- fa(matrix, nfactors = k, n.obs = n_sujets, fm = "minres", rotate = "promax")
    
    fit_k <- data.frame(
      Methode     = type,
      Facteurs    = k,
      Chi_deux    = round(modele$chi, 2),
      p_value     = round(modele$PVAL, 4),
      TLI         = round(modele$TLI, 3),
      RMSR        = round(modele$rms, 3),
      RMSEA       = round(modele$RMSEA[1], 3),     
      RMSEA_inf   = round(modele$RMSEA[2], 3),     
      RMSEA_sup   = round(modele$RMSEA[3], 3),     
      Fit_OffDiag = round(modele$fit.off, 3)
    )
    
    compare_table <- rbind(compare_table, fit_k)
  }
}

print(compare_table, row.names = FALSE)

# Loop to display factor loadings for EFA with 1 to 3 factors with Spearman correlations
for (k in 1 :3) {
  cat("Exploratory Factor Analysis")
  model_k <- fa(mat_spearman, nfactors = k, n.obs = n_sujets, fm="minres", rotation ="promax")
  print(model_k$loadings, cutoff=0.4)
}




###### Final EFA model selected based on previous explorations #####
#### EFA with 3 factors, spearman correlations, promax rotation (correlated factors) and minres estimator
efa_result <- fa(items, n.obs = nrows(items), nfactors = 2, rotate = "promax", fm = "minres", cor="spearman")
print(efa_result$loadings, cutoff = 0.4, digits = 3)

fit_metrics <- data.frame(
  Chi2     = round(efa_result$chi, 2),
  df       = round(efa_result$dof,2),
  p_value  = round(efa_result$PVAL, 4),
  TLI      = round(efa_result$TLI, 3),
  RMSR     = round(efa_result$rms, 3),
  RMSEA    = round(efa_result$RMSEA[1], 3),
  RMSEA_inf= round(efa_result$RMSEA[2], 3),
  RMSEA_sup= round(efa_result$RMSEA[3], 3),
  R2_Value = round(efa_result$R2, 3)
)
print(fit_metrics, row.names = FALSE)

loadings_matrix <- unclass(efa_result$loadings)
print(loadings_matrix, digits=3)

print(sort(efa_result$communalities))
print(efa_result$Phi)


##### Residuals analysis ##########
res <- efa_result$residual
res_vals <- res[lower.tri(res)]    

n_pairs <- length(res_vals)
cat("Résidus |r| > .05 :", sum(abs(res_vals) > .05), "/", n_pairs,
    sprintf("(%.1f%%)\n", 100 * mean(abs(res_vals) > .05)))
cat("Résidus |r| > .08 :", sum(abs(res_vals) > .08), "/", n_pairs,
    sprintf("(%.1f%%)\n", 100 * mean(abs(res_vals) > .08)))
cat("Résidus |r| > .10 :", sum(abs(res_vals) > .10), "/", n_pairs,
    sprintf("(%.1f%%)\n", 100 * mean(abs(res_vals) > .10)))

###################################################################################################################
# Adding CDS factors to data table
##################################################################################################################
df$CDS_F1 <- rowSums(df[, c("CDS1","CDS2","CDS3","CDS6","CDS8", "CDS10", "CDS13", "CDS15", "CDS23", "CDS24", "CDS26")])
df$CDS_F2 <- rowSums(df[, c("CDS5","CDS7","CDS9","CDS25","CDS28", "CDS29")])

# We only sum the items here to have the factors scores to keep a distribution with integers that allows us for
# a binomial fit of the Hurdle model. If we do a mean (divide by item number), the score becomes continuous and the
# binomial distribution would not apply anymore, thus causing discrepancy between our analysis types.
# We can create a mean factor score to compare both factors in controls using descriptive stats.

CDS <- df[, 29:57]
CDS_Fac1 <- df[, c("CDS1","CDS2","CDS3","CDS6","CDS8", "CDS10", "CDS13", "CDS15", "CDS23", "CDS24", "CDS26")]
CDS_Fac2 <- df[, c("CDS5","CDS7","CDS9","CDS25","CDS28", "CDS29")]

###### Reliability #####
# CDS Total
reliability(items=CDS, nfactors=1)
alpha(CDS)
splitHalf(CDS)

# Factor 1
reliability(items=CDS_Fac1, nfactors=1)
alpha(CDS_Fac1)
splitHalf(CDS_Fac1)

# Factor 2
reliability(items=CDS_Fac2, nfactors=1)
alpha(CDS_Fac2)
splitHalf(CDS_Fac2)



###################################################################################################################
# Linear Regression Models #
###################################################################################################################
# Filter to remove cases without complete anxiety, depression and migraine
anxdep <- c("ANXIETE_recoded", "DEPRESSION_recoded", "MIGRAINE")
df2 <- df[complete.cases(df[, anxdep]), ]
nrow(df2)

df2$AGE_centered <- df2$AGE-mean(df2$AGE)
describe(df2$AGE_centered)

#####################################################################
# Violin plots : CDS total according to anxiety/depression level (HADS)
#####################################################################
df2$ANX_grp <- factor(ifelse(df2$ANXIETE_recoded >= 8, "≥ 8", "< 8"),
                      levels = c("< 8", "≥ 8"))
df2$DEP_grp <- factor(ifelse(df2$DEPRESSION_recoded >= 8, "≥ 8", "< 8"),
                      levels = c("< 8", "≥ 8"))

table(df2$ANX_grp)
table(df2$DEP_grp)

# Convertit une p-value en étoiles (ns si non significatif)
p_to_stars_ns <- function(p) {
  case_when(p < 0.001 ~ "***",
            p < 0.01  ~ "**",
            p < 0.05  ~ "*",
            TRUE      ~ "ns")
}

make_violin <- function(data, grp, xlab, fill_cols = c("lightskyblue", "royalblue4"),
                        cutoff = 70) {
  
  # Effectifs par groupe pour l'axe x
  tab      <- table(data[[grp]])
  x_labels <- paste0(names(tab), "\n(n = ", tab, ")")
  
  # Test de Wilcoxon (distribution non normale, zéro-inflatée)
  w     <- wilcox.test(data$CDS_total_sum ~ data[[grp]])
  stars <- p_to_stars_ns(w$p.value)
  
  # Position du trait et des étoiles, au-dessus du max des données ET du seuil
  y_max   <- max(data$CDS_total_sum, na.rm = TRUE)
  y_line  <- max(y_max, cutoff) * 1.05
  tick_h  <- max(y_max, cutoff) * 0.02
  y_text  <- y_line + tick_h * 0.5
  
  ggplot(data, aes(x = .data[[grp]], y = CDS_total_sum, fill = .data[[grp]])) +
    geom_violin(alpha = 0.6, color = "grey30", trim = TRUE) +
    geom_sina(color = "grey45", alpha = 0.2, size = 0.5) +
    geom_boxplot(width = 0.12, outliers = FALSE, fill = "slategray3", alpha = 0.9) +
    # Seuil CDS
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "black", linewidth = 0.6) +
   
    # Crochet de significativité
    annotate("segment", x = 1, xend = 2, y = y_line, yend = y_line, linewidth = 0.5) +
    annotate("segment", x = 1, xend = 1, y = y_line - tick_h, yend = y_line, linewidth = 0.5) +
    annotate("segment", x = 2, xend = 2, y = y_line - tick_h, yend = y_line, linewidth = 0.5) +
    annotate("text", x = 1.5, y = y_text, label = stars,
             size = ifelse(stars == "ns", 3.5, 6), vjust = ifelse(stars == "ns", -0.2, 0.3)) +
    scale_fill_manual(values = fill_cols, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)), breaks = seq(0, 160, by = 20)) +
    labs(x = xlab, y = "CDS total score") +
    theme_classic(base_size = 12)
}

p_vio_anx <- make_violin(df2, "ANX_grp", "Anxiety (HADS-A)")
p_vio_dep <- make_violin(df2, "DEP_grp", "Depression (HADS-D)")


violin_hads <- (p_vio_anx | p_vio_dep) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.caption = element_text(size = 8, hjust = 0)))

violin_hads

ggsave(here("Figures", "violin_CDS_total_HADS.pdf"), violin_hads, width = 7, height = 4.5)
ggsave(here("Figures", "violin_CDS_total_HADS.png"), violin_hads, width = 7, height = 4.5, dpi = 300)

###################################################################################################################
# CDS Total Score
###################################################################################################################
# Description of CDS total score
summary(df2$CDS_total_sum)
describe(df2$CDS_total_sum)
shapiro.test(df2$CDS_total_sum)
# Non normal (skewed, leptokurique), overdispersion and autocorrelation of residuals in regression model see test model)

# Logarithmic transformation of CDS total score
df2$CDS_tot_log <- log1p(df2$CDS_total_sum)
summary(df2$CDS_tot_log)
describe(df2$CDS_tot_log)[, c("skew", "kurtosis")]
shapiro.test(df2$CDS_tot_log)
# Non normal (skewed, leptokurique), overdispersion and autocorrelation of residuals in regression model see test model)

# Root Square transformation of CDS total score
df2$CDS_tot_sqrt <- sqrt(df2$CDS_total_sum)
summary(df2$CDS_tot_sqrt)
describe(df2$CDS_tot_sqrt)[, c("skew", "kurtosis")]
shapiro.test(df2$CDS_tot_sqrt)


# Distribution plots for CDS total scores and its log and sqrt transformations
CDS_forms <- c("CDS_total_sum", "CDS_tot_log", "CDS_tot_sqrt")

distrib_CDS <- lapply(CDS_forms, function(x) {
  
  # Ajout de na.rm=TRUE au cas où vous auriez des données manquantes
  col_mean = mean(df2[[x]], na.rm = TRUE)
  col_median = median(df2[[x]], na.rm = TRUE)
  
  p <- ggplot(df2, aes(x = .data[[x]])) + # Correction 1 : .data[[x]] pour évaluer le texte
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "purple4", color = "white") +
    geom_density(fill = "grey", alpha = 0.5) +
    geom_vline(xintercept = col_mean, linetype = "dashed", color = "turquoise", linewidth = 1) +
    geom_vline(xintercept = col_median, linetype = "dotdash", color = "red", linewidth = 1) +
    labs(
      title = paste("Distribution of", x), # Correction 2 : paste() obligatoire ici
      x = "Score Value",
      y = "Density"
    ) + 
    theme_minimal()
  
  return(p)
})

wrap_plots(distrib_CDS, ncol = 3)

ggsave(here("figures", "distrib_CDS_transformed.png"), width=8, height=6, dpi=300)

# Zero inflation that remains throughout any transformation


########## Test of classical regression model #########
### Modèle 1 ###
Model_1 <- lm(CDS_tot_log~ AGE, na.action = na.exclude, data=df2) # or CDS_tot_log or CDS_tot_sqrt
summary(Model_1)

# Vérification des prérequis
# 1. Indépendance des résidus : pas de structure particulière dans la distribution des résidus
plot(Model_1, 1)
# Test de durbin-watson
library(lmtest)
dwtest(Model_1)
# Test de Breusch Godfrey
library(lmtest)
bgtest(Model_1)

# 2. Homoscédasticité : test de Harrison-McCabe
library(lmtest)
hmctest(Model_1)

# 3. Normalité
rstandard(Model_1)
hist(Model_1$residuals) # Histogramme des résidus standardisés
plot(Model_1, 2) # QQplot

shapiro.test(rstandard(Model_1))

library(fBasics)
dagoTest(rstandard(Model_1))

# 4. Détection des valeurs atypiques et aberrantes
# résidu std > 2 = atypique; > 3 = aberrant
plot(rstandard(Model_1))
abline(h = c(-2, 2), col="red", lty=2)
abline(h = c(-3, 3), col="forestgreen", lty=2)

# 5. Dipersion of CDS total distribution
chi2 <- sum(residuals(Model_1, "pearson")^2)
chi2 / df.residual(Model_1)
1 - pchisq(chi2, df = df.residual(Model_1))


# Assumptions are not met for classical linear regression as the model residuals are not
# normally distributed, present with heteroscedasticity and autocorrelation. 
# As the model outcome (CDS total score) is overdispersed and have a positively skewed
# distribution with many 0 values, a hurdle (negative binomial) regression would be more adapted.

#########################################################################################
# Hurdle Regression Model
#########################################################################################
# Test of fit for distribution of CDS scores > 0
scores_pos <- df2$CDS_total_sum[df2$CDS_total_sum > 0]
df_pos <- data.frame(scores_pos=scores_pos)

# Distribution of CDS total score
CDS_tot_distrib <- ggplot(data=df2, aes(x=CDS_total_sum)) +
  geom_histogram(aes(y=after_stat(density)), bins=30, fill="purple4", color="white") +
  geom_density(fill="grey", alpha = 0.5) +
  geom_vline(xintercept = mean(df2$CDS_total_sum, na.rm=TRUE), linetype="dashed", color="turquoise", linewidth=1) +
  geom_vline(xintercept = median(df2$CDS_total_sum, na.rm=TRUE), linetype="dotdash", color="red", linewidth=1) +
  theme_minimal()

# Distribution of CDS positive score (no 0 score)
CDS_pos <- ggplot(data=df_pos, aes(x=scores_pos)) +
  geom_histogram(aes(y=after_stat(density)), bins=30, fill="purple4", color="white") +
  geom_density(fill="grey", alpha = 0.5) +
  geom_vline(xintercept = mean(df_pos$scores_pos, na.rm=TRUE), linetype="dashed", color="turquoise", linewidth=1) +
  geom_vline(xintercept = median(df_pos$scores_pos, na.rm=TRUE), linetype="dotdash", color="red", linewidth=1) +
  theme_minimal()


CDS_distrib <- CDS_tot_distrib+CDS_pos
CDS_distrib


# Fitting of CDS (only scores > 0) distribution
library(fitdistrplus)

fit_lognorm <- fitdist(scores_pos, "lnorm")
fit_gamma <- fitdist(scores_pos, "gamma")
fit_binom <- fitdist(scores_pos, "nbinom")

summary(fit_lognorm)
summary(fit_gamma)
summary(fit_binom)
gofstat(list(fit_lognorm, fit_gamma, fit_binom), fitnames = c("Log-Normal", "Gamma", "Binomial"))

# QQ plot for log normal and gamma distribution compared to the distribution of CDS with positive values
plot.legend <- c("Log-normale", "Gamma", "Binomial")
denscomp(list(fit_lognorm, fit_gamma, fit_binom), legendtext = plot.legend)
qqcomp(list(fit_lognorm, fit_gamma, fit_binom), legendtext = plot.legend)

# Although a preliminary marginal fit analysis (using fitdistrplus) indicated that the positive CDS
# scores closely followed a continuous Log-Normal distribution (lowest AIC/BIC and best QQ-plot 
# alignment), this distribution failed to maintain homoscedasticity and proper residual specification
# during conditional regression modeling. Consequently, a Truncated Negative Binomial distribution
# was selected for the count component of the Hurdle model. This discrete distribution successfully
# accounted for the bounded, integer nature of the aggregated scores (1 to 10) and ensured robust,
# homoscedastic residuals across all model specifications.

# Hurdle model
library(glmmTMB)
library(DHARMa)

# Model description :
# First row = intensity model (lognormal)
# zi = binary model (log)
#family = distibution

#### Model 0 = Null Model #####
M0 <- glmmTMB(
  CDS_total_sum ~ 1,
  zi = ~ 1,
  family=truncated_nbinom2(),
  data=df2)

summary(M0)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symtpom at all)
exp(fixef(M0)$zi)
# For the insensity effect : regression on severity of reported symptoms
exp(fixef(M0)$cond)

# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M0 <- simulateResiduals(fittedModel = M0)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M0)


#### Model 1 : CDS ~AGE ####
M1 <- glmmTMB(
  CDS_total_sum ~ AGE_centered,
  zi = ~ AGE_centered,
  family=truncated_nbinom2(),
  data=df2)

summary(M1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M1)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M1)$cond)

# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M1 <- simulateResiduals(fittedModel = M1)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M1)
testDispersion(M1)



#### Model 2 : CDS ~AGE+SEX ####
M2 <- glmmTMB(
  CDS_total_sum ~ AGE_centered+SEXE,
  zi = ~ AGE_centered+SEXE,
  family=truncated_nbinom2(),
  data=df2)

summary(M2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M2)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M2)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M2 <- simulateResiduals(fittedModel = M2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M2)
testDispersion(M2)


#### Model 3 : CDS ~AGE+SEX+ANXIETY #####
M3 <- glmmTMB(
  CDS_total_sum ~ AGE_centered+SEXE+ANXIETE_recoded,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M3)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M3)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M3)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M3 <- simulateResiduals(fittedModel = M3)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M3)
testDispersion(M3)



#### Model 3bis : CDS ~AGE+SEX+ANXIETY+DEPRESSION ####
M3b <- glmmTMB(
  CDS_total_sum ~ AGE_centered+SEXE+DEPRESSION_recoded,
  zi = ~ AGE_centered+SEXE+DEPRESSION_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M3b)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M3b)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M3b)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M3b <- simulateResiduals(fittedModel = M3b)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M3b)
testDispersion(M3b)



#### Model 4 : CDS ~AGE+SEX+ANXIETY+DEPRESSION ####
M4 <- glmmTMB(
  CDS_total_sum ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M4)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M4)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M4)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M4 <- simulateResiduals(fittedModel = M4)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M4)
testDispersion(M4)



#### Model 5 : CDS ~AGE+SEX+ANXIETY+DEPRESSION+MIGRAINE ####
M5 <- glmmTMB(
  CDS_total_sum ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  family=truncated_nbinom2(),
  data=df2)

summary(M5)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M5)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M5)$cond)


# Assumptions verification
# Compute simulated residuals for hurdle model
residus_M5 <- simulateResiduals(fittedModel = M5)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M5)
testDispersion(M5)


#### Model 6 : CDS ~AGE+SEX+ANXIETY+DEPRESSION+SEX*ANXIETY ####
M6 <- glmmTMB(
  CDS_total_sum ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+SEXE*ANXIETE_recoded,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE+SEXE*ANXIETE_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M6)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M6)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M6)$cond)


# Assumptions verification
# Compute simulated residuals for hurdle model
residus_M6 <- simulateResiduals(fittedModel = M6)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M6)
testDispersion(M6)


################## Assumptions for all models #########################################
# For the QQ plot, the observed values perfectly follow the predicted values and the Kolmogorov-Smirnov test of deviation from the theoretical distribution is
# non significant, meaning that the residuals follow the expected distribution.
# The residuals vs predictions graphs indicate that homoscedasticity is respected. This is confirmed by the
# combined adjusted quantile test, that is non significant.
# The outlier test is non significant, confirming the absence of aberrant values.
# The dispersion test is non significant, indicating a correctly predicted dispersion.

######## Models comparison #############
anova(M1, M2, M3, M3b, M4, M5)

# Global summary table of all models (exportation, article ready)
# 1. On stocke le tableau dans une variable nommée "mon_tableau" (Notez le "mon_tableau <- " au début)
mon_tableau <- tab_model(M1, M2, M3, M3b, M4, M5,
                         show.aic = TRUE,
                         show.zeroinf = TRUE,
                         show.intercept = FALSE,
                         transform = "exp",
                         dv.labels = c("Model 1", "Model 2", "Model 3", "Model 3Bis", "Model 4", "Model 5"),
                         pred.labels = c("AGE_centered" = "Age (centered)",
                                         "SEXEf" = "Sex (Woman)",
                                         "ANXIETE_recoded" = "Anxiety (HADS-A)",
                                         "DEPRESSION_recoded" = "Depression (HADS-D)",
                                         "MIGRAINEoui" = "Migraine (Yes)"))


# 2. On remplace le texte d'en-tête dans le compartiment exact "page.complete" de notre variable
mon_tableau$page.complete <- gsub("R<sup>2</sup> / R<sup>2</sup> adjusted", 
                                  "Nagelkerke Pseudo-R<sup>2</sup>", 
                                  mon_tableau$page.complete)

# 3. On exporte manuellement le tableau modifié en fichier Word
writeLines(mon_tableau$page.complete, here("Figures", "tableau_modeles_CDS_tot.doc"))



### ==========================================================
### Forest plot (2 panels) - hurdle model (glmmTMB) - CDS total
### Panel A: Occurrence component (zi, logit link)  -> Odds Ratio
### Panel B: Intensity component (count, log link)  -> IRR
### + significance stars (* p<.05, ** p<.01, *** p<.001)
### ==========================================================

library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

# --------------------------------------------------------------
# 0. Helper: p-value -> stars
# --------------------------------------------------------------
p_to_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# --------------------------------------------------------------
# 1. Extract fixed effects
# --------------------------------------------------------------
tidy_intensity <- broom.mixed::tidy(
  M4, effects = "fixed", component = "cond",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

tidy_occurrence <- broom.mixed::tidy(
  M4, effects = "fixed", component = "zi",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# --------------------------------------------------------------
# 2. Relabel predictors + compute stars
# --------------------------------------------------------------
label_map <- c(
  "AGE_centered"       = "Age (centered)",
  "SEXEf"              = "Sex (female vs. male)",
  "ANXIETE_recoded"    = "Anxiety (HADS-A)",
  "DEPRESSION_recoded" = "Depression (HADS-D)"
)

prepare_data <- function(df) {
  df %>%
    mutate(
      label_text = ifelse(term %in% names(label_map), label_map[term], term),
      stars      = p_to_stars(p.value)
    ) %>%
    mutate(label = factor(label_text, levels = rev(unique(label_text))))
}

tidy_intensity  <- prepare_data(tidy_intensity)
tidy_occurrence <- prepare_data(tidy_occurrence)

# --------------------------------------------------------------
# 3. Forest plot function
# --------------------------------------------------------------
make_forest <- function(data, title, xlab, x_breaks = c(0, 1, 2), x_max = 2.5) {
  
  n_labels <- nlevels(data$label)
  bg_rects <- data.frame(y = seq_len(n_labels)) %>%
    filter(y %% 2 == 1)
  
  ggplot(data, aes(x = estimate, y = label)) +
    geom_rect(
      data = bg_rects,
      aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
      fill = "grey94",
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.15, color = "black", linewidth = 0.6
    ) +
    geom_point(size = 3, shape = 19, color = "black") +
    # Stars: fixed column at the right edge of the plotting area
    geom_text(
      aes(x = x_max - 0.05, label = stars),
      hjust = 1, vjust = 0.75, size = 5, color = "black"
    ) +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 11, hjust = 0),
      axis.text.y        = element_text(color = "black", size = 11),
      axis.text.x        = element_text(color = "black"),
      axis.title.x       = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# --------------------------------------------------------------
# 4. Build the two panels
# --------------------------------------------------------------
p_occurrence <- make_forest(
  tidy_occurrence,
  title = "A. Occurrence component (zero-inflation model)",
  xlab  = "Odds Ratio (95% CI)"
)

p_intensity <- make_forest(
  tidy_intensity,
  title = "B. Intensity component (conditional count model)",
  xlab  = "Incidence Rate Ratio (95% CI)"
)

# --------------------------------------------------------------
# 5. Stack panels, add caption, export
# --------------------------------------------------------------
combined_plot <- p_occurrence / p_intensity +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    caption = "* p < .05, ** p < .01, *** p < .001",
    theme   = theme(plot.caption = element_text(size = 8, hjust = 0))
  )

print(combined_plot)

ggsave(here("Figures", "forest_plot_hurdle_CDS_tot.pdf"), combined_plot, width = 5, height = 5.5, units = "in")
ggsave(here("Figures", "forest_plot_hurdle_CDS_tot.png"), combined_plot, width = 5, height = 5.5, units = "in", dpi = 300)

### ==========================================================
### Forest plot (3 panels) - hurdle model (glmmTMB) - CDS total
### Panel A: Absence of symptoms  P(Y = 0)   (raw zi)      -> OR
### Panel B: Presence of symptoms P(Y > 0)   (inverted zi) -> OR
### Panel C: Intensity           Y | Y > 0  (count)        -> IRR
### + significance stars (* p<.05, ** p<.01, *** p<.001)
### ==========================================================

library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

# --------------------------------------------------------------
# 0. Helper: p-value -> stars
# --------------------------------------------------------------
p_to_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# --------------------------------------------------------------
# 1. Extract fixed effects
# --------------------------------------------------------------
tidy_intensity <- broom.mixed::tidy(
  M4, effects = "fixed", component = "cond",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# Panel A: P(Y = 0), raw zi output
tidy_absence <- broom.mixed::tidy(
  M4, effects = "fixed", component = "zi",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# Panel B: P(Y > 0), multiplicative inversion (1 / OR); p-values unchanged
tidy_presence <- tidy_absence %>%
  mutate(
    conf_low_raw  = conf.low,
    conf_high_raw = conf.high,
    estimate  = 1 / estimate,
    conf.low  = 1 / conf_high_raw,
    conf.high = 1 / conf_low_raw
  )

# --------------------------------------------------------------
# 2. Relabel predictors + compute stars
# --------------------------------------------------------------
label_map <- c(
  "AGE_centered"       = "Age (centered)",
  "SEXEf"              = "Sex (female vs. male)",
  "ANXIETE_recoded"    = "Anxiety (HADS-A)",
  "DEPRESSION_recoded" = "Depression (HADS-D)"
)

prepare_data <- function(df) {
  df %>%
    mutate(
      label_text = ifelse(term %in% names(label_map), label_map[term], term),
      stars      = p_to_stars(p.value)
    ) %>%
    mutate(label = factor(label_text, levels = rev(unique(label_text))))
}

tidy_intensity <- prepare_data(tidy_intensity)
tidy_absence   <- prepare_data(tidy_absence)
tidy_presence  <- prepare_data(tidy_presence)

# --------------------------------------------------------------
# 3. Forest plot function
# --------------------------------------------------------------
make_forest <- function(data, title, xlab, x_breaks = c(0, 1, 2), x_max = 2.5) {
  
  n_labels <- nlevels(data$label)
  bg_rects <- data.frame(y = seq_len(n_labels)) %>%
    filter(y %% 2 == 1)
  
  ggplot(data, aes(x = estimate, y = label)) +
    geom_rect(
      data = bg_rects,
      aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
      fill = "grey94",
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.15, color = "black", linewidth = 0.6
    ) +
    geom_point(size = 3, shape = 19, color = "black") +
    # Stars: fixed column at the right edge of the plotting area
    geom_text(
      aes(x = x_max - 0.05, label = stars),
      hjust = 1, vjust = 0.75, size = 5, color = "black"
    ) +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 10, hjust = 0),
      axis.text.y        = element_text(color = "black", size = 10),
      axis.text.x        = element_text(color = "black", size = 9),
      axis.title.x       = element_text(size = 10),
      axis.line          = element_line(linewidth = 0.3),
      axis.ticks         = element_line(linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# --------------------------------------------------------------
# 4. Build the three panels
# --------------------------------------------------------------
p_absence <- make_forest(
  tidy_absence,
  title = "A. Occurrence component: Absence of symptoms [P(Y = 0)]",
  xlab  = "Odds Ratio (95% CI)"
)

p_presence <- make_forest(
  tidy_presence,
  title = "B. Occurrence component: Presence of symptoms [P(Y > 0)]",
  xlab  = "Odds Ratio (95% CI)"
)

p_intensity <- make_forest(
  tidy_intensity,
  title = "C. Intensity component: Symptom severity [Y | Y > 0]",
  xlab  = "Incidence Rate Ratio (95% CI)"
)

# --------------------------------------------------------------
# 5. Stack panels, add caption, export
# --------------------------------------------------------------
combined_plot <- p_absence / p_presence / p_intensity +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    caption = "* p < .05, ** p < .01, *** p < .001",
    theme   = theme(plot.caption = element_text(size = 8, hjust = 0))
  )

print(combined_plot)

ggsave(here("Figures", "forest_plot_hurdle_CDS_tot_invOR.pdf"), combined_plot, width = 6, height = 7.5, units = "in")
ggsave(here("Figures", "forest_plot_hurdle_CDS_tot_invOR.png"), combined_plot, width = 6, height = 7.5, units = "in", dpi = 300)
# La hauteur a été augmentée à 7.5 pouces pour accueillir le 3e panneau




###################################################################################################################
# CDS Factors from EFA
###################################################################################################################
# Description of CDS total score
invisible(lapply(names(df2[c("CDS_F1", "CDS_F2")]), function(nom) {
  cat("\n===", nom, "===\n")
  x <- df2[[nom]]
  print(summary(x))
  print(describe(x))
  print(shapiro.test(x))
}))
# Both factors are non normal (skewed, leptokurique), overdispersion and autocorrelation of residuals in regression model see test model)


# CDS factors distribution plots
distrib_CDS <- lapply(names(df2[c("CDS_F1", "CDS_F2")]), function(nom) {
  
  col_mean   <- mean(df2[[nom]], na.rm = TRUE)
  col_median <- median(df2[[nom]], na.rm = TRUE)
  
  p <- ggplot(df2, aes(x = .data[[nom]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "purple4", color = "white") +
    geom_density(fill = "grey", alpha = 0.5) +
    geom_vline(xintercept = col_mean, linetype = "dashed", color = "turquoise", linewidth = 1) +
    geom_vline(xintercept = col_median, linetype = "dotdash", color = "red", linewidth = 1) +
    coord_cartesian(xlim = c(0,60)) +
    labs(
      title = paste("Distribution of", nom),
      x = "Score Value",
      y = "Density"
    ) + 
    theme_minimal()
  
  return(p)
})
wrap_plots(distrib_CDS, ncol = 2)
ggsave(here("figures", "distrib_CDS_factors.png"), width = 8, height = 6, dpi = 300)



#########################################################################################
# Hurdle Regression Model for CDS factors
########################################################################################
library(fitdistrplus)

CDS_fac <- c("CDS_F1", "CDS_F2")

plots_pos <- list()
fits      <- list()

for (nom in CDS_fac) {
  
  scores_pos <- df2[[nom]][df2[[nom]] > 0]
  df_pos <- data.frame(scores_pos = scores_pos)
  
  p_pos <- ggplot(data = df_pos, aes(x = scores_pos)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "purple4", color = "white") +
    geom_density(fill = "grey", alpha = 0.5) +
    geom_vline(xintercept = mean(df_pos$scores_pos, na.rm = TRUE), linetype = "dashed", color = "turquoise", linewidth = 1) +
    geom_vline(xintercept = median(df_pos$scores_pos, na.rm = TRUE), linetype = "dotdash", color = "red", linewidth = 1) +
    labs(title = paste("Distribution of", nom, "(scores > 0)")) +
    theme_minimal()
  
  plots_pos[[nom]] <- p_pos
  print(p_pos)
  
  fit_lognorm <- fitdist(scores_pos, "lnorm")
  fit_gamma   <- fitdist(scores_pos, "gamma")
  fit_binom   <- fitdist(scores_pos, "nbinom")
  
  fits[[nom]] <- list(lognorm = fit_lognorm, gamma = fit_gamma, binom = fit_binom)
  
  cat("\n===", nom, ": summary Log-Normal ===\n")
  print(summary(fit_lognorm))
  cat("\n===", nom, ": summary Gamma ===\n")
  print(summary(fit_gamma))
  cat("\n===", nom, ": summary Binomial (nbinom) ===\n")
  print(summary(fit_binom))
  
  cat("\n===", nom, ": gofstat ===\n")
  print(gofstat(list(fit_lognorm, fit_gamma, fit_binom),
                fitnames = c("Log-Normal", "Gamma", "Binomial")))
  
  plot.legend <- c("Log-normale", "Gamma", "Binomial")
  
  # Titre passé directement à la fonction, pas via title() après coup
  denscomp(list(fit_lognorm, fit_gamma, fit_binom),
           legendtext = plot.legend,
           main = paste(nom, "- Density comparison"))
  
  qqcomp(list(fit_lognorm, fit_gamma, fit_binom),
         legendtext = plot.legend,
         main = paste(nom, "- QQ plot"))
}
# Although a preliminary marginal fit analysis (using fitdistrplus) indicated that the positive CDS FACTORS
# scores closely followed a continuous Log-Normal distribution (lowest AIC/BIC and best QQ-plot 
# alignment), this distribution failed to maintain homoscedasticity and proper residual specification
# during conditional regression modeling. Consequently, a Truncated Negative Binomial distribution
# was selected for the count component of the Hurdle model. This discrete distribution successfully
# accounted for the bounded, integer nature of the aggregated scores (1 to 10) and ensured robust,
# homoscedastic residuals across all model specifications.


########################################################################################
# Hierarchical Hurdle Models for FACTOR 1 
########################################################################################
library(glmmTMB)
#### Model 0 = Null Model #####
M0_F1 <- glmmTMB(
  CDS_F1 ~ 1,
  zi = ~ 1,
  family=truncated_nbinom2(),
  data=df2)

summary(M0_F1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symtpom at all)
exp(fixef(M0_F1)$zi)
# For the insensity effect : regression on severity of reported symptoms
exp(fixef(M0_F1)$cond)

# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M0_F1 <- simulateResiduals(fittedModel = M0)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M0_F1)


#### Model 1 : CDS ~AGE ####
M1_F1 <- glmmTMB(
  CDS_F1 ~ AGE_centered,
  zi = ~ AGE_centered,
  family=truncated_nbinom2(),
  data=df2)

summary(M1_F1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M1_F1)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M1_F1)$cond)

# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M1_F1 <- simulateResiduals(fittedModel = M1_F1)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M1_F1)
testDispersion(M1_F1)



#### Model 2 : CDS ~AGE+SEX ####
M2_F1 <- glmmTMB(
  CDS_F1 ~ AGE_centered+SEXE,
  zi = ~ AGE_centered+SEXE,
  family=truncated_nbinom2(),
  data=df2)

summary(M2_F1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M2_F1)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M2_F1)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M2_F1 <- simulateResiduals(fittedModel = M2_F1)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M2_F1)
testDispersion(M2_F1)


#### Model 3 : CDS ~AGE+SEX+ANXIETY #####
M3_F1 <- glmmTMB(
  CDS_F1 ~ AGE_centered+SEXE+ANXIETE_recoded,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M3_F1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M3_F1)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M3_F1)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M3_F1 <- simulateResiduals(fittedModel = M3_F1)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M3_F1)
testDispersion(M3_F1)



#### Model 3bis : CDS ~AGE+SEX+ANXIETY+DEPRESSION ####
M3b_F1 <- glmmTMB(
  CDS_F1 ~ AGE_centered+SEXE+DEPRESSION_recoded,
  zi = ~ AGE_centered+SEXE+DEPRESSION_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M3b_F1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M3b_F1)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M3b_F1)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M3b_F1 <- simulateResiduals(fittedModel = M3b_F1)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M3b_F1)
testDispersion(M3b_F1)



#### Model 4 : CDS ~AGE+SEX+ANXIETY+DEPRESSION ####
M4_F1 <- glmmTMB(
  CDS_F1 ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M4_F1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M4_F1)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M4_F1)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M4_F1 <- simulateResiduals(fittedModel = M4_F1)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M4_F1)
testDispersion(M4_F1)

### ==========================================================
### Forest plot (2 panels) - hurdle model (glmmTMB) - Factor 1 CDS
### Panel A: Occurrence component (zi, logit link)  -> Odds Ratio
### Panel B: Intensity component (count, log link)  -> IRR
### + significance stars (* p<.05, ** p<.01, *** p<.001)
### ==========================================================

library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

# --------------------------------------------------------------
# 0. Helper: p-value -> stars
# --------------------------------------------------------------
p_to_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# --------------------------------------------------------------
# 1. Extract fixed effects
# --------------------------------------------------------------
tidy_intensity <- broom.mixed::tidy(
  M4_F1, effects = "fixed", component = "cond",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

tidy_occurrence <- broom.mixed::tidy(
  M4_F1, effects = "fixed", component = "zi",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# --------------------------------------------------------------
# 2. Relabel predictors + compute stars
# --------------------------------------------------------------
label_map <- c(
  "AGE_centered"       = "Age (centered)",
  "SEXEf"              = "Sex (female vs. male)",
  "ANXIETE_recoded"    = "Anxiety (HADS-A)",
  "DEPRESSION_recoded" = "Depression (HADS-D)"
)

prepare_data <- function(df) {
  df %>%
    mutate(
      label_text = ifelse(term %in% names(label_map), label_map[term], term),
      stars      = p_to_stars(p.value)
    ) %>%
    mutate(label = factor(label_text, levels = rev(unique(label_text))))
}

tidy_intensity  <- prepare_data(tidy_intensity)
tidy_occurrence <- prepare_data(tidy_occurrence)

# --------------------------------------------------------------
# 3. Forest plot function
# --------------------------------------------------------------
make_forest <- function(data, title, xlab, x_breaks = c(0, 1, 2), x_max = 2.5) {
  
  n_labels <- nlevels(data$label)
  bg_rects <- data.frame(y = seq_len(n_labels)) %>%
    filter(y %% 2 == 1)
  
  ggplot(data, aes(x = estimate, y = label)) +
    geom_rect(
      data = bg_rects,
      aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
      fill = "grey94",
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.15, color = "black", linewidth = 0.6
    ) +
    geom_point(size = 3, shape = 19, color = "black") +
    # Stars: fixed column at the right edge of the plotting area
    geom_text(
      aes(x = x_max - 0.05, label = stars),
      hjust = 1, vjust = 0.75, size = 5, color = "black"
    ) +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 11, hjust = 0),
      axis.text.y        = element_text(color = "black", size = 11),
      axis.text.x        = element_text(color = "black"),
      axis.title.x       = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# --------------------------------------------------------------
# 4. Build the two panels
# --------------------------------------------------------------
p_occurrence <- make_forest(
  tidy_occurrence,
  title = "A. Occurrence component (zero-inflation model)",
  xlab  = "Odds Ratio (95% CI)"
)

p_intensity <- make_forest(
  tidy_intensity,
  title = "B. Intensity component (conditional count model)",
  xlab  = "Incidence Rate Ratio (95% CI)"
)

# --------------------------------------------------------------
# 5. Stack panels, add caption, export
# --------------------------------------------------------------
combined_plot <- p_occurrence / p_intensity +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    caption = "* p < .05, ** p < .01, *** p < .001",
    theme   = theme(plot.caption = element_text(size = 8, hjust = 0))
  )

print(combined_plot)

ggsave(here("Figures", "forest_plot_hurdle_CDS_F1.pdf"), combined_plot, width = 5, height = 5.5, units = "in")
ggsave(here("Figures", "forest_plot_hurdle_CDS_F1.png"), combined_plot, width = 5, height = 5.5, units = "in", dpi = 300)

### ==========================================================
### Forest plot (3 panels) - hurdle model (glmmTMB) - Factor 1 CDS
### Panel A: Absence of symptoms  P(Y = 0)   (raw zi)      -> OR
### Panel B: Presence of symptoms P(Y > 0)   (inverted zi) -> OR
### Panel C: Intensity           Y | Y > 0  (count)        -> IRR
### + significance stars (* p<.05, ** p<.01, *** p<.001)
### ==========================================================

library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

# --------------------------------------------------------------
# 0. Helper: p-value -> stars
# --------------------------------------------------------------
p_to_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# --------------------------------------------------------------
# 1. Extract fixed effects
# --------------------------------------------------------------
tidy_intensity <- broom.mixed::tidy(
  M4_F1, effects = "fixed", component = "cond",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# Panel A: P(Y = 0), raw zi output
tidy_absence <- broom.mixed::tidy(
  M4_F1, effects = "fixed", component = "zi",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# Panel B: P(Y > 0), multiplicative inversion (1 / OR); p-values unchanged
tidy_presence <- tidy_absence %>%
  mutate(
    conf_low_raw  = conf.low,
    conf_high_raw = conf.high,
    estimate  = 1 / estimate,
    conf.low  = 1 / conf_high_raw,
    conf.high = 1 / conf_low_raw
  )

# --------------------------------------------------------------
# 2. Relabel predictors + compute stars
# --------------------------------------------------------------
label_map <- c(
  "AGE_centered"       = "Age (centered)",
  "SEXEf"              = "Sex (female vs. male)",
  "ANXIETE_recoded"    = "Anxiety (HADS-A)",
  "DEPRESSION_recoded" = "Depression (HADS-D)"
)

prepare_data <- function(df) {
  df %>%
    mutate(
      label_text = ifelse(term %in% names(label_map), label_map[term], term),
      stars      = p_to_stars(p.value)
    ) %>%
    mutate(label = factor(label_text, levels = rev(unique(label_text))))
}

tidy_intensity <- prepare_data(tidy_intensity)
tidy_absence   <- prepare_data(tidy_absence)
tidy_presence  <- prepare_data(tidy_presence)

# --------------------------------------------------------------
# 3. Forest plot function
# --------------------------------------------------------------
make_forest <- function(data, title, xlab, x_breaks = c(0, 1, 2), x_max = 2.5) {
  
  n_labels <- nlevels(data$label)
  bg_rects <- data.frame(y = seq_len(n_labels)) %>%
    filter(y %% 2 == 1)
  
  ggplot(data, aes(x = estimate, y = label)) +
    geom_rect(
      data = bg_rects,
      aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
      fill = "grey94",
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.15, color = "black", linewidth = 0.6
    ) +
    geom_point(size = 3, shape = 19, color = "black") +
    # Stars: fixed column at the right edge of the plotting area
    geom_text(
      aes(x = x_max - 0.05, label = stars),
      hjust = 1, vjust = 0.75, size = 5, color = "black"
    ) +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 10, hjust = 0),
      axis.text.y        = element_text(color = "black", size = 10),
      axis.text.x        = element_text(color = "black", size = 9),
      axis.title.x       = element_text(size = 10),
      axis.line          = element_line(linewidth = 0.3),
      axis.ticks         = element_line(linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# --------------------------------------------------------------
# 4. Build the three panels
# --------------------------------------------------------------
p_absence <- make_forest(
  tidy_absence,
  title = "A. Occurrence component: Absence of symptoms [P(Y = 0)]",
  xlab  = "Odds Ratio (95% CI)"
)

p_presence <- make_forest(
  tidy_presence,
  title = "B. Occurrence component: Presence of symptoms [P(Y > 0)]",
  xlab  = "Odds Ratio (95% CI)"
)

p_intensity <- make_forest(
  tidy_intensity,
  title = "C. Intensity component: Symptom severity [Y | Y > 0]",
  xlab  = "Incidence Rate Ratio (95% CI)"
)

# --------------------------------------------------------------
# 5. Stack panels, add caption, export
# --------------------------------------------------------------
combined_plot <- p_absence / p_presence / p_intensity +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    caption = "* p < .05, ** p < .01, *** p < .001",
    theme   = theme(plot.caption = element_text(size = 8, hjust = 0))
  )

print(combined_plot)

ggsave(here("Figures", "forest_plot_hurdle_CDS_F1_invOR.pdf"), combined_plot, width = 6, height = 7.5, units = "in")
ggsave(here("Figures", "forest_plot_hurdle_CDS_F1_invOR.png"), combined_plot, width = 6, height = 7.5, units = "in", dpi = 300)



#### Model 5 : CDS ~AGE+SEX+ANXIETY+DEPRESSION+MIGRAINE ####
M5_F1 <- glmmTMB(
  CDS_F1 ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  family=truncated_nbinom2(),
  data=df2)

summary(M5_F1)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M5_F1)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M5_F1)$cond)


# Assumptions verification
# Compute simulated residuals for hurdle model
residus_M5_F1 <- simulateResiduals(fittedModel = M5_F1)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M5_F1)
testDispersion(M5_F1)

######## Models comparison #############
anova(M1_F1, M2_F1, M3_F1, M3b_F1, M4_F1, M5_F1)

# Global summary table of all models (exportation, article ready)
tab_model(M1_F1, M2_F1, M3_F1, M3b_F1, M4_F1, M5_F1,
          show.aic = TRUE,
          show.zeroinf = TRUE,
          show.intercept = FALSE,
          transform = "exp",
          dv.labels = c("Model 1", "Model 2", "Model 3", "Model 3Bis", "Model 4", "Model 5"),
          pred.labels = c("AGE_centered" = "Age (centered)",
                          "SEXEf" = "Sex (Woman)",
                          "ANXIETE_recoded" = "Anxiety (HADS-A)",
                          "DEPRESSION_recoded" = "Depression (HADS-D)",
                          "MIGRAINEoui" = "Migraine (Yes)"),
          file = here("Figures", "tableau_modeles_Factor1.doc"))







########################################################################################
# Hierarchical Hurdle Models for FACTOR 2
########################################################################################
#### Model 0 = Null Model #####
M0_F2 <- glmmTMB(
  CDS_F2 ~ 1,
  zi = ~ 1,
  family=truncated_nbinom2(),
  data=df2)

summary(M0_F2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symtpom at all)
exp(fixef(M0_F2)$zi)
# For the insensity effect : regression on severity of reported symptoms
exp(fixef(M0_F2)$cond)

# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M0_F2 <- simulateResiduals(fittedModel = M0_F2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M0_F2)


#### Model 1 : CDS ~AGE ####
M1_F2 <- glmmTMB(
  CDS_F2 ~ AGE_centered,
  zi = ~ AGE_centered,
  family=truncated_nbinom2(),
  data=df2)

summary(M1_F2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M1_F2)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M1_F2)$cond)

# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M1_F2 <- simulateResiduals(fittedModel = M1_F2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M1_F2)
testDispersion(M1_F2)



#### Model 2 : CDS ~AGE+SEX ####
M2_F2 <- glmmTMB(
  CDS_F2 ~ AGE_centered+SEXE,
  zi = ~ AGE_centered+SEXE,
  family=truncated_nbinom2(),
  data=df2)

summary(M2_F2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M2_F2)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M2_F2)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M2_F2 <- simulateResiduals(fittedModel = M2_F2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M2_F2)
testDispersion(M2_F2)


#### Model 3 : CDS ~AGE+SEX+ANXIETY #####
M3_F2 <- glmmTMB(
  CDS_F2 ~ AGE_centered+SEXE+ANXIETE_recoded,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M3_F2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M3_F2)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M3_F2)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M3_F2 <- simulateResiduals(fittedModel = M3_F2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M3_F2)
testDispersion(M3_F2)



#### Model 3bis : CDS ~AGE+SEX+ANXIETY+DEPRESSION ####
M3b_F2 <- glmmTMB(
  CDS_F2 ~ AGE_centered+SEXE+DEPRESSION_recoded,
  zi = ~ AGE_centered+SEXE+DEPRESSION_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M3b_F2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M3b_F2)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M3b_F2)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M3b_F2 <- simulateResiduals(fittedModel = M3b_F2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M3b_F2)
testDispersion(M3b_F2)



#### Model 4 : CDS ~AGE+SEX+ANXIETY+DEPRESSION ####
M4_F2 <- glmmTMB(
  CDS_F2 ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  family=truncated_nbinom2(),
  data=df2)

summary(M4_F2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M4_F2)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M4_F2)$cond)


# Assumptions verification
# Compute simulated residuals for Hurdle model
residus_M4_F2 <- simulateResiduals(fittedModel = M4_F2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M4_F2)
testDispersion(M4_F2)

### ==========================================================
### Forest plot (2 panels) - hurdle model (glmmTMB) - Factor 2 CDS
### Panel A: Occurrence component (zi, logit link)  -> Odds Ratio
### Panel B: Intensity component (count, log link)  -> IRR
### + significance stars (* p<.05, ** p<.01, *** p<.001)
### ==========================================================

library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

# --------------------------------------------------------------
# 0. Helper: p-value -> stars
# --------------------------------------------------------------
p_to_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# --------------------------------------------------------------
# 1. Extract fixed effects
# --------------------------------------------------------------
tidy_intensity <- broom.mixed::tidy(
  M4_F2, effects = "fixed", component = "cond",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

tidy_occurrence <- broom.mixed::tidy(
  M4_F2, effects = "fixed", component = "zi",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# --------------------------------------------------------------
# 2. Relabel predictors + compute stars
# --------------------------------------------------------------
label_map <- c(
  "AGE_centered"       = "Age (centered)",
  "SEXEf"              = "Sex (female vs. male)",
  "ANXIETE_recoded"    = "Anxiety (HADS-A)",
  "DEPRESSION_recoded" = "Depression (HADS-D)"
)

prepare_data <- function(df) {
  df %>%
    mutate(
      label_text = ifelse(term %in% names(label_map), label_map[term], term),
      stars      = p_to_stars(p.value)
    ) %>%
    mutate(label = factor(label_text, levels = rev(unique(label_text))))
}

tidy_intensity  <- prepare_data(tidy_intensity)
tidy_occurrence <- prepare_data(tidy_occurrence)

# --------------------------------------------------------------
# 3. Forest plot function
# --------------------------------------------------------------
make_forest <- function(data, title, xlab, x_breaks = c(0, 1, 2), x_max = 2.5) {
  
  n_labels <- nlevels(data$label)
  bg_rects <- data.frame(y = seq_len(n_labels)) %>%
    filter(y %% 2 == 1)
  
  ggplot(data, aes(x = estimate, y = label)) +
    geom_rect(
      data = bg_rects,
      aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
      fill = "grey94",
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.15, color = "black", linewidth = 0.6
    ) +
    geom_point(size = 3, shape = 19, color = "black") +
    # Stars: fixed column at the right edge of the plotting area
    geom_text(
      aes(x = x_max - 0.05, label = stars),
      hjust = 1, vjust = 0.75, size = 5, color = "black"
    ) +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 11, hjust = 0),
      axis.text.y        = element_text(color = "black", size = 11),
      axis.text.x        = element_text(color = "black"),
      axis.title.x       = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# --------------------------------------------------------------
# 4. Build the two panels
# --------------------------------------------------------------
p_occurrence <- make_forest(
  tidy_occurrence,
  title = "A. Occurrence component (zero-inflation model)",
  xlab  = "Odds Ratio (95% CI)"
)

p_intensity <- make_forest(
  tidy_intensity,
  title = "B. Intensity component (conditional count model)",
  xlab  = "Incidence Rate Ratio (95% CI)"
)

# --------------------------------------------------------------
# 5. Stack panels, add caption, export
# --------------------------------------------------------------
combined_plot <- p_occurrence / p_intensity +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    caption = "* p < .05, ** p < .01, *** p < .001",
    theme   = theme(plot.caption = element_text(size = 8, hjust = 0))
  )

print(combined_plot)

ggsave(here("Figures", "forest_plot_hurdle_CDS_F2.pdf"), combined_plot, width = 5, height = 5.5, units = "in")
ggsave(here("Figures", "forest_plot_hurdle_CDS_F2.png"), combined_plot, width = 5, height = 5.5, units = "in", dpi = 300)

### ==========================================================
### Forest plot (3 panels) - hurdle model (glmmTMB) - Factor 2 CDS
### Panel A: Absence of symptoms  P(Y = 0)   (raw zi)      -> OR
### Panel B: Presence of symptoms P(Y > 0)   (inverted zi) -> OR
### Panel C: Intensity           Y | Y > 0  (count)        -> IRR
### + significance stars (* p<.05, ** p<.01, *** p<.001)
### ==========================================================

library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

# --------------------------------------------------------------
# 0. Helper: p-value -> stars
# --------------------------------------------------------------
p_to_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# --------------------------------------------------------------
# 1. Extract fixed effects
# --------------------------------------------------------------
tidy_intensity <- broom.mixed::tidy(
  M4_F2, effects = "fixed", component = "cond",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# Panel A: P(Y = 0), raw zi output
tidy_absence <- broom.mixed::tidy(
  M4_F2, effects = "fixed", component = "zi",
  conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# Panel B: P(Y > 0), multiplicative inversion (1 / OR); p-values unchanged
tidy_presence <- tidy_absence %>%
  mutate(
    conf_low_raw  = conf.low,
    conf_high_raw = conf.high,
    estimate  = 1 / estimate,
    conf.low  = 1 / conf_high_raw,
    conf.high = 1 / conf_low_raw
  )

# --------------------------------------------------------------
# 2. Relabel predictors + compute stars
# --------------------------------------------------------------
label_map <- c(
  "AGE_centered"       = "Age (centered)",
  "SEXEf"              = "Sex (female vs. male)",
  "ANXIETE_recoded"    = "Anxiety (HADS-A)",
  "DEPRESSION_recoded" = "Depression (HADS-D)"
)

prepare_data <- function(df) {
  df %>%
    mutate(
      label_text = ifelse(term %in% names(label_map), label_map[term], term),
      stars      = p_to_stars(p.value)
    ) %>%
    mutate(label = factor(label_text, levels = rev(unique(label_text))))
}

tidy_intensity <- prepare_data(tidy_intensity)
tidy_absence   <- prepare_data(tidy_absence)
tidy_presence  <- prepare_data(tidy_presence)

# --------------------------------------------------------------
# 3. Forest plot function
# --------------------------------------------------------------
make_forest <- function(data, title, xlab, x_breaks = c(0, 1, 2), x_max = 2.5) {
  
  n_labels <- nlevels(data$label)
  bg_rects <- data.frame(y = seq_len(n_labels)) %>%
    filter(y %% 2 == 1)
  
  ggplot(data, aes(x = estimate, y = label)) +
    geom_rect(
      data = bg_rects,
      aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
      fill = "grey94",
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.15, color = "black", linewidth = 0.6
    ) +
    geom_point(size = 3, shape = 19, color = "black") +
    # Stars: fixed column at the right edge of the plotting area
    geom_text(
      aes(x = x_max - 0.05, label = stars),
      hjust = 1, vjust = 0.75, size = 5, color = "black"
    ) +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 10, hjust = 0),
      axis.text.y        = element_text(color = "black", size = 10),
      axis.text.x        = element_text(color = "black", size = 9),
      axis.title.x       = element_text(size = 10),
      axis.line          = element_line(linewidth = 0.3),
      axis.ticks         = element_line(linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# --------------------------------------------------------------
# 4. Build the three panels
# --------------------------------------------------------------
p_absence <- make_forest(
  tidy_absence,
  title = "A. Occurrence component: Absence of symptoms [P(Y = 0)]",
  xlab  = "Odds Ratio (95% CI)"
)

p_presence <- make_forest(
  tidy_presence,
  title = "B. Occurrence component: Presence of symptoms [P(Y > 0)]",
  xlab  = "Odds Ratio (95% CI)"
)

p_intensity <- make_forest(
  tidy_intensity,
  title = "C. Intensity component: Symptom severity [Y | Y > 0]",
  xlab  = "Incidence Rate Ratio (95% CI)"
)

# --------------------------------------------------------------
# 5. Stack panels, add caption, export
# --------------------------------------------------------------
combined_plot <- p_absence / p_presence / p_intensity +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    caption = "* p < .05, ** p < .01, *** p < .001",
    theme   = theme(plot.caption = element_text(size = 8, hjust = 0))
  )

print(combined_plot)

ggsave(here("Figures", "forest_plot_hurdle_CDS_F2_invOR.pdf"), combined_plot, width = 6, height = 7.5, units = "in")
ggsave(here("Figures", "forest_plot_hurdle_CDS_F2_invOR.png"), combined_plot, width = 6, height = 7.5, units = "in", dpi = 300)

#### Model 5 : CDS ~AGE+SEX+ANXIETY+DEPRESSION+MIGRAINE ####
M5_F2 <- glmmTMB(
  CDS_F2 ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  zi = ~ AGE_centered+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  family=truncated_nbinom2(),
  data=df2)

summary(M5_F2)

# Exponential transformation of coefficients to interpret them as OR
# For the binary effect (probability to have 0 = no symptom at all)
exp(fixef(M5_F2)$zi)
# For the intensity effect : regression on severity of reported symptoms
exp(fixef(M5_F2)$cond)


# Assumptions verification
# Compute simulated residuals for hurdle model
residus_M5_F2 <- simulateResiduals(fittedModel = M5_F2)

# Assumptions plot (distribution, homoscedasticity, dispersion)
plot(residus_M5_F2)
testDispersion(M5_F2)



######## Models comparison #############
anova(M1_F2, M2_F2, M3_F2, M3b_F2, M4_F2, M5_F2)

# Global summary table of all models (exportation, article ready)
tab_model(M1_F2, M2_F2, M3_F2, M3b_F2, M4_F2, M5_F2,
          show.aic = TRUE,
          show.zeroinf = TRUE,
          show.intercept = FALSE,
          transform = "exp",
          dv.labels = c("Model 1", "Model 2", "Model 3", "Model 3Bis", "Model 4", "Model 5"),
          pred.labels = c("AGE_centered" = "Age (centered)",
                          "SEXEf" = "Sex (Woman)",
                          "ANXIETE_recoded" = "Anxiety (HADS-A)",
                          "DEPRESSION_recoded" = "Depression (HADS-D)",
                          "MIGRAINEoui" = "Migraine (Yes)"),
          file = here("Figures", "tableau_modeles_Factor2.doc"))



#####################################################################
# Figure combinée : forest plots invOR, CDS total + Facteur 1 + Facteur 2
# Colonnes = scores ; lignes = Absence P(Y=0) / Présence P(Y>0) / Intensité
#####################################################################
library(broom.mixed)

p_to_stars <- function(p) {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

label_map <- c(
  "AGE_centered"       = "Age (centered)",
  "SEXEf"              = "Sex (female vs. male)",
  "ANXIETE_recoded"    = "Anxiety (HADS-A)",
  "DEPRESSION_recoded" = "Depression (HADS-D)"
)

comp_labels <- c(
  "A" = "A. Absence of symptoms\nP(Y = 0) - OR",
  "B" = "B. Presence of symptoms\nP(Y > 0) - OR",
  "C" = "C. Symptom severity\nY | Y > 0 - IRR"
)

# Extrait les 3 composantes d'un modèle hurdle dans un seul tableau
extract_hurdle <- function(model, score_name) {
  
  intensity <- broom.mixed::tidy(model, effects = "fixed", component = "cond",
                                 conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(component = "C")
  
  absence <- broom.mixed::tidy(model, effects = "fixed", component = "zi",
                               conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(component = "A")
  
  # Inversion multiplicative : P(Y > 0) ; p-values inchangées
  presence <- absence %>%
    mutate(component = "B",
           lo = 1 / conf.high,
           hi = 1 / conf.low,
           estimate = 1 / estimate,
           conf.low = lo,
           conf.high = hi) %>%
    dplyr::select(-lo, -hi)
  
  bind_rows(absence, presence, intensity) %>%
    filter(term != "(Intercept)") %>%
    mutate(score = score_name)
}

all_res <- bind_rows(
  extract_hurdle(M4,    "CDS total score"),
  extract_hurdle(M4_F1, "CDS Factor 1"),
  extract_hurdle(M4_F2, "CDS Factor 2")
) %>%
  mutate(
    stars     = p_to_stars(p.value),
    label     = factor(label_map[term], levels = rev(unname(label_map))),
    component = factor(component, levels = names(comp_labels), labels = comp_labels),
    score     = factor(score, levels = c("CDS total score", "CDS Factor 1", "CDS Factor 2"))
  )

x_max <- 2.5   # à augmenter si des IC dépassent (ils seraient tronqués)
# Bandes alternées : une ligne sur deux en gris (positions 1 et 3 de l'axe y)
n_labels <- nlevels(all_res$label)
bg_rects <- data.frame(y = seq_len(n_labels)) %>%
  filter(y %% 2 == 1)
forest_all <- ggplot(all_res, aes(x = estimate, y = label)) +
  # Fond alterné gris / blanc (dessiné en premier = derrière)
  geom_rect(
    data = bg_rects,
    aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
    fill = "grey94",
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.15, color = "black", linewidth = 0.6) +
  geom_point(size = 2.5, color = "black") +
  geom_text(aes(x = x_max - 0.05, label = stars),
            hjust = 1, vjust = 0.75, size = 5, color = "black") +
  facet_grid(component ~ score) +
  scale_x_continuous(breaks = c(0, 1, 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, x_max)) +
  labs(x = "Odds Ratio / Incidence Rate Ratio (95% CI)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.text.x       = element_text(face = "bold", size = 11),
    strip.text.y       = element_text(face = "bold", size = 9, angle = 270),
    strip.background   = element_rect(fill = "grey92", color = "grey60", linewidth = 0.3),
    axis.text          = element_text(color = "black", size = 9),
    panel.spacing      = unit(0.8, "lines"),
    panel.border       = element_rect(color = "grey60", fill = NA, linewidth = 0.3),
    plot.caption       = element_text(size = 8, hjust = 0)
  )

forest_all

ggsave(here("Figures", "forest_plot_hurdle_CDS_all_invOR.pdf"), forest_all, width = 10, height = 7.5, units = "in")
ggsave(here("Figures", "forest_plot_hurdle_CDS_all_invOR.png"), forest_all, width = 10, height = 7.5, units = "in", dpi = 300)



################################################################################################################
# SEM CDS factors WITH MEASUREMENT MODEL
################################################################################################################
library(lavaan)

# Check multivariate normality on CDS items
semTools::mardiaSkew(items)
semTools::mardiaKurtosis(items)
# not respected, justifies the use of MLR estimator

# Recoding of variables to use in model syntax
df2$Depression <- df2$DEPRESSION_recoded
df2$Anxiety <- df2$ANXIETE_recoded
df2$SEX <- df2$SEXE
df2$CDS_total <- df2$CDS_total_sum

#### SEM Model definition
SEM_CDS_factors <- '

Fac1 =~ CDS1+CDS2+CDS3+CDS6+CDS8+CDS10+CDS13+CDS15+CDS23+CDS24+CDS26
Fac2 =~ CDS5+CDS7+CDS9+CDS25+CDS28+CDS29

Depression ~~ Anxiety
Fac1 ~~ Fac2
CDS7 ~~ CDS25
CDS7 ~~ CDS28
CDS25 ~~ CDS28

Fac1 ~ bAD*Depression + bAA*Anxiety + AGE_centered + cA*SEX
Fac2 ~ bBD*Depression + bBA*Anxiety + AGE_centered + cB*SEX
Depression ~ a1*SEX
Anxiety ~ a2*SEX

bADbis:=bAD
bAAbis:=bAA
indAD:=bAD*a1
indAA:=bAA*a2
bBDbis:=bBD
bBAbis:=bBA
indBD:=bBD*a1
indBA:=bBA*a2
totA:=indAD+cA+indAA
totB:=indBD+cB+indBA
percMedAD:=indAD/totA
percMedAA:=indAA/totA
percMedBD:=indBD/totB
percMedBA:=indBA/totB
'

#### Results without boostrap but robust estimation
# For overview AND FIT INDICES (bootstrap does not correct for high multivariate non normality on fit indices)
fit_CDS_fac <- sem(SEM_CDS_factors, data=df2, missing="fiml", estimator="MLR")
summary(fit_CDS_fac, standardized = TRUE, fit.measures = TRUE, ci = TRUE)
fitMeasures(fit_CDS_fac, c("cfi.robust", "gfi.robust", "tli.robust", "rmsea.robust", "srmr", "nfi", "aic", "bic"))
lavInspect(fit_CDS_fac,"cor.lv")
A=standardizedSolution(fit_CDS_fac,type = "std.lv")
modindices(fit_CDS_fac, minimum.value = 10, sort = TRUE)

#### Bootsraped results (+ robust, compensate missing data) for more robust estimation of mediation CI
fit_CDS_fac_boot <- sem(SEM_CDS_factors, data=df2, missing="fiml", estimator="ML", se="bootstrap", bootstrap=5000)
summary(fit_CDS_fac_boot, standardized = TRUE, fit.measures = TRUE, ci = TRUE)


#### Export results as table
library(officer)
library(flextable)
# Extraction of standardized estimate
tab_SEM <- standardizedsolution(fit_CDS_fac_boot, type = "std.lv", ci = TRUE)

# Filter standardized loadings (std.lv)
loadings_std <- tab_SEM %>%
  mutate(across(c(est.std, se, ci.lower, ci.upper, z, pvalue), ~round(., 3)))
print(loadings_std)

# Export to Word
ft <- flextable(loadings_std)
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = here("Figures", "CDS_pop_gen_SEM_loadings_std.docx"))

# Export of fit measures table
fit_indices <- fitMeasures(fit_CDS_fac_boot, c("chisq", "df", "pvalue", 
                                          "cfi.robust", "gfi.robust", "tli.robust", 
                                          "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", 
                                          "srmr", "aic", "bic"))

fit_table <- data.frame(
  Model = "SEM_CDS_factors",
  `Chi2 (df)` = paste0(round(fit_indices["chisq"], 2), " (", fit_indices["df"], ")"),
  p = round(fit_indices["pvalue"], 3),
  CFI = round(fit_indices["cfi.robust"], 3),
  GFI = round(fit_indices["gfi.robust"], 3),
  TLI = round(fit_indices["tli.robust"], 3),
  RMSEA = paste0(round(fit_indices["rmsea.robust"], 3), " (", round(fit_indices["rmsea.ci.lower.robust"], 3), "-", round(fit_indices["rmsea.ci.upper.robust"], 3), ")"),
  SRMR = round(fit_indices["srmr"], 3),
  AIC = round(fit_indices["aic"], 1),
  BIC = round(fit_indices["bic"], 1)
)

# Transform to Word table
ft <- flextable(fit_table)
ft <- autofit(ft)

# Export to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = here("Figures", "CDS_pop_gen_SEM_Fit_Indices_boot.docx"))

#################################################################################################################
# Simple Mediation  Model CDS total score
#################################################################################################################
# Mediation model definition
med_CDS_tot <- 
        ' # direct effect
             CDS_total ~ c*SEX + AGE_centered
           # mediators
             Anxiety ~ aA*SEX
             Depression ~ aD*SEX
             CDS_total ~ bA*Anxiety + bD*Depression
           # indirect effects
             abA := aA*bA
             abD := aD*bD
           # total effect
             total := c + (aA*bA) + (aD*bD)
             # covariance
             Anxiety ~~ Depression
         '
#### Fitting of mediation model
fit_med <- sem(med_CDS_tot, data=df2, estimator="MLR")
summary(fit_med, standardized=TRUE, fit.measures=TRUE, ci=TRUE)


#### Export results as table
# Extraction of standardized estimate
tab_SEM <- standardizedsolution(fit_med, type = "std.all", ci = TRUE)

# Filter standardized loadings (std.all)
loadings_std <- tab_SEM %>%
  mutate(across(c(est.std, se, ci.lower, ci.upper, z, pvalue), ~round(., 3)))
print(loadings_std)

# Export to Word
ft <- flextable(loadings_std)
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = here("Figures", "CDS_pop_gen_Mediation_loadings_std.docx"))




### ==========================================================
### CDS item-level figure
### Left : % participants answering No / Yes for each item
### Right: mean item score (+/- SD) by group
### ==========================================================

library(tidyverse)
library(patchwork)
library(here)

# --------------------------------------------------------------
# PARAMETRES A ADAPTER
# --------------------------------------------------------------
data        <- df2                 # data.frame en format large (1 ligne = 1 participant)
items       <- paste0("CDS", 1:29) # colonnes des 29 items (score par item, 0-10)
group_var   <- "GROUP"             # nom de la colonne de groupe
lvl_ctrl    <- "Controls"          # modalite du groupe controle
lvl_dd      <- "DD patients"       # modalite du groupe patients
group_left  <- lvl_dd              # groupe utilise pour les % No/Yes (gauche)
yes_rule    <- function(x) x > 0   # definition de "Yes" (ici: score > 0)
item_order  <- NULL                # ex: c(5,1,17,13,...) ; NULL = tri par score moyen des patients

# --------------------------------------------------------------
# Libelles des items
# --------------------------------------------------------------
item_text <- c(
  "Feeling unreal or cut-off from the world",
  "Things look flat",
  "Body feels as if it did not belong to self",
  "Absence of fear in distressing situations",
  "Favourite activities no longer enjoyable",
  "Feeling of being a detached observer of oneself",
  "Flavour of meals no longer gives a feeling of pleasure or distaste",
  "Body feels very light",
  "No emotion felt when weeping or laughing",
  "Feeling of not having any thoughts at all",
  "Own voice sounds remote and unreal",
  "Feeling that hands and feet have become larger or smaller",
  "Surroundings feel detached or unreal",
  "Recently done things feel as if they had taken place long time ago",
  "Seeing self outside",
  "Feeling detached from personal memories",
  "When in a new situation, feeling as if it had happened before",
  "Unable to feel affection towards family and friends",
  "Objects look smaller or further away",
  "Unable to feel properly things touched with hands",
  "Unable to picture things in mind",
  "Feeling detached from pain",
  "Feeling of being outside the body",
  "Feeling mechanical and robotic when moving",
  "Smell of things no longer gives feeling of pleasure or dislike",
  "Thoughts seem to have a life of their own",
  "Urge to touch oneself to be reassured of body existence",
  "Unable to feel hunger or thirst",
  "Familiar places look strange and unfamiliar"
)
labels_df <- tibble(
  item_num = 1:29,
  label_txt = paste0("Q", 1:29, ". ", item_text)
)

# --------------------------------------------------------------
# Donnees en format long
# --------------------------------------------------------------
df_long <- data %>%
  select(all_of(c(group_var, items))) %>%
  rename(group = all_of(group_var)) %>%
  pivot_longer(all_of(items), names_to = "item", values_to = "score") %>%
  mutate(
    item_num = as.integer(str_remove(item, "CDS")),
    group    = factor(group, levels = c(lvl_ctrl, lvl_dd))
  ) %>%
  filter(!is.na(score), !is.na(group))

# Statistiques descriptives par item x groupe
stats <- df_long %>%
  group_by(item_num, group) %>%
  summarise(mean = mean(score), sd = sd(score), .groups = "drop")

# Pourcentage No / Yes (groupe de gauche)
pct <- df_long %>%
  filter(group == group_left) %>%
  group_by(item_num) %>%
  summarise(yes = mean(yes_rule(score)) * 100, .groups = "drop") %>%
  mutate(no = 100 - yes) %>%
  pivot_longer(c(no, yes), names_to = "answer", values_to = "pct") %>%
  mutate(
    answer     = factor(answer, levels = c("no", "yes"), labels = c("No", "Yes")),
    pct_signed = ifelse(answer == "No", -pct, pct)
  )

# Ordre des items (premier = en haut)
if (is.null(item_order)) {
  item_order <- stats %>%
    filter(group == lvl_dd) %>%
    arrange(desc(mean)) %>%
    pull(item_num)
}
lab_levels <- labels_df %>%
  slice(match(rev(item_order), item_num)) %>%
  pull(label_txt)

add_label <- function(d) {
  d %>%
    left_join(labels_df, by = "item_num") %>%
    mutate(label = factor(label_txt, levels = lab_levels))
}
pct   <- add_label(pct)
stats <- add_label(stats)

# --------------------------------------------------------------
# Panneau gauche : % participants
# --------------------------------------------------------------
col_no  <- "#8B4513"
col_yes <- "#3A8A1F"

p_left <- ggplot(pct, aes(x = pct_signed, y = label, fill = answer)) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = "#FFFFE0") +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(No = col_no, Yes = col_yes), name = NULL) +
  scale_x_continuous(
    limits = c(-100, 100), breaks = seq(-100, 100, 25),
    labels = abs, position = "top", expand = c(0, 0)
  ) +
  labs(x = "% participants", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y     = element_text(color = "black", size = 8),
    axis.text.x     = element_text(color = "black", size = 8),
    axis.title.x    = element_text(size = 9),
    legend.position = "bottom"
  )

# --------------------------------------------------------------
# Panneau droit : score moyen +/- SD
# --------------------------------------------------------------
dodge <- position_dodge(width = 0.5)

p_right <- ggplot(stats, aes(x = mean, y = label, color = group)) +
  geom_errorbarh(aes(xmin = pmax(mean - sd, 0), xmax = mean + sd),
                 height = 0, linewidth = 0.4, position = dodge) +
  geom_point(size = 1.8, position = dodge) +
  scale_color_manual(values = c("#0033FF", "#FF9500") %>% setNames(c(lvl_ctrl, lvl_dd)),
                     name = NULL) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2),
                     position = "top", expand = c(0, 0)) +
  labs(x = "Average score", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.line.y     = element_blank(),
    axis.text.x     = element_text(color = "black", size = 8),
    axis.title.x    = element_text(size = 9),
    legend.position = "bottom"
  )

# --------------------------------------------------------------
# Assemblage et export
# --------------------------------------------------------------
fig <- p_left + p_right + plot_layout(widths = c(1.7, 1))
print(fig)

ggsave(here("Figures", "CDS_items_figure.pdf"), fig, width = 10, height = 7.5, units = "in")
ggsave(here("Figures", "CDS_items_figure.png"), fig, width = 10, height = 7.5, units = "in", dpi = 300)