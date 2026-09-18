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




###################################################################################################################
# Adding CDS factors to data table
##################################################################################################################
df$CDS_F1 <- rowSums(df[, c("CDS1","CDS2","CDS3","CDS6","CDS8", "CDS10", "CDS13", "CDS15", "CDS23", "CDS24", "CDS26")])
df$CDS_F2 <- rowSums(df[, c("CDS5","CDS7","CDS9","CDS25","CDS28", "CDS29")])

# We only sum the items here to have the factors scores to keep a distribution with integers that allows us for
# a binomial fit of the Hurdle model. If we do a mean (divide by item number), the score becomes continuous and the
# binomial distribution would not apply anymore, thus causing discrepancy between our analysis types.
# We can create a mean factor score to compare both factors in controls using descriptive stats.



###################################################################################################################
# Linear Regression Models #
###################################################################################################################
# Filter to remove cases without complete anxiety, depression and migraine
anxdep <- c("ANXIETE_recoded", "DEPRESSION_recoded", "MIGRAINE")
df2 <- df[complete.cases(df[, anxdep]), ]
nrow(df2)

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
  CDS_total_sum ~ AGE,
  zi = ~ AGE,
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
  CDS_total_sum ~ AGE+SEXE,
  zi = ~ AGE+SEXE,
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
  CDS_total_sum ~ AGE+SEXE+ANXIETE_recoded,
  zi = ~ AGE+SEXE+ANXIETE_recoded,
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
  CDS_total_sum ~ AGE+SEXE+DEPRESSION_recoded,
  zi = ~ AGE+SEXE+DEPRESSION_recoded,
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
  CDS_total_sum ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  zi = ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
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
  CDS_total_sum ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  zi = ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
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
  CDS_total_sum ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+SEXE*ANXIETE_recoded,
  zi = ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE+SEXE*ANXIETE_recoded,
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
tab_model(M1, M2, M3, M3b, M4, M5,
          show.aic = TRUE,
          show.zeroinf = TRUE,
          show.intercept = FALSE,
          transform = "exp",
          dv.labels = c("Model 1", "Model 2", "Model 3", "Model 3Bis", "Model 4", "Model 5"),
          pred.labels = c("AGE" = "Age",
                          "SEXEf" = "Sex (Woman)",
                          "ANXIETE_recoded" = "Anxiety (HADS-A)",
                          "DEPRESSION_recoded" = "Depression (HADS-D)",
                          "MIGRAINEoui" = "Migraine (Yes)"),
          file = here("Figures", "tableau_modeles.doc"))




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
  CDS_F1 ~ AGE,
  zi = ~ AGE,
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
  CDS_F1 ~ AGE+SEXE,
  zi = ~ AGE+SEXE,
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
  CDS_F1 ~ AGE+SEXE+ANXIETE_recoded,
  zi = ~ AGE+SEXE+ANXIETE_recoded,
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
  CDS_F1 ~ AGE+SEXE+DEPRESSION_recoded,
  zi = ~ AGE+SEXE+DEPRESSION_recoded,
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
  CDS_F1 ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  zi = ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
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


#### Model 5 : CDS ~AGE+SEX+ANXIETY+DEPRESSION+MIGRAINE ####
M5_F1 <- glmmTMB(
  CDS_F1 ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  zi = ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
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
          pred.labels = c("AGE" = "Age",
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
  CDS_F2 ~ AGE,
  zi = ~ AGE,
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
  CDS_F2 ~ AGE+SEXE,
  zi = ~ AGE+SEXE,
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
  CDS_F2 ~ AGE+SEXE+ANXIETE_recoded,
  zi = ~ AGE+SEXE+ANXIETE_recoded,
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
  CDS_F2 ~ AGE+SEXE+DEPRESSION_recoded,
  zi = ~ AGE+SEXE+DEPRESSION_recoded,
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
  CDS_F2 ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
  zi = ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded,
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


#### Model 5 : CDS ~AGE+SEX+ANXIETY+DEPRESSION+MIGRAINE ####
M5_F2 <- glmmTMB(
  CDS_F2 ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
  zi = ~ AGE+SEXE+ANXIETE_recoded+DEPRESSION_recoded+MIGRAINE,
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
          pred.labels = c("AGE" = "Age",
                          "SEXEf" = "Sex (Woman)",
                          "ANXIETE_recoded" = "Anxiety (HADS-A)",
                          "DEPRESSION_recoded" = "Depression (HADS-D)",
                          "MIGRAINEoui" = "Migraine (Yes)"),
          file = here("Figures", "tableau_modeles_Factor2.doc"))





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

Fac1 ~ bAD*Depression + bAA*Anxiety + cA*AGE + SEX
Fac2 ~ bBD*Depression + bBA*Anxiety + cB*AGE + SEX
Depression ~ a1*SEX + AGE
Anxiety ~ a2*SEX + AGE

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

#### Bootsraped results (+ robust, compensate missing data)
fit_CDS_fac <- sem(SEM_CDS_factors, data=df2, missing="fiml", estimator="MLR", se="boostrap", boostrap=5000)

#### Results without boostrap for quick look 
fit_CDS_fac <- sem(SEM_CDS_factors, data=df2, missing="fiml", estimator="MLR")
summary(fit_CDS_fac, standardized = TRUE, fit.measures = TRUE, ci = TRUE)
fitMeasures(fit_CDS_fac, c("cfi", "gfi", "tli", "rmsea", "srmr", "nfi", "aic", "bic"))
lavInspect(fit_CDS_fac,"cor.lv")
A=standardizedSolution(fit_CDS_fac,type = "std.lv")


#### Export results as table
library(officer)
library(flextable)
# Extraction of standardized estimate
tab_SEM <- standardizedsolution(fit_CDS_fac, type = "std.lv", ci = TRUE)

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
fit_indices <- fitMeasures(fit_CDS_fac, c("chisq","df","pvalue","cfi","gfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr","aic","bic"))

fit_table <- data.frame(
  Model = "SEM_CDS_factors_general_population",
  `Chi2 (df)` = paste0(round(fit_indices["chisq"], 2), " (", fit_indices["df"], ")"),
  p = round(fit_indices["pvalue"], 3),
  CFI = round(fit_indices["cfi"], 3),
  GFI = round(fit_indices["gfi"], 3),
  TLI = round(fit_indices["tli"], 3),
  RMSEA = paste0(round(fit_indices["rmsea"], 3), " (", round(fit_indices["rmsea.ci.lower"], 3), "-", round(fit_indices["rmsea.ci.upper"], 3), ")"),
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
print(doc, target = here("Figures", "CDS_pop_gen_SEM_Fit_Indices.docx"))




#################################################################################################################
# Simple Mediation  Model CDS total score
#################################################################################################################
# Mediation model definition
med_CDS_tot <- 
        ' # direct effect
             CDS_total ~ c*SEX
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