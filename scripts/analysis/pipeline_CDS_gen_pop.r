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

######### Charge and visualize dataset ##########

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
########## Preprocessing ##########
# NaN inspection

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


####### Identify aberrant values (> total score) ######
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



############ See variables distribution ##################

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




######### Outliers identification #########
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

# Number of categories used in each itemm
n_categories <- sapply(items, function(x) length(unique(na.omit(x))))
print(n_categories)

min_nonzero_n <- sapply(items, function(x) {
  tab <- table(x)
  tab <- tab[names(tab) != "0"]
  if (length(tab) == 0) return(NA)
  min(tab)
})
print(sort(min_nonzero_n))

# Items to flag : >90% zeros AND less than 10 observations in the smallest non null category
items_flag <- names(items)[pct_zero > 90 & min_nonzero_n < 10]
cat("Items to heck :\n")
print(items_flag)

# Bivariate crossed tables for flagged items
if (length(items_flag) >= 2) {
  combn(items_flag, 2, FUN = function(pair) {
    cat("\n---", pair[1], "x", pair[2], "---\n")
    print(table(items[[pair[1]]], items[[pair[2]]]))
  })
}


##### Items binarization #####
# Justification: Almost all items exhibit a minimum sample size per active category 
# of less than 10-15 (cf. preliminary diagnostics), making multi-category polychoric 
# estimation unstable for the majority of the scale. Binarization (0 = absence, 
# 1 = presence, regardless of any score > 0) eliminates this issue. 
# The corresponding correlations will be tetrachoric (a special case of the polychoric 
# correlation for binary variables), which corrects for the marginal bias (dependence 
# on prevalence rates) that affects the phi/Pearson coefficient on binary data with 
# heterogeneous prevalence.

recode_binary <- function(x) {
  ifelse(x == 0, 0, 1)
}

# Numeric version (0/1) : for tetrachoric(), KMO(), fa.parallel()
items_recoded_num <- as.data.frame(lapply(items, recode_binary))
names(items_recoded_num) <- item_names

# Factor version : for lavaan::efa()
items_recoded <- as.data.frame(lapply(items_recoded_num, function(x) factor(x, ordered = TRUE)))
names(items_recoded) <- item_names

# Post Binarization Check : Prevalence rate per item
prevalence <- sapply(items_recoded_num, function(x) mean(x == 1, na.rm = TRUE) * 100)
prevalence_df <- data.frame(item = names(prevalence), pct_present = round(prevalence, 1)) |>
  arrange(pct_present)
print(prevalence_df)

# Items with a low prevalence (< 5%) to check
items_low_prevalence <- prevalence_df$item[prevalence_df$pct_present < 5]
cat("Low prevalence items :\n")
print(items_low_prevalence)

###### Polychoric Matrix #####
# NOTE: psych::tetrachoric() can fail on certain data structures (e.g., residual factors, 
# "haven_labelled" columns from SPSS imports, tibbles, etc.) due to its internal matrix 
# conversion via as.matrix(). lavaan::lavCor() utilizes the same pairwise estimation 
# algorithm as efa()/cfa() (which is more robust to sparse data, cf. previous discussion) 
# and is therefore preferred here for consistency and reliability.
library(lavaan)
poly_mat <- lavCor(
  items_recoded,           
  ordered = item_names,
  output = "cor"
)
poly_mat <- as.matrix(poly_mat)


# Search for extreme correlations (>0.90)
extreme_cors <- which(abs(poly_mat) > 0.90 & poly_mat != 1, arr.ind = TRUE)
if (nrow(extreme_cors) > 0) {
  cat("Paires d'items à corrélation tétrachorique extrême :\n")
  print(data.frame(
    item1 = rownames(poly_mat)[extreme_cors[,1]],
    item2 = colnames(poly_mat)[extreme_cors[,2]],
    r = poly_mat[extreme_cors]
  ))
}

# Thresholds of detection of symptom for each item
# Thresholds (expressed as Z-scores) represent the mathematical cutoff 
# on an underlying latent continuous distribution where a response 
# switches from 0 (absence) to 1 (presence). 
# Higher positive values indicate rarer, more severe symptoms.
# Using lavCor(..., output = "th") ensures perfect estimation consistency 
# with the subsequent WLSMV factor analysis.
thresholds <- lavCor(items_recoded, ordered = item_names, output = "th")
print(thresholds)

###### Assumptions checks #####
# KMO test
library(EFAtools)
KMO(items_recoded_num)

# Bartlett's test (H0 : correlation matrix = identity matrix)
cortest.bartlett(poly_mat, n = nrow(items_recoded))

# NOTE: The low KMO value (0.419) is a known mathematical artifact. 
# Inverting a sparse tetrachoric matrix inflates partial correlations, 
# which mechanically deflates the KMO index.
# Factorability is instead validated by the highly significant Bartlett's 
# test (p < 0.001) and subsequent WLSMV fit indices (CFI=0.986).

# Keiser criterion
eigen <- eigen(poly_mat)$values
sum(eigen>1)
# Suggests 7 factors

# Scree plot
df_scree <- data.frame(
  Factor = 1:length(eigen),
  Eigenvalue = eigen
)

# 3. Tracer le Scree Plot
ggplot(df_scree, aes(x = Factor, y = Eigenvalue)) +
  geom_line(color = "purple4", linewidth = 1) +
  geom_point(color = "turquoise", size = 3) +
  # Ligne repère à l'éligibilité classique de Kaiser (Eigenvalue = 1)
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  theme_minimal() +
  labs(
    title = "Scree Plot (Tetrachoric Matrix)",
    x = "Factor number",
    y = "Eigenvalues"
  )
SCREE(x = items_recoded_num, cor_method = "tetra")

# Parallel Analysis
library(EFAtools)
efa_parallel(x = items_recoded_num, N = nrow(items_recoded), eigen_type = "PCA", cor_method = "tetra")

# MAP and HULL test
library(EFAtools)

retention_tests <- EFAtools::efa_retain(
  x = items_recoded_num,                # On injecte directement la matrice
  N = nrow(items_recoded_num),
  cor_method = "tetra", # On spécifie la taille de l'échantillon (574)
  criteria = c("map", "hull")         # Sélectionne le MAP et le HULL
)

# Afficher les résultats
print(retention_tests)


# NOTE: Standard fa.parallel() over-extracted (12 factors) due to matrix 
# non-positive definiteness and smoothing artifacts on sparse binary data.
# Alternative robust methods (MAP and HULL via EFAtools) are preferred here 
# as they prevent over-factoring and confirm a parsimonious 2-factor solution.



# --------------------------------------------------------------
# 7. EFA exploratoire via lavaan avec estimateur WLSMV
# --------------------------------------------------------------

# items_recoded est déjà un data.frame de facteurs ordonnés à 2 niveaux
# (0 = absence, 1 = présence). lavaan traite automatiquement les
# variables "ordered" via des corrélations tétrachoriques sous-jacentes
# lorsque estimator = "WLSMV".
items_ord <- items_recoded

efa_wlsmv <- efa(
  data = items_ord,
  ordered = item_names,
  estimator = "WLSMV",
  nfactors = 1:5,
  rotation = "oblimin" 
)

summary(efa_wlsmv, cutoff=0.4)

# Indices de fit pour choisir le nombre de facteurs
fitMeasures(efa_wlsmv)

# --------------------------------------------------------------
# 8. Extraction de la solution retenue (exemple : 3 facteurs)
# --------------------------------------------------------------

efa_final <- efa(
  data = items_ord,
  ordered = item_names,
  estimator = "WLSMV",
  nfactors = 2,
  rotation = "oblimin"
)

summary(efa_final, nd = 3, cutoff = 0.4, dot.cutoff = 0.2)

# --------------------------------------------------------------
# 9. Vérification des cas de Heywood (signal d'estimation artefactuelle)
# --------------------------------------------------------------
# Sous asymétrie sévère, WLS/WLSMV peut surestimer les saturations
# (Marôco, 2024, Stats). Un signe direct : saturations standardisées
# >= 1, ou variances résiduelles négatives/nulles.

loadings_final <- lavInspect(efa_final, "std")$lambda
cat("Saturations standardisées (recherche de valeurs >= 1) :\n")
print(round(loadings_final, 3))
if (any(abs(loadings_final) >= 0.98)) {
  cat("ATTENTION : au moins une saturation proche ou supérieure à 1 -",
      "signe possible de cas de Heywood / surestimation WLSMV.\n")
}

resid_var <- lavInspect(efa_final, "theta")
resid_diag <- diag(resid_var)
cat("Variances résiduelles (recherche de valeurs <= 0) :\n")
print(round(resid_diag, 3))
if (any(resid_diag <= 0.01)) {
  cat("ATTENTION : au moins une variance résiduelle proche ou inférieure à 0 -",
      "cas de Heywood probable pour cet item.\n")
}

### ==========================================================
### Fin du script
### ==========================================================


########## Reliability analysis ########
rel <- reliability(items = items, nfactors=1)
rel
plot(rel)

### Split-half reliability CDS complete #####
splitHalf(items)

### Alpha reliability CDS complete #####
alpha(items)


###################################################################################################################
# Linear Regression Models #
###################################################################################################################

# We transform the CDS total score using a squared root transformation as is is more adapted to a variable with a strong positive asymetry
# with many zero values.
summary(df$CDS_total_sum)
describe(df$CDS_total_sum)

df$CDS_tot_sqrt <- sqrt(df$CDS_total_sum)
summary(df$CDS_tot_sqrt)

describe(df$CDS_tot_sqrt)[, c("skew", "kurtosis")]
shapiro.test(df$CDS_tot_sqrt)

distrib_CDS_sqrt <- ggplot(df, aes(x= CDS_tot_sqrt)) +
    geom_histogram(aes(y=after_stat(density)), bins=30, fill="purple4", color="white") +
    geom_density(fill="grey", alpha = 0.5) +
    geom_vline(xintercept = mean(df[[col]], na.rm=TRUE), linetype="dashed", color="turquoise", linewidth=1) +
    geom_vline(xintercept = median(df[[col]], na.rm=TRUE), linetype="dotdash", color="red", linewidth=1) +
    labs(title="Distribution of transformed CDS (squared root)")
    theme_minimal()

distrib_CDS_sqrt

ggsave(here("figures", "distrib_CDS_transformed.png"), width=8, height=6, dpi=300)

df$CDS_tot_log <- log1p(df$CDS_total_sum)

########## CDS TOTAL #########

##########

### Modèle 0 ###
Model_0 <- lm(CDS_tot_sqrt~1, na.action = na.exclude, data=df)
summary(Model_0)



### Modèle 1 ###
Model_1 <- lm(CDS_tot_log~ AGE, na.action = na.exclude, data=df)
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
# distribution with many 0 values, a negative binomial regression would be more adapted.

############# Negative Binomial Regression ############
###### Evaluation of model validity ######
Model_A <- glm.nb(CDS_total_sum ~ AGE, data=df, na.action = na.exclude)
summary(Model_A)
exp(coef(Model_A))

# For negative binomial regression, the coefficients are interpreted after being exponentiated
# The exp(coef) is the IRR (Incidence Ratios Rate), intepreted as a % of increase 
# (e.g. : IRR = 1.15 : 15% of increase; IRR = 0.85 : 15% of decrease)
# The theta parameter corresponds to distribution dispersion. For a Poisson distribution, it
# tends to + infinite. The closer it is to 0, the more it corresponds to a negative 
# binomial distribution. The standard error gives the precision of estimation of theta.
# The smaller it is, the best is the theta estimate.
# To confirm the advantage of negative binomial regression over a Poisson regression
# (used for count data but normally dispersed), we have to compare the two models. 
# We use AIC and BIC indicators, as well as a chi-square to compare the fitting of both
# models on the data.

poisson_A <- glm(CDS_total_sum ~ AGE, data=df, na.action = na.exclude)
summary(poisson_A)
AIC(poisson_A, Model_A)
pchisq(2 * (logLik(Model_A) - logLik(poisson_A)), df = 1, lower.tail = FALSE) / 2
# The results confirm that the negative binomial distribution is better adapted to our data
# than a poisson distribution (AIC negative binomial < AIC Poisson, significant X² test).

###### Hierarchical Negative Binomial Regression Models ######
#### Null model ###
Model_nul <- glm.nb(CDS_total_sum ~ 1, data=df, na.action = na.exclude)
summary(Model_nul)
exp(coef(Model_nul))

### Model A : CDS <- Age ###
Model_A <- glm.nb(CDS_total_sum ~ AGE, data=df, na.action = na.exclude)
summary(Model_A)
exp(coef(Model_A))

### Model B : CDS <- Age + Sex ###
Model_B <- glm.nb(CDS_total_sum ~ AGE + SEXE, data=df, na.action = na.exclude)
summary(Model_B)
exp(coef(Model_B))

### Model C : CDS <- Age + Sex + Anxiety ###
Model_C <- glm.nb(CDS_total_sum ~ AGE + SEXE + ANXIETE_recoded, data=df, na.action = na.exclude)
summary(Model_C)
exp(coef(Model_C))


### Model C' : CDS <- Age + Sex + Depression ###
Model_C2 <- glm.nb(CDS_total_sum ~ AGE + SEXE + DEPRESSION_recoded, data=df, na.action = na.exclude)
summary(Model_C2)
exp(coef(Model_C2))


### Model D : CDS <- Age + Sex + Anxiety + Depression ###
Model_D <- glm.nb(CDS_total_sum ~ AGE + SEXE + ANXIETE_recoded + DEPRESSION_recoded, data=df, na.action = na.exclude)
summary(Model_D)
exp(coef(Model_D))

### Model E : CDS <- Age + Sex + Anxiety + Depression + MIGRAINE ###
Model_E <- glm.nb(CDS_total_sum ~ AGE + SEXE + ANXIETE_recoded + DEPRESSION_recoded + MIGRAINE, data=df, na.action = na.exclude)
summary(Model_E)
exp(coef(Model_E))

# Assumptions verification
# Calculer les résidus simulés adaptés au Hurdle
residus_Model_E <- Model_E$residuals

# Tracer les graphiques de diagnostic (recherche de patterns aberrants)
plot(residus_Model_E)


# Statistical comparison of hierarchical models, (chi squared)
aov_mod <- anova(Model_nul, Model_A, Model_B, Model_C, Model_C2, Model_D, Model_E, test="ChiSq")
print(aov_mod)

# Summary table of models coefficients 
stargazer(Model_A, Model_B, Model_C, Model_C2, Model_D, Model_E,
          type = "text",
          title = "Hierarchical Regression Models (Negative Binomial Regression",
          dep.var.labels = "Score CDS total",
          column.labels = c("A", "B", "C", "C2", "D", "E"),
          apply.coef = exp,  # IRR (exp(coef))
          p.auto = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digits = 3,
          out = here("figures", "coef_hierarchical_models.txt"))

# Summary table of models metrics
models_list <- list(Model_nul = Model_nul, Model_A = Model_A, Model_B = Model_B, Model_C = Model_C, 
                    Model_C2 = Model_C2, Model_D = Model_D, Model_E = Model_E)

comparison_df <- data.frame(
  Modele = names(models_list),
  AIC = sapply(models_list, AIC),
  LogLik = sapply(models_list, function(m) as.numeric(logLik(m))),
  Theta = sapply(models_list, function(m) m$theta),
  PseudoR2_McFadden = sapply(models_list, function(m) pR2(m)["McFadden"])
)

print(comparison_df)


# Global summary table of all models (exportation, article ready)
tab_model(Model_A, Model_B, Model_C, Model_C2, Model_D, Model_E,
          show.aic = TRUE,
          transform = "exp",
          dv.labels = c("Model A", "Model B", "Model C", "Model C2", "Model D", "Model E"),
          pred.labels = c("(Intercept)" = "Intercept",
                          "AGE" = "Age",
                          "SEXEh" = "Sex (Man)",
                          "ANXIETE_recoded" = "Anxiety (HADS-A)",
                          "DEPRESSION_recoded" = "Depression (HADS-D)",
                          "MIGRAINEoui" = "Migraine (Yes)"),
          file = here("Figures", "tableau_modeles.doc"))




############ Hurdle Regression ############
anxdep <- c("ANXIETE_recoded", "DEPRESSION_recoded", "MIGRAINE")
df2 <- df[complete.cases(df[, anxdep]), ]
nrow(df2)

scores_pos <- df2$CDS_total_sum[df2$CDS_total_sum > 0]
df_pos <- data.frame(scores_pos=scores_pos)


CDS_tot_distrib <- ggplot(data=df2, aes(x=CDS_total_sum)) +
  geom_histogram(aes(y=after_stat(density)), bins=30, fill="purple4", color="white") +
  geom_density(fill="grey", alpha = 0.5) +
  geom_vline(xintercept = mean(df$CDS_total_sum, na.rm=TRUE), linetype="dashed", color="turquoise", linewidth=1) +
  geom_vline(xintercept = median(df$CDS_total_sum, na.rm=TRUE), linetype="dotdash", color="red", linewidth=1) +
  theme_minimal()

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

# The AIC and BIC differences between log normal and gamma are > 10, meaning that Log Normal has the best fit for our CDS positive value.
# This is visually confirmed by the QQ plot



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