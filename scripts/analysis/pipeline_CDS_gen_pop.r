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


describe(items)[, c("skew", "kurtosis")]


########## Reliability analysis ########
rel <- reliability(items = items, nfactors=1)
rel
plot(rel)

### Split-half reliability CDS complete #####
splitHalf(items)

### Alpha reliability CDS complete #####
alpha(items)


########## EFA ###########
##### Assumptions ######
# Polychoric Correlation Matrix
poly_cor <- polychoric(items, na.rm=TRUE, max.cat = 12)
poly_matrix <- poly_cor$rho



# Multicolinearity

cor_matrix <- cor(items, use = "pairwise.complete.obs")

corr_values <- cor_matrix[lower.tri(cor_matrix)]

range(corr_values)

mean(cor_matrix[lower.tri(cor_matrix)], na.rm = TRUE)


# Bartlett sphericity test
# Checks whether correlation matrix is significantly different from identity matrix
cor_matrix <- cor(items)                                   # correlation matrix of items
bartlett_test <- cortest.bartlett(cor_matrix, n = nrow(items))
print(bartlett_test)
# Significant so ok to EFA

# Kaiser-Meye-Oklin measure : sampling adequacy from proportion of variance among items
kmo_result <- KMO(items)
print(kmo_result$MSA) 
#Result > 0.894 = OK (méritoire)

# Kaiser criterion
eigenvalues <- eigen(cor_matrix)$values
sum(eigenvalues > 1)   # number of eigenvalues > 1
# Indicates 7 factors but sensitive tu number of items, overestimation

# Scree plot
pca <- prcomp(items, scale. = TRUE)
fviz_screeplot(pca, addlabels = TRUE, ncp = 10)  # show first 10 components
ggsave(here("figures", "Screeplot.png"), width=8, height=6, dpi=300)
# I would say 3-4 factors from the elbow on the scree plot


# Parallel analysis
fa.parallel(items, fm = "ml", fa = "fa", n.iter = 100, main = "Parallel Analysis Scree") 
# Indicates 7 factors (but sensitive to sample size and number of items, might overestimate number of factors)


#### EFA with 3 factors and promax rotation (correlated factors)
efa1_result <- fa(items, nfactors = 3, rotate = "promax", fm = "minres")
print(efa1_result$loadings, cutoff = 0.4, digits = 3)

loadings_matrix1 <- unclass(efa1_result$loadings)
print(loadings_matrix1, digits=3)

#### EFA with 3 factors and promax rotation (correlated factors) _ Lavaan
library(psych)


library(EFAtools)

# Lancer l'analyse parallèle adaptée aux données ordinales/non-normales
# type = "XG" utilise la méthode de Goerz qui est ultra-résiliente
ap_result <- PARALLEL(items, 
                      N = nrow(items), 
                      data_gender = "none", 
                      engine = "EFAtools",
                      plot = TRUE)

ap_result


library(lavaan)

fit_efa <- efa(data = items,
               nfactors = 3, 
               rotation = "promax",
               estimator = "WLSMV",
               ordered=TRUE)


summary(fit_efa, standardized=TRUE, cutoff=0.4)


fit_multi <- efa(data = items,
                 nfactors = 3:7, 
                 rotation = "promax",
                 estimator = "WLSMV",
                 ordered=TRUE)


summary(fit_multi, standardized=TRUE, cutoff=0.4, fit.measures=TRUE)

# Comparaison directe des modèles via un test de rapport de vraisemblance robuste
lavTestLRT(fit_multi)
# Extraire uniquement les indicateurs clés pour vos modèles ordonnés
fitMeasures(fit_multi, c("chisq.scaled", "df.scaled", "cfi", "tli", "rmsea"))


#### EFA with 4 factors and promax rotation (correlated factors)  
efa_result2 <- fa(items, nfactors = 4, rotate = "promax", fm = "minres")
print(efa_result2$loadings, cutoff = 0.4, digits = 3)

loadings_matrix2 <- unclass(efa_result2$loadings)
print(loadings_matrix2, digits=3)


#### EFA with 5 factors and promax rotation (correlated factors)
efa_result3 <- fa(items, nfactors = 5, rotate = "promax", fm = "minres")
print(efa_result3$loadings, cutoff = 0.4, digits = 3)

loadings_matrix3 <- unclass(efa_result3$loadings)
print(loadings_matrix3, digits=3)




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