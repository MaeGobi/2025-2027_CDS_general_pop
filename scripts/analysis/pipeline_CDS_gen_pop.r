# Pipeline Analysis General Population

library(here)
library(dplyr)
library(rstatix)
library(ggplot2)
library(patchwork)
library(factoextra)
library(psych)

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

# Visualize (boxplots)
box_plot_list <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(y= .data[[col]])) +
    geom_boxplot(outliers_color="black", outliers_size=0.5, fill="skyblue") +
    labs(y=col) +
    theme_light()
  return(p)
})
names(box_plot_list) <- num_cols

# Combine all plots
wrap_plots(box_plot_list, ncol=6)

# Save figurehttp://127.0.0.1:31017/graphics/plot_zoom_png?width=1266&height=636
ggsave(here("figures", "numeric_box_plots.png"),
       wrap_plots(box_plot_list, ncol=6),
       width = 30, height = 20, dpi = 300)

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

########## EFA ###########
items <- df %>% 
  select(all_of(CDS_items)) %>% 
  filter(complete.cases(.))
summary(items)
sum(is.na(items))


describe(items)[, c("skew", "kurtosis")]

##### Assumptions ######

# Multicolinearity

cor_matrix <- cor(items, use = "pairwise.complete.obs")

corr_values <- cor_matrix[lower.tri(cor_matrix)]

range(corr_values)

mean(cor_matrix[lower.tri(cor_matrix)], na.rm = TRUE)


# Bartlett sphericity test
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
# I would say 3-4 factors from the elbow on the scree plot

# Parallel analysis
fa.parallel(items, fm = "ml", fa = "fa", n.iter = 100, main = "Parallel Analysis Scree") 
# Indicates 7 factors (but sensitive to sample size and number of items, might overestimate number of factors)


#### EFA with 3 factors and promax rotation (correlated factors)
efa1_result <- fa(items, nfactors = 3, rotate = "promax", fm = "ml")
print(efa1_result$loadings, cutoff = 0.4, digits = 3)

loadings_matrix1 <- unclass(efa1_result$loadings)
print(loadings_matrix1, digits=3)


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



