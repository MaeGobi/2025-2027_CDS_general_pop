# Pipeline Analysis General Population

library(here)
library(dplyr)
library(rstatix)
library(ggplot2)
library(patchwork)

######### Charge and visualize dataset ##########

# Repository definition for R project. All files called from this repo.
here()

df <- read.csv(file = here("data", "2026-ALL350PATIENTS350CONTROLS_PK_15-06-26.csv"), header = TRUE, sep = ";", dec = ",")
head(df)
colnames(df)
nrow(df)

# Transform df to exclude rows from VESTICOR study
df <- df[1:665, ]

nrow(df)


########## Preprocessing ##########
# NaN inspection

# Define interest columns
cols_int <- c("SEXE", "LATERALITE", "AGE", "PROFESSION", "FAMILLE", "FUMEUR", "MIGRAINE", "ALCOOL", "aucuntroublevisuel",
              "myopie", "astigmatisme", "maladierétine", "glaucome", "lunettes", "lentilles", "ETUDES", "CDS1", "CDS2",
              "CDS3", "CDS4", "CDS5", "CDS6", "CDS7", "CDS8", "CDS9", "CDS10", "CDS11", "CDS12", "CDS13", "CDS14", 
              "CDS15","CDS16", "CDS17", "CDS18", "CDS19", "CDS20", "CDS21", "CDS22", "CDS23", "CDS24", "CDS25", "CDS26",
              "CDS27","CDS28", "CDS29", "OBE", "ANXIETE", "DEPRESSION", "CDStotal", "FREQUENCYall", "DURATIONall")

# Count number of missing/empty values for each interest column
colSums(is.na(df[cols_int]) | df[cols_int]=="")

# Display row number of missing/empty values for each variable
missing <- sapply(df[cols_int], function(x) which(is.na(x) | x == "") + 1)
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
av <- sum(df$CDStotal > 290, na.rm = TRUE)
cat ("CDStotal : Aberrant values (>290):", av, "\n")


# For HADS A and D
anxdep <- c("ANXIETE", "DEPRESSION")
for (col in anxdep) {
  av <- sum(df[[col]]  > 21, na.rm = TRUE)
  cat (col, " : Aberrant values (>21):", av, "\n")
}




######### Outliers identification #########
num_cols <- c("AGE", "CDS1", "CDS2","CDS3", "CDS4", "CDS5", "CDS6", "CDS7", "CDS8", "CDS9", "CDS10", "CDS11", "CDS12",
              "CDS13", "CDS14", "CDS15","CDS16", "CDS17", "CDS18", "CDS19", "CDS20", "CDS21", "CDS22", "CDS23", 
              "CDS24", "CDS25", "CDS26", "CDS27","CDS28", "CDS29", "ANXIETE", "DEPRESSION", "CDStotal", 
              "FREQUENCYall", "DURATIONall")

df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

summary(df[num_cols])

####### Plot numerical variables distribution #######
plot_list <- lapply(num_cols, function(col) {
  p <- ggplot(df, aes(x= .data[[col]])) +
  geom_histogram(aes(y=after_stat(density)), bins=30, fill="purple4", color="white") +
  geom_density(fill="grey", alpha = 0.5) +
  geom_vline(xintercept = mean(df[[col]], na.rm=TRUE), linetype="dashed", color="turquoise", linewidth=1) +
  geom_vline(xintercept = median(df[[col]], na.rm=TRUE), linetype="dotdash", color="red", linewidth=1) +
  theme_minimal()
})
names(plot_list) <- num_cols

# Combine all plots
wrap_plots(plot_list, ncol=6)

# Save figure
ggsave(here("figures", "numeric_distribution_plots.png"),
       wrap_plots(plot_list, ncol=6),
       width = 30, height = 20, dpi = 300)

# Visualize (boxplots), count and identify values of outliers
par(mfrow = c(3,3))

for (col in num_cols) {
  out_values <- boxplot.stats(df[[col]])$out
  
  cat("\n=== Colonne:", col, "===\n")
  cat("Nombre d'outliers:", length(out_values), "\n")
  if (length(out_values) > 0) {
    cat("Valeurs:", paste(out_values, collapse = ", "), "\n")
  }
}


