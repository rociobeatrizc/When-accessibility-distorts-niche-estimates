# report of the data in area of applicability and DI
install.packages("car")

library(CAST)
library(caret)
library(sf)
library(devtools)
library(raster)
library(viridis)
library(ggplot2)
library(tidyverse)
library(broom)
library(terra)
library(lme4)
library(scales)
library(car)
library(performance)
library(ggpredict)
library(mgcv)
library(ggeffects)

# in this folder, there are all the outputs from the step3_all_species.R script
setwd("virtualspecies_output")

### files processing 
# the first function takes all the Area of Applicability rasters and calculate summary statistics
analyze_AOA <- function(raster_stack, id) {
  
  # raster AOA
  aoa_null <- raster_stack$AOA_null
  aoa_biased <- raster_stack$AOA_biased
  aoa_all <- raster_stack$AOA_all
  
  # pixels 
  common_pixels <- sum(aoa_null[] == 1 & aoa_biased[] == 1, na.rm = TRUE)
  exclusive_null <- sum(aoa_null[] == 1 & aoa_biased[] == 0, na.rm = TRUE)
  exclusive_biased <- sum(aoa_null[] == 0 & aoa_biased[] == 1, na.rm = TRUE)
  count_all <- sum(aoa_all[] == 1, na.rm = TRUE)
  
  # data frame
  result <- data.frame(
    ID = id,
    Common_Pixels = common_pixels,
    Exclusive_Null = exclusive_null,
    Exclusive_Biased = exclusive_biased,
    Count_All = count_all
  )
  
  return(result)
}

# the second function takes all the Dissimilarity Index rasters and calculate summary statistics
analyze_DI <- function(raster_stack, species_id) {
  
  # DI null 
  # DI biased 
  # DI all
  di_null <- raster_stack$DI_null
  di_biased <- raster_stack$DI_biased
  di_all <- raster_stack$DI_all
  
  # total number of pixels in di_all
  total_pixels_di_all <- sum(!is.na(values(di_all)))
  
  # difference
  difference_DI <- di_null - di_biased
  difference_DI_all_biased <- di_all - di_biased
  difference_DI_all_null <- di_all - di_null
  
  # mean difference
  mean_difference <- mean(values(difference_DI), na.rm = TRUE)
  mean_difference_DI_all_biased <- mean(values(difference_DI_all_biased), na.rm = TRUE)
  mean_difference_DI_all_null <- mean(values(difference_DI_all_null), na.rm = TRUE)
  
  # na.remove
  total_pixels <- sum(!is.na(values(difference_DI)))
  
  # DI null > DI biased
  pixels_null_greater <- sum(values(di_null) > values(di_biased), na.rm = TRUE)
  
  # DI_null < DI_biased
  pixels_biased_greater <- sum(values(di_null) < values(di_biased), na.rm = TRUE)
  
  # DI_null == DI_biased
  unchanged_pixels <- sum(values(di_null) == values(di_biased), na.rm = TRUE)
  
  # mean and sd
  mean_null <- mean(values(di_null), na.rm = TRUE)
  mean_biased <- mean(values(di_biased), na.rm = TRUE)
  mean_all <- mean(values(di_all), na.rm = TRUE)
  sd_null <- sd(values(di_null), na.rm = TRUE)
  sd_biased <- sd(values(di_biased), na.rm = TRUE)
  sd_all <- sd(values(di_all), na.rm = TRUE)
  
  # Calculate ddt and ddb
  ddt <- di_all - di_null
  ddb <- di_all - di_biased
  
  # Count pixels where ddt and ddb are equal to zero
  pixels_ddt_zero <- sum(values(ddt) == 0, na.rm = TRUE)
  pixels_ddb_zero <- sum(values(ddb) == 0, na.rm = TRUE)
  
  # Count pixels where ddt and ddb are positive
  pixels_ddt_positive <- sum(values(ddt) >  0, na.rm = TRUE)
  pixels_ddb_positive <- sum(values(ddb) > 0, na.rm = TRUE)
  
  # ddb negative
  pixels_ddb_negative <- sum(values(ddb) < 0, na.rm = TRUE)
  
  # DI all less than null/biased
  all_great_than_null_2 <- all(values(ddt) <  0, na.rm = TRUE)
  all_great_than_biased_2 <- all(values(ddb) < 0, na.rm = TRUE)
  
  
  # df
  result <- data.frame(
    Species_ID = species_id,
    Total_pixels_DI_all = total_pixels_di_all,
    Mean_DI_null = mean_null,
    SD_DI_null = sd_null,
    Mean_DI_biased = mean_biased,
    SD_DI_biased = sd_biased,
    Mean_DI_all = mean_all,
    SD_DI_all = sd_all,
    Pixels_null_greater = pixels_null_greater,
    Pixels_biased_greater = pixels_biased_greater,
    Unchanged_pixels = unchanged_pixels,
    Mean_difference_DI = mean_difference,
    Mean_difference_DI_all_biased = mean_difference_DI_all_biased,
    Mean_difference_DI_all_null = mean_difference_DI_all_null,
    Pixels_ddt_zero = pixels_ddt_zero,
    Pixels_ddb_zero = pixels_ddb_zero,
    Pixels_ddb_positive = pixels_ddb_positive,
    Pixels_ddb_negative = pixels_ddb_negative,
    Pixels_ddt_positive = pixels_ddt_positive,
    all_great_than_null = all_great_than_null_2,
    all_great_than_biased = all_great_than_biased_2
    
  )
  
  return(result)
}


# process all the rasters in the folder with AOA and DI stats
process_rasters <- function(folder = ".") {
  # files with _aoa
  raster_files <- list.files(folder, pattern = "_aoa\\.tif$", full.names = TRUE)
  
  # df
  final_results <- data.frame()
  
  for (file in raster_files) {
    # ID species
    base_name <- tools::file_path_sans_ext(basename(file))
    id <- sub("_aoa$", "", base_name)
    
    # species prevalence and number of occurrences 
    matches <- regmatches(id, regexec("sp_prevalence_([0-9.]+).*n_occ_([0-9]+)", id))
    if (length(matches[[1]]) > 2) {
      sp_prev <- as.numeric(matches[[1]][2])
      n_occ <- as.integer(matches[[1]][3])
    } else {
      sp_prev <- NA
      n_occ <- NA
    }
    
    # upload stack
    raster_stack <- rast(file)
    
    # AOA
    aoa_results <- analyze_AOA(raster_stack, id)
    aoa_results$Species_Prevalence <- sp_prev
    aoa_results$N_Occ <- n_occ
    
    # DI
    di_results <- analyze_DI(raster_stack, id)
    di_results$Species_Prevalence <- sp_prev
    di_results$N_Occ <- n_occ
    
    # same df
    combined_results <- cbind(aoa_results, di_results[, -1])
    final_results <- rbind(final_results, combined_results)
  }
  
  # save in csv
  write.csv(final_results, file = "aoa_di.csv", row.names = FALSE)
  return(final_results)
}

# use the function
csvs <- process_rasters()
head(csvs)


# upload virtual species and get the number of occurrences (unbiased, biased) from there
# these informations needed for further analysis
merge_with_species <- function() {
  # aoa and di
  aoa_data <- read.csv("aoa_di.csv")
  
  # original csvs with the species
  csv_files <- list.files(pattern = "^species_.*\\.csv$")
  
  # unbiased and biased points 
  calculate_occurrences <- function(csv_file) {
    # ID
    file_id <- tools::file_path_sans_ext(basename(csv_file))
    
    # read file
    species_data <- read.csv(csv_file)
    
    # split between biased and unbiased
    biased <- species_data[species_data$BIASED == TRUE, ]
    unbiased_random <- species_data[species_data$UNBIASED == TRUE & species_data$BIASED == FALSE, ]
    
    # unbiased + 20%
    n_random_to_add <- ceiling(0.2 * nrow(biased))
    set.seed(123)
    random_points <- unbiased_random[sample(1:nrow(unbiased_random), n_random_to_add), ]
    
    # combine
    unbiased_20 <- rbind(biased, random_points)
    
    # return ID
    data.frame(
      ID = file_id,
      Occurrences_Biased = nrow(biased),
      Occurrences_Unbiased_20 = nrow(unbiased_20)
    )
  }
  
  # lapply for each species
  occurrence_data <- do.call(rbind, lapply(csv_files, calculate_occurrences))
  
  # merge occurrence data with aoa-di 
  final_data <- merge(aoa_data, occurrence_data, by.x = "ID", by.y = "ID", all.x = TRUE)
  
  # output
  write.csv(final_data, "output_final_csv_file.csv", row.names = FALSE)
  
  return(final_data)
}

# apply function
csvs_1 <- merge_with_species()

# check 
head(csvs_1)
names(csvs_1)

# useless columns
csvs_1$N_Occ.1 <- NULL
csvs_1$Species_Prevalence.1 <- NULL

## some more... 
# % biased occurrences
csvs_1$Occurrences_Biased_norm <- csvs_1$Occurrences_Biased / csvs_1$N_Occ

# % null occurrences
csvs_1$Occurrences_null_norm <- csvs_1$Occurrences_Unbiased_20 / csvs_1$N_Occ

# ddb zero norm
csvs_1$ddb_zero_norm <- csvs_1$Pixels_ddb_zero / csvs_1$Total_pixels_DI_all

# ddb positive norm
csvs_1$ddb_postive_norm <- csvs_1$Pixels_ddb_positive /csvs_1$Total_pixels_DI_all

# ddb negative norm
csvs_1$ddb_negative_norm <- csvs_1$Pixels_ddb_negative /csvs_1$Total_pixels_DI_all

# check
head(csvs_1)

## DI analysis
quartile_limits <- quantile(csvs_1$Occurrences_Biased_norm, probs = seq(0, 1, 0.25), na.rm = TRUE)

csvs_1 <- csvs_1 %>%
  mutate(
    Biased_Occ_Range = cut(
      Occurrences_Biased_norm,
      breaks = quartile_limits,
      include.lowest = TRUE,
      labels = paste0(
        "Range: ",
        round(head(quartile_limits, -1), 3),
        "-",
        round(tail(quartile_limits, -1), 3)
      )
    )
  )

summary_table <- csvs_1 %>%
  group_by(Species_Prevalence, Biased_Occ_Range) %>%
  summarise(
    Mean_DI_null = mean(Mean_DI_null, na.rm = TRUE),
    Mean_DI_biased = mean(Mean_DI_biased, na.rm = TRUE),
    Mean_DDB = mean(Mean_difference_DI, na.rm = TRUE),
    Pct_DDB_gt0 = mean(ddb_postive_norm, na.rm = TRUE) * 100,
    Pct_DDB_eq0 = mean(ddb_zero_norm, na.rm = TRUE) * 100,
    Pct_DDB_lt0 = mean(ddb_negative_norm, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# preliminary statistics 
head(summary_table)
write.csv(summary_table, "summaryDI.csv", row.names = FALSE)

## GLM
glm_ddb <- glm(ddb_zero_norm ~ Occurrences_Biased_norm * Species_Prevalence, 
               data = csvs_1, 
               family = quasibinomial)
summary(glm_ddb)

# residuals
plot(resid(glm_ddb) ~ predict(glm_ddb))
abline(h = 0, col = "red")

# histogram
hist(resid(glm_ddb), 
     main = "lm", 
     xlab = "residuals", 
     breaks = 30)

# Q-Q Plot
qqnorm(resid(glm_ddb), main = "Q-Q Plot (GLM)")
qqline(resid(glm_ddb), col = "red")


# predictions with ggpredict
# prediction 1 (fixed values of species prevalence)
pred1 <- ggpredict(glm_ddb, terms = c("Occurrences_Biased_norm[all]", "Species_Prevalence[0.15, 0.3, 0.4]"))

# prediction 2 (fixed values of number of biased occurrences)
pred2 <- ggpredict(glm_ddb, terms = c("Species_Prevalence", "Occurrences_Biased_norm"))

# limits (for plot purposes)
ymin <- min(c(pred1$predicted, pred1$conf.low, pred2$predicted, pred2$conf.low), na.rm = TRUE)
ymax <- max(c(pred1$predicted, pred1$conf.high, pred2$predicted, pred2$conf.high), na.rm = TRUE)

# plot 1
plot(pred1, show_ci = TRUE, ci_style = "ribbon",
     show_residuals_line = TRUE,
     colors = c("purple", "green3", "orange"),
     alpha = 0.2,
     data_labels = TRUE) + 
  xlab("Normalized biased occurrences") + 
  ylab("Pixels with ddb = 0 normalized") +
  guides(color = guide_legend(title = "Species Prevalence")) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_blank()) 

# plot 2
plot2 <- plot(pred2, show_ci = TRUE, ci_style = "ribbon",
              show_residuals_line = TRUE,
              colors = c("red", "darkblue", "orange4"),
              alpha = 0.2,
              data_labels = TRUE) + 
  xlab("Species Prevalence") + 
  ylab("Pixels with ddb = 0 normalized") +
  guides(color = guide_legend(title = "Normalized biased occurrences")) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_blank())

# GLM validation
# 1. predictions 
predictions <- predict(glm_ddb, type = "response")


# 2. deviance and Pseudo-R² (McFadden's R²)
model_deviance <- deviance(glm_ddb)
null_deviance <- deviance(glm(ddb_zero_norm ~ 1, data = csvs_1, family = quasibinomial))
pseudo_r_squared <- 1 - model_deviance / null_deviance

# 3. cross-validation (10-fold)
train_control <- trainControl(method = "cv", number = 10)
cv_model <- train(ddb_zero_norm ~ Occurrences_Biased_norm * Species_Prevalence, 
                  data = csvs_1, 
                  method = "glm", 
                  family = quasibinomial, 
                  trControl = train_control)
print(cv_model)

# 4. dispersion
dispersion <- summary(glm_ddb)$dispersion
print(paste("Dispersion parameter:", dispersion))

# residuals
residuals_glm <- residuals(glm_ddb, type = "response")

# RMSE
rmse_residuals <- sqrt(mean(residuals_glm^2))

### Relationship bewteen number of biased occurrences and species prevalence
ggplot(csvs_1, aes(x = as.factor(Species_Prevalence), y = Occurrences_Biased_norm)) +
  geom_boxplot(fill = "white", color = "black", outlier.color = "black", outlier.shape = 16) +
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "red", size = 0.8) +  # Linea rossa sulla mediana
  labs(x = "Species prevalence", 
       y = "Normalized biased occurrences"
       #      title = "Distribution of biased occurrences in relation to species prevalence"
  ) +
  theme_bw(base_size = 10) +  
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# levene test (dispersion)
leveneTest(Occurrences_Biased_norm ~ as.factor(Species_Prevalence), data = csvs_1)

# variance
var_data <- csvs_1 %>%
  group_by(Species_Prevalence) %>%
  summarise(variance = var(Occurrences_Biased_norm))

# linear model on variance
lm_var <- lm(variance ~ Species_Prevalence, data = var_data)
summary(lm_var)

ggpredict(lm_var) %>% 
  plot(., show_ci = TRUE, ci_style = "ribbon",
       show_residuals_line = TRUE,
       colors = c("darkgreen"),
       data_labels = TRUE) +
  labs(title = "variance of biased occurrences as the species prevalence increases") +
  xlab("species prevalence") + 
  ylab("variance of biased occurrences") 

dev.off()