# hypervolume as occurrences increase
# original sampling, unbiased and biased

library(hypervolume)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)

## functions 

# === 1. Z-standard ===
z_transform <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      sd_val <- sd(col, na.rm = TRUE)
      if (sd_val == 0) return(rep(0, length(col)))
      return((col - mean(col, na.rm = TRUE)) / sd_val)
    } else {
      return(col)
    }
  })
  return(df)
}

# === 2. hypervolume calculation ===
compute_global_bw <- function(env_data) {
  estimate_bandwidth(env_data)
}

hyp_calc <- function(data, bw) {
  tryCatch({
    hv <- hypervolume_gaussian(data, kde.bandwidth = bw * 0.5)
    return(hv@Volume)
  }, error = function(e) {
    return(NA)
  })
}

# === 3. accumulation curve ===
acc_curve <- function(env_data, bw, step_size = 10, init = 1) {
  all_idx <- sample(seq_len(nrow(env_data)))
  selected_idx <- all_idx[1:init]
  remaining_idx <- all_idx[-(1:init)]
  
  fx <- env_data[selected_idx, , drop = FALSE]
  ipervolumi <- c(hyp_calc(fx, bw))
  num_occurrences <- c(nrow(fx))
  
  while (length(remaining_idx) > 0) {
    to_add <- min(step_size, length(remaining_idx))
    new_idx <- remaining_idx[1:to_add]
    fx <- rbind(fx, env_data[new_idx, , drop = FALSE])
    remaining_idx <- remaining_idx[-(1:to_add)]
    
    ipervolumi <- c(ipervolumi, hyp_calc(fx, bw))
    num_occurrences <- c(num_occurrences, nrow(fx))
  }
  
  tibble(occ = num_occurrences, iperv = ipervolumi)
}

# === 4. multiple simulations for the same species & average ===
simulate_and_average <- function(env_data, num_sim = 10, step_size = 10, type = "Total") {
  bw <- compute_global_bw(env_data)
  
  all_sims <- replicate(num_sim, acc_curve(env_data, bw, step_size = step_size), simplify = FALSE)
  df_all <- bind_rows(all_sims, .id = "sim")
  
  df_mean <- df_all %>%
    group_by(occ) %>%
    summarise(mean_iperv = mean(iperv, na.rm = TRUE), .groups = "drop") %>%
    mutate(type = type)
  
  return(df_mean)
}

# === 5. example with one species ===

# load dataset (previously generated)
df <- read.csv("vs_server/species_3_sp_prevalence_0.35_sample_prev_0.9_n_occ_800.csv")

# split in unbiased and biased
unbiased_random <- df %>% filter(UNBIASED == TRUE & BIASED == FALSE)
biased <- df %>% filter(BIASED == TRUE)

# resample unbiased
set.seed(13)
n_to_add <- ceiling(0.2 * nrow(biased))
random_sample <- unbiased_random[sample(nrow(unbiased_random), n_to_add), ]
unbiased_20 <- rbind(biased, random_sample)

# select environmental variables
drops <- c("X", "Y", "distance", "ID", "ID.1", "probability", "UNBIASED", "BIASED", "suitability")
biased_env <- biased[ , !(names(biased) %in% drops)]
unbiased_env <- unbiased_20[ , !(names(unbiased_20) %in% drops)]
total_env <- df[ , !(names(df) %in% drops)]

# z-standard 
biased_env <- z_transform(biased_env)
unbiased_env <- z_transform(unbiased_env)
total_env <- z_transform(total_env)

# curves for each dataset: all the occurrences, non biased resampled and biased
curve_total    <- simulate_and_average(total_env,    num_sim = 10, step_size = 30, type = "Total")
curve_biased   <- simulate_and_average(biased_env,   num_sim = 10, step_size = 30, type = "Biased")
curve_unbiased <- simulate_and_average(unbiased_env, num_sim = 10, step_size = 30, type = "Unbiased")

# combine all
curve_all <- bind_rows(curve_total, curve_biased, curve_unbiased)

# === 6. save data & plot ===
write.csv(curve_all, "hypervolume_accumulation_curves.csv", row.names = FALSE)

ggplot(curve_all, aes(x = occ, y = mean_iperv, color = type)) +
  geom_line() +
  geom_smooth(se = FALSE, method = "loess", span = 0.3) +
  theme_minimal() +
  labs(x = "num occurrences", y = "mean hypervolume", color = "Dataset")


# plot (for the manuscript) 
ggplot(curve_all, aes(x = occ, y = mean_iperv, color = type)) +
  geom_step(linewidth = 1) +
  scale_color_manual(
    values = c("Unbiased" = "#FFA500",  # arancione brillante
               "Biased"   = "#800080",  # viola
               "Total"    = "#00C853"), # verde brillante
    labels = c("Unbiased" = "non-biased",
               "Biased"   = "biased",
               "Total"    = "total occurrences")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank()
  ) +
  labs(
    x = "occurrences",
    y = "hypervolume"
  )