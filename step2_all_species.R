# once we have generated all the virtual species, this script calculates for each of them the 3 curves 
# mean iperv and occurrences
# and saves everything into one csv, that will be the at the core of the analysis
library(tidyverse)
library(hypervolume)

# folder that contains all the csvs
setwd("virtualspecies_output")

# same functions as step2.R script
# z transform
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

# bandwidth
compute_global_bw <- function(env_data) {
  estimate_bandwidth(env_data)
}

# hypervolume calculation
hyp_calc <- function(data, bw) {
  tryCatch({
    hv <- hypervolume_gaussian(data, kde.bandwidth = bw * 0.5)
    return(hv@Volume)
  }, error = function(e) {
    return(NA)
  })
}

# accumulation curve
acc_curve <- function(x, bw, step_size = 10) {
  all_idx <- sample(1:nrow(x))
  init_idx <- all_idx[1:5]
  remaining_idx <- all_idx[-(1:5)]
  fx <- x[init_idx, , drop = FALSE]
  ipervolumi <- c(hyp_calc(fx, bw))
  num_occurrences <- c(nrow(fx))
  while (length(remaining_idx) > 0) {
    n_to_add <- min(step_size, length(remaining_idx))
    new_idx <- remaining_idx[1:n_to_add]
    fx <- bind_rows(fx, x[new_idx, , drop = FALSE])
    remaining_idx <- remaining_idx[-(1:n_to_add)]
    ipervolumi <- c(ipervolumi, hyp_calc(fx, bw))
    num_occurrences <- c(num_occurrences, nrow(fx))
  }
  tibble(occ = num_occurrences, iperv = ipervolumi)
}

# more simulations per species and average
simulate_and_average <- function(env_data, num_sim = 10, step_size = 10) {
  bw <- compute_global_bw(env_data)
  all_sims <- lapply(1:num_sim, function(i) acc_curve(env_data, bw, step_size = step_size))
  df_combined <- bind_rows(all_sims, .id = "sim") %>%
    mutate(sim = as.integer(sim))
  df_mean <- df_combined %>%
    group_by(occ) %>%
    summarise(mean_iperv = mean(iperv, na.rm = TRUE), .groups = "drop")
  return(df_mean)
}

# === file list ===
file_list <- list.files("vs_server", pattern = "^species_.*\\.csv$", full.names = TRUE)
file_list <- file_list[!grepl("_aoa", file_list)]

# === process each file ===
results <- map_dfr(seq_along(file_list), function(i) {
  file <- file_list[i]
  
  # extract info from filename
  meta <- str_match(basename(file),
                    "species_(\\d+)_sp_prevalence_([\\d.]+)_sample_prev_[\\d.]+_n_occ_(\\d+).csv"
  )
  if (any(is.na(meta))) return(NULL)
  
  sp_id <- i
  sp_prev <- as.numeric(meta[3])
  n_occ_total <- as.integer(meta[4])
  
  # upload data
  df <- read.csv(file)
  
  # split biased/unbiased
  unbiased_random <- df %>% filter(UNBIASED == TRUE & BIASED == FALSE)
  biased <- df %>% filter(BIASED == TRUE)
  
  # resample unbiased +20%
  n_add <- ceiling(0.2 * nrow(biased))
  random_points <- unbiased_random[sample(1:nrow(unbiased_random), n_add), ]
  unbiased_20 <- rbind(biased, random_points)
  
  # environmental subset
  drops <- c("X", "Y", "distance", "ID", "ID.1", "probability", "UNBIASED", "BIASED", "suitability")
  biased_env <- biased[ , !(names(biased) %in% drops)]
  unbiased_env <- unbiased_20[ , !(names(unbiased_20) %in% drops)]
  total_env <- df[ , !(names(df) %in% drops)]
  
  # Z-transform
  biased_env <- z_transform(biased_env)
  unbiased_env <- z_transform(unbiased_env)
  total_env <- z_transform(total_env)
  
  # curves
  curve_total <- simulate_and_average(total_env)
  curve_biased <- simulate_and_average(biased_env)
  curve_unbiased <- simulate_and_average(unbiased_env)
  
  tibble(
    species_id = sp_id,
    filename = basename(file),
    species_prevalence = sp_prev,
    n_occ = n_occ_total,
    n_occ_biased = nrow(biased),
    n_occ_unbiased = nrow(unbiased_20),
    ipervolume_total = iperv_total,
    ipervolume_biased = iperv_biased,
    ipervolume_unbiased = iperv_unbiased
  )
})

# === final csv ===
write.csv(results, "ipervolume_summary.csv", row.names = FALSE)