# hypervolume analysis

library(dplyr)
library(tidyr)
library(stringr)
library(ggeffects)

df <- read.csv("summary_hypervolume.csv")
df$index <- (df$ipervolume_unbiased - df$ipervolume_biased) / (df$ipervolume_biased + df$ipervolume_unbiased)
df$biased_norm <- (df$n_occ_biased / df$n_occ)

# plot
ggplot(df, aes(x = biased_norm, y = index, color = species_prevalence)) +
  geom_point(size = 2) +  
  scale_color_viridis_c(limits = c(0.05, 0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Points Biased e Index",
       x = "Points Biased",
       y = "Index",
       color = "Species Prevalence") +
  theme_minimal()


df$index_scaled <- (df$index + 1) / 2
df$index_scaled

# GLM on index scaled
glm_hyp <- glm(index_scaled ~ biased_norm + species_prevalence, 
               data = df, 
               family = quasibinomial)

# summary
summary(glm_hyp)
plot(residuals(glm_hyp))

# predictions
pred1 <- ggpredict(glm_hyp, c("biased_norm[all]", "species_prevalence[0.1, 0.3, 0.5]"))
pred1$predicted <- (pred1$predicted * 2) - 1
pred1$conf.high <- (pred1$conf.high * 2) - 1
pred1$conf.low <- (pred1$conf.low * 2) - 1

pred2 <- ggpredict(glm_hyp, terms = c("species_prevalence", "biased_norm[0.24, 0.35, 0.46]"))
pred2$predicted <- (pred2$predicted * 2) - 1
pred2$conf.high <- (pred2$conf.high * 2) - 1
pred2$conf.low <- (pred2$conf.low * 2) - 1

# plot 1
pred1 %>% 
  plot(., show_ci = TRUE, ci_style = "ribbon",
       show_residuals_line = TRUE,
       colors = c("purple", "green3", "orange"),
       data_labels = TRUE) +
  xlab("Normalized biased occurrences") +
  ylab("Hypervolume iIndex") +
  ylim(-0.3,0.3) +
  guides(color = guide_legend(title = "Species Prevalence")) +
  xlim( min(pred1$x), max(pred1$x)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  ggtitle(NULL)  

# plot 2
pred2 %>% 
  plot(., show_ci = TRUE, ci_style = "ribbon",
       show_residuals_line = TRUE,
       colors = c("red", "darkblue", "orange4"),
       data_labels = TRUE) +
  # labs(title = "ipervolume index as species prevalence increase") +
  xlab("Species Prevalence") +
  ylab("Hypervolume Index") +
  ylim(-0.3,0.3) +
  xlim( min(pred2$x), max(pred2$x)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  guides(color = guide_legend(title = "Normalized biased occurrences")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  ggtitle(NULL) 