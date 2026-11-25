#' Revision_Replot.R
#' Draws the two Supplementary figures
#'   output/figures/PBVI_Convergence.pdf            (Fig. S: PBVI convergence)
#'   output/figures/RewardSensitivity_Combined.pdf  (Fig. S: reward sensitivity)
#' from the data saved by analyses/Revision_PBVI_Convergence.R and
#' analyses/Revision_RewardSensitivity.R. No simulation.

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(patchwork)
})

palette_policy <- c(
  "Optimal MDP"     = "#213C51",
  "POMDP (PBVI)"    = "#574964",
  "k-Step (k=1)"    = "#8B6F47",
  "MCDA (TOPSIS)"   = "#5A7863",
  "No Intervention" = "darkred"
)

theme_paper <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line        = element_line(colour = "black"),
  plot.title       = element_text(hjust = 0.5, size = 20),
  axis.text        = element_text(size = 15),
  text             = element_text(size = 18),
  legend.position  = "bottom",
  legend.text      = element_text(size = 14),
  legend.title     = element_blank()
)

# =====================================================================
# PBVI convergence
# =====================================================================

df_pbvi <- readRDS(out_path("analyses", "Revision_PBVI_history.rds"))
n_iter  <- nrow(df_pbvi)
last_r  <- tail(df_pbvi$residual, 1)

fig_pbvi <- ggplot(df_pbvi, aes(iter, residual)) +
  geom_line(color = "#213C51", linewidth = 1.0) +
  geom_point(color = "#213C51", size = 1.6) +
  geom_hline(yintercept = 1e-10, linetype = "dashed",
             color = "darkred", linewidth = 0.8) +
  annotate("text", x = n_iter * 0.65, y = 1e-9,
           label = expression(epsilon == 10^{-10}),
           color = "darkred", size = 5) +
  annotate("text", x = n_iter, y = last_r,
           label = sprintf("terminated at k = %d", n_iter),
           hjust = 1.05, vjust = -0.6, size = 4.2) +
  scale_y_log10() +
  labs(x = "PBVI iteration k",
       y = expression("max Bellman residual " ~ delta[k]),
       title = "PBVI convergence (MDP warm start)") +
  theme_minimal(base_size = 14) + theme_paper +
  theme(legend.position = "none")

pdf(out_path("figures", "PBVI_Convergence.pdf"), width = 7, height = 4.5)
print(fig_pbvi); dev.off()
cat("-> output/figures/PBVI_Convergence.pdf\n")

# =====================================================================
# Reward sensitivity (3-panel)
# =====================================================================

dfA <- readRDS(out_path("analyses", "Revision_Sensitivity_dfA.rds"))
dfB <- readRDS(out_path("analyses", "Revision_Sensitivity_dfB.rds"))
dfC <- readRDS(out_path("analyses", "Revision_Sensitivity_dfC.rds"))

# k-Step label follows the k selected by the k-Step sweep (e.g. "k-Step (k=2)")
k_lab <- grep("^k-Step", unique(as.character(dfA$policy)), value = TRUE)[1]
names(palette_policy)[grep("^k-Step", names(palette_policy))] <- k_lab
pol_levels <- names(palette_policy)
dfA$policy <- factor(dfA$policy, levels = pol_levels)
dfB$policy <- factor(dfB$policy, levels = pol_levels)
dfC$policy <- factor(dfC$policy, levels = pol_levels)

# 95% confidence half-widths: A, B over the trajectories of each run (N = 200);
# C over the replicate means (n_rep = 20 replicates of 100 trajectories each).
N_AB <- 200
dfA$ci <- 1.96 * dfA$sd / sqrt(N_AB); dfB$ci <- 1.96 * dfB$sd / sqrt(N_AB)
dfC_agg <- aggregate(mean ~ policy + level, data = dfC, FUN = mean)
dfC_sd  <- aggregate(mean ~ policy + level, data = dfC, FUN = sd)
dfC_n   <- aggregate(mean ~ policy + level, data = dfC, FUN = length)
dfC_agg$ci <- 1.96 * dfC_sd$mean / sqrt(dfC_n$mean)
dfC_agg$policy <- factor(dfC_agg$policy, levels = pol_levels)

plot_panel <- function(df, xlab, tag, title) {
  ggplot(df, aes(level, mean, fill = policy)) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             colour = "black", linewidth = 0.25) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci),
                  position = position_dodge(0.8), width = 0.2,
                  colour = "black") +
    scale_fill_manual(values = palette_policy) +
    labs(x = xlab, y = "Mean discounted return",
         title = title, tag = tag) +
    theme_minimal(base_size = 14) + theme_paper +
    theme(plot.tag = element_text(face = "bold", size = 18),
          plot.tag.position = c(0.02, 0.98))
}

pA <- plot_panel(dfA,     expression(alpha), "A", "Uniform scaling")
pB <- plot_panel(dfB,     expression(beta),  "B", "Penalty-only scaling")
pC <- plot_panel(dfC_agg, expression(sigma), "C", "Gaussian perturbation")

# Shared right-hand legend via patchwork
combined <- (pA | pB | pC) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))

pdf(out_path("figures", "RewardSensitivity_Combined.pdf"),
    width = 17, height = 5.2)
print(combined)
dev.off()
cat("-> output/figures/RewardSensitivity_Combined.pdf\n")

cat("\n-- Revision_Replot complete. --\n")
