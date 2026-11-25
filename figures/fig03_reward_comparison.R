#Fig03 remake for the RESS revision (Reviewer 3, comment 6: readability).
#new layout: 95% CI of the mean in panel A, box + jitter instead of density in panel B,
#one shared legend, sized for full text width in the two-column format.
#Uses the trajectories of analyses/table5_final.R (2,000 per policy, common seeds), so Fig. 3 and
#the policy summary table agree. Requires output/table5_final/table5_runs.rds.
#Usage (repository root): Rscript figures/fig03_reward_comparison.R

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
fig_out <- out_path("figures", "fig03_reward_comparison.pdf")

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})

policies <- c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention")
palette_policy <- c("Optimal MDP" = "#213C51", "POMDP (PBVI)" = "#574964",
                    "MCDA (TOPSIS)" = "#5A7863", "No Intervention" = "darkred")

runs <- readRDS(out_path("table5_final", "table5_runs.rds"))
pick <- function(lbl) runs$res[[which(sapply(runs$res, `[[`, "label") == lbl)]]
sim <- list(t_pts = runs$t_pts,
            cr   = setNames(lapply(policies, function(p) pick(p)$curves), policies),
            frac = setNames(lapply(policies, function(p) pick(p)$fracs), policies))

#Summary numbers (for checking against Table 5 and the text)
cat("\nFinal discounted return (T = 20 s) and fraction in Nominal:\n")
for (p in policies) {
  v <- sim$cr[[p]][, length(sim$t_pts)]
  cat(sprintf("  %-16s return %.1f +/- %.1f (SD)   frac.nom %.3f\n",
              p, mean(v), sd(v), mean(sim$frac[[p]], na.rm = TRUE)))
}

band_df <- do.call(rbind, lapply(policies, function(p) {
  m <- sim$cr[[p]]
  data.frame(time = sim$t_pts, policy = p, mean = colMeans(m),
             lo = colMeans(m) - 1.96 * apply(m, 2, sd) / sqrt(nrow(m)),
             hi = colMeans(m) + 1.96 * apply(m, 2, sd) / sqrt(nrow(m)))
}))
frac_df <- do.call(rbind, lapply(policies, function(p)
  data.frame(policy = p, frac = sim$frac[[p]])))
band_df$policy <- factor(band_df$policy, levels = policies)
frac_df$policy <- factor(frac_df$policy, levels = policies)

theme_paper <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  plot.title = element_text(hjust = 0.5, size = 15),
  plot.tag = element_text(size = 16, face = "bold"),
  axis.text = element_text(size = 12, colour = "black"),
  axis.title = element_text(size = 13),
  text = element_text(size = 13),
  legend.position = "bottom",
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 12),
  legend.key.width = unit(1.2, "cm")
)

pA <- ggplot(band_df, aes(time, mean, colour = policy, fill = policy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.30, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = palette_policy) +
  scale_fill_manual(values = palette_policy) +
  labs(tag = "A", title = "Cumulative discounted reward",
       x = "Time (s)", y = "Cumulative discounted reward",
       colour = "Policy", fill = "Policy") +
  theme_paper

pB <- ggplot(frac_df, aes(policy, frac, colour = policy, fill = policy)) +
  geom_jitter(width = 0.18, height = 0, size = 0.5, alpha = 0.18, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha = 0.35, outlier.shape = NA, linewidth = 0.6,
               colour = "black", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.8,
               fill = "white", colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = palette_policy) +
  scale_fill_manual(values = palette_policy) +
  scale_x_discrete(labels = c("Optimal\nMDP", "POMDP\n(PBVI)", "MCDA\n(TOPSIS)", "No\nIntervention")) +
  scale_y_continuous(labels = percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(tag = "B", title = "Time spent in Nominal",
       x = NULL, y = "Fraction of time in Nominal") +
  theme_paper

p <- (pA | pB) + plot_layout(widths = c(1.15, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(fig_out, p, width = 10, height = 4.6, device = if (capabilities("cairo")) cairo_pdf else "pdf")
cat("saved", fig_out, "\n")
