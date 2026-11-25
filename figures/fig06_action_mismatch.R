#Fig. fig:action_mismatch from the policy-summary trajectories of analyses/table5_final.R
#(2,000 per policy, common seeds). Usage (repository root): Rscript figures/fig06_action_mismatch.R
suppressPackageStartupMessages({ library(ggplot2); library(scales) })
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
fig_out <- out_path("figures", "fig06_action_mismatch.pdf")
# per-trajectory mismatch of the policy-summary runs (analyses/table5_final.R)
runs <- readRDS(out_path("table5_final", "table5_runs.rds"))$res
r <- Filter(function(x) x$label %in% c("POMDP (PBVI)", "MCDA (TOPSIS)"), runs)
d <- data.frame(policy = sapply(r, `[[`, "label"),
                m = sapply(r, function(x) mean(x$mism)),
                se = sapply(r, function(x) sd(x$mism) / sqrt(length(x$mism))))
d$policy <- factor(d$policy, levels = c("POMDP (PBVI)", "MCDA (TOPSIS)"))
p <- ggplot(d, aes(policy, m, fill = policy)) +
  geom_col(width = 0.5, alpha = 0.9) +
  geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), width = 0.12) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * m)), vjust = -1.1, size = 4.6, fontface = "bold") +
  scale_fill_manual(values = c("POMDP (PBVI)" = "#574964", "MCDA (TOPSIS)" = "#5A7863"), guide = "none") +
  scale_y_continuous(labels = percent, limits = c(0, 0.8), expand = c(0, 0)) +
  labs(x = NULL, y = "Action mismatch rate") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13))
ggsave(fig_out, p, width = 4.5, height = 3.4, device = if (capabilities("cairo")) cairo_pdf else "pdf")
print(d)
