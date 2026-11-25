#' model_summary.R
#' ----------------------------------------------------------------------
#' Solves the MDP (value iteration) and the POMDP (PBVI) on the
#' data-driven matrices and writes the solutions:
#'   output/tables/model_summary.txt   A_NoAction, O, MDP policy / V*,
#'                                     continuous-time V, POMDP corner actions
#'   output/tables/mdp_solution.csv    policy, V* and V_ct per regime
#'   output/analyses/model_objects.rds A, O, R, Q_rate, vi, pomdp, belief points
#' ----------------------------------------------------------------------

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
load_core()

corner <- sapply(seq_len(K), function(s) {
  b <- rep(0, K); b[s] <- 1; pomdp_action(b, pomdp$Gamma, pomdp$actions)
})
mdp_tbl <- data.frame(regime = states, mdp_action = unname(vi$policy[states]),
                      V_star = unname(vi$V), V_ct = unname(V_ct),
                      pomdp_action_at_corner = corner,
                      O_diag = unname(diag(O_mat)), A_NoAction_diag = unname(diag(A_NoAction)))
write.csv(mdp_tbl, out_path("tables", "mdp_solution.csv"), row.names = FALSE)

sink(out_path("tables", "model_summary.txt"))
cat("Data-driven baseline transition matrix A_NoAction (dt = 0.02 s)\n")
print(round(A_NoAction, 6))
cat("\nData-driven observation (confusion) matrix O = mcda_confusion\n")
print(round(O_mat, 4))
cat(sprintf("\nNominal accuracy O[Nominal, Nominal] = %.4f\n", O_mat["Nominal", "Nominal"]))
cat(sprintf("\nValue iteration: %d iterations (gamma = %.2f)\n", vi$iterations, gamma_disc))
print(mdp_tbl, row.names = FALSE)
cat(sprintf("\nPBVI: %d belief points, %d iterations, %d alpha-vectors\n",
            length(belief_pts), pomdp$iterations, length(pomdp$Gamma)))
sink()
cat(readLines(out_path("tables", "model_summary.txt")), sep = "\n")

saveRDS(list(A = A, A_NoAction = A_NoAction, O_mat = O_mat, R_mat = R_mat,
             rho = rho, gamma_disc = gamma_disc, dt = dt, Q_rate = Q_rate,
             c_rate = c_rate, r_disc = r_disc, vi = vi, V_ct = V_ct,
             belief_pts = belief_pts, pomdp = pomdp),
        out_path("analyses", "model_objects.rds"))
cat("\n-> output/tables/{model_summary.txt, mdp_solution.csv}, output/analyses/model_objects.rds\n")
