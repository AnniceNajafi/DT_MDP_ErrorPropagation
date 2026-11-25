#' Revision_PBVI_Convergence.R
#' ──────────────────────────────────────────────────────────────────────
#' Reviewer Comment 3:
#'   "Convergence is defined when the maximum change in value across all
#'    belief points falls below a prescribed tolerance - but the tolerance
#'    value is not stated for PBVI. Report the actual tolerance and number
#'    of iterations to convergence."
#'
#' This script re-runs PBVI with:
#'   - tolerance  epsilon = 1e-10    (matching value iteration)
#'   - max_horizon = 2000 backups
#'   - MDP warm start via V_init = vi$V
#'
#' and records the max Bellman residual
#'   delta_k = max_b max( |V_k(b) - V_{k-1}(b)| )
#' at every iteration k. Output is a log-scale convergence curve.
#'
#' Outputs (the figure is drawn by figures/Revision_Replot.R):
#'   - output/analyses/Revision_PBVI_history.rds
#'   - output/tables/PBVI_Convergence.csv
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
load_core()   # pipeline -> MDP -> POMDP (defines belief_pts, A, O_mat, ...)

set.seed(42)

# =====================================================================
# 1. PBVI SOLVER WITH CONVERGENCE HISTORY
# =====================================================================
# Identical to pbvi_solve() in POMDPFormulation.R but records max Bellman
# residual at every iteration.

pbvi_solve_history <- function(A_list, O_mat, R_mat, gamma, belief_pts,
                               max_horizon = 500, tol = 1e-10,
                               V_init = NULL) {
  K <- ncol(O_mat)
  act <- rownames(R_mat)
  M <- length(act)
  n_pts <- length(belief_pts)
  B <- do.call(rbind, belief_pts)

  if (!is.null(V_init)) {
    Gamma_mat <- matrix(V_init, nrow = 1, ncol = K)
    Gamma_labels <- act[1]
  } else {
    Gamma_mat <- matrix(0, nrow = 1, ncol = K)
    Gamma_labels <- act[1]
  }

  history <- numeric(max_horizon)

  for (iter in seq_len(max_horizon)) {
    backed <- vector("list", M)
    for (ai in seq_len(M)) {
      a <- act[ai]
      backed[[ai]] <- vector("list", K)
      for (z in seq_len(K)) {
        weighted <- t(Gamma_mat) * O_mat[, z]
        backed[[ai]][[z]] <- gamma * (A_list[[a]] %*% weighted)
      }
    }

    Gamma_new  <- matrix(NA, n_pts, K)
    Labels_new <- character(n_pts)

    for (bi in seq_len(n_pts)) {
      b <- B[bi, ]
      best_val <- -Inf; best_alpha <- NULL; best_action <- NULL
      for (ai in seq_len(M)) {
        alpha_a <- R_mat[ai, ]
        for (z in seq_len(K)) {
          scores <- as.numeric(crossprod(backed[[ai]][[z]], b))
          best_gi <- which.max(scores)
          alpha_a <- alpha_a + backed[[ai]][[z]][, best_gi]
        }
        val <- sum(alpha_a * b)
        if (val > best_val) {
          best_val <- val; best_alpha <- alpha_a; best_action <- act[ai]
        }
      }
      Gamma_new[bi, ] <- best_alpha
      Labels_new[bi]  <- best_action
    }

    V_new <- rowSums(Gamma_new * B)
    V_old <- apply(B %*% t(Gamma_mat), 1, max)
    max_diff <- max(abs(V_new - V_old))
    history[iter] <- max_diff

    keep <- !duplicated(round(Gamma_new, 10))
    Gamma_mat    <- Gamma_new[keep, , drop = FALSE]
    Gamma_labels <- Labels_new[keep]

    if (max_diff < tol) {
      history <- history[seq_len(iter)]
      break
    }
  }

  list(history = history, iterations = length(history),
       final_residual = tail(history, 1),
       Gamma = lapply(seq_len(nrow(Gamma_mat)), function(i) Gamma_mat[i, ]),
       actions = Gamma_labels)
}

# =====================================================================
# 2. RUN PBVI WITH HISTORY
# =====================================================================

cat("\n── Running PBVI with convergence history (budget 500 backups) ──\n")
pomdp_hist <- pbvi_solve_history(A, O_mat, R_mat, gamma_disc, belief_pts,
                                 max_horizon = 500, tol = 1e-10,
                                 V_init = vi$V)

iters_to_converge <- pomdp_hist$iterations
final_resid       <- pomdp_hist$final_residual
cat(sprintf("   iterations: %d   final residual: %.3e\n",
            iters_to_converge, final_resid))

# =====================================================================
# 3. SAVE (figure: figures/Revision_Replot.R)
# =====================================================================

df_pbvi <- data.frame(iter = seq_along(pomdp_hist$history),
                      residual = pmax(pomdp_hist$history, 1e-14))
saveRDS(df_pbvi, out_path("analyses", "Revision_PBVI_history.rds"))
write.csv(df_pbvi, out_path("tables", "PBVI_Convergence.csv"), row.names = FALSE)
cat("  -> Saved output/analyses/Revision_PBVI_history.rds, output/tables/PBVI_Convergence.csv\n")

cat("\n══ Revision_PBVI_Convergence complete. ══\n")
