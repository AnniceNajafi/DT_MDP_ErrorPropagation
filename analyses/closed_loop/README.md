# Closed-loop validation on the six-module twin (Reviewer 2, comments 6 and 8)

The decision model is loaded unchanged from `R/POMDPFormulation.R` (via `load_core()`), and the physics/ARX helpers of `R/SimulationPipeline.R` are copied verbatim into `cl_sim.R`. All results are written to `output/closed_loop/`.

## Files
| file | role |
|---|---|
| `cl_setup.R` | Loads the paper's decision model unchanged (A^(a), O = in-sample Viterbi confusion, R, `vi`, PBVI with 500 iterations, seed 2024, giving 72 alpha vectors). Fits **one** causal decoder HMM (4-state Gaussian, 10 residual features) on pooled training seeds 11-20 (canonical schedule) and checks it on held-out seeds 111-120. Writes `output/closed_loop/setup.rds`. |
| `cl_sim.R` | Sample-by-sample (dt = 0.02 s) six-module simulator: interventions act on the twin; online forward-filter decoder; controllers. |
| `verify.R` | Checks that the stepwise twin with no actions reproduces `simulate_with_faults_pipe()` (max abs difference about 5e-14 on all six residuals, seeds 11 and 111), and prints calibration diagnostics. |
| `cl_repair.R` | Repair-time experiment (comment 6), giving `repair_runs*.csv` and `repair_summary*.csv`. |
| `cl_closed_loop.R` | Closed-loop controller comparison (comment 8), giving `closed_loop_{runs,summary,paired}*.csv` and `closed_loop_results*.rds`. |
| `run_closed_loop.R` | Reruns everything (called by `run_all.R`, step `closed_loop`). Total wall time is about 25 min with 8-10 parallel workers; `--quick` uses 20 seeds per experiment. |

## Twin, faults and decoder
- **Physical asset and ARX twin:** identical to `SimulationPipeline.R` (same RNG order). The ARX surrogates are fitted on the first 10 s of each run's asset data and then free-run.
- **Faults:** the injected faults are identical to the pipeline.
  - SensorNoisy: white noise with sd 0.06 on twin sensor A and sd 0.05 on twin sensor B.
  - DynamicsOff: the reported plant-surrogate position is scaled by 0.8.
  - Drift: the sensor-A offset is 0.0007·(k − k_on + 725). This equals the pipeline's 0.0007·(k−1) when k_on = 726, so the fault signature at onset matches the decoder's training data (a step of 0.5075 plus a ramp of 0.035 /s).
- **Decoder:** depmixS4 default starting values were infeasible on the pooled data, so EM is started from k-means. There are 6 restarts and the best log-likelihood is kept; its training Viterbi accuracy is 0.869.
  - Standardised features are clipped to ±6 SD. Without the clip, EM underflows on the 20-40 SD spikes in e1_diff and e4 at the DynamicsOff edges.
  - The online forward filter matches `depmixS4::forwardbackward` to 7e-16.
  - Held-out causal confusion (diagonal): Nominal 0.803, SensorNoisy 0.672, DynamicsOff 1.000, Drift 0.996.
  - **Transition smoothing:** EM on the canonical schedule learns exact zeros: Drift is absorbing and nothing returns to Nominal from DynamicsOff or Drift. A filter built on those zeros can never report a repair. The default therefore applies the paper's own A_NoAction Laplace smoothing (eps = 1e-3) to the decoder transition matrix. `CL_DEC_EPS=0` runs the raw decoder instead.

## Interventions (operations on the twin, active from the sample after they are chosen)
| action | implementation |
|---|---|
| DwSensorA | Twin fusion weight on sensor A changes from 0.6 to 0.1. |
| DwSensorB | Twin fusion weight on sensor B changes from 0.4 to 0.1 (alpha = 0.9). |
| IncFilter | An extra first-order low-pass filter with tau_f = 0.10 s is placed on both twin sensor channels. |
| ReidentPlant | Collects W_rp = 2 s (100 samples) of asset data (u_act, x_pos) while the action is held, then refits ARX(2,2) by least squares. The refit is **installed only if** its AR part is stable and |DC gain / (1/k) − 1| ≤ 5 %; otherwise a new window starts. Switching away from the action discards the window. |
| BiasCorrect | Collects W_bc = 1 s (50 samples) of the sensor-A residual e2a while the action is held. It fits an offset plus slope by OLS and adds the extrapolated line to the sensor-A correction. Switching away discards the window. |

## Repair definitions (ground-truth regime)
- **DynamicsOff** ends when an accepted refit is installed. At that point the asset-twin gain mismatch has been removed to within 5 %.
- **Drift** is present whenever the net sensor-A offset error is above tol_b = 0.05, i.e. |injected drift − installed correction| > 0.05. The nominal e2a SD is 0.0386, so 0.05 is roughly one noise SD.
  - This is evaluated every sample, so an inaccurate slope estimate can make Drift return ("relapse").
  - A spurious correction applied in Nominal can create a Drift regime.
- **SensorNoisy** is repaired when the fault-induced noise in the fused twin output has an RMS ≤ 0.02 (half the unmitigated 0.041). The RMS is taken over the last 0.5 s (25 samples). All 25 samples must lie inside the fault episode and under a single, non-nominal configuration.
  - On success the noise source is cleared, which mirrors the Markov model's jump to Nominal.
  - Steady-state ratios of fused fault noise to the unmitigated level: DwSensorA 1.10, DwSensorB 1.32, IncFilter 0.32.
  - So **only IncFilter can satisfy the criterion** with the pipeline's two-sensor fault.
  - Sensitivity check: with noise on sensor A only, DwSensorA gives a ratio of 0.17.

## Experiments
- **Repair times** (`cl_repair.R`): the fault starts at 6 s and persists. The action is held from onset, or from onset + U(0,2) s, until the fault is repaired, and NoAction is used afterwards.
  - 500 held-out seeds (400001-400500) per pair. NoAction controls use the same seeds.
  - rho_hat = dt / mean repair time. A censoring-aware geometric MLE (repairs / steps at risk) is also reported.
  - Decoder-based metrics:
    - time from action start to the first run of 25 consecutive "Nominal" labels;
    - fraction of runs where such a run begins within 1 s after the repair;
    - fraction of Nominal labels in the 1 s after the repair.
- **Closed loop** (`cl_closed_loop.R`): 1000 held-out seeds (500001-501000) per controller, with common random numbers: the same asset, noise, fault types and gaps for every controller.
  - **Schedule:** the first onset is drawn from U(1,4) s. The fault type is uniform over SensorNoisy, DynamicsOff and Drift. Each fault persists until it is repaired or reaches Dmax = 8 s (self-clear). The next onset follows after U(1,4) s.
  - **Each sample:** compute residuals, then the causal decoder label z_k, then the controller's action a_k. The reward is R(a_k, s_k), with s_k the ground truth, discounted by 0.99 per sample over T = 20 s.
  - **Controllers:**
    - Oracle MDP: `vi$policy` applied to the true regime.
    - POMDP: PBVI alpha vectors, with belief b' ∝ O[,z]·(b A^(a)) per sample, using the paper's O.
    - Myopic belief: argmax R b.
    - MCDA/TOPSIS: the decision rule of `gillespie_mcda` (EMA alpha = 0.3, 0.5 s dwell, gate), applied per sample.
    - No Intervention.
  - **Variants:** `_aonly` (SensorNoisy noise on sensor A only), `_rawdecoder` (eps = 0), `_Wbc2s` (BiasCorrect window 2 s), `_Wrp4s` (ReidentPlant window 4 s).

## Rerun
```
Rscript analyses/closed_loop/run_closed_loop.R          # full
Rscript analyses/closed_loop/run_closed_loop.R --quick  # smoke test
```
Single runs are controlled by environment variables:
- `CL_N`: closed-loop seeds.
- `CL_N_REPAIR`: repair-experiment seeds.
- `CL_VARIANT=aonly`
- `CL_DEC_EPS`: decoder transition smoothing (default 1e-3).
- `CL_CFG="W_bc=100,W_rp=200"`: overrides any entry in `cl_cfg_default`.
- `CL_TAG`: suffix added to output file names.
