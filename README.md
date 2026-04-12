Notebooks provided to generate the results in the manuscript "Optimal Sequential Decision-Making for Error Propagation Mitigation in Modular Digital Twins" by Najafi et al.

The figure notebooks set their working directory to the parent
  `ErrorPropagation_MDP/` folder and source `POMDPFormulation.R` from
  there. That source pulls in `BaseMDPFormulation.R` and
  `SimulationPipeline.R`, which run the physics to ARX to HMM pipeline
  and export the data-driven baseline transition matrix `A_NoAction` and
  the confusion matrix `mcda_confusion`.

  `MDP_POMDP_Formulation.R` in this folder is a self-contained reference
  copy that compiles the same MDP / POMDP definitions, value iteration,
  PBVI solver, and Gillespie simulators in a single file. It is not
  sourced by the notebooks but is provided for readers who want to study
  the methods without navigating the full pipeline.

  ## Reproducing the figures

  1. Clone or download the full `ErrorPropagation_MDP/` directory so that
     `POMDPFormulation.R`, `BaseMDPFormulation.R`, and
     `SimulationPipeline.R` are present in the parent of this folder.
  2. Install the R dependencies (see below).
  3. Knit any notebook from RStudio (`Knit` button) or from the command
     line:

  ```sh
  Rscript -e "rmarkdown::render('Fig01_MDP_Validation.Rmd')"

  Each notebook writes its HTML output next to the source file. The two
  heavy notebooks (Fig04_Sensitivity_VoI.Rmd and
  ModelFreeRL_Significance.Rmd) cache their training and Monte Carlo
  chunks, so the first knit takes a few minutes and subsequent knits are
  fast.

  Dependencies

  R 4.5 or newer. CRAN packages:

  ggplot2     reshape2    expm        gridExtra
  scales      knitr       depmixS4    clue
  dplyr       tibble      kableExtra

  depmixS4 and clue are pulled in transitively by
  SimulationPipeline.R to fit the HMM and Hungarian-map the Viterbi
  decoder against the ground-truth regime labels.
