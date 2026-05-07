# paper_CEE-efficient-ML

Code for the paper "Efficient and Globally Robust Causal Excursion Effect
Estimation" by Zhaoxi Cheng, Lauren Bell, Tianchen Qian.

## Folders

- `application code/` — code for the "Application: Drink Less MRT" section.
- `simulation code/` — code for the "Simulation" section and Appendix G.
  - `5.1 continuous outcome/` — Section 5.1 (continuous outcome, marginal CEE).
  - `5.2 time-varying CEE/` — Section 5.2 (time-moderated CEE).
  - `5.3 many decision points/` — Section 5.3 (large-T regime, marginal CEE).
  - `G.1 binary outcome/` — Appendix G.1 (binary outcome, marginal CEE).
  - `G.2 count outcome/` — Appendix G.2 (count outcome, marginal CEE).
- `other code/` — code for technical artifacts in the appendices.
  - `A.3 efficiency comparison counter-example/` — Appendix A.3 counter-example.
  - `D.1 wallclock benchmark/` — Appendix D.1 wall-clock benchmark (Table D.1).

## How to use the code to replicate results in the paper

For the **application** results, run each `analysis Drink Less - <outcome>.R`
script in the corresponding subfolder of `application code/`, then run
`application code/analysis Drink Less - make plot.R` to make the figure.
The data for the binary outcome is available [here](https://osf.io/mtcfa).
The data for the continuous outcome and the count outcome is not publicly
available.

For the **simulation** results, the workflow follows the same SLURM-based
template across folders: a data-generation step, an estimation step, and a
plot-making step. The figure-by-figure entry points are below.

| Figure / Table | Steps to reproduce |
| - | - |
| Figure 1 (paper §5.1) | 1. `simulation code/5.1 continuous outcome/datasets/`: run `sh batch.sh` <br> 2. `simulation code/5.1 continuous outcome/simu-control_pattern/`: run `sh batch.sh` <br> 3. `simulation code/5.1 continuous outcome/plot making/continuous - eif - control_pattern/`: run `plot-making-cont-control_pattern.R` |
| Figure 2 (paper §5.1) | 1. `simulation code/5.1 continuous outcome/datasets/`: run `sh batch.sh` <br> 2. `simulation code/5.1 continuous outcome/simu-error_var_pattern/`: run `sh batch.sh` <br> 3. `simulation code/5.1 continuous outcome/plot making/continuous - eif - error_var_pattern/`: run `plot-making-cont-error_var_pattern.R` |
| Figure 3 (paper §5.2) | In `simulation code/5.2 time-varying CEE/`: <br> 1. `bash batch_simulation.sh` (submits 360 SLURM tasks) <br> 2. `Rscript collect_results.R` <br> 3. `Rscript make_plot.R` |
| Figure 4 (paper §5.3) | In `simulation code/5.3 many decision points/`: <br> 1. `bash batch_simulation.sh` (submits 100 SLURM tasks) <br> 2. `Rscript collect_results.R` <br> 3. `Rscript make_plot.R` |
| Figure 5 (paper §6, application) | 1. In `application code/`, run `analysis Drink Less - binary outcome.R`, `analysis Drink Less - continuous outcome.R`, and `analysis Drink Less - count outcome.R` independently. <br> 2. Run `analysis Drink Less - make plot.R`. |
| Figure G.1 (paper §G.1) | 1. `simulation code/G.1 binary outcome/datasets/`: run `sh batch.sh` <br> 2. `simulation code/G.1 binary outcome/simu-control_pattern/`: run `sh batch.sh` <br> 3. `simulation code/G.1 binary outcome/plot making/`: run `plot-making-bin-control_pattern.R` |
| Figure G.2 (paper §G.2) | 1. `simulation code/G.2 count outcome/datasets/`: run `sh batch.sh` <br> 2. `simulation code/G.2 count outcome/simu-control_pattern/`: run `sh batch.sh` <br> 3. `simulation code/G.2 count outcome/plot making/`: run `plot-making-count-control_pattern.R` |
| Table D.1 (paper §D.1) | In `other code/D.1 wallclock benchmark/`: run `run_benchmark.R`. The continuous Drink Less data is required and is not publicly available. |
| Counter-example (paper §A.3) | In `other code/A.3 efficiency comparison counter-example/`: run `sim_counterexample_v2.R`. |
