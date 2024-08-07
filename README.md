# paper_CEE-efficient-ML

Code for paper "Efficient and Globally Robust Causal Excursion Effect Estimation" by Zhaoxi Cheng, Lauren Bell, Tianchen Qian

## Files

-   application code: code for results in the "Application" section.
-   simulation code: code for results in the "Simulation" section.

## How to use the code to replicate results in the paper

-   For results in the "Application" section, run each R script named "analysis Drink Less -" of each outcome type (binary, continuous, and count) in "application code" folder. Then, run "analysis Drink Less - make plot.R" file to make the figure. The data for the binary outcome is available [here](https://osf.io/mtcfa). The data for the continuous outcome and the count outcome is not publicly available.
-   For results in the "Simulation" section, run the R scripts in each folder of the outcome types with the following orders in "simulation code" folder. Then, run the R script from the "plot making" folder within each outcome type folder.

| Figure  | Steps to make the figure |
| ------------- | ------------- |
| 1  | 1\. run sh bash.sh in simulation code/continuous outcome/datasets   |
|  | 2\. run sh bash.sh in simulation code/continuous outcome/simu-control_pattern |
|  | 3\. run plot-making-cont-control_pattern.R in simulation code/continuous outcome/plot making/continuous - eif - control_pattern   |
| 2 | 1\. run sh bash.sh in simulation code/continuous outcome/datasets |
|  | 2\. run sh bash.sh in simulation code/continuous outcome/simu-error_var_pattern |
|  | 3\. run plot-making-cont-error_var_pattern.R in simulation code/continuous outcome/plot making/continuous - eif - error_var_pattern  |
| 3 | 1\. run sh bash.sh in simulation code/binary outcome/datasets  |
| | 2\. run sh bash.sh in simulation code/binary outcome/simu-control_pattern |
| | 3\. run plot-making-bin-control_pattern.R in simulation code/binary outcome/plot making |
| 4 | 1\. run sh bash.sh in simulation code/count outcome/datasets |
| | 2\. run sh bash.sh in simulation code/count outcome/simu-control_pattern   |
| | 3\. run plot-making-count-control_pattern.R in simulation code/count outcome/plot making |
| 5 | 1\. run analysis Drink Less - binary outcome.R in application code/analysis Drink Less - binary outcome, analysis Drink Less - continuous outcome.R in application code/analysis Drink Less - continuous outcome, and analysis Drink Less - count outcome.R in application code/analysis Drink Less - count outcome independently.  |
| | 2\. run analysis Drink Less - make plot.R in application code.  |
