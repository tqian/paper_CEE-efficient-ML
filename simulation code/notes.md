# Notes on the simulation code 

## Differences across functions in the "functions" folder in each outcome type folder

1. estimating the nuisance function: 
  - functions in binary outcome folder use family = binomial; functions in count outcome folder use family = poisson
  - functions in count outcome folder allow stacking as a method to estimate the nuisance function
  