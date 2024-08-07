# Notes on the application code 

## Differences compared to functions in simulation

1. estimating the nuisance function: 
  - user-defined variables in the spline terms for gam.
  - stacking method available for estimating the nuisance function.
2. more choices to model d term:
  - allow linear regression and spline with exponential in t term.
  - random forest to model d term
  - gam to model d term
3. code structure:
  - use try() in emee_eif() to re-draw a sample if an error occurs for cross-fitting.