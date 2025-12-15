# On the Graphical Rules for Recovering the Average Treatment Effect Under Selection Bias
## Yichi Zhang, Haidong Lu

R simulation code for addressing selection bias and recovering the average treatment effect (ATE) under graphical rules.

This repository contains R code to reproduce the simulation studies reported in:

Zhang Y, Lu H. *On the Graphical Rules for Recovering the Average Treatment Effect Under Selection Bias.*

## Files

- `simulation_selectionbias_v3_case1.R`  
  Simulation code for Case 1, illustrating selection bias in a randomized trial of a prophylactic medication for malaria among susceptible individuals.

- `simulation_selectionbias_v3_case2.R`  
  Simulation code for Case 2, illustrating selection bias in a randomized trial of a vaccine for an infectious disease, with selection affected by a mediator.

## How to run

```r
source("simulation_selectionbias_v3_case1.R")
source("simulation_selectionbias_v3_case2.R")

## System details

All simulations were conducted using R version 4.5.2.
