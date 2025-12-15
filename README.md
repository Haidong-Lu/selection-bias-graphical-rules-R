# selection-bias-graphical-rules-R
## Yichi Zhang, Haidong Lu

R simulation code for addressing selection bias and recovering ATE under graphical rules

This repository contains R code to reproduce the simulation studies in:

Zhang Y & Lu H. *On the Graphical Rules for Recovering the Average Treatment Effect Under Selection Bias.*

## Files

- simulation_selectionbias_v3_case1.R  
  Simulation for Case 1 regarding selection bias in a randomized trial of a prophylactic medication for malaria among susceptible individuals

- simulation_selectionbias_v3_case2.R  
  Simulation for Case 2 regarding selection bias in a randomized trial of some vaccine for an infectious disease that assigns half of the individuals to the treatment groug

## How to run

```r
source("simulation_selectionbias_v3_case1.R")
source("simulation_selectionbias_v3_case2.R")
