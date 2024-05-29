# Column generation for multistage stochastic mixed-integer nonlinear programs with discrete state variables

This repository contains code for the following case studies used in the paper *Column generation for multistage stochastic mixed-integer nonlinear programs with discrete state variables*:
- Multistage blending [`Blending/`] and 
- Mobile generator routing problem [`OPF/`] 

Each case study is an independent project with its own `Project.toml` files to manage dependencies.

## Setup

Before running the case studies, ensure you have set up the environment by downloading dependencies in the respective environment. This can be done by using the `Pkg.instantiate()` command.

We use Gurobi as the solver, which can be downloaded [here](https://www.gurobi.com/downloads/).

## Instances Overview

Each case study contains eight cases. For the blending case study, each case is a different combination of time periods/scenarios and number of input/output tanks. For the power distribution case study, each case is a different combination of time periods/scenarios and number of gensets. Five instances are solved for each of the eight cases. So instances 1-5 correspond to case 1, instances 6-10 to case 2, ..., and instances 36-40 to case 8.

## Model Types

Different models, i.e. the fullspace, CG, and CGCS are assigned the following numbers: `1`, `2`, and `3`, respectively. 

## Running a Case Study Instance

In any case study, to solve a specific instance `x` using model type `y`, follow these steps:
1. `cd` to the folder for the respective case study.
2. From the terminal, run the following command: `julia runs/main.jl y x`


For example, to solve instance 5 in the mobile generator routing case study using the CG method, run the following commands:

```
cd OPF/
julia runs/main.jl 2 5
```

