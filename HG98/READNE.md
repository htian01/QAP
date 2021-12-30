# HG98 Implementation

## Data preparation
Put the .dat of the problems to `./dataset/qaplib/problems_dat`


## Build

Simply `make` is enough. `openmp` is required

## Usage

`./main ${path_to_problems_dat} ${save path} ${problem name} ${iteration}`
for example
`./main ./dataset/qaplib/problems_dat ./dataset/qaplib/hg98_reparam_it2000 had16 2000`

## Format of reparam instance

- No incentive costs are added to the unary cost, all unary cost are 0.
- There are two values in first line of .dd files. The first one is the lower bound and the second one is the primal cost of the solution of HG98 leading problem.

## Implementation Details

- The implement is the DP variant of HG98, which only support balanced matching

- The DP would return part of bound to cost matrix, the return part follows a exp decay.  

  `float return_super_leader = super_leader * exp(-float(steps) / (float(max_steps) / 20.0));`

