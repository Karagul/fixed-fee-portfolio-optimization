# Fixed-fee Portfolio Optimization
Two separate approaches to solving portfolio optimization problems where assets have a fixed transaction fee, or a minimize investment size per asset


## Iterative_fee_reweighted_optimization.m
This is a heuristic approach to solving portfolio optimization problems where assets have a fixed fee for investing in each asset.  In practice, 
standard portfolio optimization techniques may allocate small amounts to a large number of assets when the universe of investible assets is large.
This creates problems in practice where often investing in assets can have a fixed fee.  This heuristic iterative approach approximates a solution.
The approach was inspired by the following paper:
https://faculty.washington.edu/mfazel/portfolio-final.pdf


## miqp.m
This MATLAB function solves Mixed Integer Quadratic Programming problems where in the universe of investible assets, there a specific minimum bound 
for how much can be invested in a given asset.  




