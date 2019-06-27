# Assemble

The **assemble** function aims to provide an `ODEProblem` encapsulating
  all the necessary information for modeling earthquake cycles. It will
  automatically generate caches based on given mesh sizes aiming to minimize
  allocation during solving. The full functionality of `ODEProblem`
  as well as solving interface can be viewed at documents of
  [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/index.html).

Also, please checkout [RecursiveArrayTools.jl](https://github.com/JuliaDiffEq/RecursiveArrayTools.jl) for usage of `ArrayPartition`
which is adopted by default for multiple variables in this package.

## Public Interface
```@autodocs
Modules = [JuEQ]
Pages = ["assemble.jl"]
Private = false
Order = [:type, :function, :constant, :macro]
```

## References
Rice, J. (1993). Spatio-temporal complexity of slip on a fault. Journal of Geophysical Research: Solid Earth, 98(B6), 9885–9907. https://doi.org/10.1029/93JB00191

Liu, Y., McGuire, J. J., & Behn, M. D. (2012). Frictional behavior of oceanic transform faults and its influence on earthquake characteristics. Journal of Geophysical Research: Solid Earth, 117(B4). https://doi.org/10.1029/2011JB009025

Barbot, S. (2018). Asthenosphere Flow Modulated by Megathrust Earthquake Cycles. Geophysical Research Letters, 45(12), 6018–6031. https://doi.org/10.1029/2018GL078197
