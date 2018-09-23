# JuEQ.jl
[![BuildStatus](https://travis-ci.com/shipengcheng1230/JuEQ.jl.svg?token=zsZu59CsqQTTp7wzi7zP&branch=master)](https://travis-ci.com/shipengcheng1230/JuEQ.jl)
[![codecov.io](https://codecov.io/gh/shipengcheng1230/JuEQ.jl/coverage.svg?token=ag6kv61zOW&branch=master)](https://codecov.io/gh/shipengcheng1230/JuEQ.jl?branch=master)
[![License:MIT](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://shipengcheng1230.github.io/JuEQ.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://shipengcheng1230.github.io/JuEQ.jl/latest)


This is a suite for numerically simulating earthquake sequences in [Julia](https://julialang.org/). The purpose of this package is to provide efficient Julia implementations for simulations in the field of earthquake physics.

Features of this package are listed as below:

- [x] Rate-State Friction Law
- [x] Okada's Dislocation Method
- [x] Boundary Element Method (Quasi-dynamic)


Features to be implemented:
- [ ] Viscoelastic Relaxation
- [ ] Fully Elastodynamic Effect
- [ ] Off-Fault Materials effect

## Developing
[Demos](https://github.com/shipengcheng1230/JuEQ.jl/tree/master/demos) contains the most primitive implementation of simulation algorithms.

**JuEQ** is still in alpha-stage. Breaking changes are imminent.


## Supporting
This software in this ecosystem is developed as part of academic research in
[Earthquake Physics Lab](http://weilab.uri.edu/) at
[Graduate School of Oceanography](https://web.uri.edu/gso/), University of Rhode Island.


## References

Please consider to cite the following papers if you find this package useful.

* Rice, J. (1993). Spatio-temporal complexity of slip on a fault. Journal of Geophysical Research: Solid Earth, 98(B6), 9885–9907. https://doi.org/10.1029/93JB00191

* Liu, Y., & Rice, J. R. (2005). Aseismic slip transients emerge spontaneously in three-dimensional rate and state modeling of subduction earthquake sequences. Journal of Geophysical Research: Solid Earth, 110(B8). https://doi.org/10.1029/2004JB003424
