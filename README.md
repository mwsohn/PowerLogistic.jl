# PowerLogistic.jl
Power for simple logistic regression. This is a Julia port of powerLogisticReg.R in powerMediation package authored
by Weiliang Qiu <stwxq@channing.harvard.edu>, which in turn was based on the formulae discussed in Hsieh et al. A simple method of sample size
calcuation for linear and logistic regression. _Statistic in Medicine_ 1998;17:1623-1634.

## Installation
To install this package, type

```jldoctest
julia> Pkg.clone("https://github.com/mwsohn/PowerLogistic.jl")
```

## Syntax

```jldoctest
    powerLogistic(n = 0, p1 = .0, p2 = .0, B = .0, OR = .0, alpha = 0.05, power = 0.8)
```

This function estimates power or sample size at type I error = `alpha` and power ≥ `power`.
All parameters are keyword options. If `n` is not specified (i.e., n = 0), it produces sample size.
Otherwise, it produces power of the given sample size. Other options are:

- `p1` = Pr(Y = 1) for a continous variable. `p1` = Pr(D|X=0), event rate at X = 0.
- `p2` = Pr(D|X=1), event rate at X = 1.
- `B` = Pr(X=1)
- `alpha`: Type I error (default = 0.05)
- `power`: Power (default = 0.8)

## Example
```jldoctest
julia> powerLogistic(p1=.5,OR=exp(0.405),power=0.95)
317

julia> powerLogistic(n = 317, p1=.5,OR=exp(0.405))
0.950...

julia> powerLogistic(p1 = .4, p2 = .5, B = .5, OR=exp(0.405), alpha = 0.05, power = 0.95)
1281

julia> powerLogistic(n = 1281, p1 = .4, p2 = .5, B = .5, OR=exp(0.405), alpha = 0.05)
0.950...
```
