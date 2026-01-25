"""
    MixedModelsSmallSample

Small-sample corrections for linear mixed models fitted by
[MixedModels.jl](https://juliastats.org/MixedModels.jl/stable/).

# Linear Mixed Model

A linear mixed model expresses the ``n``-vector of responses as

```math
y = X\\beta + Zu + \\varepsilon, \\quad
u \\sim \\mathcal{N}(0, G(\\theta)), \\quad
\\varepsilon \\sim \\mathcal{N}(0, \\sigma^2 I_n),
```

where ``X`` is the ``n \\times p`` fixed-effects design matrix, ``\\beta`` is the
``p``-vector of fixed-effects coefficients, ``Z`` is the random-effects design
matrix, ``u`` is the vector of random effects with covariance ``G(\\theta)``,
and ``\\varepsilon`` is the residual error.

The marginal variance of ``y`` is

```math
V(\\theta) = ZG(\\theta)Z^\\top + \\sigma^2 I_n.
```

# The Small-Sample Problem

The generalized least squares estimator ``\\hat{\\beta}`` has covariance

```math
\\Phi(\\theta) = (X^\\top V(\\theta)^{-1} X)^{-1}.
```

Standard inference replaces ``\\theta`` with its REML estimate ``\\hat{\\theta}``,
yielding ``\\hat{\\Phi} = \\Phi(\\hat{\\theta})``, and uses asymptotic approximations
(e.g., treating ``t``-statistics as normally distributed). For small samples,
this ignores:

1. **Bias** in ``\\hat{\\Phi}`` due to estimating ``\\theta``
2. **Uncertainty** in ``\\hat{\\theta}`` affecting the reference distribution

This package implements two corrections:

- **Kenward-Roger** ([`KenwardRoger`](@ref)): Applies a bias correction to
  ``\\hat{\\Phi}`` and uses a scaled F-distribution with approximate degrees
  of freedom.
- **Satterthwaite** ([`Satterthwaite`](@ref)): Uses moment-matching to
  approximate degrees of freedom for the t/F reference distribution.

# References

- Kenward, M. G. and Roger, J. H. (1997). "Small sample inference for fixed
  effects from restricted maximum likelihood." *Biometrics*, 53(3), 983-997.
- Satterthwaite, F. E. (1946). "An Approximate Distribution of Estimates of
  Variance Components." *Biometrics Bulletin*, 2(6), 110-114.
- Fai, A. H.-T. and Cornelius, P. L. (1996). "Approximate F-tests of multiple
  degree of freedom hypotheses in generalized least squares analyses of
  unbalanced split-plot experiments." *Journal of Statistical Computation
  and Simulation*, 54(4), 363-378.

# See Also

- [`small_sample_adjust`](@ref): Apply KR or SW adjustment to a fitted model
- [`small_sample_ftest`](@ref): Perform an adjusted F-test for contrasts
"""
module MixedModelsSmallSample

using MixedModels
using Distributions
using LinearAlgebra

using StatsAPI: StatsAPI, coeftable
using StatsBase: StatsBase, CoefTable
using Markdown

# Exports
export small_sample_adjust
export small_sample_ftest
export KenwardRoger
export Satterthwaite
export ObservedFIM
export ExpectedFIM
export Unstructured
export FactorAnalytic
export LinearMixedModelAdjusted
export vcov_varpar

# Include submodules in order (types first, then helpers, then main functions)
include("types.jl")
include("utils.jl")
include("dof_computation.jl")
include("parameterization.jl")
include("fisher_information.jl")
include("adjustment.jl")
include("ftest.jl")
include("show.jl")

end
