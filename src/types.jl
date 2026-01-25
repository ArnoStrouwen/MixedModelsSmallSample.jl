# Types for MixedModelsSmallSample

# Abstract Types for Options
abstract type AbstractFIM end
abstract type AbstractParameterization end
abstract type AbstractLinearMixedModelAdjustment end

# FIM Options
"""
    ObservedFIM

Use the observed (Hessian-based) Fisher Information Matrix for the variance parameters.

Let ``V``, ``V^{-1}``, ``dV_i`` and ``d^2V_{ij}`` come from [`VarianceDecomposition`](@ref). Define

```math
\\Phi = (X^\\top V^{-1} X)^{-1},\\qquad
P = V^{-1} - V^{-1}X\\,\\Phi\\,X^\\top V^{-1},\\qquad
r = y - X\\hat{\\beta}.
```

Then the observed information is

```math
I^{\\text{obs}}_{ij}(\\theta) =
-\\frac{1}{2}\\operatorname{tr}\\left(P dV_i P dV_j\\right)
+ r^{\\top} V^{-1} dV_i P dV_j V^{-1} r
+ \\frac{1}{2}\\operatorname{tr}\\left(P d^2V_{ij}\\right)
- \\frac{1}{2} r^{\\top} V^{-1} d^2V_{ij} V^{-1} r.
```


# See Also
- [`ExpectedFIM`](@ref): Alternative using expected values
- [`vcov_varpar`](@ref): Computes the inverse of the FIM
"""
struct ObservedFIM <: AbstractFIM end

"""
    ExpectedFIM

Use the expected Fisher Information Matrix for the variance parameters.

Let ``V^{-1}``, ``dV_i``, ``P_i`` and ``Q_{ij}`` come from [`VarianceDecomposition`](@ref), and define
``\\Phi = (X^\\top V^{-1} X)^{-1}``.

Then the expected information used by [`vcov_varpar`](@ref) is

```math
I^{\\text{exp}}_{ij}(\\theta) = \\frac{1}{2}\\operatorname{tr}\\left(V^{-1} dV_i V^{-1} dV_j\\right)
- \\operatorname{tr}(\\Phi\\,Q_{ij})
+ \\frac{1}{2}\\operatorname{tr}(\\Phi\\,P_i\\,\\Phi\\,P_j).
```

# See Also
- [`ObservedFIM`](@ref): Alternative using observed Hessian
- [`vcov_varpar`](@ref): Computes the inverse of the FIM
"""
struct ExpectedFIM <: AbstractFIM end

# Parameterization Options
"""
    Unstructured

Unstructured (variance component) parameterization.

In this parameterization, the marginal variance is linear in the variance parameters.
Writing the model as ``V(\\theta) = \\sigma^2 I_n + \\sum_b Z_b G_b(\\theta) Z_b^\\top``,
`Unstructured` parameterizes each block covariance ``G_b(\\theta)`` by its unique
(co)variance components.

## Gradient Structure

Because ``V`` is linear in ``\\theta``, the derivatives are trivial: each
``\\partial V/\\partial\\theta_i`` is constant and
``\\partial^2 V/(\\partial\\theta_i\\,\\partial\\theta_j) = 0``.

This implies the ``R_{ij}`` terms used by the Kenward-Roger covariance correction vanish.

# See Also
- [`FactorAnalytic`](@ref): Alternative Cholesky-based parameterization
"""
struct Unstructured <: AbstractParameterization end

"""
    FactorAnalytic

Factor analytic (Cholesky) parameterization.

For each random-effects block ``b``, decompose the block covariance as
``G_b = L_b L_b^\\top`` with ``L_b`` lower triangular. The variance parameters are
``\\sigma^2`` and the free elements of the ``L_b``.

The marginal variance is

```math
V = \\sigma^2 I_n + \\sum_b Z_b G_b Z_b^\\top = \\sigma^2 I_n + \\sum_b Z_b L_b L_b^\\top Z_b^\\top.
```

For MixedModels.jl models, the internal random-effects factor is stored as a
relative factor ``\\lambda_b``; we use the scaled factor ``L_b = \\sigma\\,\\lambda_b``.

## Derivatives of V

The derivative of ``V`` with respect to a Cholesky element ``L_{b,ij}`` uses the
chain rule through ``G_b = L_b L_b^\\top``:

```math
\\frac{\\partial V}{\\partial L_{ij}} = Z_b \\frac{\\partial G_b}{\\partial L_{ij}} Z_b^\\top,
\\qquad
\\frac{\\partial G_b}{\\partial L_{ij}} = E_{ij} L_b^\\top + L_b E_{ij}^\\top,
```

where ``E_{ij}`` is the elementary matrix with a 1 at position ``(i,j)``.

The second derivatives are non-zero. For indices within the same block,

```math
\\frac{\\partial^2 G_b}{\\partial L_{i_1 j_1}\\,\\partial L_{i_2 j_2}}
= \\delta_{j_1 j_2}\\,(E_{i_1 i_2} + E_{i_2 i_1}),
```

and hence

```math
\\frac{\\partial^2 V}{\\partial L_{i_1 j_1}\\,\\partial L_{i_2 j_2}}
= Z_b\\,\\frac{\\partial^2 G_b}{\\partial L_{i_1 j_1}\\,\\partial L_{i_2 j_2}}\\,Z_b^\\top.
```

Unlike [`Unstructured`](@ref), this parameterization has **non-zero second
derivatives** ``\\partial^2 V / \\partial L_{ij} \\partial L_{kl}``, which must
be included in the Kenward-Roger bias correction.

# See Also
- [`Unstructured`](@ref): Alternative with zero second derivatives
"""
struct FactorAnalytic <: AbstractParameterization end

# Adjustment Options
"""
    KenwardRoger(; fim=ObservedFIM(), parameterization=Unstructured())

Kenward-Roger method for [`small_sample_adjust`](@ref) and [`small_sample_ftest`](@ref).

See also [`Satterthwaite`](@ref) for an alternative method.
"""
struct KenwardRoger{F<:AbstractFIM,P<:AbstractParameterization} <:
       AbstractLinearMixedModelAdjustment
    fim::F
    parameterization::P
end

function KenwardRoger(;
    fim::AbstractFIM=ObservedFIM(),
    parameterization::AbstractParameterization=Unstructured(),
)
    KenwardRoger(fim, parameterization)
end

"""
    Satterthwaite(; fim=ObservedFIM(), parameterization=Unstructured())

Satterthwaite method for [`small_sample_adjust`](@ref) and [`small_sample_ftest`](@ref).

See also [`KenwardRoger`](@ref) for an alternative method.
"""
struct Satterthwaite{F<:AbstractFIM,P<:AbstractParameterization} <:
       AbstractLinearMixedModelAdjustment
    fim::F
    parameterization::P
end

function Satterthwaite(;
    fim::AbstractFIM=ObservedFIM(),
    parameterization::AbstractParameterization=Unstructured(),
)
    Satterthwaite(fim, parameterization)
end

"""
    VarianceDecomposition

Cached matrices derived from the marginal covariance `V(θ)` and its derivatives.

Let `θs` denote the variance parameters, and write
`dVs[i] = ∂V/∂θ_i` and `d2Vs[i,j] = ∂²V/(∂θ_i∂θ_j)`.

This struct stores `V`, `Vinv`, `θs`, `dVs`, `d2Vs`, and the derived matrices

```math
P_i = -X^\\top V^{-1} dV_i V^{-1} X,
\\qquad
Q_{ij} = X^\\top V^{-1} dV_i V^{-1} dV_j V^{-1} X,
\\qquad
R_{ij} = X^\\top V^{-1} d^2V_{ij} V^{-1} X,
```
and `R_factor`.
`R_factor` is a parameterization-dependent scalar applied to `R_{ij}` in the
Kenward-Roger covariance bias correction (see [`compute_adjusted_covariance`](@ref)).

The definition of `θ` (and thus the rest of the calculated quantities here) depends on the chosen
parameterization; see [`Unstructured`](@ref) and [`FactorAnalytic`](@ref).
"""
struct VarianceDecomposition
    V::Matrix{Float64}
    Vinv::Matrix{Float64}
    θs::Vector{Float64}
    dVs::Vector{Matrix{Float64}}
    d2Vs::Matrix{Matrix{Float64}}
    P::Vector{Matrix{Float64}}
    Q::Matrix{Matrix{Float64}}
    R::Matrix{Matrix{Float64}}
    R_factor::Float64
    names::Vector{String}
end

"""
    LinearMixedModelAdjusted{T<:LinearMixedModel, O<:AbstractLinearMixedModelAdjustment}

A wrapper around a `LinearMixedModel` that holds the small-sample adjusted results.

# Fields
- `m`: The original `LinearMixedModel`.
- `adj`: The adjustment option used ([`KenwardRoger`](@ref) or [`Satterthwaite`](@ref)).
- `varcovar_adjusted`: The adjusted variance-covariance matrix of the fixed effects (corrected for KR, unadjusted for SW).
- `W`: The asymptotic covariance matrix of the variance parameters.
- `P`: Derived matrix \\(P_i\\).
- `Q`: Derived matrix \\(Q_{ij}\\) (optional, `nothing` for Satterthwaite).
- `R`: Derived matrix \\(R_{ij}\\) (optional, `nothing` for Satterthwaite).
- `R_factor`: Scaling factor for \\(R_{ij}\\) (0.0 for Satterthwaite).
- `ν`: The denominator degrees of freedom for each fixed effect coefficient.
"""
struct LinearMixedModelAdjusted{T<:LinearMixedModel,O<:AbstractLinearMixedModelAdjustment}
    m::T
    adj::O
    varcovar_adjusted::Matrix{Float64}
    W::Matrix{Float64}
    P::Vector{Matrix{Float64}}
    Q::Union{Matrix{Matrix{Float64}},Nothing}
    R::Union{Matrix{Matrix{Float64}},Nothing}
    R_factor::Float64
    ν::Vector{Float64}
end

# Helper to update ν in adjusted models (used during construction)
function _set_dof(m::LinearMixedModelAdjusted, ν)
    LinearMixedModelAdjusted(
        m.m, m.adj, m.varcovar_adjusted, m.W, m.P, m.Q, m.R, m.R_factor, ν
    )
end
