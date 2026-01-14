"""
    needs_d2V(::AbstractLinearMixedModelAdjustment) -> Bool

Trait indicating whether the adjustment method requires second derivatives of ``V``.

Returns `true` for [`KenwardRoger`](@ref) (needed for bias correction),
`false` for [`Satterthwaite`](@ref).
"""
needs_d2V(::KenwardRoger) = true
needs_d2V(::Satterthwaite) = false

"""
    compute_adjusted_covariance(Φ, vd, W, method) -> Matrix

Compute the (possibly adjusted) covariance matrix of fixed effects.

For [`KenwardRoger`](@ref), applies the bias correction:

```math
\\Phi_A = \\hat{\\Phi} + 2\\hat{\\Phi}\\,B\\,\\hat{\\Phi},
```

where the bias matrix ``B`` accounts for the variability in estimating ``\\theta``:

```math
B = \\sum_{i,j} W_{ij} \\bigl(Q_{ij} - P_i\\hat{\\Phi}P_j + R_{\\text{factor}}\\cdot R_{ij}\\bigr).
```

Here ``W`` is the asymptotic covariance of ``\\hat{\\theta}``, and ``P_i, Q_{ij}, R_{ij}``
are defined in [`VarianceDecomposition`](@ref). The `R_factor` is `-1/4` for
[`FactorAnalytic`](@ref) parameterization and `1` for [`Unstructured`](@ref).

For [`Satterthwaite`](@ref), returns ``\\Phi`` unchanged (no bias correction).

# Arguments
- `Φ`: Unadjusted covariance matrix ``(X^\\top V^{-1} X)^{-1}``
- `vd`: [`VarianceDecomposition`](@ref) containing derivative matrices
- `W`: Asymptotic covariance of variance parameters
- `method`: Adjustment method ([`KenwardRoger`](@ref) or [`Satterthwaite`](@ref))
"""
function compute_adjusted_covariance(Φ, vd, W, ::KenwardRoger)
    p = size(Φ, 1)
    B = zeros(p, p)
    for i in eachindex(vd.θs), j in eachindex(vd.θs)
        B += W[i, j] * (vd.Q[i, j] - vd.P[i] * Φ * vd.P[j] + vd.R_factor * vd.R[i, j])
    end
    return Φ + 2 * Φ * B * Φ
end

compute_adjusted_covariance(Φ, vd, W, ::Satterthwaite) = Φ

"""
    small_sample_adjust(m::MixedModel, method) -> LinearMixedModelAdjusted

Apply small-sample corrections to a fitted linear mixed model.

This is the main entry point for computing adjusted inference. The function:

1. Validates that `m` was fitted with REML and is unweighted
2. Computes ``W = I(\\theta)^{-1}``, the asymptotic covariance of ``\\hat{\\theta}`` (see [`vcov_varpar`](@ref))
3. For [`KenwardRoger`](@ref): computes adjusted covariance ``\\Phi_A``
   (see [`compute_adjusted_covariance`](@ref))
4. Computes denominator degrees of freedom ``\\nu_k`` for each fixed effect
   (see [`compute_dof`](@ref))

# Arguments
- `m::MixedModel`: A fitted `LinearMixedModel` from MixedModels.jl (must use REML)
- `method`: Either [`KenwardRoger`](@ref) or [`Satterthwaite`](@ref)

# Returns
- [`LinearMixedModelAdjusted`](@ref)

The returned object can be displayed to show adjusted coefficient tables with
corrected standard errors (KR only) and degrees of freedom.

# Example

```julia
using MixedModels, MixedModelsSmallSample
m = fit(MixedModel, @formula(y ~ x + (1|g)), data; REML=true)
kr = small_sample_adjust(m, KenwardRoger())
sw = small_sample_adjust(m, Satterthwaite())
```

# See Also
- [`small_sample_ftest`](@ref): Perform F-tests on contrasts
"""
function small_sample_adjust(m::MixedModel, method::AbstractLinearMixedModelAdjustment)
    validation(m)

    X = m.X
    Φ = m.vcov

    # Get variance decomposition
    include_d2V = needs_d2V(method)
    vd = VarianceDecomposition(m, X, method.parameterization; include_d2V)

    # Compute W
    W = _vcov_varpar(vd, m; fim=method.fim)

    # Compute adjusted covariance
    Φ_A = compute_adjusted_covariance(Φ, vd, W, method)

    # Construct result with placeholder ν
    m_adj = _construct_result(method, m, Φ_A, W, vd, Float64[])

    # Compute DoF for each fixed effect using type dispatch
    ν = compute_per_coefficient_dof(m_adj)

    # Return with correct ν
    return _set_dof(m_adj, ν)
end

function _construct_result(method::KenwardRoger, m, Φ_A, W, vd, ν)
    return LinearMixedModelAdjusted(m, method, Φ_A, W, vd.P, vd.Q, vd.R, vd.R_factor, ν)
end

function _construct_result(method::Satterthwaite, m, Φ_A, W, vd, ν)
    return LinearMixedModelAdjusted(m, method, Φ_A, W, vd.P, nothing, nothing, 0.0, ν)
end
