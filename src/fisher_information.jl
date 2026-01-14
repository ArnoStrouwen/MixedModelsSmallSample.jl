using LinearAlgebra
using MixedModels

function _vcov_varpar(
    vd::VarianceDecomposition, m::MixedModel; fim::AbstractFIM=ObservedFIM()
)
    return inv(_compute_fim(fim, vd, m))
end

function _compute_fim(::ObservedFIM, vd, m)
    Φ = m.vcov
    X = m.X
    y = m.y
    β = m.β
    dVs = vd.dVs
    Vinv = vd.Vinv
    nparams = length(vd.θs)

    Pvcov = Vinv - Vinv * X * Φ * X' * Vinv

    FIM = Matrix{Float64}(undef, nparams, nparams)
    for i in 1:nparams, j in 1:nparams
        val =
            1/2 * tr(-Pvcov * dVs[i] * Pvcov * dVs[j]) -
            1/2 * (y - X * β)' * Vinv * (-2 * dVs[i] * Pvcov * dVs[j]) * Vinv * (y - X * β)

        # Add second derivative terms for Observed FIM if available
        # Formula for Information Matrix correction: + 1/2 tr(Vinv d2V) - 1/2 r' Vinv d2V Vinv r
        if !isnothing(vd.d2Vs) && sum(abs, vd.d2Vs[i, j]) > 0
            val +=
                1/2 * tr(Pvcov * vd.d2Vs[i, j]) -
                1/2 * (y - X * β)' * Vinv * vd.d2Vs[i, j] * Vinv * (y - X * β)
        end
        FIM[i, j] = val
    end
    return FIM
end

function _compute_fim(::ExpectedFIM, vd, m)
    Φ = m.vcov
    dVs = vd.dVs
    Vinv = vd.Vinv
    P = vd.P
    Q = vd.Q
    nparams = length(vd.θs)

    FIM = Matrix{Float64}(undef, nparams, nparams)
    for i in 1:nparams, j in 1:nparams
        val =
            1/2 * tr(Vinv * dVs[i] * Vinv * dVs[j]) - tr(Φ * Q[i, j]) +
            1/2 * tr(Φ * P[i] * Φ * P[j])
        FIM[i, j] = val
    end
    return FIM
end

"""
    vcov_varpar(m::MixedModel; fim=ObservedFIM(), parameterization=Unstructured())

Compute `W = I(θ)^{-1}`, the asymptotic covariance of the fitted variance parameters under REML.

The Fisher information `I(θ)` is computed using either [`ObservedFIM`](@ref) or [`ExpectedFIM`](@ref).
"""
function vcov_varpar(
    m::MixedModel;
    fim::AbstractFIM=ObservedFIM(),
    parameterization::AbstractParameterization=Unstructured(),
)
    validation(m)
    # Always compute d2V for ObservedFIM to ensure correct Hessian
    include_d2V = fim isa ObservedFIM
    vd = VarianceDecomposition(m, m.X, parameterization; include_d2V)
    return _vcov_varpar(vd, m; fim)
end
