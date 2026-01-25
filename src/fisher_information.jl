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

    # REML vs ML differs in the information for θ.
    # Let P = V⁻¹ - V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹.
    # - REML: trace terms use P
    # - ML:   trace terms use V⁻¹ (no log|X'V⁻¹X| term)
    # The quadratic form part uses P in both cases (β is estimated by GLS).
    P = Vinv - Vinv * X * Φ * X' * Vinv
    Ptrace = m.optsum.REML ? P : Vinv
    Pquad = P

    FIM = Matrix{Float64}(undef, nparams, nparams)
    for i in 1:nparams, j in 1:nparams
        val =
            1/2 * tr(-Ptrace * dVs[i] * Ptrace * dVs[j]) -
            1/2 * (y - X * β)' * Vinv * (-2 * dVs[i] * Pquad * dVs[j]) * Vinv * (y - X * β)

        # Add second derivative terms for Observed FIM if available
        # Correction: + 1/2 tr(P_trace d2V) - 1/2 r' Vinv d2V Vinv r
        if !isnothing(vd.d2Vs) && sum(abs, vd.d2Vs[i, j]) > 0
            val +=
                1/2 * tr(Ptrace * vd.d2Vs[i, j]) -
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

    m.optsum.REML || error(
        "ExpectedFIM is currently only supported for REML fits. " *
        "If you need ExpectedFIM for ML fits, please open an issue on GitHub.",
    )

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

Compute `W = I(θ)^{-1}`, the asymptotic covariance of the fitted variance parameters.

The Fisher information `I(θ)` is computed using either [`ObservedFIM`](@ref) or [`ExpectedFIM`](@ref).
"""
function vcov_varpar(
    m::MixedModel;
    fim::AbstractFIM=ObservedFIM(),
    parameterization::AbstractParameterization=Unstructured(),
)
    isempty(m.sqrtwts) || error("Cannot compute for weighted models")
    # Always compute d2V for ObservedFIM to ensure correct Hessian
    include_d2V = fim isa ObservedFIM
    vd = VarianceDecomposition(m, m.X, parameterization; include_d2V)
    return _vcov_varpar(vd, m; fim)
end
