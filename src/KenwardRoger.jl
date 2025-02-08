module KenwardRoger

using ForwardDiff
using MixedModels
using Distributions
using LinearAlgebra

include("struct.jl")

export kenwardroger_matrices, kenwardroger_estimates

function kenwardroger_matrices(m::MixedModel; FIM_σ²=:observed)
    β = m.β
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov
    σ2sq_gam = vcat([collect(m.sigmas[i]) .^ 2 for i in 1:length(m.sigmas)]...)
    σ2s = [m.sigma^2, σ2sq_gam...]
    Zs_gam = vcat(
        [
            [m.reterms[i][:, j:length(m.sigmas[i]):end] for j in 1:length(m.sigmas[i])] for
            i in 1:length(m.sigmas)
        ]...,
    )
    Zs = [I(nobs(m)), Zs_gam...]
    ZZs = [Z * Z' for Z in Zs]
    V(σ2s) = sum([σ2s[i] * ZZs[i] for i in eachindex(σ2s)])
    function modified_profile_likelihood(σ2s)
        Vinv = inv(V(σ2s))
        return -1 / 2 * logdet(V(σ2s)) - 1 / 2 * logdet(X' * Vinv * X) -
               1 / 2 * (y - X * β)' * Vinv * (y - X * β)
    end
    Vinv = inv(V(σ2s))
    P = [-X' * Vinv * ZZs[i] * Vinv * X for i in eachindex(ZZs)]
    Q = [
        X' * Vinv * ZZs[i] * Vinv * ZZs[j] * Vinv * X for i in eachindex(ZZs),
        j in eachindex(ZZs)
    ]

    if FIM_σ² == :observed
        FIM_σ² = -ForwardDiff.hessian(modified_profile_likelihood, σ2s)
    elseif FIM_σ² == :expected
        FIM_σ² = [
            1 / 2 * tr(Vinv * ZZs[i] * Vinv * ZZs[j]) - tr(Φ * Q[i, j]) +
            1 / 2 * tr(Φ * P[i] * Φ * P[j]) for i in eachindex(σ2s), j in eachindex(σ2s)
        ]
    else
        error("FIM_σ² needs to equal :observed or :expected")
    end
    W = inv(FIM_σ²)

    factor = zeros(size(m.vcov)...)
    for i in eachindex(ZZs)
        for j in eachindex(ZZs)
            factor += W[i, j] * (Q[i, j] - P[i] * m.vcov * P[j])
        end
    end
    varcovar_adjusted = m.vcov + 2 * m.vcov * factor * m.vcov
    adjusted_error = sqrt.([varcovar_adjusted[i, i] for i in 1:size(m.vcov, 1)])
    return KenwardRogerMatrices(
        m, σ2s, V(σ2s), W, P, Q, m.vcov, varcovar_adjusted, adjusted_error
    )
end

function kenwardroger_estimates(m::MixedModel, kr::KenwardRogerMatrices)
    fixed_effects = Vector{FixedEffect}()
    varcovar = kr.CovVar
    for (ic, coef) in enumerate(coefnames(m))
        c = 1
        p = length(coefnames(m))
        C = zeros(p, c)
        C[ic, 1] = 1
        M = C * inv(C' * varcovar * C) * C'
        A1 = 0.0
        A2 = 0.0
        for i in eachindex(kr.σ2s)
            for j in eachindex(kr.σ2s)
                Pi = kr.P[i]
                Pj = kr.P[j]
                Wij = kr.W[i, j]
                A1 +=
                    Wij *
                    tr(M * varcovar * Pi * varcovar) *
                    tr(M * varcovar * Pj * varcovar)
                A2 += Wij * tr(M * varcovar * Pi * varcovar * M * varcovar * Pj * varcovar)
            end
        end

        B = (A1 + 6A2) / (2c)
        g = ((c + 1)A1 - (c + 4)A2) / ((c + 2)A2)
        c1 = g / (3c + 2(1 - g))
        c2 = (c - g) / (3c + 2(1 - g))
        c3 = (c + 2 - g) / (3c + 2(1 - g))
        Estar = inv(1 - A2 / c)
        Vstar = (2 / c) * (1 + c1 * B) / ((1 - c2 * B)^2 * (1 - c3 * B))
        ρ = Vstar / (2 * Estar^2)
        v = 4 + (c + 2) / (c * ρ - 1)
        λ = v / (Estar * (v - 2))

        t = TDist(v)
        tstar = m.beta[ic] / kr.StdBetas[ic]
        p_value = 2 * ccdf(t, abs(tstar))
        α = 0.05
        t_α = quantile(t, 1 - α / 2)
        δ = t_α * kr.StdBetas[ic]

        push!(
            fixed_effects,
            FixedEffect(
                coef,
                m.β[ic],
                kr.StdBetas[ic],
                m.β[ic] - δ,
                m.β[ic] + δ,
                λ,
                v,
                tstar,
                p_value,
            ),
        )
    end
    return fixed_effects
end

end
