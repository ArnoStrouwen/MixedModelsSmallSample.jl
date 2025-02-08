module KenwardRoger

using ForwardDiff
using MixedModels
using Distributions
using LinearAlgebra

using StatsAPI: StatsAPI, coeftable
using StatsBase: StatsBase, CoefTable

include("struct.jl")

export kenwardroger_matrices #, kenwardroger_estimates

function kenwardroger_matrices(m::MixedModel, FIM_σ=:observed)
    β = m.β
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov
    σ2sq_gam = [m.sigmas[i][1][1]^2 for i in 1:length(m.sigmas)]
    σ2s = [m.sigma^2, σ2sq_gam...]
    Zs = [I(nobs(m)), m.reterms...]
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

    if FIM_σ == :observed
        FIM = -ForwardDiff.hessian(modified_profile_likelihood, σ2s)
    elseif FIM_σ == :expected
        FIM = [
            1 / 2 * tr(Vinv * ZZs[i] * Vinv * ZZs[j]) - tr(Φ * Q[i, j]) +
            1 / 2 * tr(Φ * P[i] * Φ * P[j]) for i in eachindex(σ2s), j in eachindex(σ2s)
        ]
    else
        error("FIM_σ needs to equal :observed or :expected")
    end
    W = inv(FIM)

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


function StatsAPI.coeftable(m::MixedModel, kr::KenwardRogerMatrices)
    varcovar = kr.CovVar
    β = Vector{Float64}()
    σ = Vector{Float64}()
    lb_ci_alpha05 = Vector{Float64}()
    ub_ci_alpha05 = Vector{Float64}()
    num_df = Vector{Float64}()
    den_df = Vector{Float64}()
    t_statistic = Vector{Float64}()
    p_values = Vector{Float64}()

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

        push!(β, m.β[ic])
        push!(σ, kr.StdBetas[ic])
        push!(lb_ci_alpha05, m.β[ic] - δ)
        push!(ub_ci_alpha05, m.β[ic] + δ)
        push!(num_df, λ)
        push!(den_df, v)
        push!(t_statistic, tstar)
        push!(p_values, p_value)

    end
    CoefTable(
        hcat(β, σ, lb_ci_alpha05, ub_ci_alpha05, num_df, den_df, t_statistic, p_values),
        ["Coef.", "Std. Error", "LB coef.", "UB coef.", "NumDF", "Denominator degrees of freedom", "t", "Pr(>|t|)"],
        coefnames(m),
        8, # pvalcol
        7, # teststatcol
    )
end

end
