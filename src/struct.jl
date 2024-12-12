using MixedModels
using ForwardDiff

struct KenwardRogerMatrices
    fit:: MixedModel
    V:: Array{Float64, 2}
    W:: Array{Float64, 2}
    P:: Vector{Matrix{Float64}}
    Q:: Matrix{Matrix{Float64}}
    CovVar:: Array{Float64, 2}
    CovBeta:: Array{Float64, 2}
    StdBetas:: Array{Float64, 1}
end


function KenwardRogerMatrices(m:: MixedModel)
    β = m.β
    y = m.y
    X = m.X
    σ2sq_gam = [m.sigmas[i][1][1]^2 for i in 1:length(m.sigmas)]
    σ2s = [m.sigma^2, σ2sq_gam...]
    Zs = [I(nobs(m)), m.reterms...]
    ZZs = [Z*Z' for Z in Zs]
    V(σ2s) = sum([σ2s[i]*ZZs[i] for i in eachindex(σ2s)])
    function modified_profile_likelihood(σ2s)
        Vinv = inv(V(σ2s))
        -1 / 2 * logdet(V(σ2s)) - 1 / 2 * logdet(X' * Vinv * X) -
        1 / 2 * (y - X * β)' * Vinv * (y - X * β)
    end

    FIM_obs = -ForwardDiff.hessian(modified_profile_likelihood, σ2s)
    W = inv(FIM_obs)

    Vinv = inv(V(σ2s))
    P = [-X' * Vinv * ZZs[i] * Vinv * X for i in eachindex(ZZs)]
    Q = [
        X' * Vinv * ZZs[i] * Vinv * ZZs[j] * Vinv * X for i in eachindex(ZZs), j in eachindex(ZZs)
    ]

    factor = zeros(size(m.vcov)...)
    for i in eachindex(ZZs)
        for j in eachindex(ZZs)
            Pi = P[i]
            Pj = P[j]
            Qij = Q[i, j]
            Wij = W[i, j]
            factor += Wij * (Qij - Pi * m.vcov * Pj)
        end
    end
    varcovar_adjusted = m.vcov + 2 * m.vcov * factor * m.vcov
    adjusted_error = sqrt.([varcovar_adjusted[i, i] for i = 1:size(m.vcov, 1)])
    adjusted_error

    KenwardRogerMatrices(
        m,
        V(σ2s), W, P, Q, m.vcov, varcovar_adjusted, adjusted_error
    )
end
