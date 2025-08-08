module MixedModelsSmallSample

using MixedModels
using Distributions
using LinearAlgebra

using StatsAPI: StatsAPI, coeftable
using StatsBase: StatsBase, CoefTable
using Markdown

export adjust_KR
export adjust_SW
export adjust

export KenwardRoger
export Satterthwaite
export AbstractLMMSS

export vcov_varpar

export ftest_KR
export ftest_SW

# Abstract type for adjustment methods
abstract type AbstractLMMSS end

# Concrete method types
struct KenwardRoger <: AbstractLMMSS end
struct Satterthwaite <: AbstractLMMSS end

struct LinearMixedModelKR{Float64} <: MixedModel{Float64}
    m::LinearMixedModel{Float64}
    varcovar_adjusted::Matrix{Float64}
    W::Matrix{Float64}
    P::Vector{Matrix{Float64}}
    Q::Matrix{Matrix{Float64}}
    v::Vector{Float64}
end
struct LinearMixedModelSW{Float64} <: MixedModel{Float64}
    m::LinearMixedModel{Float64}
    W::Matrix{Float64}
    v::Vector{Float64}
end

function validation(m)
    @assert m.optsum.REML "Restricted maximum likelihood must be used, instead of maximum likelihood."
    return nothing
end

function vcov_varpar(m::MixedModel; FIM_σ²=:observed)
    validation(m)

    β = m.β
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    nθs = [length(sigmas) for sigmas in m.sigmas]
    k(b, r, c) = (nθs[b] * (c - 1) + r - sum(1:c))
    σ²γ = [
        if r == c
            m.sigmas[b][r]^2
        else
            m.sigmarhos[b][2][k(b, r, c)] * m.sigmas[b][r] * m.sigmas[b][c]
        end for (b, r, c) in m.parmap
    ]
    σ²s = [m.sigma^2, σ²γ...]
    Zsγ = [
        if r == c
            m.reterms[b][:, r:length(m.sigmas[b]):end]
        else
            (
                m.reterms[b][:, c:length(m.sigmas[b]):end],
                m.reterms[b][:, r:length(m.sigmas[b]):end],
            )
        end for (b, r, c) in m.parmap
    ]
    Zs = [I(n), Zsγ...]
    ZZs = [Z isa Tuple ? (Z[1] * Z[2]') + (Z[2] * Z[1]') : Z * Z' for Z in Zs]
    V = sum([σ²s[i] * ZZs[i] for i in eachindex(σ²s)])
    Vinv = inv(V)
    P = [-transpose(X) * Vinv * ZZ * Vinv * X for ZZ in ZZs]
    Q = [X' * Vinv * ZZi * Vinv * ZZj * Vinv * X for ZZi in ZZs, ZZj in ZZs]
    if FIM_σ² == :observed
        Pvcov = Vinv - Vinv * X * Φ * X' * Vinv
        FIMσ² = [
            (
                1 / 2 * tr(-Pvcov * ZZs[i] * Pvcov * ZZs[j]) -
                1 / 2 *
                (y - X * β)' *
                Vinv *
                (-2 * ZZs[i] * Vinv * ZZs[j]) *
                Vinv *
                (y - X * β)
            ) for i in eachindex(σ²s), j in eachindex(σ²s)
        ]
    elseif FIM_σ² == :observed_SAS_MATCHING
        Pvcov = Vinv - Vinv * X * Φ * X' * Vinv
        FIMσ² = [
            (
                1 / 2 * tr(-Pvcov * ZZs[i] * Pvcov * ZZs[j]) -
                1 / 2 *
                (y - X * β)' *
                Vinv *
                (-2 * ZZs[i] * Pvcov * ZZs[j]) *
                Vinv *
                (y - X * β)
            ) for i in eachindex(σ²s), j in eachindex(σ²s)
        ]
    elseif FIM_σ² == :expected
        FIMσ² = [
            1 / 2 * tr(Vinv * ZZs[i] * Vinv * ZZs[j]) - tr(Φ * Q[i, j]) +
            1 / 2 * tr(Φ * P[i] * Φ * P[j]) for i in eachindex(σ²s), j in eachindex(σ²s)
        ]
    else
        error("FIM_σ² needs to equal :observed or :expected")
    end
    return W = inv(FIMσ²)
end

function adjust_KR(m::MixedModel; FIM_σ²=:observed)
    validation(m)

    β = m.β
    p = length(β)
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    nθs = [length(sigmas) for sigmas in m.sigmas]
    k(b, r, c) = (nθs[b] * (c - 1) + r - sum(1:c))
    σ²γ = [
        if r == c
            m.sigmas[b][r]^2
        else
            m.sigmarhos[b][2][k(b, r, c)] * m.sigmas[b][r] * m.sigmas[b][c]
        end for (b, r, c) in m.parmap
    ]
    σ²s = [m.sigma^2, σ²γ...]
    Zsγ = [
        if r == c
            m.reterms[b][:, r:length(m.sigmas[b]):end]
        else
            (
                m.reterms[b][:, c:length(m.sigmas[b]):end],
                m.reterms[b][:, r:length(m.sigmas[b]):end],
            )
        end for (b, r, c) in m.parmap
    ]
    Zs = [I(n), Zsγ...]
    ZZs = [Z isa Tuple ? (Z[1] * Z[2]') + (Z[2] * Z[1]') : Z * Z' for Z in Zs]
    V = sum([σ²s[i] * ZZs[i] for i in eachindex(σ²s)])
    Vinv = inv(V)
    P = [-transpose(X) * Vinv * ZZ * Vinv * X for ZZ in ZZs]
    Q = [X' * Vinv * ZZi * Vinv * ZZj * Vinv * X for ZZi in ZZs, ZZj in ZZs]

    W = vcov_varpar(m; FIM_σ²=FIM_σ²)

    factor = zeros(size(m.vcov)...)
    for i in eachindex(ZZs)
        for j in eachindex(ZZs)
            factor += W[i, j] * (Q[i, j] - P[i] * Φ * P[j])
        end
    end
    varcovar_adjusted = Φ + 2 * Φ * factor * Φ
    error_adjusted = sqrt.(diag(varcovar_adjusted))

    v = zeros(p)
    for k in eachindex(β)
        c = 1
        C = zeros(p, c)
        C[k, 1] = 1
        M = C * inv(C' * Φ * C) * C'
        A1 = 0.0
        A2 = 0.0
        for i in eachindex(σ²s)
            for j in eachindex(σ²s)
                A1 += W[i, j] * tr(M * Φ * P[i] * Φ) * tr(M * Φ * P[j] * Φ)
                A2 += W[i, j] * tr(M * Φ * P[i] * Φ * M * Φ * P[j] * Φ)
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
        v[k] = 4 + (c + 2) / (c * ρ - 1)
        #λ = v / (Estar * (v - 2))
    end
    return LinearMixedModelKR(m, varcovar_adjusted, W, P, Q, v)
end

function adjust_SW(m::MixedModel; FIM_σ²=:observed)
    validation(m)

    β = m.β
    p = length(β)
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    nθs = [length(sigmas) for sigmas in m.sigmas]
    k(b, r, c) = (nθs[b] * (c - 1) + r - sum(1:c))
    σ²γ = [
        if r == c
            m.sigmas[b][r]^2
        else
            m.sigmarhos[b][2][k(b, r, c)] * m.sigmas[b][r] * m.sigmas[b][c]
        end for (b, r, c) in m.parmap
    ]
    σ²s = [m.sigma^2, σ²γ...]
    Zsγ = [
        if r == c
            m.reterms[b][:, r:length(m.sigmas[b]):end]
        else
            (
                m.reterms[b][:, c:length(m.sigmas[b]):end],
                m.reterms[b][:, r:length(m.sigmas[b]):end],
            )
        end for (b, r, c) in m.parmap
    ]
    Zs = [I(n), Zsγ...]
    ZZs = [Z isa Tuple ? (Z[1] * Z[2]') + (Z[2] * Z[1]') : Z * Z' for Z in Zs]
    V = sum([σ²s[i] * ZZs[i] for i in eachindex(σ²s)])
    Vinv = inv(V)

    W = vcov_varpar(m; FIM_σ²=FIM_σ²)
    v = zeros(p)
    for k in eachindex(β)
        c = 1
        C = zeros(p, c)
        C[k, 1] = 1
        grad = [first(C' * Φ * X' * Vinv * ZZ * Vinv * X * Φ * C) for ZZ in ZZs]
        v[k] = 2 * (first(C' * inv(X' * inv(V) * X) * C))^2 / (grad' * W * grad)
    end
    return LinearMixedModelSW(m, W, v)
end

# Unified adjust function with method dispatch
function adjust(m::MixedModel; method::AbstractLMMSS=KenwardRoger(), FIM_σ²=:observed)
    if method isa KenwardRoger
        return adjust_KR(m; FIM_σ²=FIM_σ²)
    elseif method isa Satterthwaite
        return adjust_SW(m; FIM_σ²=FIM_σ²)
    else
        error("Unsupported method: $(typeof(method))")
    end
end

function ftest_SW(m::LinearMixedModel, L; FIM_σ²=:observed)
    validation(m)

    q = size(L, 2)
    β = m.β
    p = length(β)
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    covLβ = L' * Φ * L
    M = L * inv(L' * Φ * L) * L'
    F = eigen(Hermitian(covLβ))
    d = F.values
    P = F.vectors
    L̃ = L * P

    σ²γ = vcat([collect(sigmas) .^ 2 for sigmas in m.sigmas]...)
    σ²s = [m.sigma^2, σ²γ...]
    Zsγ = vcat(
        [
            [m.reterms[i][:, j:length(m.sigmas[i]):end] for j in 1:length(m.sigmas[i])] for
            i in 1:length(m.sigmas)
        ]...,
    )
    Zs = [I(n), Zsγ...]
    ZZs = [Z * Z' for Z in Zs]
    V = sum([σ²s[i] * ZZs[i] for i in eachindex(σ²s)])
    Vinv = inv(V)

    W = vcov_varpar(m; FIM_σ²=FIM_σ²)
    vs = zeros(q)
    for i in eachindex(vs)
        grad = [
            first(L̃[:, i]' * Φ * X' * Vinv * ZZ * Vinv * X * Φ * L̃[:, i]) for ZZ in ZZs
        ]
        vs[i] =
            2 * (first(L̃[:, i]' * inv(X' * inv(V) * X) * L̃[:, i]))^2 / (grad' * W * grad)
    end

    EQ = sum(νᵢ / (νᵢ - 2) for νᵢ in vs)
    v = 2 * EQ / (EQ - q)

    Fstar = (1 / q * β' * M * β)
    return (v, Fstar)
end
function ftest_KR(m::LinearMixedModel, L; FIM_σ²=:observed)
    validation(m)

    c = q = size(L, 2)
    β = m.β
    p = length(β)
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    σ²γ = vcat([collect(sigmas) .^ 2 for sigmas in m.sigmas]...)
    σ²s = [m.sigma^2, σ²γ...]
    Zsγ = vcat(
        [
            [m.reterms[i][:, j:length(m.sigmas[i]):end] for j in 1:length(m.sigmas[i])] for
            i in 1:length(m.sigmas)
        ]...,
    )
    Zs = [I(n), Zsγ...]
    ZZs = [Z * Z' for Z in Zs]
    V = sum([σ²s[i] * ZZs[i] for i in eachindex(σ²s)])
    Vinv = inv(V)
    P = [-transpose(X) * Vinv * ZZ * Vinv * X for ZZ in ZZs]
    Q = [X' * Vinv * ZZi * Vinv * ZZj * Vinv * X for ZZi in ZZs, ZZj in ZZs]

    W = vcov_varpar(m; FIM_σ²=FIM_σ²)

    factor = zeros(size(m.vcov)...)
    for i in eachindex(ZZs)
        for j in eachindex(ZZs)
            factor += W[i, j] * (Q[i, j] - P[i] * Φ * P[j])
        end
    end
    varcovar_adjusted = Φ + 2 * Φ * factor * Φ
    error_adjusted = sqrt.(diag(varcovar_adjusted))

    M = L * inv(L' * Φ * L) * L'
    A1 = 0.0
    A2 = 0.0
    for i in eachindex(σ²s)
        for j in eachindex(σ²s)
            A1 += W[i, j] * tr(M * Φ * P[i] * Φ) * tr(M * Φ * P[j] * Φ)
            A2 += W[i, j] * tr(M * Φ * P[i] * Φ * M * Φ * P[j] * Φ)
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
    covLβ = L' * Φ * L
    Fstar = λ * (1 / q * β' * M * β)
    return (v, Fstar)
end
include("show.jl")
end
