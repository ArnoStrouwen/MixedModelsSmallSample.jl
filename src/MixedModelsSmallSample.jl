module MixedModelsSmallSample

using MixedModels
using Distributions
using LinearAlgebra

using StatsAPI: StatsAPI, coeftable
using StatsBase: StatsBase, CoefTable
using Markdown

export adjust_KR
export adjust_SW

export vcov_varpar

export ftest_KR
export ftest_SW

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

"""
    compute_dG_dL(L, i, j)

Compute the derivative of G = L*L' with respect to L[i,j] where L is lower triangular.
Returns dG/dL[i,j] = E_ij * L' + L * E_ji (where E_ij has a 1 at position i,j and 0 elsewhere).
"""
function compute_dG_dL(L::LowerTriangular, i::Int, j::Int)
    t = size(L, 1)
    dG = zeros(t, t)
    # dG/dL[i,j] = e_i * (L e_j)' + (L e_j) * e_i'
    # = e_i * L[j,:]' + L[:,j] * e_i'
    dG[i, :] .+= L[:, j]  # row i gets column j of L
    dG[:, i] .+= L[:, j]  # column i gets column j of L (symmetric)
    return dG
end

"""
    compute_d2G_dL2(t, i1, j1, i2, j2)

Compute the second derivative d²G/d(L[i1,j1])d(L[i2,j2]) where G = L*L'.
Non-zero only when j1 == j2 (both derivatives reference the same column of L).
"""
function compute_d2G_dL2(t::Int, i1::Int, j1::Int, i2::Int, j2::Int)
    d2G = zeros(t, t)
    # d²G / dL[i1,j1] dL[i2,j2] = d/dL[i2,j2] (E_i1,j1 * L' + L * E_j1,i1)
    # = E_i1,j1 * E_j2,i2 + E_i2,j2 * E_j1,i1
    # This is non-zero only when j1 == i2 (first term) or j2 == i1 (second term)
    if j1 == j2
        # Only non-zero when both reference the same column
        d2G[i1, i2] += 1
        d2G[i2, i1] += 1  # symmetric
    end
    return d2G
end

"""
    get_cholesky_elements(m::MixedModel)

Extract the Cholesky factor L = σ * λ for each random effect term.
Returns a vector of (L matrix, Z matrices per column, block index).
"""
function get_cholesky_elements(m::MixedModel)
    n = length(m.y)
    σ = m.sigma
    result = []

    for (b, rt) in enumerate(m.reterms)
        λ = rt.λ
        L = λ * σ  # Actual Cholesky factor (not σ-scaled)
        t = size(L, 1)  # dimension of this random effect block

        # Extract Z columns for this block
        Zs_block = [rt[:, j:t:end] for j in 1:t]

        push!(result, (L=L, Zs=Zs_block, block=b, t=t))
    end

    return result
end

"""
    compute_fa0_derivatives(m::MixedModel)

Compute dV/dθ for FA0 (Cholesky) parameterization.
θ = [σ², L₁₁, L₂₁, L₂₂, ...] (residual variance + lower triangular Cholesky elements)

Returns:
- θs: parameter values
- dVs: derivatives dV/dθᵢ (as matrices)
- d2Vs: second derivatives d²V/dθᵢdθⱼ (as matrix of matrices, only needed for KR)
"""
function compute_fa0_derivatives(m::MixedModel)
    n = length(m.y)
    σ² = m.sigma^2

    # First parameter is always residual variance σ²
    θs = [σ²]
    dVs = [Matrix{Float64}(I(n))]  # dV/dσ² = I

    # Get Cholesky elements for each random effect block
    chol_blocks = get_cholesky_elements(m)

    # Parameter indices for building d2V
    param_info = [(0, 0, 0)]  # (block, row, col) - dummy for σ²

    for cb in chol_blocks
        L = cb.L
        Zs = cb.Zs
        t = cb.t
        block = cb.block

        # Precompute Z*Z' matrices for this block
        ZZs = [Zs[i] * Zs[j]' for i in 1:t, j in 1:t]

        # For each element L[i,j] of the lower triangular matrix (j ≤ i)
        for j in 1:t
            for i in j:t  # lower triangular: i >= j
                # Compute dG/dL[i,j]
                dG = compute_dG_dL(LowerTriangular(L), i, j)

                # dV/dL[i,j] = Z * (dG/dL[i,j]) * Z' = Σₖₗ dG[k,l] * Zₖ * Zₗ'
                dV = zeros(n, n)
                for k in 1:t
                    for l in 1:t
                        if dG[k, l] != 0
                            dV .+= dG[k, l] .* ZZs[k, l]
                        end
                    end
                end

                push!(θs, L[i, j])
                push!(dVs, dV)
                push!(param_info, (block, i, j))
            end
        end
    end

    # Compute second derivatives (for KR adjustment)
    nparams = length(θs)
    d2Vs = Matrix{Matrix{Float64}}(undef, nparams, nparams)

    for i in 1:nparams
        for j in 1:nparams
            d2Vs[i, j] = zeros(n, n)
        end
    end

    # d²V/dσ²dσ² = 0, d²V/dσ²dL... = 0
    # Only non-zero: d²V/dL[i1,j1]dL[i2,j2] when they're in the same block and same column

    param_idx = 2  # Start after σ²
    for cb in chol_blocks
        t = cb.t
        Zs = cb.Zs
        ZZs = [Zs[i] * Zs[j]' for i in 1:t, j in 1:t]

        # Get all L parameters for this block
        block_params = []
        for j in 1:t
            for i in j:t
                push!(block_params, (param_idx, i, j))
                param_idx += 1
            end
        end

        # Compute d²V for pairs of parameters in this block
        for (idx1, i1, j1) in block_params
            for (idx2, i2, j2) in block_params
                d2G = compute_d2G_dL2(t, i1, j1, i2, j2)

                if any(d2G .!= 0)
                    d2V = zeros(n, n)
                    for k in 1:t
                        for l in 1:t
                            if d2G[k, l] != 0
                                d2V .+= d2G[k, l] .* ZZs[k, l]
                            end
                        end
                    end
                    d2Vs[idx1, idx2] = d2V
                end
            end
        end
    end

    return θs, dVs, d2Vs
end

function vcov_varpar(m::MixedModel; FIM_σ²=:observed, parameterization=:UN)
    validation(m)

    β = m.β
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    if parameterization == :FA0
        # FA0 (Cholesky) parameterization - non-linear derivatives
        θs, dVs, d2Vs = compute_fa0_derivatives(m)
        ZZs = dVs  # Use dV/dθ directly as "ZZs" for FIM computation

        # Reconstruct V from the model
        V = zeros(n, n)
        V .+= m.sigma^2 * I(n)
        for cb in get_cholesky_elements(m)
            L = cb.L
            Zs = cb.Zs
            t = cb.t
            G = L * L'
            for i in 1:t, j in 1:t
                V .+= G[i, j] .* (Zs[i] * Zs[j]')
            end
        end
        Vinv = inv(V)

    else  # :UN parameterization (original behavior)
        nθs = [length(sigmas) for sigmas in m.sigmas]
        k(b, r, c) = (nθs[b] * (c - 1) + r - sum(1:c))
        σ²γ = [
            if r == c
                m.sigmas[b][r]^2
            else
                m.sigmarhos[b][2][k(b, r, c)] * m.sigmas[b][r] * m.sigmas[b][c]
            end for (b, r, c) in m.parmap
        ]
        θs = [m.sigma^2, σ²γ...]
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
        V = sum([θs[i] * ZZs[i] for i in eachindex(θs)])
        Vinv = inv(V)
    end

    P = [-transpose(X) * Vinv * ZZ * Vinv * X for ZZ in ZZs]
    Q = [X' * Vinv * ZZi * Vinv * ZZj * Vinv * X for ZZi in ZZs, ZZj in ZZs]

    if FIM_σ² == :observed
        # Profiled observed FIM (matching SAS) - accounts for β̂(θ) dependence
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
            ) for i in eachindex(θs), j in eachindex(θs)
        ]
    elseif FIM_σ² == :expected
        # Profiled expected FIM (matching lmertest) - accounts for β̂(θ) dependence
        FIMσ² = [
            1 / 2 * tr(Vinv * ZZs[i] * Vinv * ZZs[j]) - tr(Φ * Q[i, j]) +
            1 / 2 * tr(Φ * P[i] * Φ * P[j]) for i in eachindex(θs), j in eachindex(θs)
        ]
    else
        error("FIM_σ² needs to equal :observed or :expected")
    end
    return W = inv(FIMσ²)
end

function adjust_KR(m::MixedModel; FIM_σ²=:observed, parameterization=:UN)
    validation(m)

    β = m.β
    p = length(β)
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    if parameterization == :FA0
        # FA0 (Cholesky) parameterization - includes second-order derivatives
        θs, dVs, d2Vs = compute_fa0_derivatives(m)
        ZZs = dVs  # dV/dθ

        # Reconstruct V from the model
        V = zeros(n, n)
        V .+= m.sigma^2 * I(n)
        for cb in get_cholesky_elements(m)
            L = cb.L
            Zs = cb.Zs
            t = cb.t
            G = L * L'
            for i in 1:t, j in 1:t
                V .+= G[i, j] .* (Zs[i] * Zs[j]')
            end
        end
        Vinv = inv(V)

        P = [-transpose(X) * Vinv * ZZ * Vinv * X for ZZ in ZZs]
        Q = [
            X' * Vinv * ZZs[i] * Vinv * ZZs[j] * Vinv * X for
            i in eachindex(θs), j in eachindex(θs)
        ]

        # R term: second-order derivatives (non-zero for FA0!)
        # R[i,j] = X' * Vinv * d²V/dθidθj * Vinv * X
        R = [X' * Vinv * d2Vs[i, j] * Vinv * X for i in eachindex(θs), j in eachindex(θs)]

    else  # :UN parameterization (original behavior)
        nθs = [length(sigmas) for sigmas in m.sigmas]
        k(b, r, c) = (nθs[b] * (c - 1) + r - sum(1:c))
        σ²γ = [
            if r == c
                m.sigmas[b][r]^2
            else
                m.sigmarhos[b][2][k(b, r, c)] * m.sigmas[b][r] * m.sigmas[b][c]
            end for (b, r, c) in m.parmap
        ]
        θs = [m.sigma^2, σ²γ...]
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
        V = sum([θs[i] * ZZs[i] for i in eachindex(θs)])
        Vinv = inv(V)
        P = [-transpose(X) * Vinv * ZZ * Vinv * X for ZZ in ZZs]
        Q = [
            X' * Vinv * ZZs[i] * Vinv * ZZs[j] * Vinv * X for
            i in eachindex(θs), j in eachindex(θs)
        ]

        # R term is zero for UN (second derivatives of linear function are zero)
        R = [zeros(p, p) for i in eachindex(θs), j in eachindex(θs)]
    end

    W = vcov_varpar(m; FIM_σ²=FIM_σ², parameterization=parameterization)

    # Include R term in the adjustment factor
    # For FA0: the second-order term uses factor -1/4 = (-1/2) * (1/2):
    #   - The -1/2 comes from the bias correction direction (shrinks variance)
    #   - The 1/2 comes from the second-order Taylor expansion
    # For UN: R is zero, so the factor doesn't matter
    R_factor = (parameterization == :FA0) ? -0.25 : 1.0

    factor = zeros(size(m.vcov)...)
    for i in eachindex(ZZs)
        for j in eachindex(ZZs)
            # The full KR adjustment: Q[i,j] - P[i]*Φ*P[j] + R_factor*R[i,j]
            factor += W[i, j] * (Q[i, j] - P[i] * Φ * P[j] + R_factor * R[i, j])
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
        for i in eachindex(θs)
            for j in eachindex(θs)
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

function adjust_SW(m::MixedModel; FIM_σ²=:observed, parameterization=:UN)
    validation(m)

    β = m.β
    p = length(β)
    y = m.y
    X = m.X
    n = length(y)
    Φ = m.vcov

    if parameterization == :FA0
        # FA0 (Cholesky) parameterization
        θs, dVs, _ = compute_fa0_derivatives(m)
        ZZs = dVs  # dV/dθ

        # Reconstruct V from the model
        V = zeros(n, n)
        V .+= m.sigma^2 * I(n)
        for cb in get_cholesky_elements(m)
            L = cb.L
            Zs = cb.Zs
            t = cb.t
            G = L * L'
            for i in 1:t, j in 1:t
                V .+= G[i, j] .* (Zs[i] * Zs[j]')
            end
        end
        Vinv = inv(V)

    else  # :UN parameterization (original behavior)
        nθs = [length(sigmas) for sigmas in m.sigmas]
        k(b, r, c) = (nθs[b] * (c - 1) + r - sum(1:c))
        σ²γ = [
            if r == c
                m.sigmas[b][r]^2
            else
                m.sigmarhos[b][2][k(b, r, c)] * m.sigmas[b][r] * m.sigmas[b][c]
            end for (b, r, c) in m.parmap
        ]
        θs = [m.sigma^2, σ²γ...]
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
        V = sum([θs[i] * ZZs[i] for i in eachindex(θs)])
        Vinv = inv(V)
    end

    W = vcov_varpar(m; FIM_σ²=FIM_σ², parameterization=parameterization)
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
