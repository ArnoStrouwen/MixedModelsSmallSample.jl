# Parameterization: UN and FA0 derivative computation

#=============================================================================
# FA0 (Cholesky) Parameterization Helpers
=============================================================================#

"""
    compute_dG_dL(L, i, j)

Derivative of \\(G = LL^\\top\\) with respect to the element \\(L_{ij}\\).

```math
G = LL^\\top,\\qquad
\\frac{\\partial G}{\\partial L_{ij}} = E_{ij}L^\\top + L E_{ij}^\\top,
```

where \\(E_{ij}\\) is the elementary matrix with a 1 in entry \\((i,j)\\).
"""
function compute_dG_dL(L::LowerTriangular, i::Int, j::Int)
    t = size(L, 1)
    dG = zeros(t, t)
    dG[i, :] .+= L[:, j]  # row i gets column j of L (since L' has column j as row j)
    dG[:, i] .+= L[:, j]  # column i gets column j of L (symmetric)
    return dG
end

"""
    compute_d2G_dL2(t, i1, j1, i2, j2)

Second derivative of \\(G = LL^\\top\\) with respect to \\(L_{i_1 j_1}\\) and \\(L_{i_2 j_2}\\).

```math
\\frac{\\partial^2 G}{\\partial L_{i_1 j_1}\\,\\partial L_{i_2 j_2}}
= \\delta_{j_1 j_2}\\,(E_{i_1 i_2} + E_{i_2 i_1}).
```
"""
function compute_d2G_dL2(t::Int, i1::Int, j1::Int, i2::Int, j2::Int)
    d2G = zeros(t, t)
    if j1 == j2
        d2G[i1, i2] += 1
        d2G[i2, i1] += 1  # symmetric
    end
    return d2G
end

"""
    get_cholesky_elements(m::MixedModel)

Extract Cholesky-factor blocks used in the `:FA0` parameterization.

For each random-effects term \\(b\\), MixedModels stores a relative factor \\(\\lambda_b\\); we use

```math
L_b = \\sigma\\,\\lambda_b,\\qquad G_b = L_b L_b^\\top,
```

so that

```math
V = \\sigma^2 I_n + \\sum_b Z_b G_b Z_b^\\top.
```
"""
function get_cholesky_elements(m::MixedModel)
    σ = m.sigma
    result = []
    for (b, rt) in enumerate(m.reterms)
        λ = rt.λ
        L = λ * σ
        t = size(L, 1)
        Zs_block = [rt[:, j:t:end] for j in 1:t]
        push!(result, (L=L, Zs=Zs_block, block=b, t=t))
    end
    return result
end

"""
    _compute_PQR(X, Vinv, dVs, d2Vs)

Compute the derivative interaction matrices P, Q, and R.
"""
function _compute_PQR(X, Vinv, dVs, d2Vs)
    nparams = length(dVs)
    P = [-X' * Vinv * dV * Vinv * X for dV in dVs]
    Q = [X' * Vinv * dVs[i] * Vinv * dVs[j] * Vinv * X for i in 1:nparams, j in 1:nparams]
    R = [X' * Vinv * d2Vs[i, j] * Vinv * X for i in 1:nparams, j in 1:nparams]
    return P, Q, R
end

#=============================================================================
# Variance Decomposition: Unstructured Parameterization
=============================================================================#

"""
    VarianceDecomposition(m::MixedModel, X, ::Unstructured; include_d2V=false)

Compute variance decomposition using the [`Unstructured`](@ref) parameterization.

In this parameterization, the variance is linear in the parameters:

```math
V(\\theta) = \\sum_i \\theta_i D_i,
```

where ``\\theta_1 = \\sigma^2`` (residual variance) and the remaining ``\\theta_i`` are the
(co)variance components. The derivative matrices are simply:

```math
\\frac{\\partial V}{\\partial \\theta_i} = D_i,
```

and all second derivatives vanish: ``\\partial^2 V / \\partial\\theta_i\\partial\\theta_j = 0``.
This makes the ``R_{ij}`` matrices zero, simplifying the Kenward-Roger bias correction.

# Arguments
- `m`: fitted `MixedModel`
- `X`: design matrix
- `::Unstructured`: dispatch tag for this parameterization
- `include_d2V`: ignored (always zero for Unstructured)

# Returns
A [`VarianceDecomposition`](@ref) with `R_factor = 1.0`.

See also: [`VarianceDecomposition(m, X, ::FactorAnalytic)`](@ref)
"""
function VarianceDecomposition(
    m::MixedModel,
    X::Matrix{Float64},
    ::Unstructured;
    include_d2V::Bool=false,
    fa0_tol::Float64=1e-8,  # unused, kept for interface consistency
)
    n = length(m.y)

    # Build θs and ZZs using variance-component parameterization
    nθs = [length(sigmas) for sigmas in m.sigmas]
    k(b, r, c) = (nθs[b] * (c - 1) + r - sum(1:c))

    σ²γ = [
        if r == c
            m.sigmas[b][r]^2
        else
            m.sigmarhos[b][2][k(b, r, c)] * m.sigmas[b][r] * m.sigmas[b][c]
        end for (b, r, c) in m.parmap
    ]
    θs = Float64[m.sigma ^ 2, σ²γ...]
    names = String["sigma2", ["theta_$i" for i in 1:length(σ²γ)]...]

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
    dVs = [Z isa Tuple ? (Z[1] * Z[2]') + (Z[2] * Z[1]') : Z * Z' for Z in Zs]

    # Build V
    V = sum(θs[i] * dVs[i] for i in eachindex(θs))
    Vinv = inv(V)

    # Second derivatives are zero for UN (linear parameterization)
    nparams = length(θs)
    d2Vs = Matrix{Matrix{Float64}}(undef, nparams, nparams)
    for i in 1:nparams, j in 1:nparams
        d2Vs[i, j] = zeros(n, n)
    end

    # Compute P, Q, R matrices
    P, Q, R = _compute_PQR(X, Vinv, dVs, d2Vs)

    dVs = [Matrix{Float64}(dV) for dV in dVs]

    return VarianceDecomposition(
        Matrix(V), Matrix(Vinv), θs, dVs, d2Vs, P, Q, R, 1.0, names
    )
end

#=============================================================================
# Variance Decomposition: FactorAnalytic Parameterization
=============================================================================#

"""
    VarianceDecomposition(m::MixedModel, X, ::FactorAnalytic; include_d2V=false, fa0_tol=1e-8)

Compute variance decomposition using the [`FactorAnalytic`](@ref) (Cholesky) parameterization.

In this parameterization, the variance parameters are the Cholesky factor elements:

```math
\\theta = (\\sigma^2, \\{L_{b,ij}\\}_{b, i \\ge j}),
```

where ``L_b`` is the lower-triangular Cholesky factor of the random-effects covariance
for term ``b``. The total variance is:

```math
V = \\sigma^2 I_n + \\sum_b Z_b L_b L_b^\\top Z_b^\\top.
```

## Gradient Computation

The derivative of ``V`` with respect to ``L_{b,ij}`` uses the chain rule through the
random-effects covariance ``G_b = L_b L_b^\\top``:

```math
\\frac{\\partial V}{\\partial L_{ij}} = Z_b \\frac{\\partial G_b}{\\partial L_{ij}} Z_b^\\top,
\\qquad
\\frac{\\partial G_b}{\\partial L_{ij}} = E_{ij} L_b^\\top + L_b E_{ij}^\\top,
```

where ``E_{ij}`` is the elementary matrix with 1 at position ``(i,j)``.

Unlike [`Unstructured`](@ref), this parameterization has **non-zero second derivatives**,
which affects the ``R_{ij}`` matrices in the Kenward-Roger bias correction.

# Arguments
- `m`: fitted `MixedModel`
- `X`: design matrix
- `::FactorAnalytic`: dispatch tag for this parameterization
- `include_d2V`: whether to compute second derivatives (needed for Kenward-Roger)
- `fa0_tol`: threshold below which Cholesky elements are treated as boundary (skipped)

# Returns
A [`VarianceDecomposition`](@ref) with `R_factor = -0.25`.

See also: [`VarianceDecomposition(m, X, ::Unstructured)`](@ref), [`compute_dG_dL`](@ref)
"""
function VarianceDecomposition(
    m::MixedModel,
    X::Matrix{Float64},
    ::FactorAnalytic;
    include_d2V::Bool=false,
    fa0_tol::Float64=1e-8,
)
    n = length(m.y)
    chol_blocks = get_cholesky_elements(m)

    θs = Float64[m.sigma ^ 2]
    dVs = Matrix{Float64}[Matrix{Float64}(I(n))]
    names = String["sigma2"]
    param_info = [(0, 0, 0)]  # (block, row, col) - dummy for σ²

    # Build V matrix
    V = m.sigma^2 * Matrix{Float64}(I(n))
    for cb in chol_blocks
        L, Zs, t = cb.L, cb.Zs, cb.t
        G = L * L'
        for i in 1:t, j in 1:t
            V .+= G[i, j] .* (Zs[i] * Zs[j]')
        end
    end

    for cb in chol_blocks
        L, Zs, t, block = cb.L, cb.Zs, cb.t, cb.block
        ZsZs_block = [Zs[i] * Zs[j]' for i in 1:t, j in 1:t]

        for j in 1:t, i in j:t  # lower triangular
            abs(L[i, j]) < fa0_tol && continue  # skip boundary parameters

            dG = compute_dG_dL(LowerTriangular(L), i, j)
            dV = zeros(n, n)
            for k in 1:t, l in 1:t
                dG[k, l] != 0 && (dV .+= dG[k, l] .* ZsZs_block[k, l])
            end

            push!(θs, L[i, j])
            push!(dVs, dV)
            push!(names, "L_$(block)_$(i)_$(j)")
            push!(param_info, (block, i, j))
        end
    end

    nparams = length(θs)
    Vinv = inv(V)

    # Compute second derivatives if needed
    d2Vs = Matrix{Matrix{Float64}}(undef, nparams, nparams)
    for i in 1:nparams, j in 1:nparams
        d2Vs[i, j] = zeros(n, n)
    end

    if include_d2V
        for cb in chol_blocks
            t, Zs, block = cb.t, cb.Zs, cb.block
            ZsZs_block = [Zs[i] * Zs[j]' for i in 1:t, j in 1:t]
            block_params = [
                (idx, pi[2], pi[3]) for (idx, pi) in enumerate(param_info) if pi[1] == block
            ]

            for (idx1, i1, j1) in block_params, (idx2, i2, j2) in block_params
                d2G = compute_d2G_dL2(t, i1, j1, i2, j2)
                if any(d2G .!= 0)
                    d2V = zeros(n, n)
                    for k in 1:t, l in 1:t
                        d2G[k, l] != 0 && (d2V .+= d2G[k, l] .* ZsZs_block[k, l])
                    end
                    d2Vs[idx1, idx2] = d2V
                end
            end
        end
    end

    # Compute P, Q, R matrices
    P, Q, R = _compute_PQR(X, Vinv, dVs, d2Vs)

    return VarianceDecomposition(V, Vinv, θs, dVs, d2Vs, P, Q, R, -0.25, names)
end
