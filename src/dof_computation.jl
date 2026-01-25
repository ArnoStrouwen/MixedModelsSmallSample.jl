"""
    compute_dof(m::LinearMixedModelAdjusted, L) -> (ν, λ)

Compute degrees of freedom and F-statistic scaling factor for contrast `L`.

Dispatches to the specific method based on `m.adj`:
- For [`KenwardRoger`](@ref): computes scaled F-statistic with ``λ ≠ 1``
- For [`Satterthwaite`](@ref): returns ``λ = 1`` (no scaling)

# Arguments
- `m`: A [`LinearMixedModelAdjusted`](@ref) from [`small_sample_adjust`](@ref)
- `L`: Contrast matrix (``p × q``) where ``p`` = number of fixed effects

# Returns
- `ν`: Denominator degrees of freedom
- `λ`: F-statistic scaling factor (1.0 for Satterthwaite)
"""
function compute_dof(m::LinearMixedModelAdjusted, L::AbstractVecOrMat)
    return compute_dof(m.adj, m, L)
end

"""
    compute_dof(::KenwardRoger, m, L)

Kenward-Roger denominator degrees of freedom and F-scaling factor.

For testing ``H_0: L^\\top\\beta = 0`` with ``q`` constraints, the standard Wald
F-statistic ``F = \\beta^\\top M \\beta / q`` (where ``M = L(L^\\top\\Phi L)^{-1}L^\\top``)
is scaled to ``F^* = \\lambda F``. The pair ``(\\nu, \\lambda)`` is chosen so that
``F^*`` approximately follows an ``F(q, \\nu)`` distribution.

## Algorithm

1. Compute moment quantities:
   ```math
   A_1 = \\sum_{i,j} W_{ij} \\operatorname{tr}(M\\Phi P_i \\Phi) \\operatorname{tr}(M\\Phi P_j \\Phi)
   ```
   ```math
   A_2 = \\sum_{i,j} W_{ij} \\operatorname{tr}(M\\Phi P_i \\Phi\\, M\\Phi P_j \\Phi)
   ```

2. Compute intermediate terms:
   ```math
   B = \\frac{A_1 + 6A_2}{2q}, \\quad
   g = \\frac{(q+1)A_1 - (q+4)A_2}{(q+2)A_2}
   ```
   ```math
   c_1 = \\frac{g}{3q + 2(1-g)}, \\quad
   c_2 = \\frac{q-g}{3q + 2(1-g)}, \\quad
   c_3 = \\frac{q+2-g}{3q + 2(1-g)}
   ```

3. Compute expectation and variance of the scaled statistic:
   ```math
   E^* = \\frac{1}{1 - A_2/q}, \\quad
   V^* = \\frac{2}{q} \\cdot \\frac{1 + c_1 B}{(1-c_2 B)^2 (1-c_3 B)}
   ```

4. Match to F-distribution moments:
   ```math
   \\rho = \\frac{V^*}{2(E^*)^2}, \\quad
   \\nu = 4 + \\frac{q+2}{q\\rho - 1}, \\quad
   \\lambda = \\frac{\\nu}{E^*(\\nu - 2)}
   ```
"""
function compute_dof(::KenwardRoger, m::LinearMixedModelAdjusted, L::AbstractVecOrMat)
    Φ = m.m.vcov
    q = L isa AbstractVector ? 1 : size(L, 2)
    M = L * inv(L' * Φ * L) * L'

    # KR algorithm
    A1 = 0.0
    A2 = 0.0
    nθ = size(m.W, 1)

    for i in 1:nθ, j in 1:nθ
        t1 = tr(M * Φ * m.P[i] * Φ)
        t2 = tr(M * Φ * m.P[j] * Φ)
        A1 += m.W[i, j] * t1 * t2
        A2 += m.W[i, j] * tr(M * Φ * m.P[i] * Φ * M * Φ * m.P[j] * Φ)
    end

    B = (A1 + 6A2) / (2q)
    g = ((q + 1) * A1 - (q + 4) * A2) / ((q + 2) * A2)
    c1 = g / (3q + 2(1 - g))
    c2 = (q - g) / (3q + 2(1 - g))
    c3 = (q + 2 - g) / (3q + 2(1 - g))
    Estar = inv(1 - A2 / q)
    Vstar = (2 / q) * (1 + c1 * B) / ((1 - c2 * B)^2 * (1 - c3 * B))
    ρ = Vstar / (2 * Estar^2)
    ν = 4 + (q + 2) / (q * ρ - 1)
    λ = ν / (Estar * (ν - 2))

    return ν, λ
end

"""
    compute_dof(::Satterthwaite, m, L)

Satterthwaite denominator degrees of freedom for contrast `L`.

For a scalar contrast ``c^\\top\\beta``, the variance is ``v = c^\\top\\Phi\\,c``. The Satterthwaite
approximation matches the first two moments of ``\\hat{v}`` to a scaled chi-squared
distribution, yielding:

```math
\\nu = \\frac{2v^2}{(\\nabla_\\theta v)^\\top W (\\nabla_\\theta v)},
```

where the gradient is ``\\nabla_\\theta v = [-c^\\top\\Phi\\,P_i\\,\\Phi\\,c]_i`` and
``W`` is the asymptotic covariance of ``\\hat{\\theta}``.

## Multivariate Contrasts

For a matrix ``L`` (multiple constraints), the algorithm:

1. Computes the eigendecomposition of ``L^\\top\\Phi\\,L``
2. Transforms to orthogonal contrasts ``\\tilde{L} = L \\cdot \\text{eigenvectors}``
3. Computes ``\\nu_i`` for each orthogonal contrast
4. Combines using: ``\\nu = 2 E_Q / (E_Q - q)`` where ``E_Q = \\sum_i \\nu_i/(\\nu_i - 2)``

# Reference
Satterthwaite, F. E. (1946). "An Approximate Distribution of Estimates of
Variance Components." *Biometrics Bulletin*, 2(6), 110-114.
"""
function compute_dof(::Satterthwaite, m::LinearMixedModelAdjusted, L::AbstractVecOrMat)
    Φ = m.m.vcov

    if L isa AbstractVector
        Ltilde = reshape(L, :, 1)
        q = 1
    else
        covLβ = L' * Φ * L
        F = eigen(Hermitian(covLβ))
        Ltilde = L * F.vectors
        q = size(L, 2)
    end

    vs = zeros(q)
    for i in 1:q
        Li = Ltilde[:, i]
        variance = first(Li' * Φ * Li)
        grad = [-first(Li' * Φ * P_k * Φ * Li) for P_k in m.P]
        vs[i] = 2 * variance^2 / (grad' * m.W * grad)
    end

    if q == 1
        return vs[1], 1.0
    end

    EQ = sum(νᵢ / (νᵢ - 2) for νᵢ in vs)
    ν = 2 * EQ / (EQ - q)

    return ν, 1.0
end

# Internal helper: computes DOF for each coefficient by calling compute_dof with e_k vectors
function compute_per_coefficient_dof(m::LinearMixedModelAdjusted)
    p = length(m.m.β)
    ν = zeros(p)
    for k in 1:p
        L = zeros(p, 1)
        L[k, 1] = 1.0
        ν[k] = first(compute_dof(m, L))
    end
    return ν
end
