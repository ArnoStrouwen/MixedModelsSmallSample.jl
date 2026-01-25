"""
    small_sample_ftest(m::LinearMixedModelAdjusted, L) -> (ν, F*)

Compute an adjusted F-test for the null hypothesis ``L^\\top\\beta = 0``.

Returns ``(ν, F^*)`` where ``(ν, λ)`` is obtained from [`compute_dof`](@ref) and
``F^* = λF`` with ``F = β^\\top Mβ/q`` and ``M = L(L^\\top\\Phi L)^{-1}L^\\top``
(using the unadjusted covariance ``\\Phi``).

# Arguments
- `m`: An adjusted model ([`LinearMixedModelAdjusted`](@ref))
- `L`: Contrast matrix (``p \\times q``) where ``p`` = number of fixed effects

# Returns
- `ν`: Denominator degrees of freedom
- `F*`: Scaled F-statistic

Under ``H_0``, ``F^* \\sim F(q, \\nu)`` approximately.

# Example

```julia
kr = small_sample_adjust(model, KenwardRoger())
L = [0 1 0; 0 0 1]'  # Test β₂ = β₃ = 0
ν, Fstar = small_sample_ftest(kr, L)
pvalue = 1 - cdf(FDist(size(L,2), ν), Fstar)
```
"""
function small_sample_ftest(m::LinearMixedModelAdjusted, L)
    q = size(L, 2)
    β = m.m.β
    Φ = m.m.vcov
    M = L * inv(L' * Φ * L) * L'
    ν, λ = compute_dof(m, L)
    Fstar = λ * (β' * M * β / q)
    return (ν, Fstar)
end

"""
    small_sample_ftest(m::MixedModel, L; method=KenwardRoger()) -> (ν, F*)

Convenience wrapper that first applies [`small_sample_adjust`](@ref) then
computes the F-test.

Equivalent to:
```julia
m_adj = small_sample_adjust(m, method)
small_sample_ftest(m_adj, L)
```

# Arguments
 - `m`: A fitted `MixedModel` (must use REML for KR; ML or REML for SW)
- `L`: Contrast matrix
- `method`: [`KenwardRoger`](@ref) or [`Satterthwaite`](@ref)
"""
function small_sample_ftest(
    m::MixedModel, L; method::AbstractLinearMixedModelAdjustment=KenwardRoger()
)
    m_adj = small_sample_adjust(m, method)
    return small_sample_ftest(m_adj, L)
end
