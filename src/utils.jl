"""
    validation(m::MixedModel, method)

Validate a fitted model for a given small-sample adjustment method.

- [`KenwardRoger`](@ref) requires a REML fit.
- [`Satterthwaite`](@ref) is supported for both REML and ML fits.

All methods reject weighted fits.
"""
function validation(m::MixedModel, ::KenwardRoger)
    m.optsum.REML || error("Cannot compute REML-based statistics for non-REML fit")
    isempty(m.sqrtwts) || error("Cannot compute for weighted models")
    return nothing
end

function validation(m::MixedModel, ::Satterthwaite)
    isempty(m.sqrtwts) || error("Cannot compute for weighted models")
    return nothing
end
