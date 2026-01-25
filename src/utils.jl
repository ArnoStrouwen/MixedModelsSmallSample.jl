"""
    validation(m::MixedModel)

Checks if the model `m` is fitted with REML and is not weighted.
Throws an error if these conditions are not met.
"""
function validation(m::MixedModel)
    m.optsum.REML || error("Cannot compute REML-based statistics for non-REML fit")
    isempty(m.sqrtwts) || error("Cannot compute for weighted models")
    return nothing
end
