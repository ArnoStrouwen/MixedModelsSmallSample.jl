@kwdef struct KenwardRogerMatrices
    fit::MixedModel
    σ2s::Vector{Float64}
    V::Array{Float64,2}
    W::Array{Float64,2}
    P::Vector{Matrix{Float64}}
    Q::Matrix{Matrix{Float64}}
    CovVar::Array{Float64,2}
    CovBeta::Array{Float64,2}
    StdBetas::Array{Float64,1}
end

@kwdef struct CovMatrixVarComponents
    ids::Vector{String}
    values::Vector{String}
end

@kwdef struct FixedEffectsCovariance
    cov_matrix::Vector{Vector{Float64}}
    ids::Vector{String}
end
