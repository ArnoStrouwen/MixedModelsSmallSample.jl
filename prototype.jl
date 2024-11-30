using CSV
using DataFrames
using MixedModels
using LinearAlgebra
using ForwardDiff
df = DataFrame(CSV.File("Data Pastry Dough Experiment Chapter 7.csv"))
rename!(df, "Flow Rate" => :FR, "Moisture Content" => :MC, "Screw Speed" => :SS)
rename!(df, "Longitudinal Expansion Index" => :LEI)
function one_minus_one_coding!(x)
    minx = minimum(x)
    maxx = maximum(x)
    span = maxx-minx
    x .-= minx .+span/2
    x ./= (span/2)
    return nothing
end
one_minus_one_coding!(df[!,:FR])
one_minus_one_coding!(df[!,:MC])
one_minus_one_coding!(df[!,:SS])

fm = @formula(LEI~ 1 + (1|Day) + FR + MC + SS + FR&MC + FR&SS + MC&SS + FR&FR + MC&MC + SS&SS)
m = fit(MixedModel, fm, df; REML = true)

σsq_eps = m.sigma^2
σsq_gam = m.sigmas[1][1]^2
Z = m.reterms[1]'
V(σsq_eps,σsq_gam) = diagm(fill(σsq_eps, size(df,1))) + σsq_gam*Z'*Z
X = m.X
varcovar = inv(X'*inv(V(σsq_eps,σsq_gam))*X)
m.vcov

dVinv_dσ = [ForwardDiff.derivative(σsq_eps->inv(V(σsq_eps,σsq_gam)),σsq_eps),
            ForwardDiff.derivative(σsq_gam->inv(V(σsq_eps,σsq_gam)),σsq_gam)]

W = [0.3637156079 -0.027399087; -0.027399087 0.1169879811] # wald style var-co-var for estimate of variance components. Not yet reproduced in MixedModels.jl

factor = zeros(size(varcovar)...)
for i in eachindex(dVinv_dσ)
    for j in eachindex(dVinv_dσ)
        Pi = X'*dVinv_dσ[i]*X
        Pj = X'*dVinv_dσ[j]*X
        Qij = X'*dVinv_dσ[i]*V(σsq_eps,σsq_gam)*dVinv_dσ[j]*X
        Wij = W[i,j]
        factor += Wij*(Qij - Pi*m.vcov*Pj)
    end
end

varcovar_adjusted = m.vcov + 2*m.vcov*factor*m.vcov
adjusted_error = sqrt.([varcovar_adjusted[i,i] for i in 1:size(m.vcov,1)])
