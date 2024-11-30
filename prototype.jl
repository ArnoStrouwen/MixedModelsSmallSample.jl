using CSV
using DataFrames
using MixedModels
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

σ_eps = m.sigma
σ_gam = m.sigmas[1][1]
Z = m.reterms[1]'
V = diagm(fill(σ_eps^2, size(df,1))) + σ_gam^2*Z'*Z
X = m.X
varcovar = inv(X'*inv(V)*X)
m.vcov
