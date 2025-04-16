function StatsAPI.coeftable(m::LinearMixedModelKR)
    error_adjusted = sqrt.(diag(m.varcovar_adjusted))
    return coeftable(m, error_adjusted)
end

function StatsAPI.coeftable(m::LinearMixedModelSW)
    error = m.m.stderror
    return coeftable(m, error)
end

function StatsAPI.coeftable(m::Union{LinearMixedModelKR,LinearMixedModelSW}, error)
    t = TDist.(m.v)
    tstar = m.m.β ./ error
    p = 2 * ccdf.(t, abs.(tstar))
    α = 0.05
    t_α = quantile.(t, 1 - α / 2)
    δ = t_α .* error
    return CoefTable(
        hcat(m.m.β, error, m.v, tstar, p),
        ["Coef.", "Std. Error", "DenDF", "t", "Pr(>|t|)"],
        coefnames(m.m),
        5, # pvalcol
        4, # teststatcol
    )
end

Base.show(io::IO, m::Union{LinearMixedModelKR,LinearMixedModelSW}) = show(io, coeftable(m))
function create_table(m::Union{LinearMixedModelKR,LinearMixedModelSW})
    ct = coeftable(m)
    first_row = [vcat("", ct.colnms)]
    body_rows = [
        vcat(ct.rownms[i], [string(round(col[i]; sigdigits=5)) for col in ct.cols]) for
        i in eachindex(ct.cols[1])
    ]
    rows = vcat(first_row, body_rows)
    align = [:l, :r, :r, :r, :r, :r, :r, :r, :r]
    return table = Markdown.Table(rows, align)
end
function Base.show(
    io::IO, ::MIME"text/latex", m::Union{LinearMixedModelKR,LinearMixedModelSW}
)
    table = create_table(m)
    return print(io, Markdown.latex(table))
end
function Base.show(
    io::IO, ::MIME"text/html", m::Union{LinearMixedModelKR,LinearMixedModelSW}
)
    table = create_table(m)
    return print(io, Markdown.html(table))
end
function Base.show(
    io::IO, ::MIME"text/markdown", m::Union{LinearMixedModelKR,LinearMixedModelSW}
)
    table = create_table(m)
    return print(io, Markdown.MD(table))
end
