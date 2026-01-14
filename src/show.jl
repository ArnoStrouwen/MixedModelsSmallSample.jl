function StatsAPI.coeftable(m::LinearMixedModelAdjusted)
    if m.adj isa KenwardRoger
        error = sqrt.(diag(m.varcovar_adjusted))
    else # Satterthwaite
        error = m.m.stderror
    end
    return coeftable(m, error)
end

function StatsAPI.coeftable(m::LinearMixedModelAdjusted, error)
    t = TDist.(m.ν)

    # Avoid division by zero
    tstar = m.m.β ./ max.(error, eps())

    p = 2 * ccdf.(t, abs.(tstar))
    α = 0.05
    t_α = quantile.(t, 1 - α / 2)
    # δ = t_α .* error # This was unused in original code

    return CoefTable(
        hcat(m.m.β, error, m.ν, tstar, p),
        ["Coef.", "Std. Error", "DenDF", "t", "Pr(>|t|)"],
        coefnames(m.m),
        5, # pvalcol
        4, # teststatcol
    )
end

Base.show(io::IO, m::LinearMixedModelAdjusted) = show(io, coeftable(m))

function create_table(m::LinearMixedModelAdjusted)
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

function Base.show(io::IO, ::MIME"text/latex", m::LinearMixedModelAdjusted)
    table = create_table(m)
    return print(io, Markdown.latex(table))
end

function Base.show(io::IO, ::MIME"text/html", m::LinearMixedModelAdjusted)
    table = create_table(m)
    return print(io, Markdown.html(table))
end

function Base.show(io::IO, ::MIME"text/markdown", m::LinearMixedModelAdjusted)
    table = create_table(m)
    return print(io, Markdown.MD(table))
end
