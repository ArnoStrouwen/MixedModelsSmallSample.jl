using CSV
using DataFrames
using Statistics

function _lowercase_names!(df::DataFrame)
    rename!(df, Dict(n => Symbol(lowercase(String(n))) for n in names(df)))
    return df
end

function _covp_cols(df::DataFrame)
    return [n for n in names(df) if occursin(r"^covp[0-9]+$", lowercase(String(n)))]
end

function _detect_drop_from_asycov(asycov::DataFrame; tol::Float64=0.0)
    asycov = copy(asycov)
    _lowercase_names!(asycov)

    # DataFrames.jl returns names(::DataFrame) as Strings; use propertynames() for Symbols
    if !(:sim_id in propertynames(asycov))
        error("Expected BY variable `sim_id` to appear in AsyCov output")
    end
    if !(:covparm in propertynames(asycov))
        error("Expected column `CovParm` to appear in AsyCov output")
    end

    covp = _covp_cols(asycov)
    isempty(covp) && error("No CovP columns found in AsyCov output")

    g = groupby(asycov, :sim_id)
    out = DataFrame(; sim_id=Int[], dropped_fa22=Bool[])
    for sdf in g
        fa_row = sdf[lowercase.(String.(sdf.covparm)) .== "fa(2,2)", :]
        if nrow(fa_row) == 0
            push!(out, (sim_id=sdf.sim_id[1], dropped_fa22=true))
            continue
        end

        vals = Vector{Float64}(undef, length(covp))
        for (i, c) in enumerate(covp)
            v = fa_row[1, c]
            vals[i] = v isa Missing ? 0.0 : Float64(v)
        end
        dropped = all(abs.(vals) .<= tol)
        push!(out, (sim_id=sdf.sim_id[1], dropped_fa22=dropped))
    end
    return out
end

function analyze_fa0_drop_results(;
    sim_data_csv::AbstractString=joinpath(
        @__DIR__, "..", "data", "Data fa0 drop simulation.csv"
    ),
    mixed_asycov_csv::AbstractString=joinpath(
        @__DIR__, "..", "results", "Results fa0 drop simulation mixed asycov.csv"
    ),
    glimmix_asycov_csv::AbstractString=joinpath(
        @__DIR__, "..", "results", "Results fa0 drop simulation glimmix asycov.csv"
    ),
    tol::Float64=0.0,
)
    sim = DataFrame(CSV.File(sim_data_csv))
    _lowercase_names!(sim)
    rho_map = combine(groupby(sim, :sim_id), :rho => first => :rho)
    sort!(rho_map, :sim_id)

    mixed_asy = DataFrame(CSV.File(mixed_asycov_csv))
    glimmix_asy = DataFrame(CSV.File(glimmix_asycov_csv))

    mixed_drop = _detect_drop_from_asycov(mixed_asy; tol)
    glimmix_drop = _detect_drop_from_asycov(glimmix_asy; tol)

    mixed = leftjoin(rho_map, mixed_drop; on=:sim_id)
    glimmix = leftjoin(rho_map, glimmix_drop; on=:sim_id)

    function summarize(df)
        # dropped_fa22 can be missing if the procedure failed to produce an AsyCov row
        s = combine(
            groupby(df, :rho),
            :dropped_fa22 => (x -> mean(skipmissing(x))) => :drop_rate,
            :dropped_fa22 => (x -> count(ismissing, x)) => :n_missing,
            nrow => :n_sims,
        )
        sort!(s, :rho)
        return s
    end

    mixed_s = summarize(mixed)
    glimmix_s = summarize(glimmix)

    println("PROC MIXED drop rates (FA(2,2) row all zeros in AsyCov)")
    show(mixed_s; allrows=true, allcols=true)
    println("\n\nPROC GLIMMIX drop rates (FA(2,2) row all zeros in AsyCov)")
    show(glimmix_s; allrows=true, allcols=true)
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_fa0_drop_results()
end
