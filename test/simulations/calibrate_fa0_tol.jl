using CSV
using DataFrames
using Statistics

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels

function _lowercase_names!(df::DataFrame)
    rename!(df, Dict(n => Symbol(lowercase(String(n))) for n in names(df)))
    return df
end

function _covp_cols(df::DataFrame)
    return [n for n in propertynames(df) if occursin(r"^covp[0-9]+$", String(n))]
end

function sas_drop_from_asycov(asycov_csv::AbstractString)
    df = DataFrame(CSV.File(asycov_csv))
    _lowercase_names!(df)
    covp = _covp_cols(df)

    out = DataFrame(; sim_id=Int[], dropped_fa22=Bool[])
    for sdf in groupby(df, :sim_id)
        fa = sdf[lowercase.(String.(sdf.covparm)) .== "fa(2,2)", :]
        if nrow(fa) == 0
            push!(out, (sim_id=sdf.sim_id[1], dropped_fa22=true))
            continue
        end
        vals = [coalesce(fa[1, c], 0.0) for c in covp]
        dropped = all(abs.(Float64.(vals)) .== 0.0)
        push!(out, (sim_id=sdf.sim_id[1], dropped_fa22=dropped))
    end
    return out
end

function fit_julia_models(; sim_csv::AbstractString, out_csv::AbstractString)
    sim = DataFrame(CSV.File(sim_csv))
    _lowercase_names!(sim)

    fm = @formula(y ~ 1 + time + (1 + time | subject))

    res = DataFrame(;
        sim_id=Int[],
        rho=Float64[],
        converged=Bool[],
        sigma=Float64[],
        l11=Float64[],
        l21=Float64[],
        l22=Float64[],
    )

    for (k, sdf) in enumerate(groupby(sim, :sim_id))
        sim_id = Int(sdf.sim_id[1])
        rho = Float64(sdf.rho[1])

        converged = false
        sigma = NaN
        l11 = NaN
        l21 = NaN
        l22 = NaN

        try
            m = fit(MixedModel, fm, sdf; REML=true)
            # MixedModels reports convergence in optsum; fallback to true if missing
            converged = try
                getproperty(m.optsum, :conv)
            catch
                true
            end
            sigma = m.sigma

            L = m.reterms[1].λ * m.sigma
            l11 = L[1, 1]
            l21 = L[2, 1]
            l22 = L[2, 2]
        catch
            # leave NaNs, converged=false
        end

        push!(
            res,
            (
                sim_id=sim_id,
                rho=rho,
                converged=converged,
                sigma=sigma,
                l11=l11,
                l21=l21,
                l22=l22,
            ),
        )

        if k % 25 == 0
            println("fit ", k, "/", length(unique(sim.sim_id)))
        end
    end

    CSV.write(out_csv, res)
    return out_csv
end

function _score_threshold(
    l22::AbstractVector{<:Real}, sas_drop::AbstractVector{Bool}, tol::Real
)
    pred = abs.(l22) .< tol
    tp = sum(pred .& sas_drop)
    tn = sum(.!pred .& .!sas_drop)
    fp = sum(pred .& .!sas_drop)
    fn = sum(.!pred .& sas_drop)
    err = fp + fn
    return (err=err, tp=tp, tn=tn, fp=fp, fn=fn)
end

function calibrate_fa0_tol(;
    sim_csv::AbstractString=joinpath(
        @__DIR__, "..", "data", "Data fa0 drop simulation.csv"
    ),
    mixed_asycov_csv::AbstractString=joinpath(
        @__DIR__, "..", "results", "Results fa0 drop simulation mixed asycov.csv"
    ),
    glimmix_asycov_csv::AbstractString=joinpath(
        @__DIR__, "..", "results", "Results fa0 drop simulation glimmix asycov.csv"
    ),
    julia_fit_csv::AbstractString=joinpath(
        @__DIR__, "..", "results", "Results fa0 drop simulation julia fits.csv"
    ),
)
    if !isfile(julia_fit_csv)
        println("Fitting MixedModels.jl for all sim_id...")
        fit_julia_models(; sim_csv=sim_csv, out_csv=julia_fit_csv)
    end

    julia = DataFrame(CSV.File(julia_fit_csv))
    _lowercase_names!(julia)

    sas_mixed = sas_drop_from_asycov(mixed_asycov_csv)
    sas_glimmix = sas_drop_from_asycov(glimmix_asycov_csv)

    df_mixed = leftjoin(julia, sas_mixed; on=:sim_id)
    df_glimmix = leftjoin(julia, sas_glimmix; on=:sim_id)

    function best_tol_for(df, label)
        ok = df[
            (.!isnan.(df.l22)) .& (df.converged .== true) .& (.!ismissing.(df.dropped_fa22)),
            :,
        ]
        println(label, ": usable fits = ", nrow(ok), " / ", nrow(df))

        # Candidate tolerances (log grid)
        tols = 10.0 .^ range(-12, -1; length=200)

        best = nothing
        best_tol_val = NaN
        for tol in tols
            sc = _score_threshold(ok.l22, Bool.(ok.dropped_fa22), tol)
            if best === nothing || sc.err < best.err
                best = sc
                best_tol_val = tol
            end
        end

        println(label, ": best fa0_tol ~= ", best_tol_val)
        println(label, ": confusion (tp/tn/fp/fn) = ", (best.tp, best.tn, best.fp, best.fn))
        return best_tol_val
    end

    best_tol_for(df_mixed, "PROC MIXED")
    best_tol_for(df_glimmix, "PROC GLIMMIX")
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_fa0_tol()
end
