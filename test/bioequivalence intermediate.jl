using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

# ============================================================================
# Fit Julia Model
# ============================================================================

df = DataFrame(CSV.File(joinpath(@__DIR__, "data", "Data bioequivalence.csv")))

# Model: correlated random effects for formulation per subject
# This maps to SAS: RANDOM formulation / TYPE=FA0(2) SUBJECT=subject
fm = @formula(log_data ~ 1 + formulation + sequence + period + (0 + formulation | subject))

m = fit(
    MixedModel,
    fm,
    df;
    REML=true,
    contrasts=Dict(
        :period => DummyCoding(; base=4),
        :formulation => DummyCoding(; base="T"),
        :sequence => DummyCoding(; base="BABA"),
    ),
)

# ============================================================================
# SAS KR Results
# ============================================================================

@testset "SAS KR Results" begin
    res_pe = DataFrame(
        CSV.File(
            joinpath(@__DIR__, "results", "Results bioequivalence intermediate sas kr.csv")
        ),
    )

    kr = small_sample_adjust(m, KenwardRoger(parameterization=FactorAnalytic()))

    julia_coefs = coef(m)
    julia_se = sqrt.(diag(kr.varcovar_adjusted))
    julia_dof = kr.ν

    # Intercept
    row_int_idx = findfirst(r -> r["Effect"] == "Intercept", eachrow(res_pe))
    @test !isnothing(row_int_idx)
    row_int = res_pe[row_int_idx, :]

    @test isapprox(row_int["Estimate"], julia_coefs[1], rtol=1e-4)
    @test isapprox(row_int["StdErr"], julia_se[1], rtol=1e-4)
    @test isapprox(row_int["DF"], julia_dof[1], rtol=1e-4)

    # Formulation (R)
    row_form_idx = findfirst(
        r -> r["Effect"] == "formulation" && r["formulation"] == "R", eachrow(res_pe)
    )
    @test !isnothing(row_form_idx)
    row_form = res_pe[row_form_idx, :]

    @test isapprox(row_form["Estimate"], julia_coefs[2], rtol=1e-4)
    @test isapprox(row_form["StdErr"], julia_se[2], rtol=1e-4)
    @test isapprox(row_form["DF"], julia_dof[2], rtol=1e-4)

    # Sequence (ABAB)
    row_seq_idx = findfirst(
        r -> r["Effect"] == "sequence" && r["sequence"] == "ABAB", eachrow(res_pe)
    )
    @test !isnothing(row_seq_idx)
    row_seq = res_pe[row_seq_idx, :]

    @test isapprox(row_seq["Estimate"], julia_coefs[3], rtol=1e-4)
    @test isapprox(row_seq["StdErr"], julia_se[3], rtol=1e-4)
    @test isapprox(row_seq["DF"], julia_dof[3], rtol=1e-4)

    # Period (1, 2, 3)
    for (i, p) in enumerate([1, 2, 3])
        row_p_idx = findfirst(
            r -> r["Effect"] == "period" && parse(Int, r["period"]) == p, eachrow(res_pe)
        )
        @test !isnothing(row_p_idx)
        row_p = res_pe[row_p_idx, :]

        @test isapprox(row_p["Estimate"], julia_coefs[3 + i], rtol=1e-4)
        @test isapprox(row_p["StdErr"], julia_se[3 + i], rtol=1e-4)
        @test isapprox(row_p["DF"], julia_dof[3 + i], rtol=1e-4)
    end
end

# ============================================================================
# SAS Asymptotic Covariance (W)
# ============================================================================

@testset "SAS Asymptotic Covariance (W)" begin
    sas_asycov_df = DataFrame(
        CSV.File(
            joinpath(
                @__DIR__, "results", "Results bioequivalence intermediate sas asycov.csv"
            ),
        ),
    )

    W = vcov_varpar(m; parameterization=FactorAnalytic())
    vd = MixedModelsSmallSample.VarianceDecomposition(m, m.X, FactorAnalytic())
    names = vd.names

    # Mapping strategy:
    # Julia names: "sigma2", "L_1_1_1", "L_1_2_1" (L_1_2_2 is dropped)
    # SAS CovParm: "Residual", "FA(1,1)", "FA(2,1)", "FA(2,2)"
    # We construct a SAS matrix subset corresponding to Julia's names.

    function map_julia_to_sas(name)
        if name == "sigma2"
            return "Residual"
        elseif startswith(name, "L_")
            # Format L_block_i_j -> FA(i,j)
            parts = split(name, "_")
            if length(parts) == 4
                i = parts[3]
                j = parts[4]
                return "FA($i,$j)"
            end
        end
        return nothing
    end

    sas_indices = Int[]
    for name in names
        sas_name = map_julia_to_sas(name)
        idx = findfirst(==(sas_name), sas_asycov_df[!, "CovParm"])
        if idx !== nothing
            push!(sas_indices, idx)
        else
            error("Could not find SAS CovParm for Julia parameter: $name")
        end
    end

    sas_matrix = zeros(length(names), length(names))
    for (i, r_idx) in enumerate(sas_indices)
        for (j, c_idx) in enumerate(sas_indices)
            # SAS columns named CovP1, CovP2... corresponding to rows 1, 2...
            # So we need to map c_idx (row index of parameter) to "CovP$(c_idx)"
            col_name = "CovP$(c_idx)"
            sas_matrix[i, j] = sas_asycov_df[r_idx, col_name]
        end
    end

    @test isapprox(W, sas_matrix, rtol=1e-2)
end

# ============================================================================
# SAS SW Results
# ============================================================================

@testset "SAS SW Results" begin
    res_pe_sw = DataFrame(
        CSV.File(
            joinpath(@__DIR__, "results", "Results bioequivalence intermediate sas sw.csv")
        ),
    )

    sw = small_sample_adjust(m, Satterthwaite(parameterization=FactorAnalytic()))

    julia_dof_sw = sw.ν

    # Formulation (R)
    row_form_idx = findfirst(
        r -> r["Effect"] == "formulation" && r["formulation"] == "R", eachrow(res_pe_sw)
    )
    @test !isnothing(row_form_idx)
    row_form = res_pe_sw[row_form_idx, :]

    @test isapprox(row_form["DF"], julia_dof_sw[2], rtol=1e-2)
end
