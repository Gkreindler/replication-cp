
include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using Plots

using BSON
using InlineStrings

using CSV
using DelimitedFiles
using Glob

using DataFrames
using FreqTables
using Statistics
using StatsBase
using Printf



## All functions to read and organize GMM results
    include("functions-table.jl")



### PATHS
    sm8b_path = rootpath_base * "paper/figures/smfigure8/"
    isdir(sm8b_path) || mkdir(sm8b_path)

    # simulation results
    rootpath_sims = rootpath_base * "model/simulation_results/factors/"

    # read simulation results
all_delays = []
all_dwls = []
for i=2:9
    path = rootpath_sims * "factor_" * string(i) * "/" * "/so_eqm.bson"
    current_so_eqm = BSON.load(path)
    push!(all_delays, mean(current_so_eqm[:nash_ttimes] ./ current_so_eqm[:agents].mean_km, dims=1) |> vec)
    push!(all_dwls, (current_so_eqm[:so_logsum] - current_so_eqm[:nash_logsum])/abs(current_so_eqm[:nash_logsum]))
end

    max_delays = maximum.(all_delays)

    scatter(max_delays, all_dwls .* 100.0 ,
            xlabel="Maximum average travel delay (min/km)",
            ylabel="Deadweight Loss (% of Nash Welfare)",
            label="Fraction of All Daily Trips",
            xtickfontsize=10,ytickfontsize=10,
            xguidefontsize=10, yguidefontsize=10,
            legend=:topleft)
    vline!([max_delays[4]], label="Benchmark Specification", color=:red,)
    strs = text.( f2.((2.0:9.0) ./ 24.0), :black, :right, 10)
    myplot = annotate!(max_delays .- 0.05, all_dwls .* 100.0 .+ 0.2, strs)
    plot!(myplot, legendtitlefontsize=10, legendfontsize=10)

    myplot |> display

    savefig(myplot, sm8b_path * "panel_D_density_factor.pdf")

## Elasticity
    using GLM
    df = DataFrame(
        "max_delays" => max_delays,
        "factor" => 2.0:9.0,
        "dwl" => all_dwls .* 100.0
        )

    min_delay = minimum(all_delays[1])
    display(min_delay)

    df[!, :log_excess_delay] = log.(df.max_delays .- min_delay)
    df[!, :log_dwl] = log.(df.dwl)
    df[!, :log_factor] = log.(df.factor)

    ols_trip_volume = lm(@formula(log_dwl ~ log_factor), df)
    ols_trip_volume |> display

    ols_excess_delay = lm(@formula(log_dwl ~ log_excess_delay), df)
    ols_excess_delay |> display

    elasticity_dw_peak_delay = coef(ols_excess_delay)[2]

## Text

    # higher
    (all_dwls[5] .* 100.0) |> display
    max_delays[5] |> display

    # benchmark
    (all_dwls[4] .* 100.0) |> display
    max_delays[4] |> display

    # lower
    (all_dwls[3] .* 100.0) |> display
    max_delays[3] |> display

### Save to file
    mytext = "The elasticity of Deadweight Loss with respect to equilibrium peak excelle travel delay is = " * string(elasticity_dw_peak_delay) * "\n\n"

    mytext *= "If eqm travel delay were " * string(max_delays[5]) * " then DWL would be " * string((all_dwls[5] .* 100.0)) * "\n\n"
    mytext *= "If eqm travel delay were " * string(max_delays[3]) * " then DWL would be " * string((all_dwls[3] .* 100.0)) * "\n\n"
    mytext *= "Benchmark: eqm     delay " * string(max_delays[4]) * " and DWL =         " * string((all_dwls[4] .* 100.0)) * "\n\n"


    text_results_path = rootpath_base * "paper/text_results/sims_factor.txt"
    file = open(text_results_path, "w")
    write(file, mytext)
    close(file)

    println(mytext)