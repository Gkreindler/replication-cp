
include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using Plots

using CSV
using DelimitedFiles
using Glob
using BSON

using DataFrames
using FreqTables
using Statistics
using StatsBase
using Printf

using InlineStrings

## All functions to read and organize GMM results
    include("functions-table.jl")


    function get_ci(mymat)
        ncol = size(mymat)[2]

        cis = zeros(2, ncol)

        for icol=1:ncol
            cis[1, icol] = percentile(mymat[:,icol] |> vec, 2.5)
            cis[2, icol] = percentile(mymat[:,icol] |> vec, 97.5)
        end
        return cis
    end

### PATHS
    fig5_path = rootpath_base * "paper/figures/figure5/"
    isdir(fig5_path) || mkdir(fig5_path)

    sm8b_path = rootpath_base * "paper/figures/smfigure8/"
    isdir(sm8b_path) || mkdir(sm8b_path)

    # simulation results
    rootpath_sims = rootpath_base * "model/simulation_results/main/"
    rootpath_boot = rootpath_base * "model/simulation_results/main_boot/"

# plot settings
    Plots.scalefontsizes()
    # Plots.scalefontsizes(1.1)

## Read BSON objects -- Bootstrap

boot_nash_delays = zeros(120,85)
boot_so_delays = zeros(120,85)
for i=1:120
    # Paths
    path = rootpath_boot * "boot_" * string(i) * "/so_eqm.bson"
    so = BSON.load(path)

    # Average travel time
    nash_delays = mean(so[:nash_ttimes] ./ so[:agents].mean_km, dims=1)
    so_delays   = mean(so[:so_ttimes]   ./ so[:agents].mean_km, dims=1)

    boot_nash_delays[i,:] .= nash_delays[1:85]
    boot_so_delays[i,:] .= so_delays[1:85]
end

    ### Compute 95% CI
    boot_nash_delays_cis = get_ci(boot_nash_delays)
    boot_so_delays_cis = get_ci(boot_so_delays)

## Read main simulation
    path = rootpath_sims * "so_eqm.bson"
    so = BSON.load(path)

    nash_delays = mean(so[:nash_ttimes] ./ so[:agents].mean_km, dims=1)
    nash_delays = nash_delays[1:85]
    so_delays   = mean(so[:so_ttimes]   ./ so[:agents].mean_km, dims=1)
    so_delays   = so_delays[1:85]

# X axis
    hdgrid = 5:(1/12):14
    hdgrid = hdgrid[1:85]

 plot(hdgrid, nash_delays              , color="black", linestyle=:dash, linewidth=3, label="Nash equilibrium")
plot!(hdgrid, boot_nash_delays_cis[1,:] |> vec, color="black", linestyle=:dash, linewidth=1.5, alpha=0.75, label="Nash 95% CI")
plot!(hdgrid, boot_nash_delays_cis[2,:] |> vec, color="black", linestyle=:dash, linewidth=1.5, alpha=0.75, label="")
plot!(hdgrid, so_delays                     , color="red", linewidth=3, label="Social Optimum")
plot!(hdgrid, boot_so_delays_cis[2,:] |> vec, color="red", linewidth=1.5, alpha=0.75, label="SO 95% CI")
myplot =
plot!(hdgrid, boot_so_delays_cis[1,:] |> vec, color="red", linewidth=1.5, alpha=0.75, label="",
        xlabel="Departure Time",
        ylabel="Travel Delay (min/km)",
        xticks=5:12,
        xtickfontsize=10,ytickfontsize=10,
        xguidefontsize=10, yguidefontsize=10,
        legend=:topleft)
plot!(myplot, legendfontsize=10, legendtitlefontsize=10)

# save to file
    savefig(myplot, fig5_path * "figure5.pdf")

    myplot |> display




    
## Departure Rates
    hdgrid = 5:(1/12):14
    nash_probs = sum(so[:nash_choice_probs], dims=1) |> vec
    nash_probs = nash_probs ./ sum(nash_probs)
    # nash_probs = nash_probs[1:85]
    so_probs = sum(so[:so_choice_probs], dims=1) |> vec
    so_probs = so_probs ./ sum(so_probs)
    # so_probs = so_probs[1:85]

    plot(hdgrid, nash_probs, color="black", linestyle=:dash, linewidth=2, label="Nash equilibrium")
    # plot!(hdgrid, so_probs)
    myplot = plot!(hdgrid, so_probs, color="red", linewidth=2, label="Social Optimum",
            xlabel="Departure Time",
            ylabel="Departure Rate",
            xticks=5:14,
            xtickfontsize=10,ytickfontsize=10,
            xguidefontsize=10, yguidefontsize=10,
            legend=:topleft)
    plot!(myplot, legendfontsize=10, legendtitlefontsize=10)

    # save to file
    savefig(myplot, sm8b_path * "panel_B_volumes_nash_so.pdf")

    myplot |> display
