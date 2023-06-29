include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using Plots
using Plots.PlotMeasures

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


## Paths
    plot_save_path = rootpath_base * "paper/figures/smfigure8/"
    isdir(plot_save_path) || mkdir(plot_save_path)

    # simulation results
    rootpath_sims = rootpath_base * "model/simulation_results/main/"
    rootpath_boot = rootpath_base * "model/simulation_results/main_boot/"

# plot settings
    Plots.scalefontsizes()

## Read BSON objects -- Bootstrap

boot_nash_delays = zeros(120,85)
boot_so_delays = zeros(120,85)
for i=1:120
    # Paths
    path = rootpath_boot * "boot_" * string(i) * "/so_eqm.bson"
    main_so_eqm = BSON.load(path)

    # Average travel time
    nash_delays = mean(main_so_eqm[:nash_ttimes] ./ main_so_eqm[:agents].mean_km, dims=1)
    so_delays   = mean(main_so_eqm[:so_ttimes]   ./ main_so_eqm[:agents].mean_km, dims=1)

    boot_nash_delays[i,:] .= nash_delays[1:85]
    boot_so_delays[i,:] .= so_delays[1:85]
end

    ### Compute 95% CI
    boot_nash_delays_cis = get_ci(boot_nash_delays)
    boot_so_delays_cis = get_ci(boot_so_delays)

## Read main simulation
    path = rootpath_sims * "so_eqm.bson"
    main_so_eqm = BSON.load(path)

    hdgrid = 5:(1/12):14

    nash_delays = mean(main_so_eqm[:nash_ttimes] ./ main_so_eqm[:agents].mean_km, dims=1)
    nash_delays = nash_delays[1:85]
    so_delays   = mean(main_so_eqm[:so_ttimes]   ./ main_so_eqm[:agents].mean_km, dims=1)
    so_delays   = so_delays[1:85]

    # x axis
    hdgrid = hdgrid[1:85]

    # avg km
    # so_km = sum(main_so_eqm[:so_choice_probs] .* main_so_eqm[:agents].mean_km, dims=1) ./
    #         sum(main_so_eqm[:so_choice_probs], dims=1)
    # so_km = so_km[1:85] |> vec
    # plot(hdgrid, so_km)

    so_msc = -main_so_eqm[:so_charges]   ./ main_so_eqm[:agents].mean_km
    so_msc_ps = zeros(5, 109)
    ps = [5, 25, 50, 75, 90]
    for i=1:5
        for j=1:109
            so_msc_ps[i,j] = percentile(so_msc[:,j], ps[i])
        end
    end


    plot(hdgrid, so_msc_ps[5,1:85]  |> vec, color=:black,  linewidth=2, alpha=1  , label="p95")
    plot!(hdgrid, so_msc_ps[4,1:85]  |> vec, color=:black, linewidth=2, alpha=0.8, label="p75")
    plot!(hdgrid, so_msc_ps[3,1:85]  |> vec, color=:black, linewidth=2, alpha=0.6, label="median")
    plot!(hdgrid, so_msc_ps[2,1:85]  |> vec, color=:black, linewidth=2, alpha=0.4, label="p25")
    plot!(hdgrid, so_msc_ps[1,1:85]  |> vec, color=:black, linewidth=2, alpha=0.2, label="p5",
            xlabel="Departure Time",
            xlims=(5, 12), xticks=5:12,
            xtickfontsize=10, ytickfontsize=10,
            xguidefontsize=10, yguidefontsize=10,
            legendtitle="Optimal Charges:\n(per KM)", legendfontsize=8,
            legend_title_font_pointsize=8,
            right_margin=15mm,
            ylabel="Charge per KM (INR)", legend=:topleft)

    p = twinx()
    plot!(p, hdgrid, main_so_eqm[:so_instant_delays][1:85] |> vec, color="blue",
            linewidth=2, linestyle=:dash,
            xlims=(5, 12), xticks=5:12,
            xtickfontsize=10,
            label="Instantaneous Delay", ytickfontsize=10,
            legendfontsize=8, legend=:bottomright)
    myplot = plot!(p, hdgrid, so_delays |> vec, color="red", linewidth=2,
            xlims=(5, 12), xticks=5:12,
            xtickfontsize=10,
            ylabel="Delay (min/km)", ytickfontsize=10,
            label="Average Trip Delay",
            legendtitle="Social Optimum:",
            legend_title_font_pointsize=8,
            right_margin=15mm,
            bottom_margin=8mm,
            legend=:bottomright)
    plot!(myplot, legendfontsize=10, legendtitlefontsize=10, fg_legend = :false)

# save to file
savefig(myplot, plot_save_path * "panel_A_so_charges.pdf")

myplot |> display
