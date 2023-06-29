
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

using Loess

using InlineStrings

## All functions to read and organize GMM results
    include("functions-table.jl")


## PATHS
    sm_path = rootpath_base * "paper/figures/smfigure10/"
    isdir(sm_path) || mkdir(sm_path)

    rootpath_sims = rootpath_base * "model/simulation_results/main_2routes/v15_low_mu/"

## Read main simulation
    path = rootpath_sims * "so_eqm.bson"
    so = BSON.load(path)

    hdgrid = 5:(1/12):14
    hdgrid = hdgrid[1:85]

    so |> display

    nash_delays0 = mean(so[:nash_ttimes0] ./ so[:agents].mean_km, dims=1) |> vec
    nash_delays1 = mean(so[:nash_ttimes1] ./ so[:agents].mean_km, dims=1) |> vec
    nash_delays0 = nash_delays0[1:85]
    nash_delays1 = nash_delays1[1:85]

    temp0 = sum(so[:nash_choice_probs0], dims=1)
    temp1 = sum(so[:nash_choice_probs1], dims=1)
    nash_frac1 = temp1 ./ (temp0 .+ temp1)
    nash_frac1 = nash_frac1 |> vec
    nash_frac1 = nash_frac1[1:85]

    so_delays0   = mean(so[:so_ttimes0]   ./ so[:agents].mean_km, dims=1) |> vec
    so_delays1   = mean(so[:so_ttimes1]   ./ so[:agents].mean_km, dims=1) |> vec
    so_delays0 = so_delays0[1:85]
    so_delays1 = so_delays1[1:85]

    temp0 = sum(so[:so_choice_probs0], dims=1)
    temp1 = sum(so[:so_choice_probs1], dims=1)
    so_frac1 = temp1 ./ (temp0 .+ temp1)
    so_frac1 = so_frac1 |> vec
    so_frac1 = so_frac1[1:85]


## Graphs equilibrium and SO delays
    plot(hdgrid, nash_delays0,  color=:black, linewidth=2, label="Nash route 0")
    plot!(hdgrid, nash_delays1, color=:black, linewidth=2, label="Nash route 1", linestyle=:dash)

    plot!(hdgrid, so_delays0, color=:red, linewidth=2, label="S.O. route 0")
    myplot = plot!(hdgrid, so_delays1, color=:red, linewidth=2, label="S.O. route 1", linestyle=:dash,
            xlabel="Departure Time",
            ylabel="Travel Delay (min/km)",
            xticks=5:12, ylim=(1.5,4.2),
            legendtitlefontsize=10,
            xtickfontsize=10,ytickfontsize=10,
            xguidefontsize=10, yguidefontsize=10,
            legend=:topleft, legendfontsize=10)

    # save to file
    savefig(myplot, sm_path * "2routes_delays.pdf")

    myplot |> display

## Graph fraction taking route 1
    # plot(hdgrid, nash_frac1, ylim=(0,1))
    # plot!(hdgrid, so_frac1, ylim=(0,1))


    xs = hdgrid
    ys = nash_frac1
    model = loess(xs, ys, span=0.5)
    us = range(extrema(xs)...; step = 0.1)
    vs_nash = Loess.predict(model, us)

    ys = so_frac1
    model = loess(xs, ys, span=0.5)
    vs_so = Loess.predict(model, us)

    scatter(hdgrid, nash_frac1, markercolor=:white, markerstrokecolor=:black,
            label="Nash")
    hline!([0.5], color=:gray, linestyle=:dash, label="")
    plot!(us, vs_nash, color=:black, linewidth=2, label="Nash")
    scatter!(hdgrid, so_frac1, markercolor=:white, markerstrokecolor=:red,
            markershape=:square, markersize=2, label="Social Optimum")
    myplot = plot!(us, vs_so, color=:red, linewidth=2, linestyle=:dash,
        label="Social Optimum",
        ylim=(0.4,0.9), yticks=0.4:0.05:0.9, xticks=5:12,
        xlabel="Departure Time",
        ylabel="Fraction Choosing Route 1",
        legendtitlefontsize=10,
        legendfontsize=10,
        xtickfontsize=10, ytickfontsize=10,
        xguidefontsize=10, yguidefontsize=10)

    # save to file
    savefig(myplot, sm_path * "2routes_r1_prob.pdf")

    myplot |> display
