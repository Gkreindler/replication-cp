include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using CSV
using DelimitedFiles
using Glob

using DataFrames
using FreqTables
using Statistics
using StatsBase
using Printf

using Loess

using Plots
using Plots.PlotMeasures
# using LatexStrings

## All functions to read and organize GMM results
    include("../functions-table.jl")

## Paths
    savepath = rootpath_base * "paper/figures/smfigure4/"
    isdir(savepath) || mkdir(savepath)
    # for table SM.12 use the do file: "Table SM12 finite sample check.do"

    # num_id_type = "asy"
    num_id_type = "fin"

## Read all true parameters and GMM results -- FINITE
    all_results_fin = []
    for idx=1:120
        println("reading results ", idx)

        ### Read random parameters
        folder_path = rootpath_base * "data/coded_model/simdata_finite/sim" * string(idx) * "/"
        file_path = folder_path * "true_params.csv"
        true_params = readdlm(file_path, ',', Float64)
        true_params_df = DataFrame(transpose(true_params) |> Array, :auto)
        true_params_df[!, "type"] .= "finite"

        ### Read GMM results
        folder_path = rootpath_base * "model/estimation_results/num_identification/results_finite_sim" * string(idx) * "/"
        file_path = folder_path * "gmm-stage2-all.csv"

        gmm_df = CSV.read(file_path, DataFrame)
        gmm_df[!, "idx"] .= idx
        gmm_df = gmm_df[gmm_df.is_best_vec .== 1, :]

        push!(all_results_fin, hcat(true_params_df, gmm_df))
    end
    all_results_fin = vcat(all_results_fin...)

## Read all true parameters and GMM results -- ASYMPTOTIC
    all_results_asy = []
    for idx=1:120
        println("reading results ", idx)

        ### Read random parameters
        folder_path = rootpath_base * "data/coded_model/simdata_asymptotic/sim" * string(idx) * "/"
        file_path = folder_path * "true_params.csv"
        true_params = readdlm(file_path, ',', Float64)
        true_params_df = DataFrame(transpose(true_params) |> Array, :auto)
        true_params_df[!, "type"] .= "asymptotic"

        ### Read GMM results
        folder_path = rootpath_base * "model/estimation_results/num_identification/results_asym_sim" * string(idx) * "/"
        file_path = folder_path * "gmm-stage2-all.csv"

        gmm_df = CSV.read(file_path, DataFrame)
        gmm_df[!, "idx"] .= idx
        gmm_df = gmm_df[gmm_df.is_best_vec .== 1, :]

        push!(all_results_asy, hcat(true_params_df, gmm_df))
    end
    all_results_asy = vcat(all_results_asy...)

## Normalize alpha, beta_E and beta_L
    all_results_asy.param_1 .*= 6.0
    all_results_fin.param_1 .*= 6.0
    all_results_asy.x1 .*= 6.0
    all_results_fin.x1 .*= 6.0

    # beta_E
    all_results_asy.param_4 .*= 6.0
    all_results_fin.param_4 = min.(2000, 6.0 * all_results_fin.param_4)
    all_results_asy.x4 .*= 6.0
    all_results_fin.x4 .*= 6.0

    # beta_L
    all_results_asy.param_5 .*= 6.0
    all_results_fin.param_5 = min.(1200, 6.0 * all_results_fin.param_5)
    all_results_asy.x5 .*= 6.0
    all_results_fin.x5 .*= 6.0

    # sigma^DT
    all_results_fin.param_6 = min.(70, all_results_fin.param_6)

## Analysis

    # alpha
    maximum(all_results_asy.x1) |> display # 1158.8539955417473
    maximum(all_results_fin.param_1) |> display # 2033
    myfs = 10

    # xs = all_results_fin.x1
    # ys = all_results_fin.param_1
    # model = loess(xs, ys, span=0.5)
    # us = range(extrema(xs)...; step = 10)
    # vs = Loess.predict(model, us)

    # markerstrokewidth=0
    # markerstrokewidth=0

    plot(0:100:1200, 0:100:1200, color=:black, label="")
    scatter!(all_results_fin.x1, all_results_fin.param_1, label="finite sample",
            color=:blue, markerstrokewidth=0, markershape=:utriangle, markersize=5)
    scatter!(all_results_asy.x1, all_results_asy.param_1, label="asymptotic",
            color=:white, markerstrokewidth=2, markerstrokecolor=:red, markersize=5)
    plot!(xlim=(0, 1200), xticks=0:500:1000,
            ylim=(0, 2100), yticks=0:500:2000,
            title="VOTT \$\\alpha\$",
            # xlabel="true parameter",
            ylabel="estimated parameter",
            xtickfontsize=myfs, ytickfontsize=myfs,
            xguidefontsize=myfs, yguidefontsize=myfs,
            legend=:topleft,
            right_margin=5mm,
            legendfontsize=myfs, legendtitlefontsize=myfs)
    p1 = plot!(size=(350,350))
    # savefig(savepath * "fig_alpha.pdf")
    p1 |> display

# gamma
    maximum(all_results_asy.x2) |> display # 131
    maximum(all_results_fin.param_2) |> display #161

    plot(0:100:700, 0:100:700, color=:black, label="")
    scatter!(all_results_fin.x2, all_results_fin.param_2, label="finite sample",
            color=:blue, markershape=:utriangle, markerstrokewidth=0, markersize=5)
    scatter!(all_results_asy.x2, all_results_asy.param_2, label="asymptotic",
            color=:white, markerstrokewidth=2, markerstrokecolor=:red, markersize=5)
    plot!(xlim=(0, 150), xticks=0:50:150,
            ylim=(0, 175), yticks=0:50:150,
            title="Route Switch Cost \$\\gamma\$",
            xlabel="true parameter", ylabel="estimated parameter",
            xtickfontsize=myfs, ytickfontsize=myfs,
            xguidefontsize=myfs, yguidefontsize=myfs,
            legend=nothing,
            right_margin=5mm,
            legendfontsize=myfs, legendtitlefontsize=myfs)
    p2 = plot!(size=(350,350))
    # savefig(savepath * "fig_alpha.pdf")
    p2 |> display

# beta_E
    maximum(all_results_asy.x4) |> display # 1009
    maximum(all_results_fin.param_4) |> display #

    plot(0:100:1200, 0:100:1200, color=:black, label="")
    scatter!(all_results_fin.x4, all_results_fin.param_4, label="finite sample",
            color=:blue, markershape=:utriangle, markerstrokewidth=0, markersize=5)
    scatter!(all_results_asy.x4, all_results_asy.param_4, label="asymptotic",
            color=:white, markerstrokewidth=2, markerstrokecolor=:red, markersize=5)
    plot!(xlim=(0, 1200), xticks=0:500:1000,
            ylim=(0, 2100), yticks=0:500:2000,
            title="Schedule Cost Early \$\\beta_E\$",
            # xlabel="true parameter", ylabel="estimated parameter",
            xtickfontsize=myfs, ytickfontsize=myfs,
            xguidefontsize=myfs, yguidefontsize=myfs,
            legend=nothing,
            right_margin=5mm,
            legendfontsize=myfs, legendtitlefontsize=myfs)
    p3 = plot!(size=(350,350))
    # savefig(savepath * "fig_alpha.pdf")
    p3 |> display

# beta_L
    maximum(all_results_asy.x5) |> display # 635
    maximum(all_results_fin.param_5) |> display #

    plot(0:100:700, 0:100:700, color=:black, label="")
    scatter!(all_results_fin.x5, all_results_fin.param_5, label="finite sample",
            color=:blue, markershape=:utriangle, markerstrokewidth=0, markersize=5)
    scatter!(all_results_asy.x5, all_results_asy.param_5, label="asymptotic",
            color=:white, markerstrokewidth=2, markerstrokecolor=:red, markersize=5)
    plot!(xlim=(0, 650), xticks=0:200:600,
            ylim=(0, 1250), yticks=0:400:1200,
            title="Schedule Cost Late \$\\beta_L\$",
            # xlabel="true parameter", ylabel="estimated parameter",
            xtickfontsize=myfs, ytickfontsize=myfs,
            xguidefontsize=myfs, yguidefontsize=myfs,
            legend=nothing,
            right_margin=5mm,
            legendfontsize=myfs, legendtitlefontsize=myfs)
    p4 = plot!(size=(350,350))
    # savefig(savepath * "fig_alpha.pdf")
    p4 |> display

# sigma^DT
    maximum(all_results_asy.x6) |> display # 34
    maximum(all_results_fin.param_6) |> display #

    plot(0:10:80, 0:10:80, color=:black, label="")
    scatter!(all_results_fin.x6, all_results_fin.param_6, label="finite sample",
            color=:blue, markershape=:utriangle, markerstrokewidth=0, markersize=5)
    scatter!(all_results_asy.x6, all_results_asy.param_6, label="asymptotic",
            color=:white, markerstrokewidth=2, markerstrokecolor=:red, markersize=5)
    plot!(xlim=(0, 40), xticks=0:10:40,
            ylim=(0, 80), yticks=0:20:80,
            title="Logit Dep. Time \$\\sigma^{DT}\$",
            xlabel="true parameter",
            # ylabel="estimated parameter",
            xtickfontsize=myfs, ytickfontsize=myfs,
            xguidefontsize=myfs, yguidefontsize=myfs,
            legend=nothing,
            right_margin=5mm,
            legendfontsize=myfs, legendtitlefontsize=myfs)
    p5 = plot!(size=(350,350))
    # savefig(savepath * "fig_alpha.pdf")
    p5 |> display

# mu
    maximum(all_results_asy.x7) |> display # 102
    maximum(all_results_fin.param_7) |> display # 109

    plot(0:10:80, 0:10:80, color=:black, label="")
    scatter!(all_results_fin.x7, all_results_fin.param_7, label="finite sample",
            color=:blue, markershape=:utriangle, markerstrokewidth=0, markersize=5)
    scatter!(all_results_asy.x7, all_results_asy.param_7, label="asymptotic",
            color=:white, markerstrokewidth=2, markerstrokecolor=:red, markersize=5)
    plot!(xlim=(0, 110), xticks=0:50:100,
            ylim=(0, 110), yticks=0:50:100,
            title="Logit Route Choice \$\\sigma^{R}\$",
            xlabel="true parameter",
            # ylabel="estimated parameter",
            xtickfontsize=myfs, ytickfontsize=myfs,
            xguidefontsize=myfs, yguidefontsize=myfs,
            legend=nothing,
            right_margin=5mm,
            legendfontsize=myfs, legendtitlefontsize=myfs)
    p6 = plot!(size=(350,350))
    # savefig(savepath * "fig_alpha.pdf")
    p6 |> display

# combine!
    plot(p1, p3, p4, p2, p5, p6, layout=(2,3), left_margin=5mm, bottom_margin=1mm)
    plot!(size=(1050,700)) |> display
    savefig(savepath * "fig_num_id_all.pdf")
