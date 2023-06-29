
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

# using GLM
# using QuantileRegressions
# using QuantReg
# using RegressionTables

# using Loess

# using Plots
# using Plots.PlotMeasures
# using LatexStrings

## All functions to read and organize GMM results
    include("../functions-table.jl")

    ## ! Run both options!
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
        true_params_df = DataFrame(transpose(true_params) |> Array, :auto) # , :auto
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
    all_results_fin.param_4 = 6.0 * all_results_fin.param_4
    all_results_asy.x4 .*= 6.0
    all_results_fin.x4 .*= 6.0

    # beta_L
    all_results_asy.param_5 .*= 6.0
    all_results_fin.param_5 = 6.0 * all_results_fin.param_5
    all_results_asy.x5 .*= 6.0
    all_results_fin.x5 .*= 6.0

    # sigma^DT
    all_results_fin.param_6 = all_results_fin.param_6

## Save to later analyze in Stata
    CSV.write(rootpath_base * "data/coded_model/simdata/all_results_fin.csv", all_results_fin)
