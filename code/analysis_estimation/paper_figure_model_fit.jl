"""
	Model fit figure SM.2. (all four panels)
	
	Date: June 25 2023 (Jan 13 2022)
"""

include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using Plots

using CSV
using JSON
using DataFrames

using Optim
using LsqFit
using FiniteDifferences
using LineSearches

using StatsBase
using Statistics

using LinearAlgebra
using Random
using Distributions

using TableView
using FreqTables
using Loess
using Interpolations
using Dates
using BenchmarkTools
using NonNegLeastSquares

println("Using statements completed!...")
	include("functions-loadprep.jl")
	include("functions-graphs-full.jl")
	include("functions-model.jl")
	include("functions-model-definitions.jl")
	include("gmm_wrappers.jl")
	include("functions-gmm-full.jl")

	include("functions-table.jl")

println("Import other functions completed!...")

Random.seed!(123)


## Start analysis
    # Paths
	rootpath_input  = rootpath_base * "data/coded_model/"

	sm_folder = rootpath_base * "paper/figures/smfigure2/"
	isdir(sm_folder) || mkdir(sm_folder)

## 1. parameter values
	theta_fixed_dict = Dict{String,Float64}(
		"delta" => 90.0,
		"fe_late" => 0.0,
		"fe_w1" => 0.0,
		"fe_w2" => 0.0,
		"fe_w3" => 0.0,
		"fe_w4" => 0.0,
		"gamma_01_factor" => 1.0,  # switch cost 0->1 is this times larger than 0->1
		"prob_respond_alpha" => 100.0,
		"prop_to_wage" => 0.0
	)

	rootpath_results      =  rootpath_base * "model/estimation_results/results_delta_90/"
    rootpath_boot_results = [rootpath_base * "model/estimation_results/results_boot_delta_90/"]
    param_names = ["alpha", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)

	theta_optimum = select(main2_df, r"^param_") |> Array |> vec

	theta_optimum |> display

## 2. load data
	hdgrid_2525 = vec(Array(range(-2.5,stop=2.5,length=61)))

	# @elapsed BIGMAT, CHARGEMAT, COMPMAT, DATAMOMS = load_and_prep(rootpath_input)
	DATA = load_and_prep(rootpath_input)

	DATA["DATAMOMS"]["hdgrid_2525"] = hdgrid_2525
	# n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = DATA["DATAMOMS"]["n_resp_vec"]

	# model paramters and
	theta_fixed_dict = Dict{String,Float64}(
		"delta" => 90.0,
		"fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0, "fe_w3" => 0.0, "fe_w4" => 0.0,
		"gamma_01_factor" => 1.0,  # switch cost 0->1 is this times larger than 0->1
		"prob_respond_alpha" => 100.0,
		"prop_to_wage" => 0.0
	)
	ESTIMATION_PARAMS = model_dt_dynamic_params(theta_fixed_dict=theta_fixed_dict)

## Simulate the full model
	sim_model = solve_dynamicvot(
		theta=theta_optimum,
		param_factors=ESTIMATION_PARAMS["param_factors"],
		BIGMAT=DATA["BIGMAT"],
		CHARGEMAT=DATA["CHARGEMAT"],
		COMPMAT=DATA["COMPMAT"],
		DATAMOMS=DATA["DATAMOMS"])

## Graph: dynamic route choice

	get_params_fn(theta_optimum, ESTIMATION_PARAMS["param_factors"])

	plot_dynamic_vott_path = sm_folder * "panel_C_dynamic_vott.pdf"

	draw_moms_vot_dynamic(
		theta=theta_optimum,
		param_factors=ESTIMATION_PARAMS["param_factors"],
		DATAMOMS=DATA["DATAMOMS"],
		SIMMODEL=sim_model,
		display_plot=true,
		save_path=plot_dynamic_vott_path)


## Graph: route choice heterogeneity
	plot_vott_het_path = sm_folder * "panel_D_vott_het.pdf"
	# plot_vott_het_path = ""

	draw_moms_vot_dynamic_het(
		theta=theta_optimum,
		param_factors=ESTIMATION_PARAMS["param_factors"],
		DATAMOMS=DATA["DATAMOMS"],
		SIMMODEL=sim_model,
		display_plot=true,
		save_path=plot_vott_het_path)

## Graph: departure time fit
	plot_dt_path  =  sm_folder * "panel_A_dt.pdf"
	plot_dt_path2 =  sm_folder * "panel_B_dt_control.pdf"
	# plot_dt_path = ""
	# plot_dt_path2 = ""

	draw_moms_dt(
		theta=theta_optimum,
		param_factors=ESTIMATION_PARAMS["param_factors"],
		DATAMOMS=DATA["DATAMOMS"],
		SIMMODEL=sim_model,
		display_plot=true, display_plot2=true,
		save_path=plot_dt_path, save_path2=plot_dt_path2)
