"""

    Using estimated parameters, compute the bootstrap CI for the ratio of beta to alpha

    1) load estimated parameters
    2) load data
    3) invert ideal arrival time distributions
    4) define individual specific charge profile
	5) compute and summarize response!

"""

include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using Future

using Plots

using CSV
using DelimitedFiles
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
using NonNegLeastSquares

println("Using statements completed!...")
    include("../functions-loadprep.jl")
    include("../functions-model.jl")
	include("../functions-model-definitions.jl")
    include("../functions-table.jl")
    include("../gmm_wrappers.jl")
    include("../functions-gmm-full.jl")
    include("../functions-gmm-bootstrap.jl")
println("Import other functions completed!...")

Random.seed!(1212376)

## 0. prelims
	rootpath_input  = rootpath_base * "data/coded_model/"

	# fixed parameters
	theta_fixed_dict = Dict{String,Float64}(
		"delta" => 90.0,
		"fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0, "fe_w3" => 0.0, "fe_w4" => 0.0,
		"gamma_01_factor" => 1.0,  # switch cost 0->1 is this times larger than 0->1
		"prob_respond_alpha" => 100.0, "prop_to_wage" => 0.0
	)
	# model paramters and initial conditions
	gmm_options = Dict(
		"main_n_start_pts" => 0,
		"boot_n_start_pts" => 0
	)
	ESTIMATION_PARAMS = model_dt_dynamic_params(theta_fixed_dict=theta_fixed_dict, gmm_options=gmm_options)
	# momfn = model_dt_static

	# other model options
	ESTIMATION_PARAMS["hastar_inversion_type"] = "distribution"
	ESTIMATION_PARAMS["transition_moment_type"] = 0
	## Run entire GMM
	# run_gmm(momfn=momfn, ESTIMATION_PARAMS=ESTIMATION_PARAMS, gmm_options=gmm_options)

## load data
	DATA = load_and_prep(rootpath_input)
	BIGMAT = DATA["BIGMAT"]
	CHARGEMAT = DATA["CHARGEMAT"]
	COMPMAT = DATA["COMPMAT"]
	DATAMOMS = DATA["DATAMOMS"]


## Paths
	### Demand estimation params
    rootpath_params      =  rootpath_base * "model/estimation_results/results_delta_90/"
    rootpath_params_boot = [rootpath_base * "model/estimation_results/results_boot_delta_90/"]

## Read optimal parameters
	println("... reading estimates parameters (including bootstrap runs)")
	main_df1, boot_df1, main_df2, boot_df2, diag_df, full_df = get_all_results(
			rootpath_results=rootpath_params,
			rootpath_boot_results=rootpath_params_boot,
			boot_runs=120, pm_range=0.1, debug=false)

	# get matrix or vector of estimates ("^" anchors to start of the word)
	theta_optimum = select(main_df2, r"^param_") |> Array |> vec
	theta_boot    = select(boot_df2, r"^param_") |> Array

## ratio

	## unpack parameters
	theta_dict = get_params_fn(theta_optimum, ESTIMATION_PARAMS["param_factors"])

	beta_e_ratio = theta_dict["beta_e"] / theta_dict["alpha"]
	beta_l_ratio = theta_dict["beta_l"] / theta_dict["alpha"]

## boot
	beta_e_ratios = zeros(120)
	beta_l_ratios = zeros(120)
	for i=1:120
		## unpack parameters
		theta_dict = get_params_fn(theta_boot[i,:], ESTIMATION_PARAMS["param_factors"])

		beta_e_ratios[i] = theta_dict["beta_e"] / theta_dict["alpha"]
		beta_l_ratios[i] = theta_dict["beta_l"] / theta_dict["alpha"]
	end

	ci_levels = [2.5, 97.5]
	beta_e_ratio_ci = percentile(beta_e_ratios, ci_levels)
	beta_l_ratio_ci = percentile(beta_l_ratios, ci_levels)

## Display and save
	println("beta E ratio ", beta_e_ratio, " 95% CI ", beta_e_ratio_ci)
	println("beta L ratio ", beta_l_ratio, " 95% CI ", beta_l_ratio_ci)

	pval_equality = mean(beta_e_ratios .> beta_l_ratios)

	mytext = ""
	mytext *= "beta E ratio " * string(beta_e_ratio) * " 95% CI " * string(beta_e_ratio_ci) * "\n"
	mytext *= "beta L ratio " * string(beta_l_ratio) * " 95% CI " * string(beta_l_ratio_ci) * "\n"
	mytext *= "fraction of time when beta_E > beta_L = " * string(pval_equality) * "\n"

	text_results_path = rootpath_base * "paper/text_results/model_ratio_boot_CI.txt"
	file = open(text_results_path, "w")
	write(file, mytext)
	close(file)