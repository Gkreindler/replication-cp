
include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package

using Distributed

run_type = 1
# 1 = local            -> few cores (set n_procs below)
# 2 = server testing   -> few cores (set n_procs below)
# 3 = server full run  -> many cores (set n_procs below)

if run_type == 1
	rootpath = rootpath_local
	
	# Julia environment path (holds all package versions, etc.)
	env_path = rootpath * "code/analysis_estimation/cp_env/"
else
    rootpath = rootpath_server
	env_path = rootpath * "code/analysis_estimation/cp_env_server/"
end

if run_type < 3
	n_procs = 10
else
	n_procs = 43
end

using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

if (n_procs > 1) && (length(workers()) < n_procs)
	rmprocs(workers())
	display(workers())
	addprocs(n_procs)
	display(workers())

	# activate environment on all workers
    @everywhere using Pkg
	@everywhere Pkg.activate($env_path) # the dollar sign takes the variable from the local (main) worker
end

@everywhere begin
	using Future
	using CSV
	using DelimitedFiles
	using JSON
	using BSON

	using DataFrames
	using StatsBase
	using Statistics
	using LinearAlgebra
	using Random
	using Distributions
	using FreqTables
	using Loess
	using Interpolations
	using Dates
	using NonNegLeastSquares

	using SpecialFunctions

	println("Using statements completed!...")
	    include("sims-functions-loadprep.jl")
		include("sims-functions-model.jl")
		include("sims-functions-model-2routes.jl")
		include("sims-functions-rt.jl")
		include("functions-table.jl")

		include("sims-wrapper-v1.jl")

	println("Import other functions completed!...")

	Random.seed!(1212376)
	my_random_seed = 1212376
end

# for DEBUG
if run_type == 1
	using BenchmarkTools
	using TableView
	using Plots
end

## Paths
    ### Model input data 
    rootpath_input      = rootpath * "data/coded_model/"
    rootpath_input_tech = rootpath * "data/coded_model/road_tech/"

	### Demand estimation params
    rootpath_params      =  rootpath * "model/estimation_results/results_delta_90/"
    rootpath_params_boot = [rootpath * "model/estimation_results/results_boot_delta_90/"]

    ### Simulation output path
	rootpath_output_base = rootpath * "model/simulation_results/main_2routes/"
	isdir(rootpath_output_base) || mkdir(rootpath_output_base)

### Load requirements for simulations

	# benchmark (v1 and v2)
	theta_factors_dict = nothing

	# v3 -- only switch on for low_mu
	theta_factors_dict = Dict(
		"mu" => 0.25
	)


	agents, BIGMAT, delays0, rt_params, delta_t, agents_mean_km_small = load_all(
				rootpath_input=rootpath_input,
				rootpath_input_tech=rootpath_input_tech,
				rootpath_params=rootpath_params,
				rootpath_params_boot=rootpath_params_boot,
				rootpath_output=rootpath_output_base,
				theta_factors_dict=theta_factors_dict,
				rt_type="linear",
				random_seed=my_random_seed, # use the same random number generator in order to always get the same hastar draws
				proptowage=false,
				save_params_to_file=false)

##
	delays1 = copy(delays0)
	rt0_params = copy(rt_params)
	rt1_params = copy(rt_params)

	## v1
	rootpath_output_base *= "v15/"
	rt1_params[1] *= 1 - (0.15/2)
	rt1_params[2] *= 1.15

	## v2
	# rootpath_output_base *= "v30/"
	# rt1_params[1] *= 1 - (0.30/2)
	# rt1_params[2] *= 1.30

	## v3 -- make sure to also switch on dictionary above!
	# rootpath_output_base *= "v15_low_mu/"
	# rt1_params[1] *= 1 - (0.15/2)
	# rt1_params[2] *= 1.15

	# update for TWO routes -> density counts double
	rt0_params[3] = rt_params[3] * 2.0
	rt1_params[3] = rt_params[3] * 2.0

	BIGMAT["ttimes0"] = copy(BIGMAT["ttimes"])
	BIGMAT["ttimes1"] = copy(BIGMAT["ttimes"])

	isdir(rootpath_output_base) || mkdir(rootpath_output_base)

	nash_iteration_2routes(
				agents=agents,
				rt0_params=rt0_params,
				rt1_params=rt1_params,
				delays0=delays0,
				delays1=delays1,
				ttimes_from_delay=true,  # only first time
				dt_charges0=0.0 * BIGMAT["ttimes0"],
				dt_charges1=0.0 * BIGMAT["ttimes0"],
				BIGMAT=BIGMAT,
				compute_logsum=1,
				debug_level=2,
				error_if_fail=true,
				prop_update_frac=0.05,
				step_tol=1e-7,
				step_max=1000,
				debug_graph_period=50)

	socopt_iteration_2routes(;
				agents=agents,
				rt0_params=rt0_params,
				rt1_params=rt1_params,
				delays0=delays0,
				delays1=delays1,
				BIGMAT=BIGMAT,
				step_tol_out=1e-2,
				step_tol_in=1e-7,
				step_max=500,
				prop_update_frac=0.05,
				adaptive_adjust=true,
				adjust_factor_initial=0.5,
				msc_run_parallel=true,
				rootpath_output=rootpath_output_base,
				debug_level=2)
