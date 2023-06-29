
include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package

using Distributed

run_type = 3
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
	n_procs = 2
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
		include("sims-functions-rt.jl")
		include("functions-table.jl")

		include("sims-wrapper-v1.jl")

	println("Import other functions completed!...")

	Random.seed!(1212376)
end

# DEBUG:
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
	rootpath_output_base = rootpath * "model/simulation_results/varyparams/"
	isdir(rootpath_output_base) || mkdir(rootpath_output_base)
    

@everywhere begin
	theta_factors_dict = Dict(
		"beta_e" => 1.0
	)

	rootpath_input = $rootpath_input
	rootpath_input_tech = $rootpath_input_tech
	rootpath_params = $rootpath_params
	rootpath_params_boot = $rootpath_params_boot
	rootpath_output_base = $rootpath_output_base
end

# OLD
# pref_factors       = [0.25, 0.5, 1.0, 2.0, 4.0]
# pref_factors_names = ["0_25", "0_50", "1_00", "2_00", "4_00"]
# pref_vars         = [["beta_e", "beta_l"], ["beta_e"], ["alpha"]]
# pref_vars_names   = ["both_betas", "beta_e", "alpha"]
pref_factors       = [0.25, 4.0]
pref_factors_names = ["0_25", "4_00"]
pref_vars         = [["beta_e"], ["alpha"]]
pref_vars_names   = ["beta_e", "alpha"]


all_run_params = []

rt_type = "linear"

for idx_pref=1:2
	pref_var_name = pref_vars_names[idx_pref]
	pref_var      = pref_vars[idx_pref]

	for idx_factor=1:2
		pref_factor      = pref_factors[idx_factor]
		pref_factor_name = pref_factors_names[idx_factor]

		run_params = (pref_var, pref_var_name, pref_factor, pref_factor_name, rt_type)
		push!(all_run_params, run_params)
	end
end

display(all_run_params)

# all_run_params = all_run_params[end:end]

for run_params in all_run_params
	pref_var, pref_var_name, pref_factor, pref_factor_name, rt_type = run_params

	### prepare the dictionary
	theta_factors_dict = Dict{String, Float64}()
	for mykey in pref_var
		theta_factors_dict[mykey] = pref_factor
	end

	### prep output folder
	rootpath_output = rootpath_output_base * "sim_vary_" * rt_type * "_" * pref_var_name * "_" * pref_factor_name * "/"
	isdir(rootpath_output) || mkdir(rootpath_output)

	println("\n\n\n\n Processing ", run_params, "\n\n\n")

	### run!
	run_sim(rootpath_input=rootpath_input,
			rootpath_input_tech=rootpath_input_tech,
			rootpath_params=rootpath_params,
			rootpath_params_boot=rootpath_params_boot,
			rootpath_output=rootpath_output,
			theta_factors_dict=theta_factors_dict,
			rt_type=rt_type,
			proptowage=false,
			prop_update_frac=0.1,
			adaptive_adjust=true,
			adjust_factor_initial=0.5,
			step_tol=1e-7,
			step_max=250,
			msc_run_parallel=true,
			debug_level=1)

end
