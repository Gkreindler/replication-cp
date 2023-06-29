
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
	rootpath_output = rootpath * "model/simulation_results/main_boot/"
	isdir(rootpath_output) || mkdir(rootpath_output)
    isdir(rootpath_output * "boot_results/") || mkdir(rootpath_output * "boot_results/")

## Read boot samples
	all_bootstrap_samples = []
	all_boot_rootpath_outputs = []
	for i=1:120
		print(".")

		## Load bootstrap samples
			boot_folder = rootpath_params_boot[1] * "boot_run_" * string(i)* "/"
			bootstrap_samples = Dict{String, Vector{Int64}}()
			## Load boot samples to file
			path = string(boot_folder,"boot_sample_a.csv")
			bootstrap_samples["a"] = readdlm(path, ',', Float64) |> vec
			# display(bootstrap_samples["a"])

			path = string(boot_folder,"boot_sample_na.csv")
			bootstrap_samples["na"] = readdlm(path, ',', Float64) |> vec

		push!(all_bootstrap_samples, bootstrap_samples)
		push!(all_boot_rootpath_outputs, rootpath_output * "boot_" * string(i) * "/")
	end


## Parameter modifications
@everywhere begin
	theta_factors_dict = Dict(
		"beta_e" => 1.0
	)

	all_bootstrap_samples = $all_bootstrap_samples

	rootpath_input = $rootpath_input
	rootpath_input_tech = $rootpath_input_tech
	rootpath_params = $rootpath_params
	rootpath_params_boot = $rootpath_params_boot
	rootpath_output = $rootpath_output
end

## Run!
	# boot_idx = 11
	# bootstrap_samples = all_bootstrap_samples[boot_idx]

for boot_idx=1:120

	println("\n\n\n ====================================")
	println("BOOT IDX ", boot_idx)
	println("====================================\n\n\n")

	if isdir(all_boot_rootpath_outputs[boot_idx]) &&
	   isfile(all_boot_rootpath_outputs[boot_idx] * "so_eqm.bson")
	   println("SKIPPING")
   else
	   run_sim(rootpath_input=rootpath_input,
		rootpath_input_tech=rootpath_input_tech,
		rootpath_params=rootpath_params,
		rootpath_params_boot=rootpath_params_boot,
		rootpath_output=all_boot_rootpath_outputs[boot_idx],
		theta_factors_dict=theta_factors_dict,
		boot_idx=boot_idx,
		bootstrap_samples=all_bootstrap_samples[boot_idx],
		rt_type="linear",
		proptowage=false,   # <------------- het
        prop_update_frac=0.1,
		adaptive_adjust=true,
		adjust_factor_initial=0.5,
		step_tol=1e-7,
		step_max=10000,
		msc_run_parallel=true,
		debug_level=1)
	end
end
