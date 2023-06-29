
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
	rootpath_output = rootpath * "model/simulation_results/factors/"
	isdir(rootpath_output) || mkdir(rootpath_output)


@everywhere begin
	theta_factors_dict = Dict(
		"beta_e" => 1.0
	)

	rootpath_input = $rootpath_input
	rootpath_input_tech = $rootpath_input_tech
	rootpath_params = $rootpath_params
	rootpath_params_boot = $rootpath_params_boot
	rootpath_output = $rootpath_output
end

rt_factor_range = 2.0:10.0

rootpath_output_folders = []
for rt_factor in rt_factor_range
	rootpath_output_subfolder = rootpath_output * "factor_" * string(Int64(rt_factor)) * "/"
	push!(rootpath_output_folders, rootpath_output_subfolder)
	isdir(rootpath_output_subfolder) || mkdir(rootpath_output_subfolder)
end

for idx=1:length(rt_factor_range)

	println("\n\n\n ====================================")
	println("FACTOR IDX ", idx)
	println("====================================\n\n\n")

	if isdir(rootpath_output_folders[idx]) &&
	   isfile(rootpath_output_folders[idx] * "so_eqm.bson")
	   println("SKIPPING")
   else
	   run_sim(rootpath_input=rootpath_input,
			rootpath_input_tech=rootpath_input_tech,
			rootpath_params=rootpath_params,
			rootpath_params_boot=rootpath_params_boot,
			rootpath_output=rootpath_output_folders[idx],  ### IDX
			theta_factors_dict=theta_factors_dict,
			rt_type="linear",
			rt_factor=rt_factor_range[idx],  			   ### IDX
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
