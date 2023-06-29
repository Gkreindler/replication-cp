
if run_type == 1 # local folders
	env_path = rootpath_local * "code/analysis_estimation/cp_env/"

	else # server folders
	env_path = rootpath_server * "code/analysis_estimation/cp_env_server/"
end

# local
using Pkg
Pkg.activate(env_path) # the dollar sign takes the variable from the local (main) worker
Pkg.instantiate() # install all packages if necessary
Pkg.precompile() # install all packages if necessary

if n_procs > 1
	rmprocs(workers())
	display(workers())
	addprocs(n_procs)
	display(workers())

	@everywhere using Pkg
	@everywhere Pkg.activate($env_path) # the dollar sign takes the variable from the local (main) worker
end


@everywhere begin
	# n_start_pts = $n_start_pts
	# n_bootstrap_start_pts = $n_bootstrap_start_pts

	# if necessary, run this in julia
	# (v1.7) ] activate .
	# (cp_env) ] add Future, CSV, JSON, DataFrames, Optim, LsqFit, FiniteDifferences, LineSearches, StatsBase, Statistics, LinearAlgebra, Random, Distributions, FreqTables, Loess, Interpolations, Dates, NonNegLeastSquares
	# (cp_env) ] rm LsqFit
	# (cp_env) ] add https://github.com/JuliaNLSolvers/LsqFit.jl#master # as of June 2023 need the master branch for the maxTime option. can also do "add LsqFit#master"

	using Future

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

	# using TableView
	using FreqTables
	using Loess
	using Interpolations
	using Dates
	# using BenchmarkTools
	using NonNegLeastSquares

	println("Using statements completed!...")

	# include function to calculate the model
	include("functions-model-definitions.jl")
	include("functions-loadprep.jl")
	include("functions-model.jl")


	include("gmm_wrappers.jl")
	include("functions-gmm-full.jl")
	include("functions-gmm-bootstrap.jl")

	println("Import other functions completed!...")

	Random.seed!(1212376)
end

if run_type == 1 # local folders
	rootpath = rootpath_local
else # server folders
	rootpath = rootpath_server
end

## Define paths and create new folders (if needed)

	@everywhere function mkdir_if_needed(myrootpath)
		isdir(myrootpath) || mkdir(myrootpath)
	end

	@everywhere function prep_dirs(gmm_options, rootpath)
		rootpath_input  = rootpath * "data/coded_model/"
		rootpath_output     = rootpath * "model/estimation_results/" * gmm_options["local_subfolder"] * "/results"      * gmm_options["results_suffix"] * "/"
		rootpath_boot_output = rootpath * "model/estimation_results/" * gmm_options["local_subfolder"] * "/results_boot" * gmm_options["results_suffix"] * "/" 
		
		## create output folders (if necessary)
		gmm_options["run_main"] && mkdir_if_needed(rootpath_output)
		gmm_options["run_boot"] && mkdir_if_needed(rootpath_boot_output)
	
		## update dictionary
		gmm_options["rootpath_input"] = rootpath_input
		gmm_options["rootpath_output"] = rootpath_output
		gmm_options["rootpath_boot_output"] = rootpath_boot_output
	
	end

