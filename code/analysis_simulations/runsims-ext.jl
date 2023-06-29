
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
		include("sims-functions-model-ext.jl")
		include("sims-functions-rt.jl")
		include("functions-table.jl")

		include("sims-wrapper-v1.jl")

	println("Import other functions completed!...")

	Random.seed!(1212376)
	my_random_seed = 1212376
end

function get_elasticities(;
			theta_add_dict,
			debug_level=1)

	### Load requirements for simulations
	agents, BIGMAT, delays0, rt_params, delta_t, agents_mean_km_small = load_all(
					rootpath_input=rootpath_input,
					rootpath_input_tech=rootpath_input_tech,
					rootpath_params=rootpath_params,
					rootpath_params_boot=rootpath_params_boot,
					rootpath_output=rootpath_output_base,
					theta_add_dict=theta_add_dict,
					rt_type="linear",
					random_seed=my_random_seed, # use the same random number generator in order to always get the same hastar draws
					proptowage=false,
					save_params_to_file=false)

	##
		dt_choice, trip_prob, logsum = nash_iteration_ext(
					agents=agents,
					rt_params=rt_params,
					delays0=delays0,
					ttimes_from_delay=true,  # only first time
					dt_charges=0.0 * BIGMAT["ttimes"],
					BIGMAT=BIGMAT,
					compute_logsum=1,
					debug_level=debug_level,
					error_if_fail=true,
					prop_update_frac=0.1,
					step_tol=1e-7,
					step_max=1000,
					debug_graph_period=10)

	## Measure impact of flat charges
		charge_levels = (0.0:100.0:300.0) |> Array |> vec
		nc = length(charge_levels)

		avg_logsum = zeros(nc)
		avg_trip_prob = zeros(nc)
		for (idx_c, charge_level) in enumerate(charge_levels)
			println("charge", charge_level)
			println("idx", idx_c)
			dt_choice, trip_prob, logsum = make_choice_prbs_ext(
					agents_mean_km=agents.mean_km,
					agents_a=agents.a,
					agents_be=agents.be,
					agents_bl=agents.bl,
					agents_sig=agents.sig,
					agents_delta_ext=agents.delta_ext,
					agents_eta=agents.eta,
					ttimes=BIGMAT["ttimes"],
					dt_charges=-charge_level .+ 0.0 .* BIGMAT["ttimes"],
			        t_exact_mat=BIGMAT["t_exact_mat_flip"],
			        log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
					compute_logsum=1)

			avg_logsum[idx_c] = mean(logsum)
	        avg_trip_prob[idx_c] = mean(trip_prob)
		end

		logsum_0 = -avg_logsum[1]
		log_charge_change = log.(logsum_0 .- charge_levels) .- log.(logsum_0)
		log_avg_trip_prob = log.(avg_trip_prob)

		# elasticities
		level_ch_100 = avg_trip_prob[2] - avg_trip_prob[1]
	    level_ch_200 = avg_trip_prob[3] - avg_trip_prob[1]
	    level_ch_300 = avg_trip_prob[4] - avg_trip_prob[1]
		elastic_100 = (log_avg_trip_prob[2] - log_avg_trip_prob[1]) / (log_charge_change[2] - log_charge_change[1])
	    elastic_200 = (log_avg_trip_prob[3] - log_avg_trip_prob[1]) / (log_charge_change[3] - log_charge_change[1])
	    elastic_300 = (log_avg_trip_prob[4] - log_avg_trip_prob[1]) / (log_charge_change[4] - log_charge_change[1])

		println("elasticity @100  dlog(trip_prob) / dlog(total_cost) :", elastic_100)
	    println("elasticity @200  dlog(trip_prob) / dlog(total_cost) :", elastic_200)
	    println("elasticity @300  dlog(trip_prob) / dlog(total_cost) :", elastic_300)
	    println("pp change in trip prob:       pi(Ch=0) - pi(Ch=100) :", level_ch_100)
	    println("pp change in trip prob:       pi(Ch=0) - pi(Ch=200) :", level_ch_200)
	    println("pp change in trip prob:       pi(Ch=0) - pi(Ch=300) :", level_ch_300)

		# plot(charge_levels, avg_trip_prob)

		return logsum_0, avg_trip_prob,
			level_ch_100, level_ch_200, level_ch_300,
			elastic_100, elastic_200, elastic_300
end

# for DEBUG
if run_type == 1
	using BenchmarkTools
	using TableView
	using Plots
end

## Read estimated parameters
	### Model input data 
	rootpath_input  	= rootpath * "data/coded_model/"
	rootpath_input_tech = rootpath * "data/coded_model/road_tech/"

	### Demand estimation params
	rootpath_params      =  rootpath * "model/estimation_results/results_delta_90/"
	rootpath_params_boot = [rootpath * "model/estimation_results/results_boot_delta_90/"]

	### Simulation output path
	rootpath_output_base = rootpath * "model/simulation_results/main_ext/"
	isdir(rootpath_output_base) || mkdir(rootpath_output_base)

### prepare extensive margin parameters
	# theta_add_dict = Dict{String, Float64}(
	# 	"delta_ext" => -1100.0, # penalty for no trip
	# 	"eta" => 20.0 # nested logit parameter
	# )

	debug_level = 1


for (idx_eta, eta) in enumerate([20.0, 100.0])
	for (idx_delta, delta) in enumerate(600.0:200.0:1000.0)

		theta_add_dict = Dict{String, Float64}(
			"delta_ext" => -delta, # penalty for no trip
			"eta" => eta
		)

		# Paths
		rootpath_output = rootpath_output_base *
						 "eta_" * string(Int64(eta)) * "_" *
						"delta_" * string(Int64(delta)) * "/"
		isdir(rootpath_output) || mkdir(rootpath_output)

		### Load requirements for simulations
		agents, BIGMAT, delays0, rt_params, delta_t, agents_mean_km_small = load_all(
				rootpath_input=rootpath_input,
				rootpath_input_tech=rootpath_input_tech,
				rootpath_params=rootpath_params,
				rootpath_params_boot=rootpath_params_boot,
				rootpath_output=rootpath_output,
				theta_add_dict=theta_add_dict,
				random_seed=my_random_seed, # use the same random number generator in order to always get the same hastar draws
				rt_type="linear",
				proptowage=false,
				save_params_to_file=false)

		## Get the implied elasticities from a 100/200/300 Rs flat charge
		elasticity_results = get_elasticities(theta_add_dict=theta_add_dict, debug_level=1)
		logsum_0, avg_trip_prob,
		level_ch_100, level_ch_200, level_ch_300,
		elastic_100, elastic_200, elastic_300 = elasticity_results

		# write to file
		elasticity_results_path = rootpath_output * "elasticity_results.bson"
		BSON.@save elasticity_results_path logsum_0 avg_trip_prob level_ch_100 level_ch_200 level_ch_300 elastic_100 elastic_200 elastic_300

		## Nash
		dt_choice, trip_prob, logsum = nash_iteration_ext(
				agents=agents,
				rt_params=rt_params,
				delays0=delays0,
				ttimes_from_delay=true,  # only first time
				dt_charges=0.0 * BIGMAT["ttimes"],
				BIGMAT=BIGMAT,
				compute_logsum=1,
				debug_level=debug_level,
				error_if_fail=true,
				prop_update_frac=0.1,
				step_tol=1e-7,
				step_max=1000,
				debug_graph_period=10)

		## Social Optimum
		socopt_iteration_ext(
				agents=agents,
				rt_params=rt_params,
				delays0=delays0,
				BIGMAT=BIGMAT,
				step_tol_out=1e-2,
				step_tol_in=1e-7,
				step_max=500,
				prop_update_frac=0.1,
				adaptive_adjust=true,
				adjust_factor_initial=0.5,
				msc_run_parallel=true,
				rootpath_output=rootpath_output,
				debug_level=debug_level)

	end
end

# test - WORKS
# mypath = rootpath_output_base * "eta_20_delta_600/elasticity_results.bson"
# BSON.@load mypath logsum_0 avg_trip_prob level_ch_100 level_ch_200 level_ch_300 elastic_100 elastic_200 elastic_300
