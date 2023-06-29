
include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using Distributed

run_type = 1
n_procs = 12

if (n_procs > 1) && (length(workers()) < n_procs)
	rmprocs(workers())
	display(workers())
	addprocs(n_procs)
	display(workers())

	@everywhere using Pkg
	@everywhere Pkg.activate($env_path)
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
	using InlineStrings

	println("Using statements completed!...")
	    include("../sims-functions-loadprep.jl")
		include("../sims-functions-model.jl")
		include("../sims-functions-rt.jl")
		include("../functions-table.jl")

		include("../sims-wrapper-v1.jl")

	println("Import other functions completed!...")

	Random.seed!(1212376)
end

# DEBUG:
if run_type == 1
	using BenchmarkTools
	# using TableView
	using Plots
end


## Read estimated parameters
	# paths for input data and
	if run_type == 1
		rootpath_input = rootpath_base * "data/coded_model/"

		rootpath_input_tech = rootpath_input * "road_tech/"

		### Demand estimation params
		rootpath_params      =  rootpath_base * "model/estimation_results/results_delta_90/"
		rootpath_params_boot = [rootpath_base * "model/estimation_results/results_boot_delta_90/"]

		rootpath_output = rootpath_base * "model/simulation_results/main/"

	end

### Load requirements for simulations
	_, BIGMAT, delays0, rt_params, delta_t, agents_mean_km_small = load_all(
				rootpath_input=rootpath_input,
				rootpath_input_tech=rootpath_input_tech,
				rootpath_params=rootpath_params,
				rootpath_params_boot=rootpath_params_boot,
				rootpath_output=rootpath_output,
				rt_type="linear",
				theta_factors_dict=nothing,
				boot_idx=nothing,
				bootstrap_samples=nothing,
				proptowage=false,
				save_params_to_file=false)


### Load main simulation results
	sim_results_path = rootpath_output * "so_eqm.bson"
	sim_results = BSON.load(sim_results_path)

	ttimes = copy(sim_results[:nash_ttimes])

	agents = sim_results[:agents]
	agents_mean_km_small = agents.mean_km[1:10:3040]

	# REFRESH
	t_exact_mat 	  = agents.ha .- Matrix(transpose(BIGMAT["hdgrid"]))
	log_t_exact_mat = log.(max.(0, t_exact_mat))

	BIGMAT["t_exact_mat"] = t_exact_mat
	BIGMAT["log_t_exact_mat"] = log_t_exact_mat

	BIGMAT["t_exact_mat_flip"] = Matrix(transpose(BIGMAT["t_exact_mat"]))
	BIGMAT["log_t_exact_mat_flip"] = Matrix(transpose(BIGMAT["log_t_exact_mat"]))

	# zero_charges = sim_results[:so_charges] .* 0.0
#= PLAN
1. compute msc_partial
2. prepare list of idx1, idx2
3. for each (idx1, idx2) compute the Nash
=#

## 0. everywhere
n_agents, n_hd = size(sim_results[:so_choice_probs])
@everywhere begin

	n_agents = $n_agents
	n_hd = $n_hd
	zero_dt_charges = zeros(n_agents, n_hd)
	agents, BIGMAT, delays0, rt_params, delta_t, agents_mean_km_small =
		$agents, $BIGMAT, $delays0, $rt_params, $delta_t, $agents_mean_km_small

	step_tol = 1e-7
	step_max = 10000
end

## Compute nash again
nash_dt_choice, nash_eu = nash_iteration(
		agents=agents,
		rt_params=rt_params,
		delays0=delays0,
		dt_charges=zero_dt_charges,
		BIGMAT=BIGMAT,
		step_tol=step_tol,
		step_max=step_max,
		error_if_fail=true,
		prop_update_frac=0.1,
		ttimes_from_delay=true,
		rootpath_output="",
		compute_logsum=1,
		debug_level=1
	)

## compare
	mydiff = nash_dt_choice .- sim_results[:nash_choice_probs]

## 1. Compute PARTIAL msc
	### Initialize BIGMAT["ttimes"]
	densities, instant_delays =
	make_travel_times_density(
		choice_prbs=nash_dt_choice, #sim_results[:nash_choice_probs],  # N_agents x n_hd
		rt_params=rt_params,
		km_mean=agents.mean_km[1:10:3040], #agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes"],
		delta_t=delta_t
	)

	hdgrid = BIGMAT["hdgrid"]
	mydiff = ttimes - BIGMAT["ttimes"]
	plot(hdgrid,  ttimes[1,:])
	plot!(hdgrid, BIGMAT["ttimes"][1,:])

	msc_indiv_partial = msc_partial_all(
			agents=agents,
			rt_params=rt_params,
			BIGMAT=BIGMAT,
			choice_ne=nash_dt_choice, #sim_results[:nash_choice_probs],
			my_dt_charges=zero_dt_charges,
			delta_t=delta_t,
			run_parallel=true,
			debug=true
		)


## 2. Prepare EQUILIBRIUM msc

	## get all the i,j pairs needed (much less than n_hd^2)
	make_travel_times_density(
		choice_prbs=sim_results[:nash_choice_probs],  # N_agents x n_hd
		rt_params=rt_params,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes"],
		delta_t=delta_t
	)

	# matrix
	jk_msc = zeros(Bool, n_hd, n_hd)
	for i=1:n_agents
		for j=1:n_hd
			k = j + Int64(floor(BIGMAT["ttimes"][i,j] / delta_t))
			k = min(109, k)
			jk_msc[j,k] = true
		end
	end

	# list
	# sum(jk_msc)
	jk_msc_list = []
	jk_msc_list_start = Int64[]
	jk_msc_list_end = Int64[]
	for j=1:n_hd
		for k=j:n_hd
			if jk_msc[j,k]
				push!(jk_msc_list, (j,k))
				push!(jk_msc_list_start, j)
				push!(jk_msc_list_end, k)
			end
		end
	end

	### SAMPLE
	mysampleij = randperm(length(jk_msc_list))[1:100]
	jk_msc_list = jk_msc_list[mysampleij]
	jk_msc_list_start = jk_msc_list_start[mysampleij]
	jk_msc_list_end = jk_msc_list_end[mysampleij]

	msc_FULL_loaded =
		(myidx1, myidx2) -> nash_iteration(
				agents=agents,
				rt_params=rt_params,
				delays0=delays0,
				dt_charges=sim_results[:so_charges] .* 0.0,
				idx1=myidx1, idx2=myidx2,
				BIGMAT=BIGMAT,
				step_tol=step_tol,
				step_max=step_max,
				error_if_fail=true,
				prop_update_frac=0.1,
				ttimes_from_delay=false,
				rootpath_output="",
				compute_logsum=2,
				debug_level=1
			)


	# Run
	mscs = pmap(
		idx -> msc_FULL_loaded(
				jk_msc_list_start[idx],
				jk_msc_list_end[idx]),
		1:length(jk_msc_list_start))

	# For comparison: welfare without
	w0 = nash_iteration(
			agents=agents,
			rt_params=rt_params,
			delays0=delays0,
			dt_charges=sim_results[:so_charges] .* 0.0,
			# idx1=40, idx2=50,
			BIGMAT=BIGMAT,
			step_tol=step_tol,
			step_max=step_max,
			error_if_fail=true,
			prop_update_frac=0.1,
			ttimes_from_delay=false,
			rootpath_output="",
			compute_logsum=2,
			debug_level=1
		)

	mscs_scaled = (mscs .- w0) * (n_agents - 1) * 1000.0


## Refresh original travel times
	make_travel_times_density(
		choice_prbs=sim_results[:nash_choice_probs],  # N_agents x n_hd
		rt_params=rt_params,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes"],
		delta_t=delta_t
	)

	# msc_indiv = zeros(Float64, n_agents, n_hd)
	mscs_scaled_fulleqm = Float64[]
	mscs_scaled_partial = Float64[]
	for i=1:n_agents
		for j=1:n_hd
			k = j + Int64(floor(BIGMAT["ttimes"][i,j] / delta_t))
			k = min(109, k)

			if (j,k) in jk_msc_list
				push!(mscs_scaled_partial, msc_indiv_partial[i, j])

				idx = 1
				while jk_msc_list[idx] != (j,k)
					idx += 1
				end
				@assert jk_msc_list[idx] == (j,k)
				push!(mscs_scaled_fulleqm, mscs_scaled[idx])
			end
			# msc_indiv[i,j] = msc_square[j,k]
			# if msc_square[j,k] == -1.0
			# 	println(i,",",j,",",k)
			# end
		end
	end


## Compare
	df = DataFrame(
		"partial" => mscs_scaled_partial,
		"fulleqm" => mscs_scaled_fulleqm)

	using FixedEffectModels
	reg(df, @formula(fulleqm ~ partial)) |> display

	scatter(mscs_scaled_partial, mscs_scaled_fulleqm,
			title="Marginal social cost of add'l trip around Nash Eqm",
			xlabel="Partial-equilibrium",
			ylabel="Full equilibrium") |> display

	describe(mscs_scaled_fulleqm)
	describe(mscs_scaled_partial)
