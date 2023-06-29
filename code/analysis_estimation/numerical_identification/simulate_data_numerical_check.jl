#=
This script creates synthetic data sets from random parameter values.
We'll then run GMM on this synthetic data to see if it recovers the underlying
parameters that were used to generate the data.

This has two flavors:
1) Finite sample properties. I generate a data set of exactly the same size as the
original data. Given randomness, we cannot expect GMM to recover the underlying
parameters, or even for this to happen on average (unbiased). However, this is
informative about the size of the finite sample bias.
2) Asymptotic. Here we replace actions with average:
	- departure time pdf
	- route choice probabilities

Steps:
		1. 	Load true estimated parameters

		2. 	Load ESTIMATED parameters

		3. 	Simulate model (at ESTIMATED parameters) to recover the distribution
			of ideal arrival times

		4.   generate syntetic data (iterate 120 times)
		4.1. generate random NEW params
		4.2. simulate the model again at the NEW parameters (using stored hastar distribution from prev step)
		4.3. generate synthetic data for departure time (finite and asymptotic)
		4.4. generate synthetic data for route choice   (finite and asymptotic)
		4.5. save to files
=#

## part 0 - prep
	include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
	env_path = rootpath_base * "code/analysis_estimation/cp_env/"
	using Pkg
	Pkg.activate(env_path) # the dollar sign takes the variable from the local (main) worker
	Pkg.instantiate()

using Future
using CSV
using JSON
using BSON

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
using FreqTables
using Loess
using Interpolations
using Dates
# using BenchmarkTools
using NonNegLeastSquares

println("Using statements completed!...")

# include function to calculate the model
include("../functions-model-definitions.jl")
include("../functions-loadprep.jl")
include("../functions-model.jl")
include("../gmm_wrappers.jl")
include("../functions-gmm-full.jl")
include("../functions-gmm-bootstrap.jl")

include("../functions-table.jl")

println("Import other functions completed!...")

Random.seed!(1212376)


## part 0 - paths
	rootpath_input  = rootpath_base * "data/coded_model/"
	# rootpath_output = string(rootpath, "data/coded_model/id_check/")
	rootpath_output_fin = rootpath_base * "data/coded_model/simdata_finite/"
	isdir(rootpath_output_fin) || mkdir(rootpath_output_fin)

	rootpath_output_asy = rootpath_base * "data/coded_model/simdata_asymptotic/"
	isdir(rootpath_output_asy) || mkdir(rootpath_output_asy)

## part 1 - Load data

	# load all data from original location
    DATA = load_and_prep(rootpath_input)

	# "loaded" versions of these dataframes. Pretty sure sorted the same as the
	# originals but just to be sure, let's use these
	data_a_loaded = DATA["DATAMOMS"]["data_a"]
	data_all_loaded = DATA["DATAMOMS"]["data_all"]

	# laod raw dataframes -- area treatment participants
    filepath = string(rootpath_input, "coded_uidpn_a.csv")
    data_a = CSV.read(filepath, DataFrame)
    data_a[!, :area] .= 1

	# laod raw dataframes -- non area treatment participants
    filepath = string(rootpath_input, "coded_uidpn_na.csv")
    data_na = CSV.read(filepath, DataFrame)
    data_na[!, :area] .= 0

	# laod raw dataframes -- area treatment participants route choice
    filepath = string(rootpath_input, "area_trip_vot.csv")
    atdf = CSV.read(filepath, DataFrame)

    println(" ...All dataframes loaded!...")

	## Convert variables and remove missing vars
	atdf[!, :is_long_route] = convert.(Float64, atdf.is_long_route)

	## Remove missing and convert
	for j=1:61
		# A
		colname = string("mean_dt_pre", j)
		replace!(data_a[!, colname], missing => 0)
		data_a[!, colname] = convert.(Float64, data_a[!, colname])

		colname = string("mean_dt_pos", j)
		replace!(data_a[!, colname], missing => 0)
		data_a[!, colname] = convert.(Float64, data_a[!, colname])

		# NA
		colname = string("mean_dt_pre", j)
		replace!(data_na[!, colname], missing => 0)
		data_na[!, colname] = convert.(Float64, data_na[!, colname])

		colname = string("mean_dt_pos", j)
		replace!(data_na[!, colname], missing => 0)
		data_na[!, colname] = convert.(Float64, data_na[!, colname])
	end

	## Copies for "finite sample" versions (use original for "asymptotic")
	data_a_finite 	= copy(data_a)
	data_na_finite 	= copy(data_na)
	atdf_finite 	= copy(atdf)

	## Other coding (for sampling)
	n_hd = 79 # 61
	hdgrid = vec(Array(range(-2.5,stop=4.0,length=n_hd)))

	bin_l = -2.5
	bin_h = 2.5
	n_hd_bins = 61

	bin_delta = (bin_h-bin_l)/(n_hd_bins - 1)
	hd_bins = range(bin_l - bin_delta/2, stop=bin_h + bin_delta/2, length=n_hd_bins + 1)
	hd_bins = Array(hd_bins)

	hdgrid_2525 = vec(Array(range(-2.5,stop=2.5,length=61)))



## part 2 - load ESTIMATED parameters
    rootpath_results      =  string(rootpath_base, "model/estimation_results/results_delta_90/")
    rootpath_boot_results = [string(rootpath_base, "model/estimation_results/results_boot_delta_90/")]
    param_names = ["alpha", "delta", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    spec_diagnostics_df |> display

    theta_optimum = select(main2_df, r"^param_") |> Matrix |> vec
	println("getting optimum: ", theta_optimum)

## part 3 - simulate model (at ESTIMATED parameters) to recover the
# 				distribution of ideal arrival times

	# Simulate model at ESTIMATED theta
    # other fixed parameters
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
    gmm_options = Dict(
        "main_n_start_pts" => 2,
        "boot_n_runs" => 2,
        "boot_n_start_pts" => 2
    )

    # run this to get the param_factors matrix
    ESTIMATION_PARAMS = model_dt_dynamic_params(theta_fixed_dict=theta_fixed_dict, gmm_options=gmm_options)

    # simulate model in order to get the hastar distribution (under true parameters)
    solve_dynamicvot(
        theta=theta_optimum,
        param_factors=ESTIMATION_PARAMS["param_factors"],
        BIGMAT=DATA["BIGMAT"],
		CHARGEMAT=DATA["CHARGEMAT"],
		COMPMAT=DATA["COMPMAT"],
		DATAMOMS=DATA["DATAMOMS"],
        hastar_inversion_type="distribution"
        )

	# copy the ideal arrival time distribution
	DATA["COMPMAT"]["hastars_given"] = copy(DATA["COMPMAT"]["hastars"])

	println("done simulating model at estimated parameter. Copied hastar.")

## part 4. generate syntetic data (iterate 120 times)
# 		4.1. generate random NEW params
# 		4.2. simulate the model again at the NEW parameters (using stored hastar distribution from prev step)
# 		4.3. generate synthetic data for departure time (finite and asymptotic)
# 		4.4. generate synthetic data for route choice   (finite and asymptotic)
# 		4.5. save to files

for idx_datasim=1:120
	println("\n\n Iteration ", idx_datasim, "\n\n")

	### 4.1. generate random NEW params

		# multiplicative ok for all parameters except FE which can be negative
		rand_params = theta_optimum .* (0.25 .+ 1.5 * rand(7))

		# randomly flip sign of FE
		rand_params[3] = rand_params[3] * (1 - 2 * (rand() > 0.5))

		println("done generating new param vector: ", rand_params)

	### 4.2. simulate the model again at the NEW parameters (using stored hastar distribution from prev step)
	# reset route model initial conditions (should not matter)
		fill!(DATA["COMPMAT"]["initial_conditions_2d"], 0.0)

	    sim_model = solve_dynamicvot(
	        theta=rand_params,
	        param_factors=ESTIMATION_PARAMS["param_factors"],
	        BIGMAT=DATA["BIGMAT"],
			CHARGEMAT=DATA["CHARGEMAT"],
			COMPMAT=DATA["COMPMAT"],
			DATAMOMS=DATA["DATAMOMS"],
			hastar_inversion_type="given"   # <---------- use the same "true" distribution
	        )

		println("done simulating model for new param vector")


	### 4.3. generate synthetic data for ***departure time***

		# probability densities during and before experiment
		density_model_wcha   = hist_subgrid(sim_model["dt_wcha"], hdgrid, hd_bins)
		density_model_before = hist_subgrid(sim_model["dt_before"], hdgrid, hd_bins)

		# update for first part of the data (_a)
		for i=1:size(data_a, 1)
			uidp = data_a[i, :uidp]

			# get row number for this person -> use this to access model results for this uidp
			idx_row = findfirst(data_a_loaded.uidp .== uidp)

			# normalize
			density_before = density_model_before[idx_row, :]
			density_before = density_before ./ sum(density_before)

			density_wcha = density_model_wcha[idx_row, :]
			density_wcha = density_wcha ./ sum(density_wcha)

			# FINITE sampling
				# actual sample sizes
				n_pre = data_a[i, :n_pre]
				ismissing(n_pre) && (n_pre = 1)
				n_pos = data_a[i, :n_pos]
				ismissing(n_pos) && (n_pos = 1)

				# testing: large sample
				# n_pre = 500
				# n_pos = 500

				# draw N samples from the model distribution
				dep_times_before = sample(hdgrid_2525, Weights(density_before), n_pre)
				dep_times_pos    = sample(hdgrid_2525, Weights(density_wcha), n_pos)

				# normal distribution that approximates pre- distribution
				mymean = mean(dep_times_before)
				mysd = Statistics.std(dep_times_before)
				(n_pre == 1)  && (mysd = .4073334)   # this is the median, as in the original coding
				(mysd == 0.0) && (mysd = .4073334)  # if the random draws are exactly equal
				data_a_finite[i, :norm_mean] = mymean
				data_a_finite[i, :norm_sd] = mysd

				# frequencies pre&post
				for j=1:61
					data_a_finite[i, string("mean_dt_pre", j)] = sum( abs.(dep_times_before .- hdgrid_2525[j]) .< 0.001 )
					data_a_finite[i, string("mean_dt_pos", j)] = sum( abs.(dep_times_pos    .- hdgrid_2525[j]) .< 0.001 )
				end
				if sum([data_a_finite[i, string("mean_dt_pre", j)] for j=1:61]) == 0.0
					println(n_pre)
					println(dep_times_before)
					error("zero shares!")
				end

			# ASYMPTOTIC
				# mean and sd before
				mymean = sum(hdgrid_2525 .* density_before)
				mysd   = ((hdgrid_2525 .- mymean) .^ 2 .* density_before) |> sum |> sqrt

				data_a[i, :norm_mean] = mymean
				data_a[i, :norm_sd] = mysd

				# probability density
				for j=1:61
					data_a[i, string("mean_dt_pre", j)] = density_before[j]
					data_a[i, string("mean_dt_pos", j)] = density_wcha[j]
				end
		end

		# update for second part of the data (_na)
		for i=1:size(data_na, 1)
			uidp = data_na[i, :uidp]

			# get row number for this uidp
			idx_row = findfirst(data_all_loaded.uidp .== uidp)

			density_before = density_model_before[idx_row, :]
			density_before = density_before ./ sum(density_before)

			density_wcha = density_model_wcha[idx_row, :]
			density_wcha = density_wcha ./ sum(density_wcha)

			# FINITE
				# actual sample sizes
				n_pre = data_na[i, :n_pre]
				ismissing(n_pre) && (n_pre = 1)
				n_pos = data_na[i, :n_pos]
				ismissing(n_pos) && (n_pos = 1)

				# testing: large sample
				# n_pre = 500
				# n_pos = 500

				# draw N samples from the model distribution
				dep_times_before = sample(hdgrid_2525, Weights(density_before), n_pre)
				dep_times_pos    = sample(hdgrid_2525, Weights(density_wcha), n_pos)

				# normal distribution that approximates pre- distribution
				mymean = mean(dep_times_before)
				mysd = Statistics.std(dep_times_before)
				(n_pre == 1) && (mysd = .4073334)  # this is the median, as in the original coding
				(mysd == 0.0) && (mysd = .4073334)  # if the random draws are exactly equal
				data_na_finite[i, :norm_mean] = mymean
				data_na_finite[i, :norm_sd] = mysd

				# frequencies pre&post
				for j=1:61
					data_na_finite[i, string("mean_dt_pre", j)] = sum( abs.(dep_times_before .- hdgrid_2525[j]) .< 0.001 )
					data_na_finite[i, string("mean_dt_pos", j)] = sum( abs.(dep_times_pos    .- hdgrid_2525[j]) .< 0.001 )
				end

			# ASYMPTOTIC
				# mean and sd before
				mymean = sum(hdgrid_2525 .* density_before)
				mysd   = ((hdgrid_2525 .- mymean) .^ 2 .* density_before) |> sum |> sqrt

				data_na[i, :norm_mean] = mymean
				data_na[i, :norm_sd] = mysd

				# probability density
				for j=1:61
					data_na[i, string("mean_dt_pre", j)] = density_before[j]
					data_na[i, string("mean_dt_pos", j)] = density_wcha[j]
				end
		end

		println("done updating DT data")

	### 4.4. generate synthetic data for ***route choice***

		# n_resp_a x 5 matrix, columns are weeks 0,1,2,3,4
		probs_route0_wcha = sim_model["probs_route0_wcha"]

		for i=1:size(atdf, 1)
			uidp = atdf[i, :uidp]
			study_cycle = atdf[i, :study_cycle]

			# get row number for this uidp
			idx_row = findfirst(data_a_loaded.uidp .== uidp)

			prob_long_route = 1.0 - probs_route0_wcha[idx_row, study_cycle + 1]
			@assert 0.0 <= prob_long_route <= 1.0

			## FINITE
				atdf_finite[i, :is_long_route] = (rand() <= prob_long_route) + 0.0

			## ASYMPTOTIC
				atdf[i, :is_long_route] = prob_long_route

		end

		println("done updating route data")


	### 4.5 -- save to file(s)
	## 4.5a FINITE
		# rootpath_output = string(rootpath, "data/coded_struct/simdata_revision/finite/")
		path_output = string(rootpath_output_fin, "sim", idx_datasim)
		isdir(path_output) || mkdir(path_output)

		CSV.write(path_output * "/coded_uidpn_a.csv", data_a_finite)
		CSV.write(path_output * "/coded_uidpn_na.csv", data_na_finite)
		CSV.write(path_output * "/area_trip_vot.csv", atdf_finite)
		CSV.write(path_output * "/true_params.csv", Tables.table(rand_params), header=false)

		hastar = DATA["COMPMAT"]["hastars_given"]
		hastar_mat_path = path_output * "/hastar.bson"
		BSON.@save hastar_mat_path hastar

		println("done saving finite!")

	## 4.5b ASYMPTOTIC
	    # rootpath_output = string(rootpath, "data/coded_struct/simdata_revision/asymptotic/")
		path_output = string(rootpath_output_asy, "sim", idx_datasim)
		isdir(path_output) || mkdir(path_output)

		CSV.write(path_output * "/coded_uidpn_a.csv", data_a)
		CSV.write(path_output * "/coded_uidpn_na.csv", data_na)
		CSV.write(path_output * "/area_trip_vot.csv", atdf)
		CSV.write(path_output * "/true_params.csv", Tables.table(rand_params), header=false)

		hastar = DATA["COMPMAT"]["hastars_given"]
		hastar_mat_path = path_output * "/hastar.bson"
		BSON.@save hastar_mat_path hastar

		println("done saving asymptotic!")
end
