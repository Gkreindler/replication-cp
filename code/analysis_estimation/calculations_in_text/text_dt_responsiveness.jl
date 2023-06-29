"""

    Using estimated parameters, analyze the departure time response in minutes to
    a congestion charge that is equal to one hour's wage, and other values.

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

## Invert ideal arrival time distribution
 	## unpack parameters
	theta_dict = get_params_fn(theta_optimum, ESTIMATION_PARAMS["param_factors"])

	## DT parameters
	# convert alpha from INR per minute to INR per hour
	theta_dt = [theta_dict["alpha"], theta_dict["beta_e"], theta_dict["beta_l"],
				theta_dict["sigma_dt"], 0.0, theta_dict["prob_respond_alpha"]]

	## upack matrices
	dt_choice_free = COMPMAT["dt_choice"]
	dt_choice_free_detour = COMPMAT["dt_choice_detour"]

	dt_choice_wcha = COMPMAT["dt_choice_wcha"]
	# dt_choice_wcha_detour = copy(COMPMAT["dt_choice_detour"])


	eu_mat_free = COMPMAT["eu_mat"]
	eu_mat_free_detour = COMPMAT["eu_mat_detour"]

	eu_mat_wcha = COMPMAT["eu_mat_wcha"]
	# eu_mat_wcha_detour = copy(COMPMAT["eu_mat_detour"])


	prob_sum_mat = COMPMAT["prob_sum_mat"]
	prob_sum_mat_detour = COMPMAT["prob_sum_mat_detour"]


	dt_choice_pre = DATAMOMS["dt_choice_pre"]

	hastars = COMPMAT["hastars"]

	n_resp, n_resp_a, n_resp_a_ct23, n_resp_a_tr = DATAMOMS["n_resp_vec"]

	# sigma_adjust_detour = BIGMAT["sigma_adjust_detour"]

	sample_dynamic_vot = DATAMOMS["sample_dynamic_vot"]
	sample_dynamic_vot_full = vcat(sample_dynamic_vot, zeros(Bool, n_resp-n_resp_a))
	sample_early = DATAMOMS["sample_early"]
	sample_late  = DATAMOMS["sample_late"]

	initial_conditions = COMPMAT["initial_conditions_2d"]

	data_all = DATAMOMS["data_all"]

	n_ha = 79
	hdgrid = DATAMOMS["hdgrid"]

## steps
	### step 1: compute DT choice probs WITHOUT charges
		### Main route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat"],
				BIGMAT["EARLY"],
				BIGMAT["LATE"],
				BIGMAT["sigma_adjust"],
				CHARGEMAT["zero_dt_charges"],
				dt_choice_free, eu_mat_free, prob_sum_mat, debug=false)

		### Detour route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat_detour"],
				BIGMAT["EARLY_detour"],
				BIGMAT["LATE_detour"],
				BIGMAT["sigma_adjust_detour"],
				CHARGEMAT["zero_dt_charges_detour"],
				dt_choice_free_detour, eu_mat_free_detour, prob_sum_mat_detour,
				debug=false)

		### Main route WITH charges - used for non-detour commuters
		solve_dt(theta_dt,
				BIGMAT["tmean_mat"],
				BIGMAT["EARLY"],
				BIGMAT["LATE"],
				BIGMAT["sigma_adjust"],
				CHARGEMAT["dt_charges"],
				dt_choice_wcha, eu_mat_wcha, prob_sum_mat, debug=false)


	### step 2: invert hastar distribution
		# optimal departure time mapping, without charges
		# using ACTUAL route choice probabilities (where applicable)
		probs_route1_free = DATAMOMS["data_a"].mean_long_0
		probs_route1_free_3d = repeat(probs_route1_free, outer=(1, 1, n_ha))

		dt_choice_mapping = copy(dt_choice_free)
		dt_choice_mapping[1:n_resp_a, :, :] =
			dt_choice_mapping[1:n_resp_a, :, :] .* (1.0 .- probs_route1_free_3d) .+
			dt_choice_free_detour               .*         probs_route1_free_3d

		invert_has(dt_choice=dt_choice_mapping, dt_choice_pre=dt_choice_pre, hastars=hastars)


## 4. Define individual specific charge profile

	### median departure time for each commuter
		# ensure pdf
		dt_choice_pre = dt_choice_pre ./ sum(dt_choice_pre, dims=2)

		# compute cdf
		dt_choice_pre_cdf = cumsum(dt_choice_pre, dims=2)

		# median position = first >= 0.5
		median_positions = findfirst.( x -> (x >= 0.5), eachrow(dt_choice_pre_cdf))
		median_dep_times = hdgrid[median_positions]

		# yes, very similar
		# scatter(data_all.norm_mean, hdgrid[median_position])
		# scatter!(data_all.norm_mean, data_all.norm_mean)

	### Define charges
		println("Step 4. departure time choice arrays")
		hdgrid_mat = repeat(hdgrid, outer=[1 n_resp n_ha]);
		hdgrid_mat = permutedims(hdgrid_mat, [2 1 3])
		println(" ...size of hd mat ", size(hdgrid_mat))

		hdgrid_mat_detour = hdgrid_mat[1:n_resp_a, :, :]


	### "step" DT charges: after median pre-experiment departure time, and hourly wage
		dt_charges_step 	   = (hdgrid_mat .>= median_dep_times) .* 165.0
		dt_charges_step_detour = (hdgrid_mat_detour .>= median_dep_times[1:n_resp_a, :, :]) .* 165.0

	### "mini-ramp" DT charges: after median pre-experiment departure time, and hourly wage
		hramp = 20/60
		dt_charges_slope = hdgrid_mat .- median_dep_times
		dt_charges_slope = min.(hramp, dt_charges_slope)
		dt_charges_slope = max.(-hramp, dt_charges_slope)
		dt_charges_slope = (1 .+ dt_charges_slope ./ hramp) ./ 2.0
		dt_charges_slope = dt_charges_slope .* 165
		# plot(hdgrid, dt_charges_slope[1,:,1])

	### forver slope DT charges
		# dt_charges_slopeinf = hdgrid_mat #.* 165

## average departure time without charges
	dt_choice_free_integrated = sum(dt_choice_free .* hastars, dims=3)
	dt_choice_free_integrated = dt_choice_free_integrated[:,:, 1]
	mean_dt_free = sum(dt_choice_free_integrated .* transpose(hdgrid), dims=2) |> vec










## compute with charges

myslopes = [165, 165*2, theta_dict["beta_e"] * 95/100 ]
# myslopes = (165 .* (0.5:0.1:3.1)) |> collect
mean_delta_dt_minutes = zero(myslopes)

mytext = ""

for islope=1:length(myslopes)

	myslope = myslopes[islope]
	dt_charges_slopeinf = hdgrid_mat .* myslope

	solve_dt(theta_dt,
			BIGMAT["tmean_mat"],
			BIGMAT["EARLY"],
			BIGMAT["LATE"],
			BIGMAT["sigma_adjust"],
			dt_charges_slopeinf,
			dt_choice_wcha, eu_mat_wcha, prob_sum_mat, debug=false)

	dt_choice_wcha_integrated = sum(dt_choice_wcha .* hastars, dims=3)
	dt_choice_wcha_integrated = dt_choice_wcha_integrated[:,:, 1]
	mean_dt_wcha = sum(dt_choice_wcha_integrated .* transpose(hdgrid), dims=2) |> vec

	##
		mytext *= "\n\n\n Average change in departure time (minutes) for slope =" * string(myslope) * "\n"
		println("\n\n\n Average change in departure time (minutes) for slope =", myslope)

	## stats: how many minutes earlier would commuters leave, on average?
		delta_dt_minutes = (mean_dt_wcha .- mean_dt_free) .* 60.0
		describe(delta_dt_minutes)
		mytext *= "mean value = " * string(mean(delta_dt_minutes)) * "\n"
		mean_delta_dt_minutes[islope] = mean(delta_dt_minutes)
end

println(mytext)

text_results_path = rootpath_base * "paper/text_results/model_departure_time_stats.txt"
file = open(text_results_path, "w")
write(file, mytext)
close(file)

plot(myslopes, mean_delta_dt_minutes)
# plot(1 ./ (theta_dict["beta_e"] .- myslopes), mean_delta_dt_minutes)
