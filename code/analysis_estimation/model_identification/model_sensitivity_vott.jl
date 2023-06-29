using Distributed
using DelimitedFiles


run_type = 1
# 1 = local            -> 1 core
# 2 = server testing   -> 1 core
# 3 = server full run  -> many cores

gmm_options = Dict(
	"run_main" 					=> true, # main GMM options
		"main_Wstep1_from_moms"		=> false,
		"main_write_result_to_file" => true,
		"main_maxIter" 				=> 1000,
		"main_time_limit" 			=> 2.0 * 3600,  # stop after this time
		"main_show_trace" 			=> (run_type < 3),
		"main_debug" 				=> (run_type < 3),

	"run_boot" 					=> true,	# boot options
		"boot_round" 				=> 1, # 1,2,3 to run 40+40+40 cores
		"boot_write_result_to_file" => false,
		"boot_show_trace" 			=> false,
		"boot_maxIter" 				=> 1000,
		"boot_time_limit" 			=> 2.0 * 3600,  # stop after this time
		"boot_debug" 				=> false,
		"boot_throw_exceptions"		=> false,

	"local_subfolder"			=> ""
)

## Paths
include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
rootpath_local = rootpath_base
rootpath_gmm = rootpath_base * "model/"

savepath = rootpath_base * "paper/figures/smfigure3/"
isdir(savepath) || mkdir(savepath)

# select number of processors
	if (run_type == 1) || (run_type == 2)
		n_procs = 1
		gmm_options["main_n_start_pts"] = 2 # n start conditions main run
		gmm_options["boot_n_runs"]      = 2 # number of boostrap runs
		gmm_options["boot_n_start_pts"] = 2 # n start conditions bootstrap run
	else
		n_procs = 43
		gmm_options["main_n_start_pts"] = 120  # n start conditions main run
		gmm_options["boot_n_runs"]      = 120  # number of boostrap runs
		gmm_options["boot_n_start_pts"] = 2    # n start conditions bootstrap run
	end


## Prep
env_path = rootpath_base * "code/analysis_estimation/cp_env/"

using Pkg
Pkg.activate(env_path) # the dollar sign takes the variable from the local (main) worker
Pkg.instantiate()

# additional packages needed here
using DelimitedFiles
using Plots


include("../basic-imports.jl")
include("../functions-table.jl")
include("functions-sensitivity.jl")

## 1. Model 2: full model with delta = 90
	gmm_options["results_suffix"] = ""  # leave empty by default
	prep_dirs(gmm_options, rootpath)

	# fixed parameters
	theta_fixed_dict = Dict{String,Float64}(
		"delta" => 90.0,
		"fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0, "fe_w3" => 0.0, "fe_w4" => 0.0,
		"gamma_01_factor" => 1.0,  # switch cost 0->1 is this times larger than 0->1
		"prob_respond_alpha" => 100.0, "prop_to_wage" => 0.0
	)
	# model paramters and initial conditions
	ESTIMATION_PARAMS = model_nodt_dynamic_params(theta_fixed_dict=theta_fixed_dict, gmm_options=gmm_options)

	# other model options
	ESTIMATION_PARAMS["hastar_inversion_type"] = "distribution"
	ESTIMATION_PARAMS["transition_moment_type"] = 0
	## Run entire GMM
	# run_gmm(momfn=momfn, ESTIMATION_PARAMS=ESTIMATION_PARAMS, gmm_options=gmm_options)

	rootpath_input = gmm_options["rootpath_input"]
	DATA = load_and_prep(rootpath_input)
	BIGMAT = DATA["BIGMAT"]
	CHARGEMAT = DATA["CHARGEMAT"]
	COMPMAT = DATA["COMPMAT"]
	DATAMOMS = DATA["DATAMOMS"]

	data_a = DATAMOMS["data_a"]
	n_all = size(DATAMOMS["data_all"], 1)
	frac_w0 = sum(data_a.a_early .* (data_a.n_obs_0 .> 0)) / n_all
	frac_w1 = sum(data_a.a_early .* (data_a.n_obs_1 .> 0)) / n_all
	frac_w2 = sum(data_a.a_early .* (data_a.n_obs_2 .> 0)) / n_all
	frac_all = reshape([frac_w0, frac_w1, frac_w2], (3,1))

## Read what's necessary for Andrews et al sensitivity measure

	# without delta
	rootpath_results      =  string(rootpath_gmm,"estimation_results/results_nodt/")
	rootpath_boot_results = [string(rootpath_gmm,"estimation_results/results_boot_nodt/")]
	optimal_W_path = string(rootpath_results,"gmm-stage1-optW.csv")

	moms_names_list = ["late 0", "late 1", "late 2", "late 3", "late 4",
				 "early 0", "early 1", "early 2", "early 3", "early 4"]

 	parameter_names_list = ["alpha", "gamma", "fe_late", "mu"]

	theta_optimum, theta_boot, boot_moms, G, W, Lambda, Lambdatilde, ses =
		compute_sensitivity(
	                DATA=DATA,
					ESTIMATION_PARAMS=ESTIMATION_PARAMS,
	                mymomfunction=model_nodt_dynamic,
	                rootpath_results=rootpath_results,
	                rootpath_boot_results=rootpath_boot_results,
	                optimal_W_path=optimal_W_path,
					moms_names_list=moms_names_list,
					parameter_names_list=parameter_names_list,
					myfactor=1e16)



## Compute Jacobian for VOTT (nodt) model

	# normalize by standard deviation of estimated parameters
	# boot_param_std = std(theta_boot, dims=1)
	# myjacobian = G[6:8,[1,2,4]] .* reshape(boot_param_std[[1,2,4]], (1,3))

	# normalize by estimated parameters (alpha, gamma, mu)
	myjacobian = G[6:8,[1,2,4]] .*
			reshape(theta_optimum[[1,2,4]], (1,3)) ./
			frac_all

	# using Plots
	# plot(1:3, G[6:8,1], label="alpha")
	# plot!(1:3, G[6:8,2], label="gamma")
	# plot!(1:3, G[6:8,4], label="mu")

	# hard to interpret!
		# Lambdatilde[[1,2,4], 6:8] |> display

## Read full model and do same analysis
	rootpath_results = 		 string(rootpath_gmm,"estimation_results/results_delta_90/")
	rootpath_boot_results = [string(rootpath_gmm,"estimation_results/results_boot_delta_90/")]

	# parameters
	println("... reading estimates parameters (including bootstrap runs)")
	main_df1, boot_df1, main_df2, boot_df2, diag_df, full_df = get_all_results(
			rootpath_results=rootpath_results,
			rootpath_boot_results=rootpath_boot_results,
			boot_runs=120, pm_range=0.1, debug=false)

	# get matrix or vector of estimates ("^" anchors to start of the word)
	full_theta_optimum = select(main_df2, r"^param_") |> Array |> vec
	full_theta_boot    = select(boot_df2, r"^param_") |> Array

	full_boot_param_std = std(full_theta_boot, dims=1)

	# path = string(rootpath_results, "boot_moms.csv")
	# boot_moms = readdlm(path, ',', Float64)

	path = string(rootpath_results, "jac.csv")
	myjac_df = CSV.read(path, DataFrame)
	myjac_small_df = select(myjac_df[55:57,:], [:alpha, :gamma, :mu])

	# # normalize by standard deviation of estimated parameters
	# full_myjacobian = Matrix(myjac_small_df) .* reshape(full_boot_param_std[[1,2,7]], (1,3))

	# normalize by estimated parameters (alpha, gamma, mu)
	full_myjacobian = Matrix(myjac_small_df) .*
			reshape(full_theta_optimum[[1,2,7]], (1,3)) ./
			frac_all

## rows = moments, columns = alpha/gamma/mu
	myjacobian |> display
	full_myjacobian |> display

# save
	myjacobian_df = DataFrame(myjacobian, :auto)
	push!(myjacobian_df, theta_optimum[[1,2,4]])
	myjacobian_df[!, :moment] = ["w0", "w1", "w2", "optimal"]
	select!(myjacobian_df, [:moment, :x1, :x2, :x3])
	rename!(myjacobian_df, [:moment, :alpha, :gamma, :mu])
	CSV.write(savepath * "jacobian_nodt.csv", myjacobian_df)

	myjacobian_df = DataFrame(full_myjacobian, :auto)
	push!(myjacobian_df, full_theta_optimum[[1,2,7]])
	myjacobian_df[!, :moment] = ["w0", "w1", "w2", "optimal"]
	select!(myjacobian_df, [:moment, :x1, :x2, :x3])
	rename!(myjacobian_df, [:moment, :alpha, :gamma, :mu])
	CSV.write(savepath * "jacobian_full.csv", myjacobian_df)
