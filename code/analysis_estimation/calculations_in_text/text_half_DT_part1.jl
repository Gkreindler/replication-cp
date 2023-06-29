using Distributed

run_type = 1
n_procs = 1

gmm_options = Dict(
	"run_main" 					=> true, # main GMM options
		"main_Wstep1_from_moms"		=> false,
		"main_write_result_to_file" => true,
		"main_maxIter" 				=> 1000,
		"main_time_limit" 			=> Int64(1.0 * 3600),  # stop after this time
		"main_show_trace" 			=> false,
		"main_debug" 				=> false,

	"run_boot" 					=> false,	# boot options
		"boot_round" 				=> 1, # 1,2,3 to run 40+40+40 cores
		"boot_write_result_to_file" => false,
		"boot_show_trace" 			=> false,
		"boot_maxIter" 				=> 1000,
		"boot_time_limit" 			=> Int64(1.0 * 3600),  # stop after this time
		"boot_debug" 				=> false,
		"boot_throw_exceptions" 	=> true,

	"local_subfolder"			=> ""
)
	gmm_options["main_n_start_pts"] = 1 # n start conditions main run
	gmm_options["boot_n_runs"]      = 2 # number of boostrap runs
	gmm_options["boot_n_start_pts"] = 2 # n start conditions bootstrap run

# modules, other .jl files, define paths functions
include("../basic-imports.jl")
include("../functions-table.jl")


## 1. Model 2: full model with delta = 90
	gmm_options["results_suffix"] = "_DD_05x"  # leave empty by default
	prep_dirs(gmm_options, rootpath_base)

	# fixed parameters
	theta_fixed_dict = Dict{String,Float64}(
		"delta" => 90.0,
		"fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0, "fe_w3" => 0.0, "fe_w4" => 0.0,
		"gamma_01_factor" => 1.0,  # switch cost 0->1 is this times larger than 0->1
		"prob_respond_alpha" => 100.0, "prop_to_wage" => 0.0
	)
	# model paramters and initial conditions
	ESTIMATION_PARAMS = model_dt_dynamic_params(theta_fixed_dict=theta_fixed_dict, gmm_options=gmm_options)
	# function momfn_2x(theta, DATA, PARAMS)
	# 	return model_dt_dynamic(theta, DATA, PARAMS; dt_factor=2.0)
	# end
	function momfn_05x(theta, DATA, PARAMS)
		return model_dt_dynamic(theta, DATA, PARAMS; dt_factor=0.5)
	end

	# other model options
	ESTIMATION_PARAMS["hastar_inversion_type"] = "distribution"
	ESTIMATION_PARAMS["transition_moment_type"] = 0

	## Run entire GMM
	run_gmm(momfn=momfn_05x, ESTIMATION_PARAMS=ESTIMATION_PARAMS, gmm_options=gmm_options)


# ## Read baseline results
# 	rootpath_params = rootpath_base * "model/estimation_results/results_delta_90/"

# 	main_df = read_results_file(path=rootpath_params)
# 	theta_optimum = select(main_df[(main_df.stage .== 2) .& (main_df.is_best_vec .== 1),:], r"^param_") |> Array |> vec

# ## Read new results
# 	rootpath_base = "D:/Dropbox (Personal)/projects/bang_cp_paper/"
# 	rootpath_gmm = rootpath_base * "analysis/new_2021_gmm/"
# 	rootpath_results      = string(rootpath_gmm,"gmm_diff_DT_DD/results_DD_05x/")

# 	main_df = read_results_file(path=rootpath_results)
# 	theta_05x = select(main_df[(main_df.stage .== 2) .& (main_df.is_best_vec .== 1),:], r"^param_") |> Array |> vec

# ## Compare
# 	beta_e, beta_l = (theta_optimum[4], theta_optimum[5]) .* 6.0
# 	beta_e_x05, beta_l_x05 = (theta_05x[4], theta_05x[5]) .* 6.0

# 	println("With x0.5 DT response we get beta_e, beta_l: ", beta_e_x05, ", ", beta_l_x05)
# 	println("With x0.5 DT response we get ratios: ", beta_e_x05/beta_e , ", ", beta_l_x05/beta_l)
