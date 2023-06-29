using Distributed

run_type = 3
# 1 = local            -> 1 core
# 2 = server testing   -> 1 core
# 3 = server full run  -> many cores

include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package

gmm_options = Dict(
	"run_main" 					=> true, # main GMM options
		"main_Wstep1_from_moms"		=> false,
		"main_write_result_to_file" => true,
		"main_maxIter" 				=> 1000,
		"main_time_limit" 			=> 2.0 * 3600,  # stop after this time
		"main_show_trace" 			=> false,
		"main_debug" 				=> false,

	"run_boot" 					=> true,	# boot options
		"boot_round" 				=> 1, # 1,2,3 to run 40+40+40 cores
		"boot_write_result_to_file" => false,
		"boot_show_trace" 			=> false,
		"boot_maxIter" 				=> 1000,
		"boot_time_limit" 			=> 2.0 * 3600,  # stop after this time
		"boot_debug" 				=> false,
		"boot_throw_exceptions" 	=> true,

	"local_subfolder"			=> ""
)


# select number of processors
	if (run_type == 1) || (run_type == 2)
		n_procs = 2
		gmm_options["main_n_start_pts"] = 2 # n start conditions main run
		gmm_options["boot_n_runs"]      = 2 # number of boostrap runs
		gmm_options["boot_n_start_pts"] = 2 # n start conditions bootstrap run
	else
		n_procs = 43
		gmm_options["main_n_start_pts"] = 120  # n start conditions main run
		gmm_options["boot_n_runs"]      = 120  # number of boostrap runs
		gmm_options["boot_n_start_pts"] = 2    # n start conditions bootstrap run
	end

# modules, other .jl files, define paths functions
include("basic-imports.jl")

Random.seed!(123)

## 1. Model 2: full model with delta = 90
	gmm_options["results_suffix"] = "_with_delta"  # leave empty by default
	prep_dirs(gmm_options, rootpath)

	# fixed parameters
	theta_fixed_dict = Dict{String,Float64}(
		# "delta" => 90.0,
		"fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0, "fe_w3" => 0.0, "fe_w4" => 0.0,
		"gamma_01_factor" => 1.0,  # switch cost 0->1 is this times larger than 0->1
		"prob_respond_alpha" => 100.0, "prop_to_wage" => 0.0
	)
	# model paramters and initial conditions
	ESTIMATION_PARAMS = model_dt_dynamic_params(theta_fixed_dict=theta_fixed_dict, gmm_options=gmm_options)
	momfn = model_dt_dynamic

	# other model options
	ESTIMATION_PARAMS["hastar_inversion_type"] = "distribution"
	ESTIMATION_PARAMS["transition_moment_type"] = 3

	## Run entire GMM
	run_gmm(momfn=momfn, ESTIMATION_PARAMS=ESTIMATION_PARAMS, gmm_options=gmm_options)
