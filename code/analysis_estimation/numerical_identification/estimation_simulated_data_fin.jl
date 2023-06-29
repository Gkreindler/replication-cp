using Distributed

# "fin" for finite or "asy" for "asymptotic" synthetic data
num_id_type = "fin" # "asy"

run_type = 3
# 1 = local            -> 1 core
# 2 = server testing   -> 1 core
# 3 = server full run  -> many cores

include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
rootpath_local = rootpath_base
# rootpath_server = "/n/holystore01/LABS/kreindler_lab/Lab/cp-replication/"

gmm_options = Dict(
	"run_main" 					=> true, # main GMM options
		"main_Wstep1_from_moms"		=> true,
		"main_write_result_to_file" => false,
		"main_maxIter" 				=> 1000,
		"main_time_limit" 			=> 1.0 * 3600,  # stop after this time
		"main_show_trace" 			=> (run_type < 3),
		"main_debug" 				=> (run_type < 3),

	"run_boot" 					=> false,	# boot options
		"boot_round" 				=> 1, # 1,2,3 to run 40+40+40 cores
		"boot_write_result_to_file" => false,
		"boot_show_trace" 			=> false,
		"boot_maxIter" 				=> 1000,
		"boot_time_limit" 			=> 1.0 * 3600,  # stop after this time
		"boot_debug" 				=> false,
		"boot_throw_exceptions" 	=> true,

	"local_subfolder"			=> "num_identification"
)


# select number of processors
	if (run_type == 1) || (run_type == 2)
		n_procs = 2
		gmm_options["main_n_start_pts"] = 2 # n start conditions main run
		gmm_options["boot_n_runs"]      = 2 # number of boostrap runs
		gmm_options["boot_n_start_pts"] = 2 # n start conditions bootstrap run
	else
		n_procs = 42
		gmm_options["main_n_start_pts"] = 4  # n start conditions main run
		gmm_options["boot_n_runs"]      = 120  # number of boostrap runs
		gmm_options["boot_n_start_pts"] = 2    # n start conditions bootstrap run
	end

## Prep
# if n_procs > 1
# 	rmprocs(workers())
# 	display(workers())
# 	addprocs(n_procs)
# 	display(workers())
# end
# modules, other .jl files, define paths functions
include("../basic-imports.jl")

@everywhere rootpath = $rootpath

@everywhere gmm_options = $gmm_options
@everywhere run_type = $run_type
@everywhere num_id_type = $num_id_type

## 1. Model 2: full model with delta = 90
	# gmm_options["results_suffix"] = "_finite_sim1"  # leave empty by default
	# prep_dirs(gmm_options, run_type)

@everywhere begin
	# fixed parameters
	theta_fixed_dict = Dict{String,Float64}(
		"delta" => 90.0,
		"fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0, "fe_w3" => 0.0, "fe_w4" => 0.0,
		"gamma_01_factor" => 1.0,  # switch cost 0->1 is this times larger than 0->1
		"prob_respond_alpha" => 100.0, "prop_to_wage" => 0.0
	)
	# model paramters and initial conditions
	ESTIMATION_PARAMS = model_dt_dynamic_params(theta_fixed_dict=theta_fixed_dict, gmm_options=gmm_options)
	# momfn = model_dt_dynamic

	# other model options
	ESTIMATION_PARAMS["hastar_inversion_type"] = "distribution"
	ESTIMATION_PARAMS["transition_moment_type"] = 0
end

	## Run entire GMM
	# run_gmm(momfn=momfn, ESTIMATION_PARAMS=ESTIMATION_PARAMS, gmm_options=gmm_options)


@everywhere function run_gmm_synthetic_data(idx)
		println("\n\n launching for idx ", idx)

	# prepare dirs
		(num_id_type == "fin") && (gmm_options["results_suffix"] = "_finite_sim" * string(idx))
		(num_id_type == "asy") && (gmm_options["results_suffix"] = "_asym_sim" * string(idx))
		prep_dirs(gmm_options, rootpath)

		## asymptotic
		(num_id_type == "asy") && (gmm_options["rootpath_input"] = rootpath  * "data/coded_model/simdata_asymptotic/sim" * string(idx) * "/")

		## finite sample properties
		(num_id_type == "fin") && (gmm_options["rootpath_input"] = rootpath  * "data/coded_model/simdata_finite/sim" * string(idx) * "/")

	# run individual GMM serial -- otherwise not enough workers
		# TODO: implement pmap with specific worker pools
		gmm_options["run_parallel"] = false

	# run
		run_gmm(momfn=model_dt_dynamic, ESTIMATION_PARAMS=ESTIMATION_PARAMS, gmm_options=gmm_options)
end

# run!
pmap(run_gmm_synthetic_data, 1:120)

