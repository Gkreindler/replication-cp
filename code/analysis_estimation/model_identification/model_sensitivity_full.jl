#=
This script computes the Jacobian and Andrews et al (2017) measure for the full model (to be used in Table SM.XI), and generates Figure SM.3
=#
using Distributed


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

# paths
	# rootpath = "D:/Dropbox (Personal)/projects/bang_cp_paper/"
	# rootpath_input  = string(rootpath, "data/coded_sandbox/")
	# rootpath_output = string(rootpath, "analysis/new_2021_gmm/testing_full23/")

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
	ESTIMATION_PARAMS = model_dt_dynamic_params(theta_fixed_dict=theta_fixed_dict, gmm_options=gmm_options)
	# momfn = model_dt_static

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


## Read what's necessary for Andrews et al sensitivity measure
	

	# without delta
	# rootpath_results      = string(rootpath_gmm,"gmm37/results/")
	# rootpath_boot_results = [string(rootpath_gmm,"gmm37/results_boot1/")]
	# optimal_W_path = string(rootpath_gmm,"gmm37/results/gmm-stage1-optW.csv")

	rootpath_results      =  string(rootpath_gmm, "estimation_results/results_delta_90/")
	rootpath_boot_results = [string(rootpath_gmm, "estimation_results/results_boot_delta_90/")]
	optimal_W_path = string(rootpath_results,"gmm-stage1-optW.csv")


	moms_dt = "dt" .* string.(1:49)
	moms_area = ["late 0", "late 1", "late 2", "late 3", "late 4",
				 "early 0", "early 1", "early 2", "early 3", "early 4"]
	 moms_names_list = vcat(moms_dt, moms_area)

	 parameter_names_list = ["alpha", "gamma", "fe_late", "beta_e", "beta_l", "sigma_dt", "mu"]

	theta_optimum, theta_boot, boot_moms, G, W, Lambda, Lambdatilde, ses =
		compute_sensitivity(
	                DATA=DATA,
					ESTIMATION_PARAMS=ESTIMATION_PARAMS,
	                mymomfunction=model_dt_dynamic,
	                rootpath_results=rootpath_results,
	                rootpath_boot_results=rootpath_boot_results,
	                optimal_W_path=optimal_W_path,
					moms_names_list=moms_names_list,
					parameter_names_list=parameter_names_list,
					myfactor=1e16)


## Ploting
	hdgrid = DATAMOMS["hdgrid"]
	hdgrid_2525 = vec(Array(range(-2.5,stop=2.5,length=61)))
	n_hd_mom_l = 7 # corresponds to -2h
	n_hd_mom_h = 55 # corresponds to 2h
	hdgrid_dtmom = hdgrid_2525[n_hd_mom_l:n_hd_mom_h]

	dt_charges = max.(0.0, min.(1.0, min.(   # peak
						(hdgrid_dtmom .+ 1.5),  # up-ramp
						(1.5 .- hdgrid_dtmom)))) # down-ramp

### Sensitivity measure
	plot(hdgrid_dtmom, 50 .* dt_charges, color=:blue, alpha=0.5, linewidth=2, label="Departure Time Rate")

	xs = hdgrid_dtmom
	ys = 6 .* Lambdatilde[4,1:49]
	model = loess(xs, ys, span=0.15)
    us = range(extrema(xs)...; step = 0.1)
    vs = Loess.predict(model, us)

	scatter!(hdgrid_dtmom, ys, markercolor=:green, markeralpha=0.5,
				markerstrokewidth=0, markersize=4, label=" ")
	plot!(us, vs, color=:green, linestyle=:dash, linewidth=4,
			label="Schedule Cost Early \$\\beta_E\$")

	ys = 6 .* Lambdatilde[5,1:49]
	model = loess(xs, ys, span=0.20)
    us = range(extrema(xs)...; step = 0.1)
    vs = Loess.predict(model, us)

	scatter!(xs, ys, markershape=:square,
			color=:gray, markerstrokewidth=0, markeralpha=0.5,
			markerstrokecolor=:gray, markersize=3, label=" ")
	plot!(us, vs, color=:black, linewidth=4,
			label="Schedule Cost Late \$\\beta_L\$",
			xlabel="Departure time moment (hours relative to peak)",
			ylabel="Parameter Change (INR/h)",
			xticks=[-2.0,-1.5,-0.5,0.5,1.5,2.0],
			title="Sensitivity (\$d\\theta/d m\$)",
			legend=:topright) |> display

	savefig(savepath * "deptime_sensitivity.pdf")

### Jacobian
	ys_beta_e = G[1:49, 4] |> vec
	ys_beta_l = G[1:49, 5] |> vec

	scale_factor = maximum(abs.(ys_beta_e)) / 2

	plot(hdgrid_dtmom, scale_factor * dt_charges, color=:blue, alpha=0.5, linewidth=2, label="Departure Time Rate")

	plot!(xs, ys_beta_e, color=:green, linestyle=:dash, linewidth=3,
			label="Schedule Cost Early \$\\beta_E\$")


	plot!(xs, ys_beta_l, color=:black, alpha=0.75, linewidth=3,
			label="Schedule Cost Late \$\\beta_L\$",
			xlabel="Departure time moment (hours relative to peak)",
			ylabel="Moment Change",
			ylim=(-0.00025, 0.00025),
			xticks=[-2.0,-1.5,-0.5,0.5,1.5,2.0],
			title="Jacobian (\$dm/d\\theta\$)",
			legend=:topright) |> display

	savefig(savepath * "deptime_jacobian.pdf")
