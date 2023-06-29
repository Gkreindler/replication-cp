

function curve_fit_wrapper(idx, myobjfunction, Whalf,
				n_moms,
				theta_initial_vec, theta_lower, theta_upper;
				write_result_to_file=false, results_dir_path="",
				jacobian=nothing,
				my_show_trace=true, my_maxIter=1000, my_time_limit=Inf)

	println("starting iteration ", idx)

	# 	moments are already differenced out so we target zero:
	ymoms = zeros(n_moms)

	if isnothing(jacobian)
		timeittook = @elapsed result = curve_fit(myobjfunction, Whalf, ymoms,
					theta_initial_vec, lower=theta_lower, upper=theta_upper,
					show_trace=my_show_trace, maxIter=my_maxIter, maxTime=float(my_time_limit))
	else
		timeittook = @elapsed result = curve_fit(myobjfunction, Whalf, ymoms,
						theta_initial_vec, lower=theta_lower, upper=theta_upper,
						jacobian=jacobian,
						show_trace=my_show_trace, maxIter=my_maxIter, maxTime=float(my_time_limit))
	end
	println(">>> iter    ", idx, " time: ", timeittook, " converged? ", result.converged)
	println(">>> optimum ", result.param)
	println(">>> obj val ", norm(result.resid))

	# save each round of results to csv file
	# this option is useful if some initial conditions are taking too long
	# the corresponding files will not be created -> you can debug for the
	# corresponding initial conditions
	if write_result_to_file
		results_df = DataFrame(
			"obj_vals" => norm(result.resid),
			"opt_results" => result.converged .+ 0,
			"opt_runtime" => timeittook
		)
		for i=1:length(result.param)
			results_df[:, string("param_", i)] = [result.param[i]]
		end

	    outputfile = string(results_dir_path,"gmm-stage1-run",idx,".csv")
		CSV.write(outputfile, results_df)
	end

	return Dict(
		"result" => result,
		"timeittook" => timeittook
	)
end

"""
	momfn = the moment function
		- should take a single vector argument
		- data etc is already "loaded" in this function
	theta_initials = initial condition (vector)
	theta_lower, theta_upper = bounds (vectors, can include +/-Inf)
	n_runs = number of initial conditions
	n_moms = size of moment function (number of moments) TODO: get this auto
	results_dir_path = where to write results
	Wstep1 = weighting matrix (default = identity matrix).
			Should be Hermitian (this will be enforced).
			Will be Cholesky decomposed
	jacobian = provide jacobian function
	write_result_to_file = write individual run results to files (default=no)
	run_parallel  = individual runs in parallel vs in serial (default=parallel)
	my_show_trace = show trace of curve_fit from LsqFit?
	my_maxIter    = maximum iterations for curve_fit from LsqFit
	debug 		  = pass on to objective function
"""
function gmm_parallel_2step(;
			momfn,
			theta_initials, theta_lower, theta_upper,
			n_runs,
			results_dir_path,
			n_moms=nothing,
			Wstep1=I,
			Wstep1_from_moms=false,
			jacobian=nothing,
			write_result_to_file=false,
			run_parallel=true,
			my_show_trace=false,
			my_maxIter=1000, my_time_limit=Inf,
			debug=false)

## If not provided, compute size of moments by running moment function once
	if isnothing(n_moms)
		n_moms = size(momfn(theta_initials[1,:]), 2)
	end

## save starting conditions
## print to csv file
	start_conditions_df = DataFrame(
		"iteration" => 1:n_runs
	)
	for i=1:size(theta_initials)[2]
		start_conditions_df[!, string("param_", i)] = vec(theta_initials[:, i])
	end

	outputfile = string(results_dir_path,"gmm-initial-conditions.csv")
	CSV.write(outputfile, start_conditions_df)

## Compute initial weighting matrix based on moments covariances at initial "best guess" theta

## initial weighting matrix -> use "optimal" around plausible parameters

if Wstep1_from_moms

	# initial guess = median along all runs
	theta_initial = median(theta_initials, dims=1) |> vec

	# evaluable moments
	mom_matrix = momfn(theta_initial)

	# compute optimal W matrix -> ensure it's Hermitian
	nmomsize = size(mom_matrix, 1)
	Wstep1 = Hermitian(transpose(mom_matrix) * mom_matrix / nmomsize)

	# if super small determinant:
	if det(Wstep1) < 1e-100
		println(" Matrix determinant very low: 1e-100. Adding 0.001 * I.")
		Wstep1 = Wstep1 + 0.001 * I
	end

	# invert (still Hermitian)
	Wstep1 = inv(Wstep1)
	Wstep1 = Wstep1 * det(Wstep1)^(-1/size(Wstep1,1))

	println("Step 1 weighting matrix: determinant=", det(Wstep1))
else
	Wstep1 = diagm(ones(n_moms))
end
	## Save matrix to file
	outputfile = string(results_dir_path,"gmm-stage0-optW.csv")
	CSV.write(outputfile, Tables.table(Wstep1), header=false)

## FIRST stage
    println("GMM => 1. Launching FIRST stage, number of parallel runs: ", n_runs)

	# curve_fit in LsqFit requires to give the objective function data (x)
	#   here, we hack this to give the GMM weighting matrix (Cholesky half) to
	# identityWhalf = diagm(ones(n_moms))
	initialWhalf = Matrix(cholesky(Hermitian(Wstep1)).L)
	@assert norm(initialWhalf * transpose(initialWhalf) - Wstep1) < 1e-10

	# subdirectory for individual run results
	if write_result_to_file
		print("GMM => 1. creating subdirectory to save individual run results...")
		results_subdir_path = string(results_dir_path, "stage1")
		isdir(results_subdir_path) || mkdir(results_subdir_path)
	end

	# objective function with moment function "loaded"
	gmm_obj_fn = (anyWhalf, any_theta) ->
			gmm_obj(theta=any_theta, Whalf=anyWhalf, momfn=momfn, debug=debug)

	if run_parallel && n_runs > 1

	    all_runs_results = pmap(
	        idx -> curve_fit_wrapper(idx, gmm_obj_fn, initialWhalf,
							n_moms,
	                        theta_initials[idx,:], theta_lower, theta_upper,
							write_result_to_file=write_result_to_file,
							jacobian=jacobian,
							results_dir_path=string(results_dir_path, "stage1/"),
	                        my_show_trace=my_show_trace,
							my_maxIter=my_maxIter, my_time_limit=my_time_limit),
	        1:n_runs
	    )
	else
		all_runs_results = Vector{Any}(undef, n_runs)
		for idx=1:n_runs
	        all_runs_results[idx] = curve_fit_wrapper(idx, gmm_obj_fn, initialWhalf,
							n_moms,
	                        theta_initials[idx,:], theta_lower, theta_upper,
							write_result_to_file=write_result_to_file,
							jacobian=jacobian,
							results_dir_path=string(results_dir_path, "stage1/"),
	                        my_show_trace=my_show_trace,
							my_maxIter=my_maxIter, my_time_limit=my_time_limit)
		end
	end

    # collect results
    obj_vals = [norm(all_runs_results[idx]["result"].resid) for idx=1:n_runs]
    opt_thetas = [all_runs_results[idx]["result"].param for idx=1:n_runs]
	opt_results = [all_runs_results[idx]["result"].converged for idx=1:n_runs]
    opt_runtime = [all_runs_results[idx]["timeittook"] for idx=1:n_runs]


	# pick best
	idx_best = argmin(obj_vals)
	is_best_vec = ((1:n_runs) .== idx_best) .+ 0
	theta_stage1 = opt_thetas[idx_best]
	obj_val_stage1 = obj_vals[idx_best]

	println("GMM => 1. FIRST STAGE optimal theta   ", theta_stage1)
	println("GMM => 1. FIRST STAGE optimal obj val ", obj_val_stage1)

    ## print to csv file
	stage1_df = DataFrame(
		"run" => 1:n_runs,
		"obj_vals" => obj_vals,
		"opt_converged" => opt_results .+ 0,
		"opt_runtime" => opt_runtime,
		"is_best_vec" => is_best_vec
	)
	for i=1:length(all_runs_results[1]["result"].param)
		stage1_df[!, string("param_", i)] = [all_runs_results[idx]["result"].param[i] for idx=1:n_runs]
	end

    outputfile = string(results_dir_path,"gmm-stage1-all.csv")
	CSV.write(outputfile, stage1_df)

	# error("STOP HERE")

## Optimal Weighting Matrix
	println("GMM => 2. Computing optimal weighting matrix")

	mom_matrix = momfn(theta_stage1)

 	# compute optimal W matrix -> ensure it's Hermitian
	nmomsize = size(mom_matrix, 1)
	Wstep2 = Hermitian(transpose(mom_matrix) * mom_matrix / nmomsize)

	# if super small determinant:
	if det(Wstep2) < 1e-100
		debug && println(" Matrix determinant very low: 1e-100. Adding 0.001 * I.")
		Wstep2 = Wstep2 + 0.001 * I
	end

	# invert (still Hermitian)
	Wstep2 = inv(Wstep2)

	# normalize?
	debug && println("original Determinant and matrix")
	debug && display(det(Wstep2))
	debug && display(Wstep2)

	Wstep2 = Wstep2 * det(Wstep2)^(-1/size(Wstep2,1))

	debug && println("Determinant and matrix:")
	debug && display(det(Wstep2))  # close to 1
	debug && display(Wstep2)

	## Save matrix to file
	outputfile = string(results_dir_path,"gmm-stage1-optW.csv")
	CSV.write(outputfile, Tables.table(Wstep2), header=false)


## SECOND stage
    println("GMM => 3. Launching SECOND stage, number of parallel runs: ", n_runs)

	# curve_fit in LsqFit requires to give the objective function data (x)
	#   here, we hack this to give the GMM weighting matrix (Cholesky half) to
	optimalWhalf = Matrix(cholesky(Wstep2).L)

	debug && display(optimalWhalf)

	# check that this works
	debug && println("checking if Cholesky decomposition worked: (1) difference matrix, (2) norm")
	# display(optimalWhalf * transpose(optimalWhalf) - Wstep2)
	# display(norm(optimalWhalf * transpose(optimalWhalf) - Wstep2))
	debug && display(optimalWhalf * transpose(optimalWhalf) - Wstep2)
	debug && display(norm(optimalWhalf * transpose(optimalWhalf) - Wstep2))
	@assert norm(optimalWhalf * transpose(optimalWhalf) - Wstep2) < 1e-10

	# subdirectory for individual run results
	if write_result_to_file
		print("GMM => 3. creating subdirectory to save individual run results...")
		results_subdir_path = string(results_dir_path, "stage2")
		isdir(results_subdir_path) || mkdir(results_subdir_path)
	end

	# run GMM

	if run_parallel && n_runs > 1
	    all_runs_results = pmap(
	        idx -> curve_fit_wrapper(idx, gmm_obj_fn, optimalWhalf,
							n_moms,
	                        theta_initials[idx,:], theta_lower, theta_upper,
							write_result_to_file=write_result_to_file,
							results_dir_path=string(results_dir_path, "stage2/"),
	                        my_show_trace=my_show_trace,
							my_maxIter=my_maxIter, my_time_limit=my_time_limit),
	        1:n_runs
	    )
	else
		all_runs_results = Vector{Any}(undef, n_runs)
		for idx=1:n_runs
			all_runs_results[idx]=curve_fit_wrapper(idx, gmm_obj_fn, optimalWhalf,
							n_moms,
	                        theta_initials[idx,:], theta_lower, theta_upper,
							write_result_to_file=write_result_to_file,
							results_dir_path=string(results_dir_path, "stage2/"),
	                        my_show_trace=my_show_trace,
							my_maxIter=my_maxIter, my_time_limit=my_time_limit)
		end
	end

    # collect results
    obj_vals = [norm(all_runs_results[idx]["result"].resid) for idx=1:n_runs]
    opt_thetas = [all_runs_results[idx]["result"].param for idx=1:n_runs]
	opt_results = [all_runs_results[idx]["result"].converged for idx=1:n_runs]
    opt_runtime = [all_runs_results[idx]["timeittook"] for idx=1:n_runs]

	# pick best
	idx_best = argmin(obj_vals)
	is_best_vec = ((1:n_runs) .== idx_best) .+ 0
	theta_stage2 = opt_thetas[idx_best]

	obj_val_stage2 = obj_vals[idx_best]

	println("GMM => 3. SECOND STAGE optimal theta   ", theta_stage2)
	println("GMM => 3. SECOND STAGE optimal obj val ", obj_val_stage2)


    ## print to csv file
	stage1_df = DataFrame(
		"run" => 1:n_runs,
		"obj_vals" => obj_vals,
		"opt_converged" => opt_results .+ 0,
		"opt_runtime" => opt_runtime,
		"is_best_vec" => is_best_vec
	)
	for i=1:length(all_runs_results[1]["result"].param)
		stage1_df[!, string("param_", i)] = [all_runs_results[idx]["result"].param[i] for idx=1:n_runs]
	end

    outputfile = string(results_dir_path,"gmm-stage2-all.csv")
	CSV.write(outputfile, stage1_df)

end

function gmm_obj(;theta, Whalf, momfn, debug=false)

	# write parameter vector at current step
	debug && print("theta ", theta, " ")

	# compute moments
	mymoms = momfn(theta)

	# multiply by (Cholesky) half matrice and take means
	mean_mymoms = vec(mean(mymoms, dims=1) * Whalf)

	# write the value of the objective function
	debug && println(transpose(mean_mymoms) * mean_mymoms)

	return mean_mymoms
end

"""
This little function allows us to:
	(1) define function f
	(2) within a function load data (MYDATA)
	(3) within the same function, make it available to all workers
"""
curry(f, MYDATA, MYPARAMS) = theta -> f(theta, MYDATA, MYPARAMS)


function run_gmm(;
		momfn,
		ESTIMATION_PARAMS,
		gmm_options=nothing, ## BETTER
		rootpath_input="",
		rootpath_output="",
		rootpath_boot_output="",
		run_main=true,  ## main GMM options
		main_n_start_pts=1,
		main_Wstep1_from_moms=true,
		main_write_result_to_file=true,
		main_show_trace=false,
		main_maxIter=1000,
		main_time_limit=-1,
		main_debug=false,
		run_boot=false,       	## Bootstrap options
		boot_round=1,  			# rare cases we want to run several bootstraps in separate
		boot_n_runs=1, 			# number of bootstrap runs
		boot_n_start_pts=1,		# number of bootstrap initial conditions
		boot_write_result_to_file=true,
		boot_show_trace=false,
		boot_maxIter=1000,
		boot_time_limit=-1,
		boot_debug=false
	)

## 0. Prep

	if ~isnothing(gmm_options)
		run_main				= gmm_options["run_main"]
		main_n_start_pts		= gmm_options["main_n_start_pts"]
		main_Wstep1_from_moms	= gmm_options["main_Wstep1_from_moms"]
		main_write_result_to_file	= gmm_options["main_write_result_to_file"]
		main_show_trace			= gmm_options["main_show_trace"]
		main_maxIter			= gmm_options["main_maxIter"]
		main_time_limit			= gmm_options["main_time_limit"]
		main_debug				= gmm_options["main_debug"]
		run_boot				= gmm_options["run_boot"]
 		boot_round				= gmm_options["boot_round"]
		boot_n_runs				= gmm_options["boot_n_runs"]
		boot_n_start_pts		= gmm_options["boot_n_start_pts"]
		boot_write_result_to_file	= gmm_options["boot_write_result_to_file"]
		boot_show_trace			= gmm_options["boot_show_trace"]
		boot_maxIter			= gmm_options["boot_maxIter"]
		boot_time_limit			= gmm_options["boot_time_limit"]
		boot_debug				= gmm_options["boot_debug"]
		boot_throw_exceptions   = gmm_options["boot_throw_exceptions"]

		rootpath_input			= gmm_options["rootpath_input"]
		rootpath_output			= gmm_options["rootpath_output"]
		rootpath_boot_output	= gmm_options["rootpath_boot_output"]
	end

	if ~haskey(gmm_options, "run_parallel")
		gmm_options["run_parallel"] = true
	end

## 1. Load data (including computation matrices)
	DATA = load_and_prep(rootpath_input)

## 2. Load DATA into function that is available everywhere
	mymomfunction_loaded = curry(momfn, DATA, ESTIMATION_PARAMS)

## 4. Run main GMM

if run_main
	println("Starting GMM 2 step")

	# save params to file
	params_df = DataFrame(ESTIMATION_PARAMS["param_factors"])
	rename!(params_df, [:param_name, :factor, :fixed_value])
	params_df.fixed_value = replace(params_df.fixed_value, nothing => missing)
	CSV.write(string(rootpath_output, "gmm-param-names.csv"), params_df)

	gmm_parallel_2step(
			momfn=mymomfunction_loaded,
			theta_initials=ESTIMATION_PARAMS["theta_initials"],
			theta_lower=ESTIMATION_PARAMS["theta_lower"],
			theta_upper=ESTIMATION_PARAMS["theta_upper"],
			Wstep1_from_moms=main_Wstep1_from_moms,
			n_runs=main_n_start_pts,
			run_parallel=gmm_options["run_parallel"],
			results_dir_path=rootpath_output,
			write_result_to_file=main_write_result_to_file,
			my_show_trace=main_show_trace,
			my_maxIter=main_maxIter,
			my_time_limit=main_time_limit,
			debug=main_debug)
end




## Bootstrap
if run_boot
	println("Starting Bootstrap")
	# Random number generators (being extra careful) one per bootstrap run
	master_rng = MersenneTwister(123);
	boot_rngs = Vector{Any}(undef, boot_n_runs)

	for i=1:boot_n_runs
		println("creating random number generator for boot run ", i)

		# each bootstrap run gets a different random seed
		# as we run the bootrap in separate rounds, large initial skip
		boostrap_skip = (boot_round-1)*boot_n_runs + i

		boot_rngs[i] = Future.randjump(master_rng, big(10)^20 * boostrap_skip)
	end

	# Folders
	boot_folders = Vector{String}(undef, boot_n_runs)
	for i=1:boot_n_runs
		boot_folders[i] = string(rootpath_boot_output, "boot_run_", i, "/")
	end

	# Run bootstrap
	pmap(
	idx -> bootstrap_2step(
		n_runs=boot_n_start_pts,
		momfn=momfn,
		DATA=DATA,
		ESTIMATION_PARAMS=ESTIMATION_PARAMS,
		# rootpath_input=rootpath_input,
		rootpath_boot_output=boot_folders[idx],
		boot_rng=boot_rngs[idx],
		write_result_to_file=boot_write_result_to_file,
		my_maxIter=boot_maxIter,
		my_time_limit=boot_time_limit,
		throw_exceptions=boot_throw_exceptions
		),
	1:boot_n_runs)
end


end


# function test_fn(x)
# 	a = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# 	return [a, a]
# end
#
#
# theta_initials = zeros(1,2)
#
# theta_lower = [-Inf, -Inf]
# theta_upper = [Inf, Inf]
#
# testpath = "D:/Dropbox (Personal)/projects/bang_cp_paper/analysis/new_2021_gmm/testing_full200/"
#
# gmm_parallel_2step(
# 			momfn=test_fn,
# 			theta_initials=theta_initials,
# 			theta_lower=theta_lower,
# 			theta_upper=theta_upper,
# 			n_runs=1,
# 			results_dir_path=testpath,
# 			debug=true)
