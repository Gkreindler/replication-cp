
function load_all(;
			rootpath_input,
			rootpath_input_tech,
			rootpath_params,
			rootpath_params_boot,
			rootpath_output,
			rt_type::String,
			rt_factor::Float64=5.0, # in hours
			rt_power=nothing,
			theta_add_dict=nothing,
			theta_factors_dict=nothing,
            random_seed=nothing,   # random number generator
			boot_idx=nothing,
			bootstrap_samples=nothing,
			proptowage=false,          # preferences propto individual wage
			save_params_to_file=false,
			debug=false)

	if isnothing(random_seed)
		rng=nothing
	else
		rng = MersenneTwister(random_seed)
	end

	##
	if rootpath_output != ""
		isdir(rootpath_output) || mkdir(rootpath_output)
	end

		if ~isnothing(boot_idx)
			println("\n\n\n >>> Running BOOT simulation ", boot_idx, "\n\n\n")
		end

	## Read estimated parameters
		println("\n\n\n >>> Reading estimated parameters")
		param_names = ["alpha", "delta", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

		main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
			get_all_results(
				rootpath_results=rootpath_params,
				rootpath_boot_results=rootpath_params_boot,
				pm_range=0.1, boot_runs=120, debug=false)

		debug && spec_diagnostics_df |> display

		theta_optimum = select(main2_df, r"^param_") |> Matrix |> vec

		thetas_boot = select(boot2_df, r"^param_") |> Matrix

		## for bootstrap!
		if ~isnothing(bootstrap_samples)
			theta_optimum = thetas_boot[boot_idx, :] |> vec
		end

		println("... done getting optimal parameters and bootstrap versions.")

		## define other fixed parameters
		theta_fixed_dict = Dict{String,Float64}(
			"delta" => 90.0, "fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0,
			"fe_w3" => 0.0, "fe_w4" => 0.0, "gamma_01_factor" => 1.0,
			"prob_respond_alpha" => 100.0, "prop_to_wage" => 0.0
		)

		theta_dict = make_theta_dict(theta_optimum, theta_fixed_dict)

		## ADD preference parameters?
		if ~isnothing(theta_add_dict)
			for mykey in keys(theta_add_dict)
				println("... ADDING value of ", mykey, " = ", theta_add_dict[mykey])
				theta_dict[mykey] = theta_add_dict[mykey]
			end
		end

		## Modify preference parameters?
		if ~isnothing(theta_factors_dict)
			for mykey in keys(theta_factors_dict)
				println("... changing value of ", mykey, " from ", theta_dict[mykey], " by a factor of ", theta_factors_dict[mykey])
				theta_dict[mykey] *= theta_factors_dict[mykey]
			end
		end

		println("... final parameters to use:")
		display(theta_dict)

		if save_params_to_file
			theta_optimum_path = rootpath_output * "params.json"
			open(theta_optimum_path,"w") do f
			    JSON.print(f, theta_dict, 4)
			end
			println("... parameters saved to file")
		end

	## load road technology parameters and uncertainty (vcov matrix)

		println("\n\n\n >>> Reading road technology parameters and variance-covariance matrices")
		### LINEAR
			rt_linear_path = rootpath_input_tech * "vcov_density_linear.csv"
			rt_linear_df = CSV.read(rt_linear_path, DataFrame)
			@assert names(rt_linear_df) == ["c_density", "c_const", "vcov_density", "vcov_const"]

			# coefficients: constant, slope
			rt_linear_params = [rt_linear_df[1, :c_const], rt_linear_df[1, :c_density]]
			rt_linear_vcov = Matrix(select(rt_linear_df, [:vcov_density, :vcov_const]))
			rt_linear_vcov = rt_linear_vcov[[2, 1], [2, 1]]

			# sample
			rt_linear_mvnormal = MvNormal(zeros(2), rt_linear_vcov)

			# For bootstrap samples: draw parameters from v-cov matrix
			if ~isnothing(bootstrap_samples)
				rt_linear_params .+= rand(rt_linear_mvnormal)
			end

		### POWER
			rt_pow_path = rootpath_input_tech * "vcov_density_pow.csv"
			rt_pow_df = CSV.read(rt_pow_path, DataFrame)
			@assert names(rt_pow_df) == ["c_const", "c_density", "c_gamma", "vcov_const", "vcov_density", "vcov_gamma"]

			# coefficients: constant, slope,
			rt_pow_params = [rt_pow_df[1, :c_const], rt_pow_df[1, :c_density], rt_pow_df[1, :c_gamma]]
			rt_pow_vcov = Matrix(select(rt_pow_df, [:vcov_const, :vcov_density, :vcov_gamma]))

			rt_pow_mvnormal = MvNormal(zeros(3), rt_pow_vcov)

			# For bootstrap samples: draw parameters from v-cov matrix
			if ~isnothing(bootstrap_samples)
				rt_pow_params .+= rand(rt_pow_mvnormal)
			end


	## initialize
		println("\n\n\n >>> Init agents")
		if save_params_to_file
			rootpath_output_agents = rootpath_output
		else
			rootpath_output_agents = ""
		end
	    agents, delays0, BIGMAT = init_agents(
				rootpath_input=rootpath_input,
				rootpath_output=rootpath_output_agents,
				bootstrap_samples=bootstrap_samples,
				rng=rng,
	    		theta_dict=theta_dict,
				df_expansion_factor=10,
				proptowage=proptowage)

		n_agents = size(agents, 1)
		hdgrid = BIGMAT["hdgrid"]
		n_hd = length(hdgrid)
		delta_t = (hdgrid[2] - hdgrid[1]) * 60

		BIGMAT["ttimes"] = zeros(n_agents, n_hd)
		BIGMAT["km_left_small"] = zeros(304, n_hd)
		BIGMAT["ttimes_small"] = zeros(304, n_hd)
		agents_mean_km_small = agents.mean_km[1:10:3040]

	## road tech
		println("\n\n\n >>> Final road tech parameters")

		### this factor helps bring density in the same units as in the empirical
		# 	section where we estimated the (instant_delay,density) relationship
		# 	the normalization there is density * 24 / 0.5 / n_trips
		#	where 0.5 is half hour = approx avg trip duration
		rt_factor = rt_factor * 2.0 / n_agents

	if rt_type == "linear"
		# define linear function
		rt_params = [rt_linear_params[1], rt_linear_params[2], rt_factor]

	elseif rt_type == "pow"
		# define power function
		rt_params = [rt_pow_params[1], rt_pow_params[2], rt_factor, rt_pow_params[3]]

	else
		 @assert rt_type == "pow2"
		# define power function
		@assert ~isnothing(rt_power)

		rt_params = [rt_linear_params[1], rt_linear_params[2], rt_factor, rt_power]
	end

		println()
		println("... road technology parameters:")
		println(rt_params)
		println()

		if save_params_to_file
			# save to file
			rt_path = rootpath_output * "params_rt.csv"
			CSV.write(rt_path, Tables.table(rt_params), header=false)
		end


	return agents, BIGMAT, delays0, rt_params, delta_t, agents_mean_km_small
end



"""
	rt_factor - default value = 5.0, meaning the model sample represents 5/24 of
				all daily traffic.
"""
function run_sim(;
		rootpath_input,
		rootpath_input_tech,
		rootpath_params,
		rootpath_params_boot,
		rootpath_output,
		rt_type::String,
		rt_factor::Float64=5.0, # in hours
		rt_power=nothing,  # the density exponent for rt_type=="pow2"
		theta_factors_dict=nothing,
		boot_idx=nothing,
		bootstrap_samples=nothing,
		proptowage=false,          # preferences propto individual wage
		prop_update_frac=0.5,      # nash: update choices 0.1 towards new optimal choices
		adaptive_adjust=true,      # socopt: adaptive adjust charges
		adjust_factor_initial=0.3, # socopt: update charge 0.3 towards smc
		step_tol=1e-7,
		step_max=10000,
		msc_run_parallel::Bool=true,
		debug_level=0)

## Load all
	agents, BIGMAT, delays0, rt_params, delta_t, agents_mean_km_small = load_all(
			rootpath_input=rootpath_input,
			rootpath_input_tech=rootpath_input_tech,
			rootpath_params=rootpath_params,
			rootpath_params_boot=rootpath_params_boot,
			rootpath_output=rootpath_output,
			rt_type=rt_type,
			rt_factor=rt_factor,
			rt_power=rt_power,
			theta_factors_dict=theta_factors_dict,
			boot_idx=boot_idx,
			bootstrap_samples=bootstrap_samples,
			proptowage=proptowage,
			save_params_to_file=true)

##
	isdir(rootpath_output) || mkdir(rootpath_output)

	if ~isnothing(boot_idx)
		println("\n\n\n >>> Running BOOT simulation ", boot_idx, "\n\n\n")
	end

## Read estimated parameters
	println("\n\n\n >>> Reading estimated parameters")
	param_names = ["alpha", "delta", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

	main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
		get_all_results(
			rootpath_results=rootpath_params,
			rootpath_boot_results=rootpath_params_boot,
			pm_range=0.1, boot_runs=120, debug=false)

	(debug_level == 2) && spec_diagnostics_df |> display

	theta_optimum = select(main2_df, r"^param_") |> Matrix |> vec

	thetas_boot = select(boot2_df, r"^param_") |> Matrix

	## for bootstrap!
	if ~isnothing(bootstrap_samples)
		theta_optimum = thetas_boot[boot_idx, :] |> vec
	end

	println("... done getting optimal parameters and bootstrap versions.")

	## define other fixed parameters
	theta_fixed_dict = Dict{String,Float64}(
		"delta" => 90.0, "fe_late" => 0.0, "fe_w1" => 0.0, "fe_w2" => 0.0,
		"fe_w3" => 0.0, "fe_w4" => 0.0, "gamma_01_factor" => 1.0,
		"prob_respond_alpha" => 100.0, "prop_to_wage" => 0.0
	)

	theta_dict = make_theta_dict(theta_optimum, theta_fixed_dict)

	## Modify preference parameters?
	if ~isnothing(theta_factors_dict)
		for mykey in keys(theta_factors_dict)
			println("... changing value of ", mykey, " from ", theta_dict[mykey], " by a factor of ", theta_factors_dict[mykey])
			theta_dict[mykey] *= theta_factors_dict[mykey]
		end
	end

	println("... final parameters to use:")
	display(theta_dict)

	theta_optimum_path = rootpath_output * "params.json"
	open(theta_optimum_path,"w") do f
	    JSON.print(f, theta_dict, 4)
	end
	println("... parameters saved to file")

## load road technology parameters and uncertainty (vcov matrix)

	println("\n\n\n >>> Reading road technology parameters and variance-covariance matrices")
	### LINEAR
		rt_linear_path = rootpath_input_tech * "vcov_density_linear.csv"
		rt_linear_df = CSV.read(rt_linear_path, DataFrame)
		@assert names(rt_linear_df) == ["c_density", "c_const", "vcov_density", "vcov_const"]

		# coefficients: constant, slope
		rt_linear_params = [rt_linear_df[1, :c_const], rt_linear_df[1, :c_density]]
		rt_linear_vcov = Matrix(select(rt_linear_df, [:vcov_density, :vcov_const]))
		rt_linear_vcov = rt_linear_vcov[[2, 1], [2, 1]]

		# sample
		rt_linear_mvnormal = MvNormal(zeros(2), rt_linear_vcov)

		# For bootstrap samples: draw parameters from v-cov matrix
		if ~isnothing(bootstrap_samples)
			rt_linear_params .+= rand(rt_linear_mvnormal)
		end

	### POWER
		rt_pow_path = rootpath_input_tech * "vcov_density_pow.csv"
		rt_pow_df = CSV.read(rt_pow_path, DataFrame)
		@assert names(rt_pow_df) == ["c_const", "c_density", "c_gamma", "vcov_const", "vcov_density", "vcov_gamma"]

		# coefficients: constant, slope,
		rt_pow_params = [rt_pow_df[1, :c_const], rt_pow_df[1, :c_density], rt_pow_df[1, :c_gamma]]
		rt_pow_vcov = Matrix(select(rt_pow_df, [:vcov_const, :vcov_density, :vcov_gamma]))

		rt_pow_mvnormal = MvNormal(zeros(3), rt_pow_vcov)

		# For bootstrap samples: draw parameters from v-cov matrix
		if ~isnothing(bootstrap_samples)
			rt_pow_params .+= rand(rt_pow_mvnormal)
		end


## initialize
	println("\n\n\n >>> Init agents")
    agents, delays0, BIGMAT = init_agents(
			rootpath_input=rootpath_input,
			rootpath_output=rootpath_output,
			bootstrap_samples=bootstrap_samples,
    		theta_dict=theta_dict,
			df_expansion_factor=10,
			proptowage=proptowage)

	n_agents = size(agents, 1)
	hdgrid = BIGMAT["hdgrid"]
	n_hd = length(hdgrid)
	delta_t = (hdgrid[2] - hdgrid[1]) * 60

	BIGMAT["ttimes"] = zeros(n_agents, n_hd)
	BIGMAT["km_left_small"] = zeros(304, n_hd)
	BIGMAT["ttimes_small"] = zeros(304, n_hd)
	agents_mean_km_small = agents.mean_km[1:10:3040]

## road tech
	println("\n\n\n >>> Final road tech parameters")

	### this factor helps bring density in the same units as in the empirical
	# 	section where we estimated the (instant_delay,density) relationship
	# 	the normalization there is density * 24 / 0.5 / n_trips
	#	where 0.5 is half hour = approx avg trip duration
	rt_factor = rt_factor * 2.0 / n_agents

if rt_type == "linear"
	# define linear function
	rt_params = [rt_linear_params[1], rt_linear_params[2], rt_factor]

elseif rt_type == "pow"
	# define power function
	rt_params = [rt_pow_params[1], rt_pow_params[2], rt_factor, rt_pow_params[3]]

else
	 @assert rt_type == "pow2"
	# define power function

	# rt_params = [2.246323, .6955761, rt_factor, 1.5]  # pow = 1.5
	# rt_params = [rt_linear_params[1], rt_linear_params[2], rt_factor, 1.5]  # pow = 1.5

	# rt_params = [2.511608, .2247403, rt_factor, 3.0]  # pow = 3.0
	rt_params = [rt_linear_params[1], rt_linear_params[2], rt_factor, rt_power]  # pow = 3.0
end

	println()
	println("... road technology parameters:")
	println(rt_params)
	println()

	# save to file
	rt_path = rootpath_output * "params_rt.csv"
	CSV.write(rt_path, Tables.table(rt_params), header=false)


## Find Nash equilibrium
	println("\n\n\n >>> Starting Nash")
	zero_charges = zeros(Float64, 1, n_hd)

	choice_ne, logsum_ne = nash_iteration(
			agents=agents,
			rt_params=rt_params,
			delays0=delays0,
			dt_charges=zero_charges,
			BIGMAT=BIGMAT,
			compute_logsum=1,
			debug_level=debug_level,
			prop_update_frac=prop_update_frac,
			step_tol=step_tol,
			step_max=step_max,
			ttimes_from_delay=true,
			error_if_fail=true,
			rootpath_output=rootpath_output
		)

## Find social optimum
	println("\n\n\n >>> Starting Social Optimum ")
	socopt_iteration(
			agents=agents,
			rt_params=rt_params,
			delays0=delays0,
			BIGMAT=BIGMAT,
			step_tol_out=1e-2, # 5e-1 for testing
			step_tol_in=step_tol,
			step_max=step_max,
			prop_update_frac=prop_update_frac,
			adaptive_adjust=adaptive_adjust,
			adjust_factor_initial=adjust_factor_initial,
			rootpath_output=rootpath_output,
			msc_run_parallel=msc_run_parallel,
			debug_level=debug_level)
end
