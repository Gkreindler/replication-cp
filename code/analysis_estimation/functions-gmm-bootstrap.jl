
function sample_bigmats(
		BIGMAT, CHARGEMAT, COMPMAT, DATAMOMS,
		boot_sample_a, boot_sample_na; debug=false)

	boot_sample = vcat(
			boot_sample_a,
			boot_sample_na .+ length(boot_sample_a))
	@assert length(boot_sample_a) != length(boot_sample_na)

## BIGMAT -- all 3
	# keys(BIGMAT)
	BIGMAT_BOOT = Dict{String, Array{Float64, 3}}()
	for mykey in ["EARLY",
				  "LATE",
				  "sigma_adjust",
				  "sigma_adjust_detour",
				  "tmean_mat_detour",
				  "tmean_mat",
				  "LATE_detour",
				  "EARLY_detour"]

		debug && println("processing BIGMAT key ", mykey)
		myobj = BIGMAT[mykey]
		if size(myobj, 1) == length(boot_sample)
			BIGMAT_BOOT[mykey] = myobj[boot_sample, :, :] |> copy
		else
			@assert size(myobj, 1) == length(boot_sample_a)
			BIGMAT_BOOT[mykey] = myobj[boot_sample_a, :, :] |> copy
		end

	end

## CHARGEMAT -- all 3
	# keys(CHARGEMAT)
	CHARGEMAT_BOOT = Dict{String, Array{Float64, 3}}()
	for mykey in ["dt_charges_detour",
				  "dt_charges",
				  "area_charge_mean",
				  "area_charge_mean_detour"
				  ]
		debug && println("processing CHARGEMAT key ", mykey)
		myobj = CHARGEMAT[mykey]
		if size(myobj, 1) == length(boot_sample)
			CHARGEMAT_BOOT[mykey] = myobj[boot_sample, :, :] |> copy
		else
			@assert size(myobj, 1) == length(boot_sample_a)
			CHARGEMAT_BOOT[mykey] = myobj[boot_sample_a, :, :] |> copy
		end
	end

	for mykey in [
		"zero_dt_charges",
		"zero_area_charge_mean",
		"zero_dt_charges_detour",
		]
		CHARGEMAT_BOOT[mykey] = CHARGEMAT[mykey] |> copy
	end

## CHARGEMAT -- all 3
	# keys(COMPMAT)
	COMPMAT_BOOT = Dict{String, Array{Float64, N} where N}()
	for mykey in [
			# "initial_conditions_free_flat",
			"prob_sum_mat",
			"prob_route0_free",
			"hastars",
			"eu_mat_detour",
			"prob_route0_wcha",
			"dt_choice",
			"dt_choice_detour",
			"eu_mat_wcha",
			"eu_mat",
			"prob_sum_mat_detour",
			"dt_choice_wcha"
		]

		myobj = COMPMAT[mykey]
		if size(myobj, 1) == length(boot_sample)
			COMPMAT_BOOT[mykey] = myobj[boot_sample, :, :] |> copy
		else
			@assert size(myobj, 1) == length(boot_sample_a)
			COMPMAT_BOOT[mykey] = myobj[boot_sample_a, :, :] |> copy
		end
	end

	COMPMAT_BOOT["initial_conditions_2d"] = COMPMAT["initial_conditions_2d"][boot_sample_a, :] |> copy

## DATAMOMS -- all 3
	# keys(DATAMOMS)
	DATAMOMS_BOOT = Dict{String, Any}()

	# vectors of size
	for mykey in [
		  "sample_dynamic_vot"
		  "cha_ct23"
		  "sample_dt_pre"
		  "moment_het_w1"
		  "sample_dt_contr"
		  "dt_choice_pre"
		  "sample_tr"
		  "cha_tr"
		  "moment_het_w0"
		  "mean_long_mat"
		  "nct23"
		  "mean_dt_pos"
		  "sample_dt_treat"
		  "ntr"
		  "sample_late"
		  "sample_dt_pos"
		  "mean_dt_pre"
		  "sample_ct23"
		  "sample_early"
		]

		debug && println("processing DATAMOMS key ", mykey)
		myobj = DATAMOMS[mykey]

		if size(myobj, 1) == length(boot_sample)
			my_sample = boot_sample
		else
			@assert size(myobj, 1) == length(boot_sample_a)
			my_sample = boot_sample_a
		end

		if length(size(myobj)) == 3
			DATAMOMS_BOOT[mykey] = myobj[my_sample, :, :] |> copy
		elseif length(size(myobj)) == 2
			DATAMOMS_BOOT[mykey] = myobj[my_sample, :] |> copy
		elseif length(size(myobj)) == 1
			DATAMOMS_BOOT[mykey] = myobj[my_sample] |> copy
		else
			error("DATAMOMS error key", mykey)
		end
	end

	### Data Frame
	DATAMOMS_BOOT["data_a"] = copy(DATAMOMS["data_a"][boot_sample_a, :])
	DATAMOMS_BOOT["data_all"] = copy(DATAMOMS["data_all"][boot_sample, :])

	for mykey in ["n_hd_bins", "hdgrid", "n_resp_vec", "hd_bins_hl", "shadow_rates"]
		DATAMOMS_BOOT[mykey] = copy(DATAMOMS[mykey])
	end

	return BIGMAT_BOOT, CHARGEMAT_BOOT, COMPMAT_BOOT, DATAMOMS_BOOT
end




## Serial (not parallel) bootstrap with multiple initial conditions

function bootstrap_2step(;n_runs,
					momfn,
					DATA, ESTIMATION_PARAMS,
					rootpath_boot_output,
					boot_rng=nothing,
					bootstrap_samples=nothing,
					Wstep1_from_moms=true,
					run_parallel=false,          # currently, always should be false
					write_result_to_file=false,
					my_show_trace=false,
					my_maxIter=100,
					my_time_limit=-1,
					throw_exceptions=true
					)

try
## Announce
	println("\nBOOT => Starting Bootstrap ", rootpath_boot_output)

## folder
	isdir(rootpath_boot_output) || mkdir(rootpath_boot_output)

## load data and prepare

	if isnothing(boot_rng)
		if isnothing(bootstrap_samples)
			error("must provide boot_rng or boot_sample")
		end
	else
		if ~isnothing(bootstrap_samples)
			error("can only provide boot_rng or boot_sample, not both")
		end
		bootstrap_samples = Dict(
			"a" => StatsBase.sample(boot_rng, 1:205, 205),
			"na" => StatsBase.sample(boot_rng, 1:99, 99)
		)
	end
	# bootstrap_samples = nothing

	## Save boot samples to file
	outputfile = string(rootpath_boot_output, "boot_sample_a.csv")
	CSV.write(outputfile, Tables.table(bootstrap_samples["a"]), header=false)
	outputfile = string(rootpath_boot_output, "boot_sample_na.csv")
	CSV.write(outputfile, Tables.table(bootstrap_samples["na"]), header=false)

	# save params to file
	params_df = DataFrame(ESTIMATION_PARAMS["param_factors"])
	rename!(params_df, [:param_name, :factor, :fixed_value])
	params_df.fixed_value = replace(params_df.fixed_value, nothing => missing)
	CSV.write(string(rootpath_boot_output, "gmm-param-names.csv"), params_df)


## Sample existing loaded matrices:
	BIGMAT_BOOT, CHARGEMAT_BOOT, COMPMAT_BOOT, DATAMOMS_BOOT =
		sample_bigmats(
		DATA["BIGMAT"], DATA["CHARGEMAT"], DATA["COMPMAT"], DATA["DATAMOMS"],
		bootstrap_samples["a"], bootstrap_samples["na"])

	DATA_BOOT=Dict(
		"BIGMAT" => BIGMAT_BOOT,
		"CHARGEMAT" => CHARGEMAT_BOOT,
		"COMPMAT" => COMPMAT_BOOT,
		"DATAMOMS" => DATAMOMS_BOOT
	)

## define the moment function with Boostrap Data
	mymomfunction_loaded = curry(momfn, DATA_BOOT, ESTIMATION_PARAMS)

## run 2-step GMM and save results
	gmm_parallel_2step(
			momfn=mymomfunction_loaded,
			theta_initials=ESTIMATION_PARAMS["theta_initials_boot"],
			theta_lower=ESTIMATION_PARAMS["theta_lower"],
			theta_upper=ESTIMATION_PARAMS["theta_upper"],
			Wstep1_from_moms=Wstep1_from_moms,
			n_runs=n_runs,
			results_dir_path=rootpath_boot_output,
			write_result_to_file=write_result_to_file,  # Don't write each step
			run_parallel=run_parallel,          # Will run the entire boot function in parallel so this is serial
			my_show_trace=my_show_trace,
			my_maxIter=my_maxIter, my_time_limit=my_time_limit,
			debug=false)

catch e
	println("BOOTSTRAP_EXCEPTION for ", rootpath_boot_output)
	# println("BOOTSTRAP_EXCEPTION for ", rootpath_boot_output)
	# println("BOOTSTRAP_EXCEPTION for ", rootpath_boot_output)

	bt = catch_backtrace()
	msg = sprint(showerror, e, bt)
	println(msg)

	if throw_exceptions
		throw(e)
	end
end

end


function boot_cleanup(;
			rootpath_input, boot_folder,
			idx_batch, idx_run, idx_boot_file,
			mymomfunction,
			my_show_trace=false,
			my_maxIter=100,
			my_time_limit=-1
			)
	#=
		1. read optW and check not empty
		2. load data and big mats
		3. collect initial condtions to run
		4. define function
		5. run wrapper and write to file
		6. write all results to 2nd stage file
	=#
	println("\n\n\n... Processing batch ", idx_batch, " run ", idx_run, "\n\n\n")

	## Step 1. Read optW and check it's ok
	mypath = string(boot_folder, "gmm-stage1-optW.csv")
	optW = readdlm(mypath, ',', Float64)

	# define everywhere
	optimalWhalf = Matrix(cholesky(Hermitian(optW)).L)
	@assert size(optW) == (60, 60)

	## Read initial conditions
	mypath = string(boot_folder, "gmm-initial-conditions.csv")
	theta_initials_boot_df = CSV.read(mypath, DataFrame)
	select!(theta_initials_boot_df, r"param_")
	theta_initials_boot = Matrix(theta_initials_boot_df)

	# display(theta_initials_boot[1,:])

	## Step 2. Load data
	begin
		## initial params
		n_ha = 79
		n_hd = 79 # one per 5 minutes

		hdgrid = vec(Array(range(-2.5,stop=4.0,length=n_hd)))
		hagrid = vec(Array(range(-2.0,stop=4.5,length=n_ha)))

		bin_l = -2.5
		bin_h = 2.5
		n_hd_bins = 61

		n_hd_mom_l = 7 # corresponds to -2h
		n_hd_mom_h = 55 # corresponds to 2h

		n_moms = (n_hd_mom_h - n_hd_mom_l + 1) + 2 + 12

		## load boot samples:
		bootstrap_samples = Dict{String, Vector{Int64}}()

		## Load boot samples to file
		path = string(boot_folder,"boot_sample_a.csv")
		bootstrap_samples["a"] = readdlm(path, ',', Float64) |> vec
		# display(bootstrap_samples["a"])

		path = string(boot_folder,"boot_sample_na.csv")
		bootstrap_samples["na"] = readdlm(path, ',', Float64) |> vec

		## load data and prepare
		BIGMAT_BOOT, CHARGEMAT_BOOT, COMPMAT_BOOT, DATAMOMS_BOOT =
			load_and_prep(rootpath_input, hdgrid=hdgrid, hagrid=hagrid, bootstrap_samples=bootstrap_samples)
		DATAMOMS_BOOT["n_hd_bins"] = n_hd_bins
		DATAMOMS_BOOT["hd_bins_hl"] = [bin_h, bin_l, n_hd_mom_l, n_hd_mom_h]
		DATAMOMS_BOOT["hdgrid"] = hdgrid

		n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = DATAMOMS_BOOT["n_resp_vec"]
	end


	## Step 4. Define function (with boot data loaded)
		my_momfn = theta -> mymomfunction(theta,
							BIGMAT=BIGMAT_BOOT,
							CHARGEMAT=CHARGEMAT_BOOT,
							COMPMAT=COMPMAT_BOOT,
							DATAMOMS=DATAMOMS_BOOT)

	## Step 5. Run wrapper and write to file

		# objective function with moment function "loaded"
		gmm_obj_fn = (anyWhalf, any_theta) ->
			gmm_obj(theta=any_theta, Whalf=anyWhalf, momfn=my_momfn, debug=false)

		# subdirectory for individual run results
		print("GMM => 3. creating subdirectory to save individual run results...")
		results_subdir_path = string(boot_folder, "stage2")
		isdir(results_subdir_path) || mkdir(results_subdir_path)

		# run GMM
		# n_runs = length(step2_runs)
		n_moms = size(optW)[1]

		curve_fit_wrapper(idx_boot_file,
					gmm_obj_fn, optimalWhalf,
					n_moms,
					theta_initials_boot[idx_boot_file,:],
					theta_lower,
					theta_upper,
					write_result_to_file=true,
					results_dir_path=string(boot_folder, "stage2/"),
					my_show_trace=my_show_trace,
					my_maxIter=my_maxIter,
					my_time_limit=my_time_limit)

end
