"""
	Put vector of parameters in a dictionary
		- combine fixed parametrs with parametrs to estimate
		- fixed parameters are included in the 3rd entry of param_factors

"""
function make_theta_dict(mytheta, theta_fixed_dict)

	# full parameter list with expansion factors
	# this is the order that will be used in optimization
	# will next drop from optimization any parameters held fixed
	param_factors = [
		("alpha", 				6.0),   ## x 6
		("delta",      			0.01),  ## /100
		("gamma", 				1.0),
		("gamma_01_factor", 	1.0),
		("fe_early", 			1.0),
		("fe_late", 			1.0),
		("fe_w1", 				1.0),
		("fe_w2", 				1.0),
		("fe_w3", 				1.0),
		("fe_w4", 				1.0),
		("beta_e", 				6.0),  ## x 6
		("beta_l", 				6.0),  ## x 6
		("sigma_dt", 			1.0),
		("mu", 					1.0),
		("prob_respond_alpha", 0.01), ## /100
		("prop_to_wage", 		1.0)
	]

	theta_dict = Dict{String, Float64}()
	idx = 1
	for el in param_factors
		param_name = el[1]
		if haskey(theta_fixed_dict, param_name)
			theta_dict[param_name] = theta_fixed_dict[param_name]
		else
			theta_dict[param_name] = mytheta[idx]
			idx = idx + 1
		end

		# expansion factor
		theta_dict[param_name] *= el[2]
	end
	@assert idx == length(mytheta) + 1

	return theta_dict
end

"""
	Store results:
	dt_choice 	 -> choice probabilities
	eu_mat 		 -> expected utility
	# prob_sum_mat -> sum_j exp(beta*u_j)
"""
function solve_dt(theta::Vector{Float64},
				tmean_mat::Array{Float64, 3},
				EARLY::Array{Float64, 3},
				LATE::Array{Float64, 3},
				sigma_adjust::Array{Float64, 3},
				dt_charges::Array{Float64, 3},
				dt_choice::Array{Float64, 3},
				eu_mat::Array{Float64, 3},
				prob_sum_mat::Array{Float64, 3};
				MAXDIFF=300.0, debug=false)

	# Note: We never charge area for route=1, so area_charges_detour::Array{Float64, 3} missing

	# unpack params
	alpha, beta_e, beta_l, sigma, mu, prob_respond_alpha = theta

	### FULL
	@. dt_choice = (-1) * alpha * tmean_mat +
				   (-1) * beta_e * EARLY +
				   (-1) * beta_l * LATE +
				   (-1) * dt_charges

	@. dt_choice = dt_choice / sigma_adjust / sigma

	# at most -300.0 relative to max -> otherwise we get Inf when we exponentiate
	# for detour cases -> choose the max among both routes

	# max entry along hd dimension (use the shape of eu_mat, which is n_resp x 1 x n_ha)
	maximum!(eu_mat, dt_choice)

	# main route utility -> bound from below relative to max
	dt_choice .= max.(dt_choice, eu_mat .- MAXDIFF)

	# take out mean level to make it easier to exponentiate
	# save this -> this is the offset
	mean!(eu_mat, dt_choice)

	# exponentiate
    @. dt_choice = exp(dt_choice - eu_mat)

	# compute choice probabilities
	sum!(prob_sum_mat, dt_choice)
	@. dt_choice = dt_choice / prob_sum_mat

	### Expected utility!
	# add back the offset from eu_mat
	@. eu_mat = (eu_mat + log(prob_sum_mat)) * sigma_adjust * sigma

end


"""
	invert
"""
function invert_has(;dt_choice::Array{Float64, 3},
					dt_choice_pre::Matrix{Float64},
					hastars::Array{Float64, 3},
					pdf_weight=5.0, debug=false)
	### for each person, invert the observed distribution of departure times
	#   to get the distribution of ideal arrival times
	# 	using non-negative least squares
	#   also include a constraint that PDF sums up to 1 (has associated weight)

	n_ha::Int64 = size(hastars)[3]

	temp = zeros(n_ha)
	# a line of ones to impose density sums up to 1
	# pdf_weight is the weight we put on this identity
	myones = reshape(ones(n_ha), 1, n_ha) .* pdf_weight

	for i=1:size(dt_choice)[1]
		debug && println("processing ", i)
		 temp = nonneg_lsq(
		 			vcat(dt_choice[i, :, :], myones),
					vcat(dt_choice_pre[i, :], [pdf_weight]),
					alg=:nnls)

		 # check that we don't deviate too much from
		 debug && if abs(sum(temp) - 1.0) > 0.01
			 # display(temp)
			 println("hastar inversion in line ", i, " PDF total mass not 1, =", sum(temp))
			 # @assert abs(sum(temp) - 1.0) < 0.01
		 end
		 hastars[i, 1, :] = temp ./ sum(temp)
	end
end

"""
Smooth Google Travel times across time of day, using loess.
Separately for each respondent
	(function is identical to estimation)
"""
function travel_times_smooth_sims(data_all; myspan=0.25, debug=false)

	# 0-143 are travel times corresponding to 0:00 - 23:50 in 10 minute intervals
	travel_times = select(data_all, r"dur_gm_*")
	n_resp = size(travel_times)[1]

	dep_times = vec(0:143) ./ 6

	ttimes_smooth = zeros(n_resp, length(dep_times))

	for i=1:n_resp
		debug && println("travel time smoothing agent ", i, " of ", n_resp)

		ttimes = collect(travel_times[i, :])

		# smoothing
		# span = 0.25 gives good looking graphs with just the right amount of smoothing
		model = loess(dep_times, ttimes, span=myspan)

		# us = range(extrema(dep_times)...; step = 0.1)
		ttime_smooth = Loess.predict(model, dep_times)
		ttimes_smooth[i, :] = ttime_smooth

		if debug
			plot(dep_times, mytt)
			display(plot!(dep_times, mytt_smooth))
		end
	end

	return ttimes_smooth
end

"""
Interpolate (previously smoothed) Google Travel times at the departure time
grid points

unlike estimation, grid points are real time (not relative to the morning peak-hour for each commuter)
"""
function travel_times_interpolate_sims(ttimes_smooth, hdgrid; debug=false)

	# 0-143 are travel times corresponding to 0:00 - 23:50 in 10 minute intervals
	n_resp = size(ttimes_smooth)[1]
	dep_times = vec(0:143) ./ 6

	hdgrid_ttimes = zeros(n_resp, length(hdgrid))

	for i=1:n_resp
		debug && println("travel time interpolation agent ", i, " of ", n_resp)

		ttime_smooth = ttimes_smooth[i, :]

		# interpolation
		itp = interpolate((dep_times,), ttime_smooth, Gridded(Linear()))
		hdgrid_ttimes[i, :] = itp(hdgrid);

		if debug
			display(plot!(hdgrid, hdgrid_ttimes[i, :]))
			fdsfd
		end
	end

	return hdgrid_ttimes
end

"""
travel time uncertainty values hard-coded
"""
function SD_fun(my_delays::Matrix{Float64})
    # my_delays is a vector (or matrix) or delays (roughly 1-7 min/km)
    SD_factors::Vector{Float64} = [0.2393907, -0.049169, 0.0405207]

    return 	SD_factors[1] .+
			SD_factors[2] .* my_delays .+
			SD_factors[3] .* my_delays .* my_delays
end


"""
# for log normal variable, goes from (T,S) namely mean and SD
# to mu and sigma2 parameters of the log normal (mean and VAR of log(X)
"""
function ts2musigma(
		T::Matrix{Float64},
		S::Matrix{Float64})

	sigma2 = @. log(1.0 + S * S / (T * T))
	mu = @. log(T) - sigma2 * 0.5
	@. sigma2 = sqrt(sigma2)

	return mu, sigma2
end


"""
Load data for dynamic VOTT moments
- start from trip level data for home-work trips
- collapse at study_cycle (week=0,1,2,3,4) x commuter level
- reshape
"""
function load_area_byweek(rootpath_input; mydebug=false)
	# read raw data
        filepath = string(rootpath_input, "area_trip_vot.csv")
        atdf = CSV.read(filepath, DataFrame)
        println(" ...Area trip dataframe loaded!...")

	# collapse at study_cycle (= week=0,1,2,3,4 where 0 is all pre data) x commuter
		# atdf_mean = combine(groupby(atdf, [:uidp, :study_cycle]),
		# 	:adetourtime => mean => :adetourtime,
		# 	:a_early => mean => :a_early,
		# 	:a_late  => mean => :a_late,
		# 	:high_dttime => mean => :high_dttime,
		# 	:is_long_route => mean => :is_long_route)

	# commuter level
		u_level = select(atdf,
				[:uidp, :adetourtime, :a_early, :a_late, :high_dttime, :acharge_weekmean])
		u_level = unique(u_level)
		sort!(u_level, :uidp)

		atdf[!, :n_obs] .= 1

	# study_cycle x commuter level
	# study_cycle is week=0,1,2,3,4 where 0 is all pre data
		uw_level = combine(groupby(atdf, [:uidp, :study_cycle]),
				:is_long_route => mean => :is_long_route,
				:n_obs => sum => :n_obs
			)

	# full panel "frame" (only individual level data)
		all_weeks = DataFrame(:study_cycle => [0,1,2,3,4])
		uw_cross = crossjoin(u_level, all_weeks)

	# add data to full panel frame
		uw_full = outerjoin(uw_level, uw_cross, on=[:uidp, :study_cycle],
			indicator="_m", validate=(true, true))
		mydebug && display(freqtable(uw_full._m))
		@assert all(uw_full._m .!= "left_only")
		uw_full.has_data = (uw_full._m .== "both")
		select!(uw_full, Not(:_m))

		uw_full[uw_full.has_data .== 0, :is_long_route] .= 0.0
		uw_full[uw_full.has_data .== 0, :n_obs] .= 0

	# TODO: what does this do? drop?
		# uw_full.is_long_route = (uw_full.is_long_route)
		# uw_full.adetourtime = (uw_full.adetourtime)

		# normalized adetourtime
		uw_full.adetourtime_normalized = (uw_full.adetourtime ./ mean(uw_full.adetourtime))

		# sort by week and uidp
		sort!(uw_full, [:study_cycle, :uidp])
		mydebug && display(uw_full)

	## Unstack (reshape wide)
		idvars = [:uidp, :adetourtime, :a_early, :a_late, :high_dttime,
				  :acharge_weekmean, :adetourtime_normalized]
		uw_full_wide1 = unstack(uw_full, idvars, :study_cycle, :is_long_route,
					renamecols=(idx -> string("mean_long_", idx) ))
		uw_full_wide2 = unstack(uw_full, :uidp, :study_cycle, :n_obs        ,
					renamecols=(idx -> string("n_obs_", idx) ))
		uw_full_wide = outerjoin(uw_full_wide1, uw_full_wide2, on="uidp",
					makeunique=true, validate=(true,true))

	return uw_full_wide
end





"""
Load all data and pre-compute large arrays used in model simulation
Step 1. Load data and merge/append
Step 2. Smooth travel times
Step 3. compute the time matrices
Step 4. departure time choice arrays
Step 5. sigma adjustment: to ensure that commuters who travel farther do not
Step 6. Arrays where we store computations
Step 7. Congestion pricing arrays
Step 8. Data Moments
"""
function init_agents(;
		rootpath_input,
		rootpath_output="",
		theta_dict,
		bootstrap_samples=nothing,
		rng=nothing, # random number generator
		proptowage=false,
		# n_ha=85, n_hd=85,
		n_ha=109, n_hd=109,
		df_expansion_factor=10,
		mydebug=false)
		println(" Starting load and prep... ")

## Step 0. Grids -- every 5 minutes
	# hdgrid = vec(Array(range(5,stop=12.0,length=n_hd)))
	# hagrid = vec(Array(range(5.5,stop=12.5,length=n_ha)))

	hdgrid = vec(Array(range(5,stop=14.0,length=n_hd)))
	hagrid = vec(Array(range(5.5,stop=14.5,length=n_ha)))

## Step 1. Load data and merge/append
	# load main data (includes departure time choices, travel times, etc.)
		println("Step 1. Load data and merge/append")
		filepath = string(rootpath_input, "coded_uidpn_a.csv")
		data_a = CSV.read(filepath, DataFrame)
		data_a[!, :area] .= 1

		filepath = string(rootpath_input, "coded_uidpn_na.csv")
		data_na = CSV.read(filepath, DataFrame)
		data_na[!, :area] .= 0

	# Load area data for dynamic VOT: uidp x study cycle level (wide)
		dynamic_area_df_wide = load_area_byweek(rootpath_input, mydebug=mydebug)
		select!(dynamic_area_df_wide, Not(:adetourtime))

		# join to main data
		data_a = outerjoin(data_a, dynamic_area_df_wide, on=:uidp, indicator=:source)
		mydebug && display(freqtable(data_a.source))
		@assert all(data_a.source .!= "right_only")
		data_a.sample_dynamic_vot = data_a.source .== "both"

		# 11 observations without h-w trips post
		mydebug && display(data_a[data_a.source .== "left_only", :uidp])

		select!(data_a, Not(:source))

	# bootstrap sample
		if ~isnothing(bootstrap_samples)
			println("\n... Bootstrap Sample! ...\n")
			boot_sample_a  = bootstrap_samples["a"]
			boot_sample_na = bootstrap_samples["na"]
			data_a  = data_a[boot_sample_a, :]
			data_na = data_na[boot_sample_na, :]
		end

	# replace missing with zeros
		for mycol in ["adetourtime_normalized",
					  "mean_long_0", "mean_long_1", "mean_long_2", "mean_long_3", "mean_long_4",
					  "n_obs_0", "n_obs_1", "n_obs_2", "n_obs_3", "n_obs_4",
					  "a_early", "a_late"]
			data_a[!, mycol][.~data_a.sample_dynamic_vot] .= 0.0
			data_a[!, mycol] = convert.(Float64, data_a[:, mycol])
		end

	# append data for with detour option ("a") and without detour option ("na")
		data_all = vcat(data_a, data_na, cols = :union)
		n_resp = size(data_all)[1]
		n_resp_a = size(data_a)[1]

	# array sizes
		# n_hd = length(hdgrid)
		# n_ha = length(hagrid)

## Step 2. Smooth and interpolate travel times
		println("Step 2. Smooth and interpolate travel times")

		### 1. smooth travel times
		ttimes_smooth = travel_times_smooth_sims(data_all)

		### 2. compute (smooth) detour (when applicable)
		time_gm_9am = ttimes_smooth[:, 55]  # seconds
	    area_detour_inflate_factor = 60 * data_all.adetourtime ./ time_gm_9am .+ 1 # seconds/seconds
		ttimes_smooth_detour = ttimes_smooth .* area_detour_inflate_factor
		ttimes_smooth_detour = ttimes_smooth_detour[1:n_resp_a, :]
		ttimes_smooth_detour = convert.(Float64, ttimes_smooth_detour)

		# in KM
		area_detour_inflate_factor_km = data_all.adetourdist ./ data_all.dist_mode .+ 1.0
		# note: some are < 1 (shorter but slower)
			# scatter(area_detour_inflate_factor_km, area_detour_inflate_factor)

		### 3. sample at hdgrid level
		hdgrid_ttimes        = travel_times_interpolate_sims(ttimes_smooth, hdgrid)
		hdgrid_ttimes_detour = travel_times_interpolate_sims(ttimes_smooth_detour, hdgrid)

## Step 3. compute the time matrices
		println("Step 3. compute the time matrices")

		# parameters for travel time distribution standard deviation
		SD_factor = [0.2393907, -0.049169, 0.0405207]
		SD_f0, SD_f1, SD_f2 = SD_factor

		# in minutes per KM (= "delay")
		hdgrid_delays        = hdgrid_ttimes 		./ 60 ./ data_all.dist_mode
		hdgrid_delays_detour = hdgrid_ttimes_detour ./ 60 ./ (data_a.dist_mode .+ data_a.adetourdist)

		mydebug && describe(hdgrid_delays[:])
		mydebug && describe(hdgrid_delays_detour[:])  # max = 9.42
		@assert all( 1 .< hdgrid_delays[:] .< 10)
		@assert all( 1 .< hdgrid_delays_detour[:] .< 10)

		delay_sd 	    = SD_f0 .+ SD_f1 .* hdgrid_delays        .+ SD_f2 .* hdgrid_delays        .* hdgrid_delays;
		delay_sd_detour = SD_f0 .+ SD_f1 .* hdgrid_delays_detour .+ SD_f2 .* hdgrid_delays_detour .* hdgrid_delays_detour;

		@assert all(maximum(delay_sd[:]) .< 3)
		@assert all(maximum(delay_sd_detour[:]) .< 3.4)

		# convert back to hours
		hdgrid_ttimes_sd		= delay_sd        .* data_all.dist_mode ./ 60;
		hdgrid_ttimes_sd_detour = delay_sd_detour .* data_a.dist_mode   ./ 60;

		hdgrid_ttimes 		 = hdgrid_ttimes ./ 3600
		hdgrid_ttimes_detour = hdgrid_ttimes_detour ./ 3600

		mydebug && describe((hdgrid_ttimes_sd        ./ hdgrid_ttimes)[:])
		mydebug && describe((hdgrid_ttimes_sd_detour ./ hdgrid_ttimes_detour)[:])

		# lognormal parameters (mean and SD of log(Y))
		println(" ... from mean and sd to lognormal parameters")
		tmu       , tsi 	   = ts2musigma(hdgrid_ttimes, hdgrid_ttimes_sd)
		tmu_detour, tsi_detour = ts2musigma(hdgrid_ttimes_detour, hdgrid_ttimes_sd_detour)


## Step 4. departure time choice arrays
		println("Step 4. departure time choice arrays")
		hdgrid_mat = repeat(hdgrid, outer=[1 n_resp n_ha]);
		hdgrid_mat = permutedims(hdgrid_mat, [2 1 3])
		println(" ...size of hd mat ", size(hdgrid_mat))

		hdgrid_mat_detour = hdgrid_mat[1:n_resp_a, :, :]

		hagrid_mat = repeat(hagrid, outer=[1 n_resp n_hd]);
		hagrid_mat = permutedims(hagrid_mat, [2 3 1]);
		println(" ...size of ha mat ", size(hdgrid_mat))

		hagrid_mat_detour = hagrid_mat[1:n_resp_a, :, :]

	### aux matrices
		println(" ...Building large matrices... ")
		# travel time matrices
		println(" ...travel time...")
		tmu_mat        = repeat(tmu  	  , outer=[1 1 n_ha])
		tmu_mat_detour = repeat(tmu_detour, outer=[1 1 n_ha])

		tsi_mat        = repeat(tsi       , outer=[1 1 n_ha])
		tsi_mat_detour = repeat(tsi_detour, outer=[1 1 n_ha])

		tmean_mat 	     = exp.(tmu_mat		   .+ 0.5 .* tsi_mat .^ 2)
		tmean_mat_detour = exp.(tmu_mat_detour .+ 0.5 .* tsi_mat_detour .^ 2)

		# check that the lognormal formulas "work"
		@assert norm(tmean_mat[:,:,1] - hdgrid_ttimes) < 1.0e-12

		# aux
		tmu_sig2_mat 		= tmu_mat        .+ tsi_mat .^ 2
		tmu_sig2_mat_detour = tmu_mat_detour .+ tsi_mat_detour .^ 2

		# travel time to arrive exactly at ideal h_a^star
		t_exact_mat 		= hagrid_mat 		 - hdgrid_mat
		t_exact_mat_detour  = hagrid_mat_detour  - hdgrid_mat_detour
		log_t_exact_mat 	   = log.(max.(0, t_exact_mat))
		log_t_exact_mat_detour = log.(max.(0, t_exact_mat_detour))

		# CDF for arrival loss - these are specific to log normal distribution fn
		println(" ...cdfs...")
			# CDF: part of the formula (missing betas and sign)
			# INT_CDF: integrated probability times penalty
		CDF  	   = cdf.(Normal.(tmu_mat, tsi_mat)              , log_t_exact_mat       )
		CDF_detour = cdf.(Normal.(tmu_mat_detour, tsi_mat_detour), log_t_exact_mat_detour)

		INT_CDF 	   = cdf.(Normal.(tmu_sig2_mat, tsi_mat)			  , log_t_exact_mat)
		INT_CDF_detour = cdf.(Normal.(tmu_sig2_mat_detour, tsi_mat_detour), log_t_exact_mat_detour)

		#
		ha_b4_hd 		= hagrid_mat 		.<= hdgrid_mat
		ha_b4_hd_detour = hagrid_mat_detour .<= hdgrid_mat_detour

	### early expected cost
		println(" ...early and late mats with conditional expectations")
		EARLY 		 = t_exact_mat 		  .* CDF 	    .- tmean_mat        .* INT_CDF
		EARLY_detour = t_exact_mat_detour .* CDF_detour .- tmean_mat_detour .* INT_CDF_detour

		# for the unusual case when ha < hd (super late departure)
		EARLY[ha_b4_hd] .= 0.0
		EARLY_detour[ha_b4_hd_detour] .= 0.0

	### Late expected cost
		LATE        = (-1) .* t_exact_mat         .* (1 .- CDF       ) .+ tmean_mat        .* (1 .- INT_CDF)
		LATE_detour = (-1) .* t_exact_mat_detour  .* (1 .- CDF_detour) .+ tmean_mat_detour .* (1 .- INT_CDF_detour)

		LATE[ha_b4_hd]                = tmean_mat[ha_b4_hd]               .- t_exact_mat[ha_b4_hd]
		LATE_detour[ha_b4_hd_detour]  = tmean_mat_detour[ha_b4_hd_detour] .- t_exact_mat_detour[ha_b4_hd_detour]

	# travel time to arrive exactly at ideal h_a^star
		t_exact_mat = hagrid_mat - hdgrid_mat
		log_t_exact_mat = log.(max.(0, t_exact_mat))


## Step 5. sigma adjustment: to ensure that commuters who travel farther do not
		# mechanically have more precise choices
		println("Step 5. sigma adjustment: to ensure that commuters who travel farther do not")
		mean_mean_km = mean(data_all.mean_km)
		sigma_adjust_vec = data_all.mean_km ./ mean_mean_km
		sigma_adjust_vec_detour = sigma_adjust_vec[1:n_resp_a]

		sigma_adjust = reshape(sigma_adjust_vec, n_resp, 1, 1)
		sigma_adjust_detour = reshape(sigma_adjust_vec_detour, n_resp_a, 1, 1)

## Store all BIG matrices to use in utility fn
		BIGMAT = Dict{String, Array{Float64}}(
			"hdgrid" => hdgrid,
			"tmean_mat" => tmean_mat,
			"tmean_mat_detour" => tmean_mat_detour,
			"EARLY" => EARLY,
			"EARLY_detour" => EARLY_detour,
			"LATE" => LATE,
			"LATE_detour" => LATE_detour,
			"sigma_adjust" => sigma_adjust,
			"sigma_adjust_detour" => sigma_adjust_detour
		)


## Step 6. Arrays where we store computations
	println("Step 6. Arrays where we store computations")

	# store utility
		dt_choice 		 = zeros(n_resp, n_hd, n_ha)
		dt_choice_detour = zeros(n_resp_a, n_hd, n_ha)

		dt_choice_wcha	 = zeros(n_resp, n_hd, n_ha)

	# store expected utility across departure times
		eu_mat 			   = zeros(n_resp, 1, n_ha)
		eu_mat_detour 	   = zeros(n_resp_a, 1, n_ha)

		eu_mat_wcha		   = zeros(n_resp, 1, n_ha)

	# store sum(exp) across departure times
		prob_sum_mat 			   = zeros(n_resp, 1, n_ha)
		prob_sum_mat_detour 	   = zeros(n_resp_a, 1, n_ha)

	# still used for the nested logit model
		prob_route0_free      = ones(n_resp_a, 1, n_ha)
		prob_route0_wcha      = ones(n_resp_a, 1, n_ha)
		# prob_route0_2d = zeros(n_resp_a, n_ha)

	# store (inverted) ideal arrival time distribution
		hastars = zeros(n_resp, 1, n_ha)

	# NOT USED: individual-level route choice random effects
		i_shocks_temp = randn(n_resp_a)

	# put together
		COMPMAT = Dict(
			"dt_choice" => dt_choice,
			"dt_choice_detour" => dt_choice_detour,
			"dt_choice_wcha" => dt_choice_wcha,
			"eu_mat" => eu_mat,
			"eu_mat_detour" => eu_mat_detour,
			"eu_mat_wcha" => eu_mat_wcha,
			"prob_sum_mat" 		  => prob_sum_mat,
			"prob_sum_mat_detour" => prob_sum_mat_detour,
			"prob_route0_free" => prob_route0_free,
			"prob_route0_wcha" => prob_route0_wcha,
			"hastars" => hastars,
			"i_shocks" => repeat(i_shocks_temp, outer=[1,1,n_ha]),
			"i_shocks_flat" => reshape(i_shocks_temp, n_resp_a, 1,1),
			"initial_conditions_2d" => zeros(n_resp_a, 2),
            "initial_conditions_wcha" => zeros(n_resp_a, 1, n_ha),
			"initial_conditions_free" => zeros(n_resp_a, 1, n_ha),
			"initial_conditions_wcha_flat" => zeros(n_resp_a, 1, 1),
			"initial_conditions_free_flat" => zeros(n_resp_a, 1, 1)
		)

## Step 7. Congestion pricing arrays
	println("Step 7. Congestion pricing arrays")

	### DT charges
		dt_charges = max.(0.0, min.(1.0, min.(   # peak
							(hdgrid_mat .+ 1.5),  # up-ramp
							(1.5 .- hdgrid_mat)))) # down-ramp
		dt_charges .*= repeat(data_all.mean_km, outer=[1, n_hd, n_ha])
		dt_charges .*= repeat(data_all.peak_rate, outer=[1, n_hd, n_ha])

	### DT charges for detour (adjust distance!)
		dt_charges_detour = max.(0.0, min.(1.0, min.(   # peak
								(hdgrid_mat_detour .+ 1.5),  # up-ramp
								(1.5 .- hdgrid_mat_detour)))) # down-ramp
		dt_charges_detour .*= repeat(data_a.mean_km .+ data_a.adetourdist, outer=[1, n_hd, n_ha])
		dt_charges_detour .*= repeat(data_a.peak_rate, outer=[1, n_hd, n_ha])

	### AREA charges
		temp = data_all.area_charge_mean
		temp[ismissing.(temp)] .= 0.0
		temp = convert.(Float64, temp)
		area_charge_mean = repeat(temp, outer=[1 1 n_ha])
		area_charge_mean_detour = area_charge_mean[1:n_resp_a, :, :]
		area_charge_mean_detour_flat = reshape(area_charge_mean[1:n_resp_a, 1, 1], n_resp_a, 1, 1)

	### ZERO charges
		zero_dt_charges        = repeat([0.0], outer=[1,1,1])
		zero_dt_charges_detour = repeat([0.0], outer=[1,1,1])

		zero_area_charge_mean        = repeat([0.0], outer=[1,1,1])

	### put together
		CHARGEMAT = Dict(
			"dt_charges" => dt_charges,
			"dt_charges_detour" => dt_charges_detour,

			"zero_dt_charges" => zero_dt_charges,
			"zero_dt_charges_detour" => zero_dt_charges_detour,

			"area_charge_mean" => area_charge_mean,
			"area_charge_mean_detour" => area_charge_mean_detour,
			"area_charge_mean_detour_flat" => area_charge_mean_detour_flat,
			"zero_area_charge_mean" => zero_area_charge_mean,
		)

## Step 8. Invert distribution of ideal departure times
		println("Step 8. Data Moments")

	## DT parameters
	theta_dt = [theta_dict["alpha"], theta_dict["beta_e"], theta_dict["beta_l"],
				theta_dict["sigma_dt"], 0.0, theta_dict["prob_respond_alpha"]]


	### build density of pre departures (based on normal)
		hdgrid2 = repeat(hdgrid, outer=[1, n_resp])
		hdgrid2 = permutedims(hdgrid2, [2 1])
		dt_choice_pre = pdf.(Normal.(data_all.norm_mean .+ data_all.pk_am, data_all.norm_sd), hdgrid2)
		dt_choice_pre = dt_choice_pre ./ sum(dt_choice_pre, dims=2)

	### Invert ideal departure time distribution
		### step 1: compute DT choice probs WITHOUT charges
		### Main route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat"],
				BIGMAT["EARLY"],
				BIGMAT["LATE"],
				BIGMAT["sigma_adjust"],
				CHARGEMAT["zero_dt_charges"],
				dt_choice, eu_mat, prob_sum_mat, debug=false)

		### Detour route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat_detour"],
				BIGMAT["EARLY_detour"],
				BIGMAT["LATE_detour"],
				BIGMAT["sigma_adjust_detour"],
				CHARGEMAT["zero_dt_charges_detour"],
				dt_choice_detour, eu_mat_detour, prob_sum_mat_detour,
				debug=false)

		### Main route WITH charges - used for non-detour commuters
		# solve_dt(theta_dt,
		# 		BIGMAT["tmean_mat"],
		# 		BIGMAT["EARLY"],
		# 		BIGMAT["LATE"],
		# 		BIGMAT["sigma_adjust"],
		# 		CHARGEMAT["dt_charges"],
		# 		dt_choice_wcha, eu_mat_wcha, prob_sum_mat, debug=(debug > 0))


	### step 2: invert hastar distribution
		# optimal departure time mapping, without charges
		# using ACTUAL route choice probabilities (where applicable)
		probs_route1_free = data_a.mean_long_0
		probs_route1_free_3d = repeat(probs_route1_free, outer=(1, 1, n_ha))

		dt_choice_mapping = copy(dt_choice)
		dt_choice_mapping[1:n_resp_a, :, :] =
			dt_choice_mapping[1:n_resp_a, :, :] .* (1.0 .- probs_route1_free_3d) .+
			dt_choice_detour                    .*         probs_route1_free_3d

		invert_has(dt_choice=dt_choice_mapping, dt_choice_pre=dt_choice_pre, hastars=hastars)

	## put together in "agents" dataframe
		agents = copy(data_all)

		select!(agents, [:uidp, :mean_km, :dist_mode, :dist_mean,
				:norm_mean, :norm_sd, :pk_am, :inverse_income_with_pred])
		# :inverse_income_with_pred

		mydebug && describe(agents.inverse_income_with_pred)

		if proptowage
			prop_vec = 1.0 ./ agents.inverse_income_with_pred
		else
			prop_vec = zeros(size(agents, 1)) .+ 1.0
		end

		agents[!, "a"]   .= theta_dict["alpha"]    .* prop_vec
		agents[!, "be"]  .= theta_dict["beta_e"]   .* prop_vec
		agents[!, "bl"]  .= theta_dict["beta_l"]   .* prop_vec
		agents[!, "sig"] .= theta_dict["sigma_dt"] .* sigma_adjust_vec
		agents[!, "mu"]  .= theta_dict["mu"]       .* sigma_adjust_vec

		# extensive margin decisions
        if haskey(theta_dict, "eta")
			agents[!, "delta_ext"].= theta_dict["delta_ext"] .* sigma_adjust_vec
			agents[!, "eta"]      .= theta_dict["eta"]       .* sigma_adjust_vec
			agents.eta   = convert.(Float64, agents.eta)
			agents.delta_ext = convert.(Float64, agents.delta_ext)
		end

		agents[!, "idx"]  = 1:size(agents, 1)

		delays0 = mean(hdgrid_ttimes .* 60.0 ./ agents.mean_km, dims=1)

		## expand
		repeat!(agents, inner=df_expansion_factor)

		## TODO: parameters for extensive margin

		## draw hastar's in long format (simulate)
		agents[!, "ha"] .= 0.0
		for i=1:size(agents, 1)
			if isnothing(rng)
				agents[i, "ha"] = sample(hagrid, Weights(hastars[agents[i, "idx"], 1, :]))
			else
				agents[i, "ha"] = sample(rng, hagrid, Weights(hastars[agents[i, "idx"], 1, :]))
			end
		end

		# bigmats that we reuse
		t_exact_mat 	  = agents.ha .- Matrix(transpose(hdgrid))
		log_t_exact_mat = log.(max.(0, t_exact_mat))

		BIGMAT["t_exact_mat"] = t_exact_mat
		BIGMAT["log_t_exact_mat"] = log_t_exact_mat

		# for faster computation, flip them
		BIGMAT["t_exact_mat_flip"] = Matrix(transpose(BIGMAT["t_exact_mat"]))
		BIGMAT["log_t_exact_mat_flip"] = Matrix(transpose(BIGMAT["log_t_exact_mat"]))
		# t_exact_mat 	= hagrid_mat 		 - hdgrid_mat
		# log_t_exact_mat = log.(max.(0, t_exact_mat))

	## Final cleaning
		agents.mean_km = convert.(Float64, agents.mean_km)
		agents.a   = convert.(Float64, agents.a)
		agents.be  = convert.(Float64, agents.be)
		agents.bl  = convert.(Float64, agents.bl)
		agents.sig = convert.(Float64, agents.sig)
		agents.mu  = convert.(Float64, agents.mu)


	## Save agents to csv file
	if rootpath_output != ""
		agents_path = rootpath_output * "agents_df.csv"
		CSV.write(agents_path, agents)
	end

 	return agents, delays0, BIGMAT
end
