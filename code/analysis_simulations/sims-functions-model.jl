
function ts2musigma_v2(
		T::Matrix{Float64},
		S::Matrix{Float64})

	sigma2 = @. log(1.0 + S * S / (T * T))
	mu = @. log(T) - sigma2 * 0.5
	@. sigma2 = sqrt(sigma2)

	tmean_mat = @. exp(mu + 0.5 * sigma2 * sigma2)

	return mu, sigma2, tmean_mat
end

# check this works!
# mu = rand(5, 5)
# sigma = 2 * rand(5, 5) .^ 2
# x = randn(5, 5) .* 2
#
# t1 = cdf.(Normal.(mu,sigma), x)
# t2 = normcdf(x .- mu, sigma)

function normcdf(x_μ::Matrix{Float64}, σ::Matrix{Float64})
	return @. 0.5 * (1.0 + erf( (x_μ) / σ / sqrt(2) ))
end

function normcdf(x_μ::Float64, σ::Float64)
	return 0.5 * (1.0 + erf( (x_μ) / σ / sqrt(2) ))
end

function tt_lognormal_integrals(;
		t_exact_mat::Array{Float64, 2},
		log_t_exact_mat::Array{Float64, 2},
		tmu::Matrix{Float64},
		tsi::Matrix{Float64},
		tmean_mat::Matrix{Float64}
	)
	temp1 = transpose(t_exact_mat)
	temp2 = transpose(log_t_exact_mat)

	# CDF     = @. cdf(Normal(tmu, tsi)        , temp2)
    # INT_CDF = @. cdf(Normal(tmu + tsi^2, tsi), temp2)

	CDF     = @. normcdf(temp2 .- tmu, tsi)
    INT_CDF = @. normcdf(temp2 - tmu - tsi^2, tsi)

    EARLY   = @.        temp1 * CDF       - tmean_mat * INT_CDF
    LATE    = @. (-1) * temp1 * (1 - CDF) + tmean_mat * (1 - INT_CDF)

	return EARLY, LATE
end

"""
 Selfish choice probabilities of everyone
 Charges enter positively (use negative charges for tax)
 :param compute_logsum: 0 = no, 1 = all, 2 = only logsum (attention: mean)
"""
function make_choice_prbs(;
		agents_mean_km::Vector{Float64},
		agents_a::Vector{Float64},
		agents_be::Vector{Float64},
		agents_bl::Vector{Float64},
		agents_sig::Vector{Float64},
		ttimes::Matrix{Float64},
		dt_charges::Array{Float64, 2},
        t_exact_mat::Array{Float64, 2},
        log_t_exact_mat::Array{Float64, 2},
        compute_logsum::Int64=0,
		MAXDIFF::Float64=300.0,
		return_utility::Bool=false,
		debug=false)

    # individual specific travel time profiles
	times_sd = SD_fun(ttimes ./ agents_mean_km) .*
				agents_mean_km / 60.0  # in HOURS

    # convert to lognormal params
    tmu, tsi, tmean_mat = ts2musigma_v2(ttimes ./ 60.0, times_sd)

    """ Init BIG MATRICES tmean, EARLY, LATE, CHARGES """
	EARLY, LATE = tt_lognormal_integrals(
			t_exact_mat=t_exact_mat,
			log_t_exact_mat=log_t_exact_mat,
			tmu=tmu,
			tsi=tsi,
			tmean_mat=tmean_mat
		)

    """ compute utility """
    debug && print("choice - utility")
    # assuming single route, utility = - alpha ET - beta_e .* E|| - beta_l .* E||
    dt_choice = @. (-1) * agents_a * tmean_mat +
                   (-1) * agents_be * EARLY +
                   (-1) * agents_bl * LATE +
                   dt_charges

	debug && all(dt_choice .<= 0.0)
	debug && describe(dt_choice[:])

	if return_utility
		return dt_choice
	end

    @. dt_choice = dt_choice / agents_sig

	# at most -300.0 relative to max -> otherwise we get Inf when we exponentiate
	# for detour cases -> choose the max among both routes

	# max entry along hd dimension (use the shape of eu_mat, which is n_resp x 1 x n_ha)
	# maximum!(eu_mat, dt_choice)
	eu_mat = maximum(dt_choice, dims=2)

	# main route utility -> bound from below relative to max
	@. dt_choice = max(dt_choice, eu_mat - MAXDIFF)

	# take out mean level to make it easier to exponentiate
	# save this -> this is the offset
	# mean!(eu_mat, dt_choice)
	eu_mat = mean(dt_choice, dims=2)

	if debug
		Plots.plot(hdgrid, dt_choice[129, :, 1]) |> display
	end

	# exponentiate
    @. dt_choice = exp(dt_choice - eu_mat)

	# compute choice probabilities
	prob_sum_mat = sum(dt_choice, dims=2)
	@. dt_choice = dt_choice / prob_sum_mat

	# average paid charges
	avg_charges = sum(dt_charges .* dt_choice, dims=2)

	### Expected utility!
	# add back the offset from eu_mat
	# add back charges
	@. eu_mat = (eu_mat + log(prob_sum_mat)) *
				agents_sig +
				(-1) * avg_charges

	if compute_logsum == 0
		return dt_choice, 0

	elseif compute_logsum == 1
		return dt_choice, eu_mat

	else
		@assert compute_logsum == 2
        return mean(eu_mat)
	end

end



function nash_iteration(;
		agents::DataFrame,
		rt_params::Vector{Float64},
		delays0::Matrix{Float64},
		dt_charges::Matrix{Float64},
		idx1::Int64=-1, idx2::Int64=-1,
		BIGMAT::Dict{String, Array{Float64, N} where N},
		step_tol::Float64=1e-8,
		step_max::Int64=100,
		error_if_fail::Bool=false,
		prop_update_frac::Float64=1.0,
		ttimes_from_delay::Bool=true,
		rootpath_output::String="",
		compute_logsum::Int64=0,
		debug_level::Int64=0,
		debug_graph_period::Int64=50,
		debug_norm_period::Int64=50
	)

	t0 = time()

	# big mats
    n_agents = size(agents, 1)
	hdgrid =  BIGMAT["hdgrid"]
	ttimes = BIGMAT["ttimes"]
	km_left_small = BIGMAT["km_left_small"]
	ttimes_small = BIGMAT["ttimes_small"]
	km_mean_small = agents.mean_km[1:10:3040]

	# initialize travel times from "common" delay function
	if ttimes_from_delay
		# minutes
		@. ttimes = delays0 .* agents.mean_km
	end

	# time interval in minutes
	delta_t = (hdgrid[2] - hdgrid[1]) * 60.0

	choice_prbs_current = nothing
	choice_prbs_first = nothing

    step_norm = 1.0
    step_idx = 1
	all_norms = []
    while (step_norm > step_tol) && (step_idx < step_max)

		# compute optimal choices given current travel times
        choice_prbs, _ = make_choice_prbs(
					agents_mean_km=agents.mean_km,
					agents_a=agents.a,
					agents_be=agents.be,
					agents_bl=agents.bl,
					agents_sig=agents.sig,
					ttimes=ttimes,
					t_exact_mat=BIGMAT["t_exact_mat_flip"],
			        log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
                  	dt_charges=dt_charges)

		# rowwise norm in choice probs
		if (step_idx  == 1)
			step_norm = 1.0
			choice_prbs_current = copy(choice_prbs)
			choice_prbs_first = copy(choice_prbs)
		else

			step_norms = mapslices(norm, choice_prbs .- choice_prbs_current; dims=2)
	        step_norm = mean(step_norms)
			push!(all_norms, step_norm)

	        (debug_level > 0) && (step_idx % debug_norm_period == 0) && println("...", step_idx, " norm: ", step_norm)
		end

		# for long iterations, analyze norm
		# if step_idx == 500
		# 	plot(1:length(all_norms), log.(all_norms)) |> display
		# end

        if (debug_level == 2) && (step_idx % debug_graph_period == 0)
			# draw
			plot(hdgrid, sum(choice_prbs_first, dims=1) |> vec, color=:black)
			plot!(hdgrid, sum(choice_prbs_current, dims=1) |> vec, color=:red, right_margin=15mm)

			# plot!(twinx(), hdgrid, delays |> vec, color=:blue) |> display
			avg_delay = mean(ttimes ./ agents.mean_km, dims=1) |> vec
			plot!(twinx(), hdgrid, avg_delay, color=:blue, right_margin=15mm) |> display
		end

		# update choices towards new optimal choices
		@. choice_prbs_current = prop_update_frac  * choice_prbs +
					      (1.0 - prop_update_frac) * choice_prbs_current

        # travel times from current (updated) choice probabilities
		make_travel_times_density(
            choice_prbs=choice_prbs_current,
			rt_params=rt_params,
			extra_idx1=idx1,
			extra_idx2=idx2,
			km_mean=km_mean_small,
			km_left=km_left_small,
			ttimes_small=ttimes_small,
			ttimes=ttimes,
			delta_t=delta_t)

        (debug_level > 0) && print(".")

        step_idx += 1
	end # end while


	if (debug_level > 0)
		t1 = time()
		println("Time for Nash: ", floor((t1 - t0) * 100)/100, " seconds")
	end

	if (step_idx >= step_max)
		print("~~~ Failed - exceeded maximum iterations ", step_max)
		if error_if_fail
			error("~~~ Failed - exceeded maximum iterations ")
		end
	end

	# compute welfare
	logsum = nothing
	if compute_logsum >= 1
		myoutput = make_choice_prbs(
					agents_mean_km=agents.mean_km,
					agents_a=agents.a,
					agents_be=agents.be,
					agents_bl=agents.bl,
					agents_sig=agents.sig,
					ttimes=ttimes,
					t_exact_mat=BIGMAT["t_exact_mat_flip"],
					log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
					dt_charges=dt_charges,
					compute_logsum=compute_logsum)

		if compute_logsum == 1
			_, logsum = myoutput
		else
			@assert compute_logsum == 2
			logsum = myoutput
		end
 	end

	# save to file
	if rootpath_output != ""
		nash_output_path = rootpath_output * "nash_eqm.bson"
		BSON.bson(nash_output_path, Dict(
			:agents => agents,
			:rt_params => rt_params,
			:nash_choice_probs => choice_prbs_current,
			:nash_ttimes => ttimes,
			:nash_logsum => logsum
		))
	end

	if compute_logsum == 0
		return choice_prbs_current

	elseif compute_logsum == 1
		return choice_prbs_current, logsum

	else
		@assert compute_logsum == 2
		return logsum
	end

end

"""
	Compute partial-equilibrium marginal social cost of adding a trip between
	(the start of) departure time idx1 and (the end of) arrival time idx2
	- partial equilibrium
	1. compute new travel time profile (including old choices and this additional trip)
	2. compute expected utility (over alpha, beta_E, beta_L) for each departure time
	3. (do #2 at the old and new travel time profile)
	4. take difference weighted by baseline choice probabilities

	- note, dt_charges don't matter here (they drop out) because we use the
		baseline choice probabilities already

"""
function msc_partial(;
			choice_ne::Matrix{Float64},
			rt_params::Vector{Float64},
			idx1s::Vector{Int64}, idx2s::Vector{Int64},
			BIGMAT::Dict{String, Array{Float64}},
			delta_t::Float64,
			agents_mean_km::Vector{Float64},
			agents_mean_km_small::Vector{Float64},
			agents_a::Vector{Float64},
			agents_be::Vector{Float64},
			agents_bl::Vector{Float64},
			agents_sig::Vector{Float64},
			my_dt_charges::Matrix{Float64},
			debug::Bool=false
		)

msc_list = zeros(length(idx1s))

for i=1:length(idx1s)
	idx1 = idx1s[i]
	idx2 = idx2s[i]
	debug && println("processing ", idx1, ",", idx2)

	n_agents, n_hd = size(choice_ne)

	make_travel_times_density(
		choice_prbs=choice_ne,
		rt_params=rt_params,
		extra_idx1=idx1, extra_idx2=idx2,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes"],
		delta_t=delta_t
	)

	exp_utility_new = make_choice_prbs(
		agents_mean_km=agents_mean_km,
		agents_a=agents_a,
		agents_be=agents_be,
		agents_bl=agents_bl,
		agents_sig=agents_sig,
		ttimes=BIGMAT["ttimes"],
		t_exact_mat=BIGMAT["t_exact_mat_flip"],
		log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
		dt_charges=my_dt_charges,
		return_utility=true
		)

	debug && println("processing ", idx1, ",", idx2, " ... DONE")

	msc_list[i] = mean(sum(exp_utility_new .* choice_ne, dims=2))
end

return msc_list

end

"""
	Compute partial-equilibrium marginal social cost of adding a trip between
	(the start of) departure time idx1 and (the end of) arrival time

	- at the end of the trip, weight MSC(idx2-1) and MSC(idx2) based on fraction
		of that time period that's covered
	- note, dt_charges don't matter here (they drop out) because we use the
		baseline choice probabilities already

"""
function msc_partial_all(;
		agents::DataFrame,
		rt_params::Vector{Float64},
		BIGMAT::Dict{String, Array{Float64}},
		choice_ne::Matrix{Float64},
		my_dt_charges::Matrix{Float64},
		delta_t::Float64,
		run_parallel::Bool=true,
		debug::Bool=false,
	)

	n_agents, n_hd = size(choice_ne)
	agents_mean_km_small = agents.mean_km[1:10:n_agents]

	# benchmark utility in current equilibrium
	exp_utility = make_choice_prbs(
		agents_mean_km=agents.mean_km,
		agents_a=agents.a,
		agents_be=agents.be,
		agents_bl=agents.bl,
		agents_sig=agents.sig,
		ttimes=BIGMAT["ttimes"],
		t_exact_mat=BIGMAT["t_exact_mat_flip"],
		log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
		compute_logsum=1,
		dt_charges=my_dt_charges,
		return_utility=true
		)
	debug && describe(exp_utility[:])

## get all the i,j pairs needed (much less than n_hd^2)

	# @btime
	make_travel_times_density(
		choice_prbs=choice_ne,  # N_agents x n_hd
		rt_params=rt_params,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes"],
		delta_t=delta_t
	)

	# matrix
	jk_msc = zeros(Bool, n_hd, n_hd)
	for i=1:n_agents
		for j=1:n_hd
			k = j + Int64(floor(BIGMAT["ttimes"][i,j] / delta_t))
			k = min(109, k)
			jk_msc[j,k] = true

			# also add previous (for later interpolation)
			jk_msc[j, max(1, k-1)] = true
		end
	end

	# list
	# sum(jk_msc)
	jk_msc_list = []
	jk_msc_list_start = Int64[]
	jk_msc_list_end = Int64[]
	for j=1:n_hd
		for k=j:n_hd
			if jk_msc[j,k]
				push!(jk_msc_list, (j,k))
				push!(jk_msc_list_start, j)
				push!(jk_msc_list_end, k)
			end
		end
	end

	msc_partial_loaded = (my_choice_ne, myidx1s, myidx2s) -> msc_partial(
				choice_ne=my_choice_ne,
				rt_params=rt_params,
				idx1s=myidx1s, idx2s=myidx2s,
				BIGMAT=BIGMAT, delta_t=delta_t,
				agents_mean_km=agents.mean_km,
				agents_mean_km_small=agents_mean_km_small,
				agents_a=agents.a,
				agents_be=agents.be,
				agents_bl=agents.bl,
				agents_sig=agents.sig,
				my_dt_charges=my_dt_charges,
				debug=false
			)

	# partition jm_msc_list into batches
	chunk(arr, batch_size) = [arr[i:min(i + batch_size - 1, end)] for i in 1:batch_size:length(arr)]
	batch_size = 50
	jk_msc_list_start_batches = chunk(jk_msc_list_start, batch_size)
	jk_msc_list_end_batches = chunk(jk_msc_list_end, batch_size)

	# Run (in batches)
	if run_parallel
		mscs = pmap(
			idx -> msc_partial_loaded(
							choice_ne,
							jk_msc_list_start_batches[idx],
							jk_msc_list_end_batches[idx]),
		 	1:length(jk_msc_list_start_batches))
	else
		mscs = []
		for idx=1:length(jk_msc_list_start_batches)
			msc_batch = msc_partial_loaded(choice_ne,
							jk_msc_list_start_batches[idx],
							jk_msc_list_end_batches[idx])
			push!(mscs, msc_batch)
		end
	end

	# put together in a single list
	mscs = vcat(mscs...)
	debug && describe(mscs[:])

	# subtract expected utility without additional traffic
	mscs_scaled = (mscs .- mean(sum(exp_utility .* choice_ne, dims=2)))  * (n_agents - 1) * 1000.0
	debug && describe(mscs_scaled[:])

## Createa n_agents x n_hd matrix with MSC
	msc_square = zeros(Float64, n_hd, n_hd) .- 1.0
	# populate with existing
	for idx=1:length(jk_msc_list)
		j,k = jk_msc_list[idx]
		msc_square[j,k] = mscs_scaled[idx]
	end

	msc_indiv = zeros(Float64, n_agents, n_hd)
	for i=1:n_agents
		for j=1:n_hd
			k = j + Int64(floor(BIGMAT["ttimes"][i,j] / delta_t))
			k = min(109, k)

			# old (simple)
			# msc_indiv[i,j] = msc_square[j,k]

			# interpolate
			if k == 109
				msc_indiv[i,j] = msc_square[j,k]
			elseif k > j
				frac_density_k = BIGMAT["ttimes"][i,j] / delta_t - (k-j)
				if ~( 0.0 - eps() <= frac_density_k <= 1.0 + eps() )
					println(frac_density_k)
					println(BIGMAT["ttimes"][i,j] / delta_t)
					println(k)
					println(j)
					error("what?")
				end
				msc_indiv[i,j] = frac_density_k * msc_square[j,k] +
						(1.0 - frac_density_k) *  msc_square[j, max(1, k-1)]

			else
				frac_density_k = BIGMAT["ttimes"][i,j] / delta_t - (k-j)
				if ~( 0.0 - eps() <= frac_density_k <= 1.0 + eps() )
					println(frac_density_k)
					println(BIGMAT["ttimes"][i,j] / delta_t)
					println(k)
					println(j)
					error("what?")
				end
				msc_indiv[i,j] = frac_density_k * msc_square[j,k]
			end

			if msc_square[j,k] == -1.0
				println(i,",",j,",",k)
			end
		end
	end

	return msc_indiv
end

"""
	utility optimization function to choose how much to adjust charges in the
		 direction of the partial equilibrium msc
"""
function my_line_search(;
		nash_fn_loaded::Function,
		dt_charges::Matrix{Float64},
		old_dt_charges::Matrix{Float64},
		msc::Matrix{Float64},
		logsum_current::Float64,
		adjust_factor0::Float64)

	if adjust_factor0 <= 1.1e-4
		# success flag, factor, next initial factor
		return false, 0.0, 0.0
	end

	adjust_grid = (0.1:0.1:0.9) .* adjust_factor0
	println("... optimal factor trying: ", adjust_factor0)
	logsum_nexts = zeros(length(adjust_grid))
	for idx=1:length(adjust_grid)
		adjust_factor = adjust_grid[idx]

		@. dt_charges = adjust_factor * msc + (1 - adjust_factor) * old_dt_charges
		dt_charges = dt_charges .- dt_charges[:, 1]

		logsum_nexts[idx] = nash_fn_loaded(dt_charges)
	end

	idx = argmax(logsum_nexts)

	if (logsum_nexts[idx] > logsum_current)
		# pick one of the factors we tried
		adjust_factor = adjust_grid[idx]

		# if we picked the largest, try larger factors next time
		adjust_factor_next = adjust_factor0
		if idx == length(logsum_nexts)
			# go up
			adjust_factor_next = min(1.0, 10.0 * adjust_factor0)
		end
		# success flag, factor, next initial factor
		return true, adjust_factor, adjust_factor_next

	else
		# if none of the options was better, try a finer search
		# recursive
		return my_line_search(
				nash_fn_loaded=nash_fn_loaded,
				dt_charges=dt_charges,
				old_dt_charges=old_dt_charges,
				msc=msc,
				logsum_current=logsum_current,
				adjust_factor0=0.1 * adjust_factor0)
	end
end

"""
"""
function socopt_iteration(;
			agents::DataFrame,
			rt_params::Vector{Float64},
			delays0::Matrix{Float64},
			BIGMAT::Dict{String, Array{Float64, N} where N},
			step_tol_out=1e-2,
			step_tol_in=1e-6,
			step_max=5000,
			prop_update_frac=1.0,
			adaptive_adjust=false,
			adjust_factor_initial=0.3,
			msc_run_parallel=true,
			rootpath_output::String="",
			debug_level::Int64=0)

    # Social Optimum Algorithm

	# prep
	n_agents, n_hd = size(BIGMAT["ttimes"])
	agents_mean_km_small = agents.mean_km[1:10:n_agents]
	hdgrid = BIGMAT["hdgrid"]
    delta_t = (hdgrid[2] - hdgrid[1]) * 60.0

	dt_choice_current = zeros(Float64, n_agents, n_hd)
	logsum_current = nothing
	logsum_current_vec = zeros(Float64, n_agents, 1)

	dt_choice_prev = zeros(Float64, n_agents, n_hd)
	logsum_prev = 0.0
	adjust_factor_prev = 1.0

	dt_choice_initial = zeros(Float64, n_agents, n_hd)
	ttimes_initial = zeros(Float64, n_agents, n_hd)
	logsum_initial = 0.0
	logsum_initial_vec = zeros(Float64, n_agents, 1)

	dt_charges = 0.0 .* copy(BIGMAT["ttimes"])

	ls_norm = 1000.0
	msc_norm = 1000.0
	found_optimal_logsum = false
    step_idx = 1

	progress_df = DataFrame(
		"index"		=> [],
		"sim_idx" 	=> [],
        "hdgrid" 	=> [],
        "avg_msc" 		=> [],
        "avg_delays" 	=> [],
		"instant_delays" 	=> [],
        "vol_dep" 		=> [],
		"density" 		=> [],
        "logsum" 	=> [],
        "time" 		=> [],
	)

	## Step 0. Nash equilibrium under NO congestion
	println("first computing Nash without congestion:")
	rt_params_nocon = copy(rt_params)
	rt_params_nocon[2] = 0.0
	dt_choice_nocon, logsum_nocon = nash_iteration(
				agents=agents,
				rt_params=rt_params_nocon,
				delays0=delays0,
				ttimes_from_delay=true,
				dt_charges=dt_charges,
				BIGMAT=BIGMAT,
				compute_logsum=1,
				debug_level=1,
				error_if_fail=true,
				prop_update_frac=prop_update_frac,
				step_tol=step_tol_in,
				step_max=step_max)

    # while ((ls_norm > step_tol_out) || (msc_norm > 1.0)) && ~found_optimal_logsum && (step_idx < step_max)
	while (ls_norm > step_tol_out) && ~found_optimal_logsum && (step_idx < step_max)

        println("\nNEW SIM: Iteration ", step_idx)

		# Step 1: find Nash equilibrium with current charges
        output = nash_iteration(
					agents=agents,
					rt_params=rt_params,
					delays0=delays0,
					ttimes_from_delay=(step_idx == 1),  # only first time
					dt_charges=dt_charges,
					BIGMAT=BIGMAT,
					compute_logsum=1,
					debug_level=1,
					error_if_fail=true,
					prop_update_frac=prop_update_frac,
					step_tol=step_tol_in,
					step_max=step_max)

		# Step 2. print stats
			dt_choice_current .= output[1]
			logsum_current_vec .= output[2]
	        logsum_current = mean(output[2])

	        avg_time = sum(dt_choice_current .* BIGMAT["ttimes"], dims=2)
	        avg_time = mean(avg_time)
	        (debug_level > 0) && println("...average trip duration (new delays): ", avg_time)

	        avg_char = sum(dt_choice_current .* dt_charges, dims=2)
	        avg_char = mean(avg_char)
	        (debug_level > 0) && println("...average charges paid : ", avg_char)

		# Step 3. logsum norm
		if step_idx > 1
			ch_norm = norm(
						mean(dt_choice_prev, dims=1) .-
                        mean(dt_choice_current, dims=1))

			ls_norm = abs(logsum_prev - logsum_current)

			(debug_level > 0) && println("...EU initial: ", logsum_initial)
			# (debug_level > 0) && println("...EU current: ", logsum_current)
            (debug_level > 0) && println("...EU improve: ", logsum_current - logsum_initial)
			(debug_level > 0) && println("...logsum norm: ", ls_norm)
			(debug_level > 0) && println("...choice norm: ", ch_norm)
        else
            @. dt_choice_initial = dt_choice_current
			@. ttimes_initial = copy(BIGMAT["ttimes"])
            logsum_initial = logsum_current
			logsum_initial_vec .= logsum_current_vec
			(debug_level > 0) && println("...EU initial: ", logsum_initial)
		end

        # Step 4. compute marginal social cost (partial equilibrium)
        (debug_level > 0) && println("...computing the new social cost and updating taxes")
        t0 = time()
		msc_indiv = msc_partial_all(
							agents=agents,
							rt_params=rt_params,
							BIGMAT=BIGMAT,
							choice_ne=dt_choice_current,
							my_dt_charges=dt_charges,
							delta_t=delta_t,
							debug=false, run_parallel=msc_run_parallel
						)
        println("... marginal sw time: ", 0.1 * floor((time() - t0) * 10) )

		### MSC norm
		msc_norm = abs.(mean(msc_indiv, dims=1) .-
					mean(dt_charges, dims=1)) |> maximum
		(debug_level > 0) && println("...msc    norm: ", msc_norm)

	# only update charges if haven't already converged (so we maintain consistent dt_charges)
	# if ((ls_norm > step_tol_out) || (msc_norm > 1.0))
	if (ls_norm > step_tol_out)
        # Step 5. adjust the charges based on msc
        if adaptive_adjust
			t0 = time()

			old_dt_charges = copy(dt_charges)

			nash_fn_loaded = my_dt_charges -> nash_iteration(
					agents=agents, rt_params=rt_params, delays0=delays0, ttimes_from_delay=false,
					dt_charges=my_dt_charges, # <----------------------------------------------
					BIGMAT=BIGMAT, compute_logsum=2, debug_level=0, error_if_fail=true,
					step_tol=step_tol_in, step_max=step_max, prop_update_frac=prop_update_frac)

			success_flag, adjust_factor, adjust_factor_prev = my_line_search(
					nash_fn_loaded=nash_fn_loaded,
					dt_charges=dt_charges,
					old_dt_charges=old_dt_charges,
					msc=msc_indiv,
					logsum_current=logsum_current,
					adjust_factor0=adjust_factor_prev)

			@. dt_charges = adjust_factor * msc_indiv + (1 - adjust_factor) * old_dt_charges
			dt_charges = dt_charges .- dt_charges[:, 1]

			if ~success_flag
				found_optimal_logsum = true
			end

			println("... optimal factor: ", adjust_factor ," time ", 0.1 * floor((time() - t0) * 10) )
		else
			### always use (large) step
	        adjust_factor = adjust_factor_initial
	        old_dt_charges = copy(dt_charges)
	        @. dt_charges = adjust_factor * msc_indiv + (1 - adjust_factor) * old_dt_charges
			# normalize for each agent
			dt_charges = dt_charges .- dt_charges[:, 1]
		end

		# Plot current choices and charges
        if (debug_level == 2)
			 plot(hdgrid, mean(dt_choice_initial, dims=1) |> vec, color=:black, label="initial")
			p1 = plot!(hdgrid, mean(dt_choice_current, dims=1) |> vec, color=:red, label="current", legend=:topleft,
				title="social optimum step " * string(step_idx))

			# new plot
			mean_charges = mean(-dt_charges, dims=1) |> vec
			mean_msc = mean(-msc_indiv, dims=1) |> vec
			xrng = (5,14)
			yrng = (0, maximum(vcat(mean_charges, mean_msc)))
			plot!(twinx(p1), hdgrid, mean_charges, xlims = xrng, ylims = yrng,
					color=:blue, label="mean charges")

			plot!(twinx(p1), hdgrid, mean_msc, xlims = xrng, ylims = yrng,
					color=:purple, linestyle=:dash, label="msc",
					legend=:topright) |> display

			# sleep(1)
		end
	end

        # Step 6. saving results to df
			densities, instant_delays =
			make_travel_times_density(
				choice_prbs=dt_choice_current,  # N_agents x n_hd
	            rt_params=rt_params,
				km_mean=agents_mean_km_small,
				km_left=BIGMAT["km_left_small"],
				ttimes_small=BIGMAT["ttimes_small"],
				ttimes=BIGMAT["ttimes"],
				delta_t=delta_t
			)

			rows_to_add = Dict(
				"index"		=> 1:length(hdgrid),
				"sim_idx" 	=> step_idx,  	# constant
				"hdgrid"	=> hdgrid |> vec,

				"avg_msc" 		=> mean(msc_indiv, dims=1) |> vec,
		        "avg_delays" 	=> mean(BIGMAT["ttimes"] ./ agents.mean_km, dims=1) |> vec,
				"instant_delays" 	=> instant_delays,
		        "vol_dep" 		=> sum(dt_choice_current, dims=1) |> vec,
				"density" 		=> densities,
				"logsum" 	=> logsum_current, # constant
				"time" 		=> avg_time, 	# constant
			)

			rows_to_add = DataFrame(rows_to_add)
			progress_df = vcat(progress_df, rows_to_add)

			# write at each step!
			progress_df_path = rootpath_output * "soc_opt_progress.csv"
			CSV.write(progress_df_path, progress_df)

        # Step 7. updates and changes
	        @. dt_choice_prev = dt_choice_current
	        logsum_prev = logsum_current
	        step_idx += 1
	end

    # """ DONE! """
    if step_idx >= step_max
        println("~~~ Failed - exceeded maximum iterations ", step_max)
        return nothing, nothing
	end

    println("~~~ Success Social Optimum ~~~")
    if debug_level == 2
		    plot(hdgrid, mean(dt_choice_initial, dims=1) |> vec, color=:black)
		scatter!(hdgrid, mean(dt_choice_initial, dims=1) |> vec, color=:black)
		   plot!(hdgrid, mean(dt_choice_current, dims=1) |> vec, color=:red)
		scatter!(hdgrid, mean(dt_choice_current, dims=1) |> vec,
				color=:red, title="social optimum (final)") |> display
	end

	if rootpath_output != ""
		progress_df_path = rootpath_output * "soc_opt_progress.csv"
		CSV.write(progress_df_path, progress_df)

		densities, instant_delays =
		make_travel_times_density(
			choice_prbs=dt_choice_current,  # N_agents x n_hd
            rt_params=rt_params,
			km_mean=agents_mean_km_small,
			km_left=BIGMAT["km_left_small"],
			ttimes_small=BIGMAT["ttimes_small"],
			ttimes=BIGMAT["ttimes"],
			delta_t=delta_t
		)

		so_output_path = rootpath_output * "so_eqm.bson"
		BSON.bson(so_output_path, Dict(
			:agents => agents,
			:rt_params => rt_params,

			:so_choice_probs => dt_choice_current,
			:so_ttimes => BIGMAT["ttimes"],
			:so_logsum => logsum_current,
			:so_logsum_vec => logsum_current_vec,
			:so_charges => dt_charges,
			:so_densities => densities,
			:so_instant_delays => instant_delays,

			:nash_choice_probs => dt_choice_initial,
			:nash_ttimes => ttimes_initial,
			:nash_logsum => logsum_initial,
			:nash_logsum_vec => logsum_initial_vec,

			:dt_choice_noconm  => dt_choice_nocon,
			:logsum_nocon => logsum_nocon

		))
	end

	# 1. optimal choices
	# 2. charges under social optimum
	# 3. logsum_1: social welfare under nash (no charges)
	# 4. logsum: social welfare under social optimum (with eqm optimal charges)
    return dt_choice_current, dt_charges, logsum_initial, logsum_current
end
