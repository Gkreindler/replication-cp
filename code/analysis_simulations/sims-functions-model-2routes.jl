

"""
 Selfish choice probabilities of everyone
 Charges enter positively (use negative charges for tax)
 :param compute_logsum: 0 = no, 1 = all, 2 = only logsum (attention: mean)
"""
function make_choice_prbs_2routes(;
		agents_mean_km::Vector{Float64},
		agents_a::Vector{Float64},
		agents_be::Vector{Float64},
		agents_bl::Vector{Float64},
		agents_sig::Vector{Float64},
		agents_mu::Vector{Float64},
		ttimes0::Matrix{Float64},
		ttimes1::Matrix{Float64},
		dt_charges0::Array{Float64, 2},
		dt_charges1::Array{Float64, 2},
        t_exact_mat::Array{Float64, 2},
        log_t_exact_mat::Array{Float64, 2},
        compute_logsum::Int64=0,
		MINU::Float64=-300.0,
		return_utility::Bool=false,
		debug=false)

    # individual specific travel time profiles
	# convert to lognormal params
	times0_sd = SD_fun(ttimes0 ./ agents_mean_km) .* agents_mean_km / 60.0  # in HOURS
    tmu0, tsi0, tmean_mat0 = ts2musigma_v2(ttimes0 ./ 60.0, times0_sd)

	times1_sd = SD_fun(ttimes1 ./ agents_mean_km) .* agents_mean_km / 60.0  # in HOURS
    tmu1, tsi1, tmean_mat1 = ts2musigma_v2(ttimes1 ./ 60.0, times1_sd)

    """ Init BIG MATRICES tmean, EARLY, LATE, CHARGES """
	EARLY0, LATE0 = tt_lognormal_integrals(
			t_exact_mat=t_exact_mat,
			log_t_exact_mat=log_t_exact_mat,
			tmu=tmu0,
			tsi=tsi0,
			tmean_mat=tmean_mat0)
	EARLY1, LATE1 = tt_lognormal_integrals(
			t_exact_mat=t_exact_mat,
			log_t_exact_mat=log_t_exact_mat,
			tmu=tmu1,
			tsi=tsi1,
			tmean_mat=tmean_mat1)

    """ compute utility """
    debug && print("choice - utility")
    # assuming single route, utility = - alpha ET - beta_e .* E|| - beta_l .* E||
    dt_choice0 = @. (-1) * agents_a * tmean_mat0 +
                   (-1) * agents_be * EARLY0 +
                   (-1) * agents_bl * LATE0 +
                   dt_charges0
	dt_choice1 = @. (-1) * agents_a * tmean_mat1 +
                   (-1) * agents_be * EARLY1 +
                   (-1) * agents_bl * LATE1 +
                   dt_charges1

	if return_utility
		return dt_choice0, dt_choice1
	end

    @. dt_choice0 = dt_choice0 / agents_sig
	@. dt_choice1 = dt_choice1 / agents_sig

	# at most -300.0 relative to max -> otherwise we get Inf when we exponentiate
	# for detour cases -> choose the max among both routes

	# max entry along hd dimension (use the shape of eu_mat, which is n_resp x 1 x n_ha)
	# maximum!(eu_mat, dt_choice)
	eu_mat = max.(
				maximum(dt_choice0, dims=2),
				maximum(dt_choice1, dims=2),
				)

	# exponentiate
	@. dt_choice0 = exp( max(dt_choice0 - eu_mat, MINU))
	@. dt_choice1 = exp( max(dt_choice1 - eu_mat, MINU))

	# compute choice probabilities
	prob_sum_mat0 = sum(dt_choice0, dims=2)
	prob_sum_mat1 = sum(dt_choice1, dims=2)
	@. dt_choice0 = dt_choice0 / prob_sum_mat0
	@. dt_choice1 = dt_choice1 / prob_sum_mat1

	### Route choice (upper level)
	# add back the offset from eu_mat
	# the sigma_adjut cancels out!
	v0 = @. (eu_mat + log(prob_sum_mat0)) * agents_sig / agents_mu
	v1 = @. (eu_mat + log(prob_sum_mat1)) * agents_sig / agents_mu

	# route 1 minus route 0 utility difference
	@. v1 = exp(min(max(v1 - v0, MINU), -MINU))

	# route 1 choice probability
	rt_choice1 = @. v1 / (1.0 + v1)

	# update choice probabilities by route choice
	@. dt_choice0 = dt_choice0 * (1.0 - rt_choice1)
	@. dt_choice1 = dt_choice1 * rt_choice1

	# average paid charges
	avg_charges = sum(dt_charges0 .* dt_choice0, dims=2) .+
				  sum(dt_charges1 .* dt_choice1, dims=2)

	### Expected utility (overall, excluding charges)
	v_mat = @. (v0 + log(1.0 + v1)) * agents_mu - avg_charges

	if compute_logsum == 0
		return dt_choice0, dt_choice1, 0

	elseif compute_logsum == 1
		return dt_choice0, dt_choice1, v_mat
	else
		@assert compute_logsum == 2
        return mean(v_mat)
	end

end



function nash_iteration_2routes(;
		agents::DataFrame,
		rt0_params::Vector{Float64},
		rt1_params::Vector{Float64},
		delays0::Matrix{Float64},
		delays1::Matrix{Float64},
		dt_charges0::Matrix{Float64},
		dt_charges1::Matrix{Float64},
		idx1_r0::Int64=-1, idx2_r0::Int64=-1,
		idx1_r1::Int64=-1, idx2_r1::Int64=-1,
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

	ttimes0 = BIGMAT["ttimes0"]
	ttimes1 = BIGMAT["ttimes1"]

	km_left_small = BIGMAT["km_left_small"]
	ttimes_small = BIGMAT["ttimes_small"]
	km_mean_small = agents.mean_km[1:10:3040]

	# initialize travel times from "common" delay function
	if ttimes_from_delay
		# minutes
		@. ttimes0 = delays0 .* agents.mean_km
		@. ttimes1 = delays1 .* agents.mean_km
	end

	# time interval in minutes
	delta_t = (hdgrid[2] - hdgrid[1]) * 60.0

	dt_choice0_current = zeros(size(ttimes0))
	dt_choice1_current = zeros(size(ttimes0))
	dt_choice0_first = zeros(size(ttimes0))
	dt_choice1_first = zeros(size(ttimes0))

    step_norm = 1.0
    step_idx = 1
    while (step_norm > step_tol) && (step_idx < step_max)

		# compute optimal choices given current travel times
        dt_choice0, dt_choice1, _ = make_choice_prbs_2routes(
					agents_mean_km=agents.mean_km,
					agents_a=agents.a,
					agents_be=agents.be,
					agents_bl=agents.bl,
					agents_sig=agents.sig,
					agents_mu=agents.mu,
					ttimes0=ttimes0,
					ttimes1=ttimes1,
					t_exact_mat=BIGMAT["t_exact_mat_flip"],
			        log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
                  	dt_charges0=dt_charges0,
					dt_charges1=dt_charges1)

		# rowwise norm in choice probs
		if (step_idx  == 1)
			step_norm = 1.0
			dt_choice0_current .= dt_choice0
			dt_choice1_current .= dt_choice1
			dt_choice0_first .= dt_choice0
			dt_choice1_first .= dt_choice1
		else

			step_norms =
				mapslices(norm, dt_choice0 .- dt_choice0_current; dims=2) .+
				mapslices(norm, dt_choice1 .- dt_choice1_current; dims=2)
	        step_norm = mean(step_norms ./ 2.0)

	        (debug_level > 0) && (step_idx % debug_norm_period == 0) && println("...", step_idx, " norm: ", step_norm)
		end

        if (debug_level == 2) && (step_idx % debug_graph_period == 0)
			# draw
			 plot(hdgrid, sum(dt_choice0_first, dims=1) |> vec, color=:black, label="r0")
			plot!(hdgrid, sum(dt_choice1_first, dims=1) |> vec, color=:black, linestyle=:dash, label="r1")

			plot!(hdgrid, sum(dt_choice0, dims=1) |> vec, color=:red, label="r0")
			plot!(hdgrid, sum(dt_choice1, dims=1) |> vec, color=:red, linestyle=:dash, label="r1") |> display

			# plot!(twinx(), hdgrid, delays |> vec, color=:blue) |> display
			avg_delay0 = mean(ttimes0 ./ agents.mean_km, dims=1) |> vec
			avg_delay1 = mean(ttimes1 ./ agents.mean_km, dims=1) |> vec
			plot(hdgrid, avg_delay0, color=:blue, label="r0")
			plot!(hdgrid, avg_delay1, color=:blue, linestyle=:dash, label="r1") |> display
		end

		# update choices towards new optimal choices
		@. dt_choice0_current = prop_update_frac  * dt_choice0 +
					     (1.0 - prop_update_frac) * dt_choice0_current
		@. dt_choice1_current = prop_update_frac  * dt_choice1 +
					     (1.0 - prop_update_frac) * dt_choice1_current

        # travel times from current (updated) choice probabilities
		make_travel_times_density(
            choice_prbs=dt_choice0_current,
			rt_params=rt0_params,
			extra_idx1=idx1_r0,
			extra_idx2=idx2_r0,
			km_mean=km_mean_small,
			km_left=km_left_small,
			ttimes_small=ttimes_small,
			ttimes=ttimes0,
			delta_t=delta_t)
		make_travel_times_density(
            choice_prbs=dt_choice1_current,
			rt_params=rt1_params,
			extra_idx1=idx1_r1,
			extra_idx2=idx2_r1,
			km_mean=km_mean_small,
			km_left=km_left_small,
			ttimes_small=ttimes_small,
			ttimes=ttimes1,
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
		myoutput = make_choice_prbs_2routes(
					agents_mean_km=agents.mean_km,
					agents_a=agents.a,
					agents_be=agents.be,
					agents_bl=agents.bl,
					agents_sig=agents.sig,
					agents_mu=agents.mu,
					ttimes0=ttimes0,
					ttimes1=ttimes1,
					t_exact_mat=BIGMAT["t_exact_mat_flip"],
					log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
					dt_charges0=dt_charges0,
					dt_charges1=dt_charges1,
					compute_logsum=compute_logsum)

		if compute_logsum == 1
			_, _, logsum = myoutput
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
			:rt0_params => rt0_params,
			:rt1_params => rt1_params,
			:nash_dt_choice0 => dt_choice0_current,
			:nash_dt_choice1 => dt_choice1_current,
			:nash_ttimes0 => ttimes0,
			:nash_ttimes1 => ttimes1,
			:nash_logsum => logsum
		))
	end

	if compute_logsum == 0
		return dt_choice0_current, dt_choice1_current

	elseif compute_logsum == 1
		return dt_choice0_current, dt_choice1_current, logsum

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
function msc_partial_2routes(;
			dt_choice0_ne::Matrix{Float64},
			dt_choice1_ne::Matrix{Float64},
			rt0_params::Vector{Float64},
			rt1_params::Vector{Float64},
			# idx1s_r0::Vector{Int64}, idx2s_r0::Vector{Int64},
			# idx1s_r1::Vector{Int64}, idx2s_r1::Vector{Int64},
			idxs::Vector{Any},
			BIGMAT::Dict{String, Array{Float64}},
			delta_t::Float64,
			agents_mean_km::Vector{Float64},
			agents_mean_km_small::Vector{Float64},
			agents_a::Vector{Float64},
			agents_be::Vector{Float64},
			agents_bl::Vector{Float64},
			agents_sig::Vector{Float64},
			agents_mu::Vector{Float64},
			my_dt_charges0::Matrix{Float64},
			my_dt_charges1::Matrix{Float64},
			debug::Bool=false
		)

 # return zeros(length(idxs)), zeros(length(idxs))

msc0_list = zeros(length(idxs))
msc1_list = zeros(length(idxs))

for i=1:length(idxs)
	# idx1_r0 = idx1s_r0[i]
	# idx2_r0 = idx2s_r0[i]
	# idx1_r1 = idx1s_r1[i]
	# idx2_r1 = idx2s_r1[i]
	idx1_r0 = idxs[i][1]
	idx2_r0 = idxs[i][2]
	idx1_r1 = idxs[i][3]
	idx2_r1 = idxs[i][4]
	debug && println("processing ", idx1_r0, ",", idx2_r0,",", idx1_r1, ",", idx2_r1)

	n_agents, n_hd = size(dt_choice0_ne)

	# update both routes -- only one will have extra volume
	make_travel_times_density(
		choice_prbs=dt_choice0_ne,
		rt_params=rt0_params,
		extra_idx1=idx1_r0, extra_idx2=idx2_r0,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes0"],
		delta_t=delta_t
	)
	make_travel_times_density(
		choice_prbs=dt_choice1_ne,
		rt_params=rt1_params,
		extra_idx1=idx1_r1, extra_idx2=idx2_r1,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes1"],
		delta_t=delta_t
	)

	# get utility
	eu0_new, eu1_new = make_choice_prbs_2routes(
		agents_mean_km=agents_mean_km,
		agents_a=agents_a,
		agents_be=agents_be,
		agents_bl=agents_bl,
		agents_sig=agents_sig,
		agents_mu=agents_mu,
		ttimes0=BIGMAT["ttimes0"],
		ttimes1=BIGMAT["ttimes1"],
		t_exact_mat=BIGMAT["t_exact_mat_flip"],
		log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
		dt_charges0=my_dt_charges0,
		dt_charges1=my_dt_charges1,
		return_utility=true
		)

	debug && println("processing ", idx1_r0, ",", idx2_r0,",", idx1_r1, ",", idx2_r1, " DONE")

	msc0_list[i] = mean(sum(eu0_new .* dt_choice0_ne, dims=2))
	msc1_list[i] = mean(sum(eu1_new .* dt_choice1_ne, dims=2))
end

return msc0_list, msc1_list

end

"""
	Compute partial-equilibrium marginal social cost of adding a trip between
	(the start of) departure time idx1 and (the end of) arrival time

	- at the end of the trip, weight MSC(idx2-1) and MSC(idx2) based on fraction
		of that time period that's covered
	- note, dt_charges don't matter here (they drop out) because we use the
		baseline choice probabilities already

"""
function msc_partial_all_2routes(;
		agents::DataFrame,
		dt_choice0_ne::Matrix{Float64},
		dt_choice1_ne::Matrix{Float64},
		rt0_params::Vector{Float64},
		rt1_params::Vector{Float64},
		BIGMAT::Dict{String, Array{Float64}},
		my_dt_charges0::Matrix{Float64},
		my_dt_charges1::Matrix{Float64},
		delta_t::Float64,
		run_parallel::Bool=true,
		debug::Bool=false,
	)

	n_agents, n_hd = size(dt_choice0_ne)
	agents_mean_km_small = agents.mean_km[1:10:n_agents]

	# benchmark utility in current equilibrium
	eu0, eu1  = make_choice_prbs_2routes(
		agents_mean_km=agents.mean_km,
		agents_a=agents.a,
		agents_be=agents.be,
		agents_bl=agents.bl,
		agents_sig=agents.sig,
		agents_mu=agents.mu,
		ttimes0=BIGMAT["ttimes0"],
		ttimes1=BIGMAT["ttimes1"],
		t_exact_mat=BIGMAT["t_exact_mat_flip"],
		log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
		compute_logsum=1,
		dt_charges0=my_dt_charges0,
		dt_charges1=my_dt_charges1,
		return_utility=true
		)

## get all the i,j pairs needed (much less than n_hd^2)

	make_travel_times_density(
		choice_prbs=dt_choice0_ne,
		rt_params=rt0_params,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes0"],
		delta_t=delta_t
	)
	make_travel_times_density(
		choice_prbs=dt_choice1_ne,
		rt_params=rt1_params,
		km_mean=agents_mean_km_small,
		km_left=BIGMAT["km_left_small"],
		ttimes_small=BIGMAT["ttimes_small"],
		ttimes=BIGMAT["ttimes1"],
		delta_t=delta_t
	)

	# matrix
	jk_msc0 = zeros(Bool, n_hd, n_hd)
	jk_msc1 = zeros(Bool, n_hd, n_hd)
	for i=1:n_agents
		for j=1:n_hd
			k = j + Int64(floor(BIGMAT["ttimes0"][i,j] / delta_t))
			k = min(109, k)
			jk_msc0[j,k] = true
			jk_msc0[j, max(1, k-1)] = true # also add previous (for later interpolation)

			k = j + Int64(floor(BIGMAT["ttimes1"][i,j] / delta_t))
			k = min(109, k)
			jk_msc1[j,k] = true
			jk_msc1[j, max(1, k-1)] = true # also add previous (for later interpolation)
		end
	end

	# list
	# sum(jk_msc)
	jk_msc_list = []
	jk_msc_list_start = Int64[]
	jk_msc_list_end = Int64[]
	for j=1:n_hd
		for k=j:n_hd
			if jk_msc0[j,k]
				push!(jk_msc_list, [j,k, -1, -1])
			end
			if jk_msc1[j,k]
				push!(jk_msc_list, [-1, -1, j,k])
			end
		end
	end

	msc_partial_loaded =
		(my_dt_choice0_ne, my_dt_choice1_ne, myidxs_batch) ->
		msc_partial_2routes(
				dt_choice0_ne=my_dt_choice0_ne,
				dt_choice1_ne=my_dt_choice1_ne,
				rt0_params=rt0_params,
				rt1_params=rt1_params,
				idxs=myidxs_batch,
				BIGMAT=BIGMAT, delta_t=delta_t,
				agents_mean_km=agents.mean_km,
				agents_mean_km_small=agents_mean_km_small,
				agents_a=agents.a,
				agents_be=agents.be,
				agents_bl=agents.bl,
				agents_sig=agents.sig,
				agents_mu=agents.mu,
				my_dt_charges0=my_dt_charges0,
				my_dt_charges1=my_dt_charges1,
				debug=debug
			)

	# partition jm_msc_list into batches
	chunk(arr, batch_size) = [arr[i:min(i + batch_size - 1, end)] for i in 1:batch_size:length(arr)]
	batch_size = 100
	jk_msc_batches = chunk(jk_msc_list, batch_size)
	# jk_msc_list_start_batches = chunk(jk_msc_list_start, batch_size)
	# jk_msc_list_end_batches = chunk(jk_msc_list_end, batch_size)

	# Run (in batches)
	if run_parallel
		mscs = pmap(
			idx -> msc_partial_loaded(dt_choice0_ne, dt_choice1_ne,
									  jk_msc_batches[idx]),
		 	1:length(jk_msc_batches))
	else
		mscs = []
		for idx=1:length(jk_msc_batches)
			msc_batch = msc_partial_loaded(dt_choice0_ne, dt_choice1_ne,
							jk_msc_batches[idx]...)
			push!(mscs, msc_batch)
		end
	end

	# put together in a single list
	mscs0 = vcat([msc[1] for msc in mscs]...)
	mscs1 = vcat([msc[2] for msc in mscs]...)

	# subtract expected utility without additional traffic
	mscs0_scaled = (mscs0 .- mean(sum(eu0 .* dt_choice0_ne, dims=2))) * (n_agents - 1) * 1000.0
	mscs1_scaled = (mscs1 .- mean(sum(eu1 .* dt_choice1_ne, dims=2))) * (n_agents - 1) * 1000.0

## Createa n_agents x n_hd matrix with MSC
	msc0_square = zeros(Float64, n_hd, n_hd)
	msc1_square = zeros(Float64, n_hd, n_hd)
	msc0_square .= msc0_square .- 1.0
	msc1_square .= msc1_square .- 1.0

	# populate with existing
	for idx=1:length(jk_msc_list)
		j0, k0, j1, k1 = jk_msc_list[idx]
		if j0 > 0
			msc0_square[j0,k0] = mscs0_scaled[idx]
		else
			msc1_square[j1,k1] = mscs1_scaled[idx]
		end
	end

	# populate MSC for each individual x depature time x route
	msc0_indiv = zeros(Float64, n_agents, n_hd)
	msc1_indiv = zeros(Float64, n_agents, n_hd)
	for i=1:n_agents
		for j=1:n_hd

			### Route 0
			k = j + Int64(floor(BIGMAT["ttimes0"][i,j] / delta_t))
			k = min(109, k)

			# interpolate
			if k == 109
				msc0_indiv[i,j] = msc0_square[j,k]
			elseif k > j
				frac_density_k = BIGMAT["ttimes0"][i,j] / delta_t - (k-j)
				@assert 0.0 - eps() <= frac_density_k <= 1.0 + eps()
				msc0_indiv[i,j] = frac_density_k  * msc0_square[j,k] +
						   (1.0 - frac_density_k) * msc0_square[j, max(1, k-1)]

			else
				frac_density_k = BIGMAT["ttimes0"][i,j] / delta_t - (k-j)
				@assert 0.0 - eps() <= frac_density_k <= 1.0 + eps()
				msc0_indiv[i,j] = frac_density_k * msc0_square[j,k]
			end
			(msc0_square[j,k] == -1) &&  println(i,",",j,",",k)

			### ROUTE 1
			k = j + Int64(floor(BIGMAT["ttimes1"][i,j] / delta_t))
			k = min(109, k)

			# interpolate
			if k == 109
				msc1_indiv[i,j] = msc1_square[j,k]
			elseif k > j
				frac_density_k = BIGMAT["ttimes1"][i,j] / delta_t - (k-j)
				@assert 0.0 - eps() <= frac_density_k <= 1.0 + eps()
				msc1_indiv[i,j] = frac_density_k  * msc1_square[j,k] +
						   (1.0 - frac_density_k) * msc1_square[j, max(1, k-1)]

			else
				frac_density_k = BIGMAT["ttimes1"][i,j] / delta_t - (k-j)
				@assert 0.0 - eps() <= frac_density_k <= 1.0 + eps()
				msc1_indiv[i,j] = frac_density_k * msc1_square[j,k]
			end
			(msc1_square[j,k] == -1) &&  println(i,",",j,",",k)
		end
	end

	return msc0_indiv, msc1_indiv
end

"""
	utility optimization function to choose how much to adjust charges in the
		 direction of the partial equilibrium msc
"""
function my_line_search_2routes(;
		nash_fn_loaded::Function,
		dt_charges0::Matrix{Float64},
		dt_charges1::Matrix{Float64},
		old_dt_charges0::Matrix{Float64},
		old_dt_charges1::Matrix{Float64},
		msc0::Matrix{Float64},
		msc1::Matrix{Float64},
		logsum_current::Float64,
		adjust_factor0::Float64)

	if adjust_factor0 <= 1e-4
		# success flag, factor, next initial factor
		return false, 0.0, 0.0
	end

	adjust_grid = (0.1:0.1:0.9) .* adjust_factor0
	println("... optimal factor trying: ", adjust_factor0)
	logsum_nexts = zeros(length(adjust_grid))
	for idx=1:length(adjust_grid)
		adjust_factor = adjust_grid[idx]

		@. dt_charges0 = adjust_factor * msc0 + (1 - adjust_factor) * old_dt_charges0
		@. dt_charges1 = adjust_factor * msc1 + (1 - adjust_factor) * old_dt_charges1
		dt_charges0 = dt_charges0 .- dt_charges0[:, 1]
		dt_charges1 = dt_charges1 .- dt_charges1[:, 1]

		logsum_nexts[idx] = nash_fn_loaded(dt_charges0, dt_charges1)
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
		return my_line_search_2routes(
				nash_fn_loaded=nash_fn_loaded,
				dt_charges0=dt_charges0,
				dt_charges1=dt_charges1,
				old_dt_charges0=old_dt_charges0,
				old_dt_charges1=old_dt_charges1,
				msc0=msc0,
				msc1=msc1,
				logsum_current=logsum_current,
				adjust_factor0=0.1 * adjust_factor0)
	end
end

"""
"""
function socopt_iteration_2routes(;
			agents::DataFrame,
			rt0_params::Vector{Float64},
			rt1_params::Vector{Float64},
			delays0::Matrix{Float64},
			delays1::Matrix{Float64},
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

	dt_choice0_current = zeros(Float64, n_agents, n_hd)
	dt_choice1_current = zeros(Float64, n_agents, n_hd)
	logsum_current = nothing

	dt_choice0_prev = zeros(Float64, n_agents, n_hd)
	dt_choice1_prev = zeros(Float64, n_agents, n_hd)
	logsum_prev = 0.0
	adjust_factor_prev = 1.0

	dt_choice0_initial = zeros(Float64, n_agents, n_hd)
	dt_choice1_initial = zeros(Float64, n_agents, n_hd)
	ttimes0_initial = zeros(Float64, n_agents, n_hd)
	ttimes1_initial = zeros(Float64, n_agents, n_hd)
	logsum_initial = 0.0

	dt_charges0 = 0.0 .* copy(BIGMAT["ttimes0"])
	dt_charges1 = 0.0 .* copy(BIGMAT["ttimes0"])

	ls_norm = 1000.0
	msc_norm = 1000.0
	found_optimal_logsum = false
    step_idx = 1

	progress_df = DataFrame(
		"index"		=> [],
		"sim_idx" 	=> [],
        "hdgrid" 	=> [],
        "avg_msc0" 		=> [],
		"avg_msc1" 		=> [],
        "avg_delays0" 	=> [],
		"avg_delays1" 	=> [],
		"instant_delays0" 	=> [],
		"instant_delays1" 	=> [],
        "vol0_dep" 		=> [],
		"vol1_dep" 		=> [],
		"density0" 		=> [],
		"density1" 		=> [],
        "logsum" 	=> [],
        "time" 		=> []
	)

	## Step 0. Nash equilibrium under NO congestion
	println("first computing Nash without congestion:")
	rt_params_nocon = copy(rt0_params)
	rt_params_nocon[2] = 0.0
	dt_choice0_nocon, dt_choice1_nocon, logsum_nocon = nash_iteration_2routes(
				agents=agents,
				rt0_params=rt_params_nocon,
				rt1_params=rt_params_nocon,
				delays0=delays0,
				delays1=delays1,
				ttimes_from_delay=true,
				dt_charges0=dt_charges0,
				dt_charges1=dt_charges1,
				BIGMAT=BIGMAT,
				compute_logsum=1,
				debug_level=1,
				error_if_fail=true,
				prop_update_frac=prop_update_frac,
				step_tol=step_tol_in,
				step_max=step_max)

    while ((ls_norm > step_tol_out) || (msc_norm > 1.0)) && ~found_optimal_logsum && (step_idx < step_max)

        println("\nNEW SIM: Iteration ", step_idx)

		# Step 1: find Nash equilibrium with current charges
        output = nash_iteration_2routes(
					agents=agents,
					rt0_params=rt0_params,
					rt1_params=rt1_params,
					delays0=delays0,
					delays1=delays1,
					ttimes_from_delay=(step_idx == 1),  # only first time
					dt_charges0=dt_charges0,
					dt_charges1=dt_charges1,
					BIGMAT=BIGMAT,
					compute_logsum=1,
					debug_level=1,
					error_if_fail=true,
					prop_update_frac=prop_update_frac,
					step_tol=step_tol_in,
					step_max=step_max)

		# Step 2. print stats
			dt_choice0_current .= output[1]
			dt_choice1_current .= output[2]
	        logsum_current = mean(output[3])

	        avg_time = sum((
				dt_choice0_current .* BIGMAT["ttimes0"] .+
				dt_choice1_current .* BIGMAT["ttimes1"]), dims=2)
	        avg_time = mean(avg_time)
	        (debug_level > 0) && println("...average trip duration (new delays): ", avg_time)

	        avg_char = sum((
				dt_choice0_current .* dt_charges0 .+
				dt_choice1_current .* dt_charges1), dims=2)
	        avg_char = mean(avg_char)
	        (debug_level > 0) && println("...average charges paid : ", avg_char)

		# Step 3. logsum norm
		if step_idx > 1
			ch_norm = norm(mean(dt_choice0_prev, dims=1) .- mean(dt_choice0_current, dims=1)) .+
					  norm(mean(dt_choice1_prev, dims=1) .- mean(dt_choice1_current, dims=1))

			ls_norm = abs(logsum_prev - logsum_current)

			(debug_level > 0) && println("...EU initial: ", logsum_initial)
            (debug_level > 0) && println("...EU improve: ", logsum_current - logsum_initial)
			(debug_level > 0) && println("...logsum norm: ", ls_norm)
			(debug_level > 0) && println("...choice norm: ", ch_norm)
        else
            @. dt_choice0_initial = dt_choice0_current
			@. dt_choice1_initial = dt_choice1_current
			@. ttimes0_initial = copy(BIGMAT["ttimes0"])
			@. ttimes1_initial = copy(BIGMAT["ttimes1"])
            logsum_initial = logsum_current
			(debug_level > 0) && println("...EU initial: ", logsum_initial)
		end

        # Step 4. compute marginal social cost (partial equilibrium)
        (debug_level > 0) && println("...computing the new social cost and updating taxes")
        t0 = time()
		msc0_indiv, msc1_indiv = msc_partial_all_2routes(
							agents=agents,
							rt0_params=rt0_params,
							rt1_params=rt1_params,
							BIGMAT=BIGMAT,
							dt_choice0_ne=dt_choice0_current,
							dt_choice1_ne=dt_choice1_current,
							my_dt_charges0=dt_charges0,
							my_dt_charges1=dt_charges1,
							delta_t=delta_t,
							debug=false, run_parallel=msc_run_parallel
						)
        println("... marginal sw time: ", 0.1 * floor((time() - t0) * 10) )

		### MSC norm
		msc_norm = maximum(abs.(mean(msc0_indiv, dims=1) .- mean(dt_charges0, dims=1))) .+
				   maximum(abs.(mean(msc1_indiv, dims=1) .- mean(dt_charges1, dims=1)))
		(debug_level > 0) && println("...msc L2 norm: ", msc_norm)

	# only update charges if haven't already converged (so we maintain consistent dt_charges)
	if ((ls_norm > step_tol_out) || (msc_norm > 1.0))
        # Step 5. adjust the charges based on msc
        if adaptive_adjust
			t0 = time()

			old_dt_charges0 = copy(dt_charges0)
			old_dt_charges1 = copy(dt_charges1)

			nash_fn_loaded_2routes = (my_dt_charges0, my_dt_charges1) -> nash_iteration_2routes(
					agents=agents, rt0_params=rt0_params, rt1_params=rt1_params,
					delays0=delays0, delays1=delays1, ttimes_from_delay=false,
					dt_charges0=my_dt_charges0, # <----------------------------------------------
					dt_charges1=my_dt_charges1, # <----------------------------------------------
					BIGMAT=BIGMAT, compute_logsum=2, debug_level=0, error_if_fail=true,
					step_tol=step_tol_in, step_max=step_max, prop_update_frac=prop_update_frac)

			success_flag, adjust_factor, adjust_factor_prev = my_line_search_2routes(
					nash_fn_loaded=nash_fn_loaded_2routes,
					dt_charges0=dt_charges0,
					dt_charges1=dt_charges1,
					old_dt_charges0=old_dt_charges0,
					old_dt_charges1=old_dt_charges1,
					msc0=msc0_indiv,
					msc1=msc1_indiv,
					logsum_current=logsum_current,
					adjust_factor0=adjust_factor_prev)

			@. dt_charges0 = adjust_factor * msc0_indiv + (1 - adjust_factor) * old_dt_charges0
			@. dt_charges1 = adjust_factor * msc1_indiv + (1 - adjust_factor) * old_dt_charges1
			dt_charges0 = dt_charges0 .- dt_charges0[:, 1]
			dt_charges1 = dt_charges1 .- dt_charges1[:, 1]

			if ~success_flag
				found_optimal_logsum = true
			end

			println("... optimal factor: ", adjust_factor ," time ", 0.1 * floor((time() - t0) * 10) )
		else
			### always use (large) step
	        adjust_factor = adjust_factor_initial
			old_dt_charges0 = copy(dt_charges0)
 		    old_dt_charges1 = copy(dt_charges1)
			@. dt_charges0 = adjust_factor * msc0_indiv + (1 - adjust_factor) * old_dt_charges0
			@. dt_charges1 = adjust_factor * msc1_indiv + (1 - adjust_factor) * old_dt_charges1
			dt_charges0 = dt_charges0 .- dt_charges0[:, 1]
			dt_charges1 = dt_charges1 .- dt_charges1[:, 1]
		end

		# Plot current choices and charges
        if (debug_level == 2)
			plot(hdgrid, mean(dt_choice0_initial, dims=1) |> vec, color=:black, label="initial r0")
			plot!(hdgrid, mean(dt_choice1_initial, dims=1) |> vec, color=:black, linestyle=:dash, label="initial r1")

			plot!(hdgrid, mean(dt_choice0_current, dims=1) |> vec, color=:red, label="current r0")
			plot!(hdgrid, mean(dt_choice1_current, dims=1) |> vec, color=:red, linestyle=:dash, label="current r1") |> display

			# new plot
			plot(hdgrid, mean(-dt_charges0, dims=1) |> vec, color=:blue, label="charges r0")
			plot!(hdgrid, mean(-dt_charges1, dims=1) |> vec, color=:blue, linestyle=:dash, label="charges r1")

		    plot!(hdgrid, mean(-msc0_indiv, dims=1) |> vec, color=:purple, label="msc r0")
		    plot!(hdgrid, mean(-msc1_indiv, dims=1) |> vec, color=:purple, linestyle=:dash, label="msc r1") |> display
		end
	end

        # Step 6. saving results to df
			densities0, instant_delays0 = make_travel_times_density(
				choice_prbs=dt_choice0_current,  # N_agents x n_hd
	            rt_params=rt0_params,
				km_mean=agents_mean_km_small,
				km_left=BIGMAT["km_left_small"],
				ttimes_small=BIGMAT["ttimes_small"],
				ttimes=BIGMAT["ttimes0"],
				delta_t=delta_t)
			densities1, instant_delays1 = make_travel_times_density(
				choice_prbs=dt_choice1_current,  # N_agents x n_hd
	            rt_params=rt1_params,
				km_mean=agents_mean_km_small,
				km_left=BIGMAT["km_left_small"],
				ttimes_small=BIGMAT["ttimes_small"],
				ttimes=BIGMAT["ttimes1"],
				delta_t=delta_t)

			rows_to_add = Dict(
				"index"		=> 1:length(hdgrid),
				"sim_idx" 	=> step_idx,  	# constant
				"hdgrid"	=> hdgrid |> vec,

				"avg_msc0" 		=> mean(msc0_indiv, dims=1) |> vec,
				"avg_msc1" 		=> mean(msc1_indiv, dims=1) |> vec,
		        "avg_delays0" 	=> mean(BIGMAT["ttimes0"] ./ agents.mean_km, dims=1) |> vec,
				"avg_delays1" 	=> mean(BIGMAT["ttimes1"] ./ agents.mean_km, dims=1) |> vec,
				"instant_delays0" 	=> instant_delays0,
				"instant_delays1" 	=> instant_delays1,
		        "vol0_dep" 		=> sum(dt_choice0_current, dims=1) |> vec,
				"vol1_dep" 		=> sum(dt_choice1_current, dims=1) |> vec,
				"density0" 		=> densities0,
				"density1" 		=> densities1,
				"logsum" 	=> logsum_current, # constant
				"time" 		=> avg_time, 	# constant
			)

			rows_to_add = DataFrame(rows_to_add)
			progress_df = vcat(progress_df, rows_to_add)

			# write at each step!
			progress_df_path = rootpath_output * "soc_opt_progress.csv"
			CSV.write(progress_df_path, progress_df)

        # Step 7. updates and changes
	        @. dt_choice0_prev = dt_choice0_current
			@. dt_choice1_prev = dt_choice1_current
	        logsum_prev = logsum_current
	        step_idx += 1
	end

    # """ DONE! """
    if step_idx >= step_max
        println("~~~ Failed - exceeded maximum iterations ", step_max)
        return nothing, nothing
	end

    println("~~~ Success Social Optimum ~~~")

	if rootpath_output != ""
		progress_df_path = rootpath_output * "soc_opt_progress.csv"
		CSV.write(progress_df_path, progress_df)

		densities0, instant_delays0 = make_travel_times_density(
			choice_prbs=dt_choice0_current,  # N_agents x n_hd
			rt_params=rt0_params,
			km_mean=agents_mean_km_small,
			km_left=BIGMAT["km_left_small"],
			ttimes_small=BIGMAT["ttimes_small"],
			ttimes=BIGMAT["ttimes0"],
			delta_t=delta_t)
		densities1, instant_delays1 = make_travel_times_density(
			choice_prbs=dt_choice1_current,  # N_agents x n_hd
			rt_params=rt1_params,
			km_mean=agents_mean_km_small,
			km_left=BIGMAT["km_left_small"],
			ttimes_small=BIGMAT["ttimes_small"],
			ttimes=BIGMAT["ttimes1"],
			delta_t=delta_t)

		so_output_path = rootpath_output * "so_eqm.bson"
		BSON.bson(so_output_path, Dict(
			:agents => agents,
			:rt0_params => rt0_params,
			:rt1_params => rt1_params,

			:so_choice_probs0 => dt_choice0_current,
			:so_choice_probs1 => dt_choice1_current,
			:so_ttimes0 => BIGMAT["ttimes0"],
			:so_ttimes1 => BIGMAT["ttimes1"],
			:so_logsum => logsum_current,
			:so_charges0 => dt_charges0,
			:so_charges1 => dt_charges1,
			:so_densities0 => densities0,
			:so_densities1 => densities1,
			:so_instant_delays0 => instant_delays0,
			:so_instant_delays1 => instant_delays1,

			:nash_choice_probs0 => dt_choice0_initial,
			:nash_choice_probs1 => dt_choice1_initial,
			:nash_ttimes0 => ttimes0_initial,
			:nash_ttimes1 => ttimes1_initial,
			:nash_logsum => logsum_initial,

			:dt_choice0_noconm  => dt_choice0_nocon,
			:dt_choice1_noconm  => dt_choice1_nocon,
			:logsum_nocon => logsum_nocon

		))
	end

	# 1. optimal choices
	# 2. charges under social optimum
	# 3. logsum_1: social welfare under nash (no charges)
	# 4. logsum: social welfare under social optimum (with eqm optimal charges)
    # return dt_choice_current, dt_charges, logsum_initial, logsum_current
end
