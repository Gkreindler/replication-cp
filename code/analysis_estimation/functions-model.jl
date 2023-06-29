



## compute choice probabilities

function checkmat(theta::Vector{Float64}, mymat::Array{Float64})
	# check that there are nan and infinite
	if (sum(isnan.(mymat[:])) > 0) | (sum(isinf.(mymat[:])) > 0)
		error("choice probs have nan or inf for theta=", theta)
	end
end


## invert
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
					alg=:nnls, max_iter=1000)

		 # check that we don't deviate too much from
		 debug && if abs(sum(temp) - 1.0) > 0.01
			 # display(temp)
			 println("hastar inversion in line ", i, " PDF total mass not 1, =", sum(temp))
			 # @assert abs(sum(temp) - 1.0) < 0.01
		 end
		 hastars[i, 1, :] = temp ./ sum(temp)
	end
end

## invert to a single hastar
function invert_single_has_per_person(;dt_choice::Array{Float64, 3},
					dt_choice_pre::Matrix{Float64},
					hastars::Array{Float64, 3},
					pdf_weight=5.0, debug=false)
	### for each person, invert the observed distribution of departure times
	#   to get A SINGLE ideal arrival time
	# 	compute loss for each ideal arrival time and pick best
	# 	dt_choice[i, :, j] for the ith commuter, assuming h^A = hagrid[j]
	#

	n_ha::Int64 = size(hastars)[3]

	# reset all values to zero
	fill!(hastars, 0.0)

	for i=1:size(dt_choice)[1]
		debug && println("processing ", i)

		# compare
		# (1) actual departure time distribution at baseline (dt_choice_pre)
		# with
		# (2) predicted departure time distribution given a single h^A
		# 	(given by dt_choice_free[i, :, idx_ha])
		loss = dt_choice[i, :, :] .- dt_choice_pre[i, :]

		# for each h^A, the loss is sum of square deviation in pdf
		loss_norm = sum(loss .^ 2, dims=1) |> vec

		# get best h^A
		idx_ha = argmin(loss_norm)

		# set h^A to full probability on this value
		hastars[i, 1, idx_ha] = 1.0
	end
end

## invert to a single hastar
function invert_single_has(;dt_choice::Array{Float64, 3},
					dt_choice_pre::Matrix{Float64},
					hastars::Array{Float64, 3},
					pdf_weight=5.0, debug=false)
	### invert the observed distribution of departure times to a unique number
	#   to get A SINGLE ideal arrival time
	# 	compute loss for each ideal arrival time and pick best
	# 	dt_choice[i, :, j] for the ith commuter, assuming h^A = hagrid[j]
	#

	n_ha::Int64 = size(hastars)[3]

	# reset all values to zero
	fill!(hastars, 0.0)

	# compare
	# (1) avg actual departure time distribution at baseline (dt_choice_pre)
	# term is n_hd vec  -> broadcasted as n_hd x 1 matrix
	# with
	# (2) avg predicted departure time distribution given a single h^A
	# 	(given by mean(dt_choice_free[:, :, idx_ha], dims=1) )
	# term is n_hd x n_ha matrix
	# loss is n_hd x n_ha
	loss = mean(dt_choice, dims=1)[1, :, :] .- mean(dt_choice_pre, dims=1)[1, :]

	# for each h^A, the loss is sum of square deviation in pdf
	# this is 1 x n_ha
	loss_norm = sum(loss .^ 2, dims=1) |> vec

	# get best h^A
	idx_ha = argmin(loss_norm)

	# set h^A to full probability on this value
	hastars[:, 1, idx_ha] .= 1.0
end


## Model with static value of time choice
section_hasdt_staticvot = 1



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
				MINU::Float64=-300.0, debug::Bool=false)

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
	# dt_choice .= max.(dt_choice, eu_mat .- MAXDIFF)
	#
	# # take out mean level to make it easier to exponentiate
	# # save this -> this is the offset
	# mean!(eu_mat, dt_choice)
	#
	# if debug
	# 	Plots.plot(hdgrid, dt_choice[129, :, 1]) |> display
	# end

	# exponentiate
    @. dt_choice = exp( max(dt_choice - eu_mat, MINU) )

	# checkmat(theta, dt_choice)

	# compute choice probabilities
	sum!(prob_sum_mat, dt_choice)
	@. dt_choice = dt_choice / prob_sum_mat

	### Expected utility!
	# add back the offset from eu_mat
	@. eu_mat = (eu_mat + log(prob_sum_mat)) * sigma_adjust * sigma

	# if debug
	# 	Plots.plot(hagrid, prob_sum_mat[229, 1, 1], yaxis=:log) |> display
	# end

	# checkmat(theta, dt_choice)
end


function nestedlogit_r0(
				theta::Vector{Float64},
				eu_mat::Array{Float64, 3}, eu_mat_detour::Array{Float64, 3},
				sigma_adjust_detour::Array{Float64, 3},
				area_charges::Array{Float64,3},
				prob_route0::Array{Float64,3};
				n_resp_a::Int64, MAXDIFF=300.0)
	# detour (r=1) divided by typical route (r=0)

	alpha, beta_e, beta_l, sigma, mu, prob_respond_alpha = theta

	# prob_sum_mat's are alreayd expected utility (including offsets put back in)
	# area charge (positive) enters with "+" because this is detour - main
	@. prob_route0 = eu_mat_detour +
	 				 (-1) * eu_mat[1:n_resp_a, :, :] +
					 area_charges

	### remove "sigma_adjust" -> equivalent to saying mu_i is also proportional to KM
	@. prob_route0 = prob_route0 / sigma_adjust_detour / mu

	# AFTER power, censor at +/-300 to avoid +/- infinity when exp-ing
	@. prob_route0 = min(max(prob_route0, (-1) * MAXDIFF), MAXDIFF)
	@. prob_route0 = exp(prob_route0)

	@. prob_route0 = 1.0 / (1.0 + prob_route0)

end



## final departure time probabilities (including integrating over hastar)
function nestedlogit_integrate_hastar_dt(
			dt_choice::Array{Float64, 3}, dt_choice_detour::Array{Float64, 3},
			prob_route0::Array{Float64,3},
			hastars::Array{Float64, 3}, n_resp_a::Int64
			)
	# sum up over hastar
	return vcat(
		dropdims(sum(
			(dt_choice[1:n_resp_a, :, :] .* prob_route0 .+
			 dt_choice_detour .* (1.0 .- prob_route0)) .*
		hastars[1:n_resp_a, :, :], dims=3), dims=3),
		dropdims(sum(
			(dt_choice[(n_resp_a+1):end, :, :]) .*
		hastars[(n_resp_a+1):end, :, :], dims=3), dims=3)
	)
end


## final departurroute choice probabilities (integrating over hastar)
function nestedlogit_integrate_hastar_r0(
	prob_route0::Array{Float64,3},
	hastars::Array{Float64, 3}, n_resp_a::Int64)

	return vec(sum(prob_route0[:, 1, :] .* hastars[1:n_resp_a, 1, :], dims=2))
end


## compute with and without charges -- NESTED LOGIT MODEL

function solve_nestedlogit(;
			theta::Vector{Float64},
			param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}},
			BIGMAT::Dict{String, Array{Float64, 3}},
			CHARGEMAT::Dict{String, Array{Float64, 3}},
			COMPMAT::Dict{String, Array{Float64, N} where N},
			DATAMOMS::Dict{String, Any},
			debug::Bool=false
	)

	## unpack all parameters
	theta_dict = get_params_fn(theta, param_factors)

	## DT parameters
	alpha, beta_e, beta_l, sigma, mu, prob_respond_alpha = theta_dt =
		[theta_dict["alpha"], theta_dict["beta_e"], theta_dict["beta_l"],
				theta_dict["sigma_dt"], theta_dict["mu"], theta_dict["prob_respond_alpha"]]

	## unpack
	dt_choice_free = COMPMAT["dt_choice"]
	dt_choice_free_detour = COMPMAT["dt_choice_detour"]

	# TODO: come back to this and pre-allocate
	dt_choice_wcha = copy(COMPMAT["dt_choice"])
	dt_choice_wcha_detour = copy(COMPMAT["dt_choice_detour"])

	eu_mat_free = COMPMAT["eu_mat"]
	eu_mat_free_detour = COMPMAT["eu_mat_detour"]

	# TODO: come back to this and pre-allocate
	eu_mat_wcha = copy(COMPMAT["eu_mat"])
	eu_mat_wcha_detour = copy(COMPMAT["eu_mat_detour"])

	prob_sum_mat = COMPMAT["prob_sum_mat"]
	prob_sum_mat_detour = COMPMAT["prob_sum_mat_detour"]

	hastars = COMPMAT["hastars"]
	prob_route0_free = COMPMAT["prob_route0_free"]
	prob_route0_wcha = COMPMAT["prob_route0_wcha"]

	dt_choice_pre = DATAMOMS["dt_choice_pre"]

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])

	### step 1: compute DT choice probs WITHOUT charges
		zero_dt_charges        = repeat([0.0], outer=[1,1,1])
		zero_dt_charges_detour = repeat([0.0], outer=[1,1,1])

		zero_area_charge_mean_mat        = repeat([0.0], outer=[1,1,1])
		zero_area_charge_mean_mat_detour = repeat([0.0], outer=[1,1,1])

		### Compute dt choices for main route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat"],
				BIGMAT["EARLY"],
				BIGMAT["LATE"],
				BIGMAT["sigma_adjust"],
				CHARGEMAT["zero_dt_charges"],
				dt_choice_free, eu_mat_free, prob_sum_mat, debug=(debug > 0))
		### Compute dt choices for DETOUR route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat_detour"],
				BIGMAT["EARLY_detour"],
				BIGMAT["LATE_detour"],
				BIGMAT["sigma_adjust_detour"],
				CHARGEMAT["zero_dt_charges_detour"],
				dt_choice_free_detour, eu_mat_free_detour, prob_sum_mat_detour,
				debug=(debug > 0))

	### step 2: compute DT choice probs WITH charges
		### Main route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat"],
				BIGMAT["EARLY"],
				BIGMAT["LATE"],
				BIGMAT["sigma_adjust"],
				CHARGEMAT["dt_charges"],
				dt_choice_wcha, eu_mat_wcha, prob_sum_mat, debug=(debug > 0))
		### Detour route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat_detour"],
				BIGMAT["EARLY_detour"],
				BIGMAT["LATE_detour"],
				BIGMAT["sigma_adjust_detour"],
				CHARGEMAT["dt_charges_detour"],
				dt_choice_wcha_detour, eu_mat_wcha_detour, prob_sum_mat_detour,
				debug=(debug > 0))

	### step 3: route choice probabilities WITHOUT/WITH charges

		### Without charges
		nestedlogit_r0(theta_dt,
			eu_mat_free, eu_mat_free_detour,
			BIGMAT["sigma_adjust_detour"],
			CHARGEMAT["zero_area_charge_mean"],
			prob_route0_free, n_resp_a=n_resp_a)

		### WITH charges
		nestedlogit_r0(theta_dt,
			eu_mat_free, eu_mat_free_detour,
			BIGMAT["sigma_adjust_detour"],
			CHARGEMAT["area_charge_mean_detour"],
			prob_route0_wcha, n_resp_a=n_resp_a)


	### step 4: invert hastar distribution
		# in hastars we store the distribution of h_A*

		# optimal departure time mapping, without charges, including route choice probabilities (where applicable)

		dt_choice_mapping = copy(dt_choice_free)
		dt_choice_mapping[1:n_resp_a, :, :] =
			dt_choice_mapping[1:n_resp_a, :, :] .*         prob_route0_free .+
			dt_choice_free_detour               .* (1.0 .- prob_route0_free)

		invert_has(dt_choice=dt_choice_mapping, dt_choice_pre=dt_choice_pre, hastars=hastars)

		# DEBUG:
		# return hastars
		# 	Plots.plot(hagrid, hastars[139, 1, :]) |> display
		# 	error("boo")

	### step 5: compute DT choice probs integrating over hastar
		dt_choice_free_integrated = nestedlogit_integrate_hastar_dt(
				dt_choice_free, dt_choice_free_detour,
				prob_route0_free,
				hastars, n_resp_a)

		dt_choice_wcha_integrated = nestedlogit_integrate_hastar_dt(
				dt_choice_wcha, dt_choice_wcha_detour,
				prob_route0_free,
				hastars, n_resp_a)

		# also compute single-route version
		dt_choice_free_1route_integrated = sum(dt_choice_free .* hastars, dims=3)
		dt_choice_wcha_1route_integrated = sum(dt_choice_wcha .* hastars, dims=3)

	### step 6: AREA integrate along hastar dimension
		probs_route0_free_integrated = nestedlogit_integrate_hastar_r0(
				prob_route0_free,
				hastars, n_resp_a)

		probs_route0_wcha_integrated = nestedlogit_integrate_hastar_r0(
				prob_route0_wcha,
				hastars, n_resp_a)


	return Dict{String, Any}(
		"dt_free" => dt_choice_free_integrated,
		"dt_wcha" => dt_choice_wcha_integrated,
		"dt_free_1route" => dt_choice_free_1route_integrated,
		"dt_wcha_1route" => dt_choice_wcha_1route_integrated,
		"probs_route0_free" => probs_route0_free_integrated,
		"probs_route0_wcha" => probs_route0_wcha_integrated
	)
end

##
section_hasdt_dynamicvot = 1

# Dynamic VOTT+DT model: switching cost + nested logit

function checkmat(theta::Dict{String, Float64}, mymat::Array{Float64})
	# check that there are nan and infinite
	if (sum(isnan.(mymat[:])) > 0) | (sum(isinf.(mymat[:])) > 0)
		display(theta)
		error("choice probs have nan or inf for theta above")
	end
end

function nestedlogit_all3d(;
		theta_dict::Dict{String, Float64},
		e0::Array{Float64, 3}, e1::Array{Float64, 3},
		eu_mat_0::Array{Float64}, eu_mat_1::Array{Float64, 3},
		prob_sum_0::Array{Float64, 3}, prob_sum_1::Array{Float64, 3},
		v0::Array{Float64, 3}, v1::Array{Float64, 3},
		sigma_adjust::Array{Float64, 3},
		compute_exp_diff::Bool=true,
		MINU::Float64=-300.0
	)

## unpack
	sigma_dt = theta_dict["sigma_dt"]
	mu 		 = theta_dict["mu"]

## new version
	@. e0 /= sigma_adjust * sigma_dt
	@. e1 /= sigma_adjust * sigma_dt

	# overall maximum (over both routes and all dep times)
	maximum!(eu_mat_0, e0)
	maximum!(eu_mat_1, e1)
	@. eu_mat_0 = max(eu_mat_0, eu_mat_1)

	# minimum utility diff = MINU
	# @. e0 = max(e0 - eu_mat_0, MINU)
	# @. e1 = max(e1 - eu_mat_0, MINU)

	# DEBUG
	# if ~all( e0 .>= MINU)
	# 	println("STATS:")
	# 	println(sum(isnan.(e0)))
	# 	println(sum(isinf.(e0)))
	# 	println(sum(ismissing.(e0)))
	# 	println(minimum(e0))
	#
	# 	error("hm")
	# end

	# @. e0 = exp(e0)
	# @. e1 = exp(e1)

	@. e0 = exp( max(e0 - eu_mat_0, MINU) )
	@. e1 = exp( max(e1 - eu_mat_0, MINU) )

	# DEBUG:
	# @assert all(e0 .>= exp(MINU))

	# sum exponential
	sum!(prob_sum_0, e0)
	sum!(prob_sum_1, e1)

	# DEBUG:
	# @assert all(prob_sum_0 .> 0.0)
	# @assert all(prob_sum_1 .> 0.0)

## upper nest - route level
	### Expected utility (1 x n_ha)
	# add back the offset from eu_mat
	# the sigma_adjut cancels out!
	@. v0 = (eu_mat_0 + log(prob_sum_0)) * sigma_dt / mu
	@. v1 = (eu_mat_0 + log(prob_sum_1)) * sigma_dt / mu

	if compute_exp_diff
		# this is exp(e1-e0)
		@. v1 = exp(min(max(v1 - v0, MINU), -MINU))
	end

end

"""
Small utility function to adjust e0 and exp(e1-e0) from
	- route 0 (where -gamma is in e1)
	- route 1 (where -gamma is in e0)
"""
function update_from_route1(
		v0::Array{Float64}, v1::Array{Float64},
		gamma01::Float64,
		gamma10::Float64,
		mu::Float64,
		sigma_adjust::Array{Float64},
		MINU::Float64=-300.0)

	@. v0 -= gamma01 / sigma_adjust / mu
	@. v1 *= exp( min(max( (gamma01 + gamma10) / sigma_adjust / mu, MINU ), -MINU) )

end

"""
	Given assumed value functions for next period (two numbers):
		- compute utility from each h,r combination | h_it^A
		- compute expected utility | h_it
		- compute expected utility = v0 and v1 at the start of this period

	exp_utility_0 includes time and expected schedule costs

	exp_utility_0 = n_hd x n_ha
"""
function onestep_eu(;
		theta_dict::Dict{String, Float64},
		v0s_assumed::Vector{Float64},
		v1s_assumed::Vector{Float64},
		e0::Array{Float64}, e1::Array{Float64},
		eu_mat_0::Array{Float64}, eu_mat_1::Array{Float64},
		prob_sum_0::Array{Float64}, prob_sum_1::Array{Float64},
		v0::Array{Float64}, v1::Array{Float64},
		v_mat::Array{Float64},
		prob0::Array{Float64},
		sigma_adjust::Array{Float64},
		hastar::Array{Float64},
		exp_utility_0::Array{Float64},
		exp_utility_1::Array{Float64},
		return_probs=false
		)

## unpack
	gamma01 = theta_dict["gamma"] * theta_dict["gamma_01_factor"] # switch cost r_{t-1}=0
	gamma10 = theta_dict["gamma"] # switch cost r_{t-1}=1
	delta 	 = theta_dict["delta"]
	# sigma_dt = theta_dict["sigma_dt"]
	mu 		 = theta_dict["mu"]

## Starting with r_{t-1} = 0

	# net present value of each h x r combination
	@. e0 = exp_utility_0 + delta * v0s_assumed
	@. e1 = exp_utility_1 + delta * v1s_assumed - gamma01

	# only run once!
	nestedlogit_all3d(
		theta_dict=theta_dict,
		e0=e0, e1=e1,
		eu_mat_0=eu_mat_0, eu_mat_1=eu_mat_1,
		prob_sum_0=prob_sum_0, prob_sum_1=prob_sum_1,
		v0=v0, v1=v1,
		sigma_adjust=sigma_adjust
	)

## compute all "from route 0"

	# part 1 expected utility
	@. v_mat = (v0 + log(1 + v1)) * mu * sigma_adjust
	# integrate
	v0_integrated = sum(v_mat .* hastar, dims=3)[:,1,1]

	if return_probs
		# route 0 probability
		@. prob0 = 1.0 / (1.0 + v1)

		# pi(h | r , hA)
		@. e0 = e0 / prob_sum_0
		@. e1 = e1 / prob_sum_1

		# part 2 -- departure time probs (integrated over route and hastar)
		phd0 = sum((e0 .* prob0 .+ e1 .* (1.0 .- prob0)) .* hastar, dims=3)[:,:,1]

		# part 3
		prob_route0_from0 = sum(prob0 .* hastar, dims=3)[:,1,1]
		@. prob_route0_from0 = min(1.0, prob_route0_from0)
		@. prob_route0_from0 = max(0.0, prob_route0_from0)

		# DEBUG
		# @assert all(prob0 .<= 1.0)
		#
		# idx = findfirst(prob_route0_from0 .> 1.0)
		# println(hastar[idx, 1, :])
		# println(sum(hastar[idx, 1, :]))
		# hastar[idx, 1, :] |> describe
		# println(prob0[idx,:,:])
		# prob0[idx,1,:] |> describe
		#
		# println(prob_route0_from0[idx])
		#
		# @assert all(prob_route0_from0 .<= 1.0)

	end

## compute all "from route 0"

	## update v0 to v0 - (gamma / mu / sigma_adjut)
	## update v1 to v1 * exp(gamma / mu / sigma_adjut)
	update_from_route1(v0, v1, gamma01, gamma10, mu, sigma_adjust)

	# DEBUG
	# @assert all(v1 .>= 0.0)
	# @assert all(v1 .<= 1.0)

	# part 1 expected utility
	@. v_mat = (v0 + log(1 + v1)) * mu * sigma_adjust
	# integrate
	v1_integrated = sum(v_mat .* hastar, dims=3)[:,1,1]

	if return_probs
		# route 0 probability
		@. prob0 = 1.0 / (1.0 + v1)

		# DEBUG
		# @assert all(prob0 .>= 0.0)

		# already done
		# pi(h | r , hA)
		# @. e0 = e0 / prob_sum_0
		# @. e1 = e1 / prob_sum_1

		# part 2 -- departure time probs (integrated over route and hastar)
		phd1 = sum((e0 .* prob0 .+ e1 .* (1.0 .- prob0)) .* hastar, dims=3)[:,:,1]

		# part 3
		prob_route0_from1 = sum(prob0 .* hastar, dims=3)[:,1,1]

		return (
			v0s=v0_integrated,
			v1s=v1_integrated,
			phd0=phd0,
			phd1=phd1,
			p01=1 .- prob_route0_from0,
			p10=prob_route0_from1,
			stat0=nothing
			)
	end

	return (v0_integrated, v1_integrated)

end


function onestep_eu_route(;
		mu_adj::Matrix{Float64},
		delta_mu_adj::Matrix{Float64},
		exp_minus_gamma01_mu_adj::Matrix{Float64},
		exp_minus_gamma10_mu_adj::Matrix{Float64},
		v0_current::Vector{Float64}, v1_current::Vector{Float64},
		v0_new::Vector{Float64}, v1_new::Vector{Float64},
		ev0::Matrix{Float64}, ev1::Matrix{Float64},
		expu0::Matrix{Float64}, expu1::Matrix{Float64},
		hastar2d::Matrix{Float64},
		MINU::Float64=-300.0
		)

	@. ev0 = exp( min(max(delta_mu_adj * v0_current, MINU), -MINU) )
	@. ev1 = exp( min(max(delta_mu_adj * v1_current, MINU), -MINU) )

	# @. ev0 = exp(ev0)
	# @. ev1 = exp(ev1)

	# @. ev0 = exp(delta_mu_adj * v0_current)
	# @. ev1 = exp(delta_mu_adj * v1_current)

	v0_new .= sum(mu_adj .*
			log.(expu0 .* ev0 .+
			 	 expu1 .* ev1 .* exp_minus_gamma01_mu_adj) .* hastar2d, dims=2) |> vec

	v1_new .= sum(mu_adj .*
		 	log.(expu0 .* ev0 .* exp_minus_gamma10_mu_adj .+
			 	 expu1 .* ev1 ) .* hastar2d, dims=2) |> vec

	# # DEBUG
	# @assert sum(isnan.(hastar2d)) == 0
	# @assert sum(isnan.(exp_minus_gamma_mu_adj)) == 0
	# @assert sum(isnan.(delta_mu_adj)) == 0
	# @assert sum(isnan.(v0_current)) == 0
	# @assert sum(isnan.(ev0)) == 0
	# @assert sum(isnan.(ev1)) == 0
	# @assert sum(isnan.(v0_new)) == 0
	# @assert sum(isnan.(v1_new)) == 0

end

"""
	Take advantage of recursive property of nested logit
	1. compute two nr x nha matrices of eu across dep time (given r and hA)
	2. iterate route-level Bellman equation -- integrating only over hA
"""
function dynamic_vott_steady_state_fast(;
		theta_dict::Dict{String, Float64},
		hastar::Array{Float64},
		exp_utility_0::Array{Float64},
		exp_utility_1::Array{Float64},
		initial_conditions::Matrix{Float64},
		sigma_adjust::Array{Float64},
		norm_epsilon::Float64=1e-12, lazy_param::Float64=0.1,
		max_steps::Int64=-1, mydebug::Int64=0,
		MINU::Float64=-300.0
		)

## unpack
	gamma01 = theta_dict["gamma"] * theta_dict["gamma_01_factor"] # switch cost r_{t-1}=0
	gamma10 = theta_dict["gamma"] # switch cost r_{t-1}=1
	delta 	 = theta_dict["delta"]
	# sigma_dt = theta_dict["sigma_dt"]
	mu 		 = theta_dict["mu"]

## computation matrices
	n_resp_a = size(exp_utility_1,1)
	n_hd = size(exp_utility_1,2)
	n_ha = size(exp_utility_1,3)

	e0			= zeros(n_resp_a, n_hd, n_ha)
	e1			= zeros(n_resp_a, n_hd, n_ha)
	eu_mat_0	= zeros(n_resp_a, 1, n_ha)
	eu_mat_1	= zeros(n_resp_a, 1, n_ha)
	prob_sum_0	= zeros(n_resp_a, 1, n_ha)
	prob_sum_1	= zeros(n_resp_a, 1, n_ha)
	v0			= zeros(n_resp_a, 1, n_ha)
	v1			= zeros(n_resp_a, 1, n_ha)
	v_mat		= zeros(n_resp_a, 1, n_ha)
	prob0		= zeros(n_resp_a, 1, n_ha)

## Start with departure time expected utilities
	@. e0 = exp_utility_0
	@. e1 = exp_utility_1
	nestedlogit_all3d(;
		theta_dict=theta_dict,
		e0=e0, e1=e1,
		eu_mat_0=eu_mat_0, eu_mat_1=eu_mat_1,
		prob_sum_0=prob_sum_0, prob_sum_1=prob_sum_1,
		v0=v0, v1=v1,
		sigma_adjust=sigma_adjust,
		compute_exp_diff=false
	)

	# # # DEBUG
	# @assert sum(isnan.(v0)) == 0
	# @assert sum(isnan.(v1)) == 0

## Belman equation
#=

	V(0) = \int_h^A sigma_adjust * mu * log( exp(u0 + delta * V(0) /mu/sigma_adjust) +
											 exp(u1 + delta * V(1) /mu/sigma_adjust - gamma/mu/sigma_adjust) ) dF(h^A)

	V(1) = \int_h^A sigma_adjust * mu * log( exp(u0 + delta * V(0) /mu/sigma_adjust - gamma/mu/sigma_adjust) +
										   	 exp(u1 + delta * V(1) /mu/sigma_adjust) ) dF(h^A)

	where u0 is given by the variable v0 etc.
=#

	# initial conditions for value function
	v0_current = initial_conditions[:, 1] |> copy
	v1_current = initial_conditions[:, 2] |> copy

	norm_change = 1.0

	# 2d matrices
	expu0 = @. exp( min(max(v0[:, 1, :], MINU), -MINU) )
	expu1 = @. exp( min(max(v1[:, 1, :], MINU), -MINU) )

	# factors
	mu_adj = mu .* sigma_adjust[:,1,:]
	delta_mu_adj = (delta ./ mu ./ sigma_adjust[:,1,:])

	exp_minus_gamma01_mu_adj = (-1) .* gamma01 ./ mu ./ sigma_adjust[:,1,:]
	exp_minus_gamma01_mu_adj = min.(max.(exp_minus_gamma01_mu_adj, MINU), -MINU)
	exp_minus_gamma01_mu_adj = exp.(exp_minus_gamma01_mu_adj)

	exp_minus_gamma10_mu_adj = (-1) .* gamma10 ./ mu ./ sigma_adjust[:,1,:]
	exp_minus_gamma10_mu_adj = min.(max.(exp_minus_gamma10_mu_adj, MINU), -MINU)
	exp_minus_gamma10_mu_adj = exp.(exp_minus_gamma10_mu_adj)


	hastar2d = hastar[:, 1, :]
	ev0 = zeros(n_resp_a, n_ha)
	ev1 = zeros(n_resp_a, n_ha)

	v0_new = zeros(n_resp_a)
	v1_new = zeros(n_resp_a)

	# iterate value function until we converge
	i = 0
	while (norm_change > norm_epsilon) && ((max_steps < 0) || (i < max_steps))
		# (mydebug > 0) && println("iteration ", i, " norm ", norm_change)

		# @. ev0 = exp(delta_mu_adj * v0_current)
		# @. ev1 = exp(delta_mu_adj * v1_current)
		#
		# v0_new = sum(mu_adj .*
		# 		log(expu0 .* ev0 .+
		# 		 	expu1 .* ev1 .* exp_minus_gamma_mu_adj) .* hastar2d, dims=2) |> vec
		#
		# v1_new = sum(mu_adj .*
		# 	 	log(expu0 .* ev0 .* exp_minus_gamma_mu_adj .+
		# 		 	expu1 .* ev1 ) .* hastar2d, dims=2) |> vec

		onestep_eu_route(
				mu_adj=mu_adj,
				delta_mu_adj=delta_mu_adj,
				exp_minus_gamma01_mu_adj=exp_minus_gamma01_mu_adj,
				exp_minus_gamma10_mu_adj=exp_minus_gamma10_mu_adj,
				v0_current=v0_current, v1_current=v1_current,
				v0_new=v0_new, v1_new=v1_new,
				ev0=ev0, ev1=ev1,
				expu0=expu0, expu1=expu1,
				hastar2d=hastar2d
				)

		norm_change = maximum(
				abs.(v0_new .- v0_current) .+
				abs.(v1_new .- v1_current)
				 )
		# temp = abs.(v0_new .- v0_current) .+ abs.(v1_new .- v1_current)
		# n_c = sum(temp .< norm_epsilon)
		# println(i, " n converged ", n_c, " out of ", length(temp))

		v0_current .= lazy_param * v0_new + (1.0 - lazy_param) * v0_current
		v1_current .= lazy_param * v1_new + (1.0 - lazy_param) * v1_current

		i += 1
	end

	# # DEBUG
	# # (mydebug > 0) &&
	# println("DONE after iterations: ", i)
	# println("DONE after iterations: ", i)
	# println("DONE after iterations: ", i)

	# error if went over the max number of steps
	if (max_steps > 0) && (i >= max_steps)
		# display(plot(1:max_steps, allv0s))
		println("error for vi=", myvi)
		println("v0_current=", v0_current)
		println("v1_current=", v1_current)
		error("above max number of steps")
	end

	# # DEBUG
	# println("v0_current=", v0_current[46])
	# println("v1_current=", v1_current[46])

	## compute probabilities
	# phd0 = Pr(h| r_{it-1} = 0)
	# phd1 = Pr(h| r_{it-1} = 1)
	# p01 = Pr(r=1|r_{it-1} = 0)
	# p10 = Pr(r=0|r_{it-1} = 1)
	_, _, phd0, phd1, p01, p10, _ =
	onestep_eu(
			theta_dict=theta_dict,
			v0s_assumed=v0_current, v1s_assumed=v1_current,
			e0=e0, e1=e1, eu_mat_0=eu_mat_0, eu_mat_1=eu_mat_1,
			prob_sum_0=prob_sum_0, prob_sum_1=prob_sum_1,
			v0=v0, v1=v1, v_mat=v_mat, prob0=prob0,
			sigma_adjust=sigma_adjust,
			hastar=hastar,
			exp_utility_0=exp_utility_0,
			exp_utility_1=exp_utility_1,
			return_probs=true
			)

	# stationary distribution
	stat_probs_0 = p10 ./ (p01 .+ p10)

	# # DEBUG:
	# @assert all(p01 .>= 0.0)
	# @assert all(p10 .>= 0.0)
	# println(stat_probs_0[46])
	# println(v0_current[46])
	# println(v1_current[46])

	# v0s 			= n_resp_a vector
	# probs_01 		= n_resp_a vector
	# stat_probs_0 	= n_resp_a vector

	# phd0 = Pr(h| r_{it-1} = 0)
	# phd1 = Pr(h| r_{it-1} = 1)
	# p01 = Pr(r=1|r_{it-1} = 0)
	# p10 = Pr(r=0|r_{it-1} = 1)
	# stat0 = stationary distribution of route0
	return (
		v0s=v0_current, v1s=v1_current,
		p01=p01, p10=p10, stat0=stat_probs_0,
		phd0=phd0, phd1=phd1
		)
end


function route0_transition(;p0, prob_01, prob_10)
	# updates the probability to use route 0 (main route)

	return @. p0 * (1 - prob_01) + (1 - p0) * prob_10
end

"""
Simulate model: steady state + all weeks during experiment
"""
function dynamic_vott_experiment(;
	ss=nothing,
	theta_dict::Dict{String, Float64},
	BIGMAT::Dict{String, Array{Float64, 3}},
	CHARGEMAT::Dict{String, Array{Float64, 3}},
	DATAMOMS::Dict{String,Any},
	hastar::Array{Float64},
	initial_conditions::Array{Float64},
	response_factor=1.0, ### SET THIS TO ZERO FOR NO CHARGES EVER
	return_full_solutions=false,
	my_lazy_param=1.0, norm_epsilon=1e-12,
	max_steps=50000, # maximum number of steps to find solution
	mydebug=0
	)

## unpack
	n_hd = size(BIGMAT["EARLY"],2)

	myalpha = theta_dict["alpha"]
	beta_e = theta_dict["beta_e"]
	beta_l = theta_dict["beta_l"]

	# charges positive, so should enter utility of route 0 negatively
	area_charges = CHARGEMAT["area_charge_mean_detour"]

	sigma_adjust_detour = BIGMAT["sigma_adjust_detour"]

	sample_early = DATAMOMS["sample_early"]
	sample_late  = DATAMOMS["sample_late"]

	n_resp, n_resp_a, n_resp_a_ct23, n_resp_a_tr = DATAMOMS["n_resp_vec"]

	n_ha = size(hastar, 3)

## Make all parameters proportional to wages
	# If parameters are proportional to wages
	# equivalent to charges INVERSELY proportional to wages
	if haskey(theta_dict, "prop_to_wage") && theta_dict["prop_to_wage"] .> 0.0
		@assert theta_dict["prop_to_wage"] == 1.0
		response_factor_all = response_factor .* DATAMOMS["data_all"].inverse_income_with_pred
		response_factor_a   = response_factor_all[1:n_resp_a]
	else
		response_factor_all = response_factor
		response_factor_a   = response_factor
	end

## Step 0. Expected utilities
	exp_utility_0_free = @. (-1) * myalpha * BIGMAT["tmean_mat"] +
				   (-1) * beta_e * BIGMAT["EARLY"] +
				   (-1) * beta_l * BIGMAT["LATE"]

	exp_utility_1_free = @. (-1) * myalpha * BIGMAT["tmean_mat_detour"] +
				   (-1) * beta_e * BIGMAT["EARLY_detour"] +
				   (-1) * beta_l * BIGMAT["LATE_detour"] +
				   theta_dict["fe_late"] * sample_late +
				   theta_dict["fe_early"] * sample_early

	exp_utility_0_wcha = @. (-1) * myalpha * BIGMAT["tmean_mat"] +
				   (-1) * beta_e * BIGMAT["EARLY"] +
				   (-1) * beta_l * BIGMAT["LATE"] +
				   (-1) * CHARGEMAT["dt_charges"] * response_factor_all

	exp_utility_1_wcha = @. (-1) * myalpha * BIGMAT["tmean_mat_detour"] +
				   (-1) * beta_e * BIGMAT["EARLY_detour"] +
				   (-1) * beta_l * BIGMAT["LATE_detour"] +
				   (-1) * CHARGEMAT["dt_charges_detour"] * response_factor_a +
				   theta_dict["fe_late"] * sample_late +
				   theta_dict["fe_early"] * sample_early


## Step 1: solve for the steady state
	# return format: named tuple:
	# (v0s, v1s, p01, p10, stat0, phd0, phd1)

	# phd0 = Pr(h| r_{it-1} = 0)
	# phd1 = Pr(h| r_{it-1} = 1)
	# p01 = Pr(r=1|r_{it-1} = 0)
	# p10 = Pr(r=0|r_{it-1} = 1)
	# stat0 = stationary distribution of route0

	if isnothing(ss)
		ss = dynamic_vott_steady_state_fast(
			theta_dict=theta_dict,
			hastar=hastar,
			exp_utility_0=exp_utility_0_free[1:n_resp_a, :, :],
			exp_utility_1=exp_utility_1_free,
			initial_conditions=initial_conditions,
			sigma_adjust=sigma_adjust_detour,
			norm_epsilon=norm_epsilon, lazy_param=my_lazy_param,
			max_steps=max_steps, mydebug=max(0, mydebug-1)
			)

		# update initial conditions (value fn of route 0 and of route 1)
		initial_conditions[:, 1] = ss.v0s
		initial_conditions[:, 2] = ss.v1s
	end

	# DEBUG:
	# ss.v0s[178] |> display
	# v0s=v0_current, v1s=v1_current,
	# p01=p01, p10=p10, stat0=stat_probs_0,
	# phd0=phd0, phd1=phd1

## Step 2: Solve for value fn in weeks 1-4, working backward, starting from week 4
	# separately for early/late
	e0			= zeros(n_resp_a, n_hd, n_ha)
	e1			= zeros(n_resp_a, n_hd, n_ha)
	eu_mat_0	= zeros(n_resp_a, 1, n_ha)
	eu_mat_1	= zeros(n_resp_a, 1, n_ha)
	prob_sum_0	= zeros(n_resp_a, 1, n_ha)
	prob_sum_1	= zeros(n_resp_a, 1, n_ha)
	v0			= zeros(n_resp_a, 1, n_ha)
	v1			= zeros(n_resp_a, 1, n_ha)
	v_mat		= zeros(n_resp_a, 1, n_ha)
	prob0		= zeros(n_resp_a, 1, n_ha)

	u0 			= zeros(n_resp_a, n_hd, n_ha)
	u1 			= zeros(n_resp_a, n_hd, n_ha)

## prepare utility by week

## week 4
    # Week 4 -- LATE -- no DT charges
	u0[sample_late, :, :] .= (exp_utility_0_free[1:n_resp_a, :, :] .-
					   theta_dict["fe_w4"] .-
					   area_charges .* response_factor_a)[sample_late, :, :]
	u1[sample_late, :, :]  .= exp_utility_1_free[sample_late, :, :]

	# Week 4 -- EARLY -- yes DT charges
	u0[sample_early, :, :] .= (exp_utility_0_wcha[1:n_resp_a, :, :] .-
						theta_dict["fe_w4"])[sample_early, :, :]
	u1[sample_early, :, :] .= exp_utility_1_wcha[sample_early, :, :]

	w4_soln = onestep_eu(
				theta_dict=theta_dict,
				v0s_assumed=ss.v0s,   # <---- start from steady state in week 5
				v1s_assumed=ss.v1s,
				e0=e0, e1=e1, eu_mat_0=eu_mat_0, eu_mat_1=eu_mat_1,
				prob_sum_0=prob_sum_0, prob_sum_1=prob_sum_1,
				v0=v0, v1=v1, v_mat=v_mat, prob0=prob0,
				sigma_adjust=sigma_adjust_detour,
				hastar=hastar[1:n_resp_a, :, :],
				exp_utility_0=u0,
				exp_utility_1=u1,
				return_probs=true
				)

	# # DEBUG:
	# @assert sum(isnan.(w4_soln.v0s)) == 0
	# @assert sum(isnan.(w4_soln.v1s)) == 0
	# @assert sum(isnan.(w4_soln.phd0)) == 0
	# @assert sum(isnan.(w4_soln.phd1)) == 0

## week 3
	w3_soln = onestep_eu(
				theta_dict=theta_dict,
				v0s_assumed=w4_soln.v0s,   # <---- start from week 4
				v1s_assumed=w4_soln.v1s,
				e0=e0, e1=e1, eu_mat_0=eu_mat_0, eu_mat_1=eu_mat_1,
				prob_sum_0=prob_sum_0, prob_sum_1=prob_sum_1,
				v0=v0, v1=v1, v_mat=v_mat, prob0=prob0,
				sigma_adjust=sigma_adjust_detour,
				hastar=hastar[1:n_resp_a, :, :],
				exp_utility_0=exp_utility_0_wcha[1:n_resp_a, :, :] .- theta_dict["fe_w3"],
				exp_utility_1=exp_utility_1_wcha,
				return_probs=true
				)
## week 2
	w2_soln = onestep_eu(
				theta_dict=theta_dict,
				v0s_assumed=w3_soln.v0s,   # <---- start from week 4
				v1s_assumed=w3_soln.v1s,
				e0=e0, e1=e1, eu_mat_0=eu_mat_0, eu_mat_1=eu_mat_1,
				prob_sum_0=prob_sum_0, prob_sum_1=prob_sum_1,
				v0=v0, v1=v1, v_mat=v_mat, prob0=prob0,
				sigma_adjust=sigma_adjust_detour,
				hastar=hastar[1:n_resp_a, :, :],
				exp_utility_0=exp_utility_0_wcha[1:n_resp_a, :, :] .- theta_dict["fe_w2"],
				exp_utility_1=exp_utility_1_wcha,
				return_probs=true
				)

## week 1
	# Week 1 -- LATE -- yes DT charges
	u0[sample_late, :, :] .= (exp_utility_0_wcha[1:n_resp_a, :, :] .-
					   theta_dict["fe_w1"])[sample_late, :, :]
	u1[sample_late, :, :]  .= exp_utility_1_wcha[sample_late, :, :]

	# Week 1 -- EARLY -- no DT charges
	u0[sample_early, :, :] .= (exp_utility_0_free[1:n_resp_a, :, :] .-
						theta_dict["fe_w1"] .-
 					   area_charges .* response_factor_a)[sample_early, :, :]
	u1[sample_early, :, :] .= exp_utility_1_free[sample_early, :, :]

	w1_soln = onestep_eu(
				theta_dict=theta_dict,
				v0s_assumed=w2_soln.v0s,   # <---- start from steady state in week 5
				v1s_assumed=w2_soln.v1s,
				e0=e0, e1=e1, eu_mat_0=eu_mat_0, eu_mat_1=eu_mat_1,
				prob_sum_0=prob_sum_0, prob_sum_1=prob_sum_1,
				v0=v0, v1=v1, v_mat=v_mat, prob0=prob0,
				sigma_adjust=sigma_adjust_detour,
				hastar=hastar[1:n_resp_a, :, :],
				exp_utility_0=u0,
				exp_utility_1=u1,
				return_probs=true
				)

## Step 3: Compute fraction in each state each week using transition probabilities

	## Late: stationary probability to use route 0 (main route)
		late_route0_w0 = ss.stat0[sample_late]

		# push forward using transition probabilities
		late_route0_w1 = route0_transition(p0=late_route0_w0, prob_01=w1_soln.p01[sample_late], prob_10=w1_soln.p10[sample_late])
		late_route0_w2 = route0_transition(p0=late_route0_w1, prob_01=w2_soln.p01[sample_late], prob_10=w2_soln.p10[sample_late])
		late_route0_w3 = route0_transition(p0=late_route0_w2, prob_01=w3_soln.p01[sample_late], prob_10=w3_soln.p10[sample_late])
		late_route0_w4 = route0_transition(p0=late_route0_w3, prob_01=w4_soln.p01[sample_late], prob_10=w4_soln.p10[sample_late])

	## Early: stationary probability to use route 0 (main route)
		early_route0_w0 = ss.stat0[sample_early]

		# push forward using transition probabilities
		early_route0_w1 = route0_transition(p0=early_route0_w0, prob_01=w1_soln.p01[sample_early], prob_10=w1_soln.p10[sample_early])
		early_route0_w2 = route0_transition(p0=early_route0_w1, prob_01=w2_soln.p01[sample_early], prob_10=w2_soln.p10[sample_early])
		early_route0_w3 = route0_transition(p0=early_route0_w2, prob_01=w3_soln.p01[sample_early], prob_10=w3_soln.p10[sample_early])
		early_route0_w4 = route0_transition(p0=early_route0_w3, prob_01=w4_soln.p01[sample_early], prob_10=w4_soln.p10[sample_early])

	### Combine back route probabilities -- by week
		probs_route0 = zeros(n_resp_a, 5)
		probs_route0[sample_late, 1] = late_route0_w0
		probs_route0[sample_late, 2] = late_route0_w1
		probs_route0[sample_late, 3] = late_route0_w2
		probs_route0[sample_late, 4] = late_route0_w3
		probs_route0[sample_late, 5] = late_route0_w4

		probs_route0[sample_early, 1] = early_route0_w0
		probs_route0[sample_early, 2] = early_route0_w1
		probs_route0[sample_early, 3] = early_route0_w2
		probs_route0[sample_early, 4] = early_route0_w3
		probs_route0[sample_early, 5] = early_route0_w4

	### Departure time probabilities
		# departure time probabilities (BEFORE the experiment)
		probs_hd_before = (@. ss.stat0 * ss.phd0 + (1 - ss.stat0) * ss.phd1)

		# Pr(h | week T) = Pr(r_{T-1} == 0) * Pr(h | r_{T-1} == 0) + Pr(r_{T-1} == 1) * Pr(h | r_{T-1} == 1)
		late_hd_w0 = probs_hd_before[sample_late, :]
		late_hd_w1 = @. late_route0_w0 * w1_soln.phd0[sample_late, :] + (1 - late_route0_w0) * w1_soln.phd1[sample_late, :]
		late_hd_w2 = @. late_route0_w1 * w2_soln.phd0[sample_late, :] + (1 - late_route0_w1) * w2_soln.phd1[sample_late, :]
		late_hd_w3 = @. late_route0_w2 * w3_soln.phd0[sample_late, :] + (1 - late_route0_w2) * w3_soln.phd1[sample_late, :]
		late_hd_w4 = @. late_route0_w3 * w4_soln.phd0[sample_late, :] + (1 - late_route0_w3) * w4_soln.phd1[sample_late, :]

		early_hd_w0 = probs_hd_before[sample_early, :]
		early_hd_w1 = @. early_route0_w0 * w1_soln.phd0[sample_early, :] + (1 - early_route0_w0) * w1_soln.phd1[sample_early, :]
		early_hd_w2 = @. early_route0_w1 * w2_soln.phd0[sample_early, :] + (1 - early_route0_w1) * w2_soln.phd1[sample_early, :]
		early_hd_w3 = @. early_route0_w2 * w3_soln.phd0[sample_early, :] + (1 - early_route0_w2) * w3_soln.phd1[sample_early, :]
		early_hd_w4 = @. early_route0_w3 * w4_soln.phd0[sample_early, :] + (1 - early_route0_w3) * w4_soln.phd1[sample_early, :]

	### Combine departure time probabilities (during the experiment)
		probs_hd_during = zeros(n_resp_a, n_hd)

		probs_hd_during[sample_late,  :] = @. (late_hd_w1  + late_hd_w2  + late_hd_w3)/3
		probs_hd_during[sample_early, :] = @. (early_hd_w2 + early_hd_w3 + early_hd_w4)/3

	##
	return probs_route0, probs_hd_before, probs_hd_during,
		ss, w1_soln, w2_soln, w3_soln, w4_soln
end


"""
Main model: dynamic route choice model with departure time choice.
- uses nested logit
- performs inversion using pre-experiment route choice frequencies in data
-
"""
function solve_dynamicvot(;
		theta::Vector{Float64},
		param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}},
		BIGMAT::Dict{String, Array{Float64, 3}},
		CHARGEMAT::Dict{String, Array{Float64, 3}},
		COMPMAT::Dict{String, Array{Float64, N} where N},
		DATAMOMS::Dict{String, Any},
		hastar_inversion_type="distribution",
		debug::Int64=0
	)

	if ~in(hastar_inversion_type, ["single", "single_per_person", "distribution", "given"])
		error("argument hastar_inversion_type should be `single`, `single_per_person`, `distribution` or `given`")
	end

## unpack parameters
	theta_dict = get_params_fn(theta, param_factors)

	## DT parameters
	# convert alpha from INR per minute to INR per hour
	theta_dt = [theta_dict["alpha"], theta_dict["beta_e"], theta_dict["beta_l"],
				theta_dict["sigma_dt"], 0.0, theta_dict["prob_respond_alpha"]]

## upack matrices

	dt_choice_free = COMPMAT["dt_choice"]
	dt_choice_free_detour = COMPMAT["dt_choice_detour"]

	dt_choice_wcha = COMPMAT["dt_choice_wcha"]
	# dt_choice_wcha_detour = copy(COMPMAT["dt_choice_detour"])


	eu_mat_free = COMPMAT["eu_mat"]
	eu_mat_free_detour = COMPMAT["eu_mat_detour"]

	eu_mat_wcha = COMPMAT["eu_mat_wcha"]
	# eu_mat_wcha_detour = copy(COMPMAT["eu_mat_detour"])


	prob_sum_mat = COMPMAT["prob_sum_mat"]
	prob_sum_mat_detour = COMPMAT["prob_sum_mat_detour"]


	dt_choice_pre = DATAMOMS["dt_choice_pre"]

	hastars = COMPMAT["hastars"]

	n_resp, n_resp_a, n_resp_a_ct23, n_resp_a_tr = DATAMOMS["n_resp_vec"]

    n_ha = size(dt_choice_free,3)
	n_hd = size(dt_choice_free,2)

	# sigma_adjust_detour = BIGMAT["sigma_adjust_detour"]

	sample_dynamic_vot = DATAMOMS["sample_dynamic_vot"]
	sample_dynamic_vot_full = vcat(sample_dynamic_vot, zeros(Bool, n_resp-n_resp_a))
	sample_early = DATAMOMS["sample_early"]
	sample_late  = DATAMOMS["sample_late"]

	initial_conditions = COMPMAT["initial_conditions_2d"]

## steps
	### step 1: compute DT choice probs WITHOUT charges
		### Main route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat"],
				BIGMAT["EARLY"],
				BIGMAT["LATE"],
				BIGMAT["sigma_adjust"],
				CHARGEMAT["zero_dt_charges"],
				dt_choice_free, eu_mat_free, prob_sum_mat, debug=(debug > 0))

		### Detour route
		solve_dt(theta_dt,
				BIGMAT["tmean_mat_detour"],
				BIGMAT["EARLY_detour"],
				BIGMAT["LATE_detour"],
				BIGMAT["sigma_adjust_detour"],
				CHARGEMAT["zero_dt_charges_detour"],
				dt_choice_free_detour, eu_mat_free_detour, prob_sum_mat_detour,
				debug=(debug > 0))

		### Main route WITH charges - used for non-detour commuters
		solve_dt(theta_dt,
				BIGMAT["tmean_mat"],
				BIGMAT["EARLY"],
				BIGMAT["LATE"],
				BIGMAT["sigma_adjust"],
				CHARGEMAT["dt_charges"],
				dt_choice_wcha, eu_mat_wcha, prob_sum_mat, debug=(debug > 0))


	### step 2: invert hastar distribution
		# optimal departure time mapping, without charges
		# using ACTUAL route choice probabilities (where applicable)
		probs_route1_free = DATAMOMS["data_a"].mean_long_0
		probs_route1_free_3d = repeat(probs_route1_free, outer=(1, 1, n_ha))

		dt_choice_mapping = copy(dt_choice_free)
		dt_choice_mapping[1:n_resp_a, :, :] =
			dt_choice_mapping[1:n_resp_a, :, :] .* (1.0 .- probs_route1_free_3d) .+
			dt_choice_free_detour               .*         probs_route1_free_3d

	if hastar_inversion_type == "single_per_person"

		# invert to a single hastar value (per person)
		invert_single_has_per_person(dt_choice=dt_choice_mapping, dt_choice_pre=dt_choice_pre, hastars=hastars)
	elseif hastar_inversion_type == "single"

		# invert to a single hastar value (one value for everyone)
		invert_single_has(dt_choice=dt_choice_mapping, dt_choice_pre=dt_choice_pre, hastars=hastars)
	elseif hastar_inversion_type == "given"
		hastars = COMPMAT["hastars_given"]
		# println("USING GIVEN HASTAR MATRIX")
	else
		# invert to an entire distribution of hastar
		@assert hastar_inversion_type == "distribution"
		invert_has(dt_choice=dt_choice_mapping, dt_choice_pre=dt_choice_pre, hastars=hastars)
	end

	### step 3: define utility matrices (only detour length)

	### step 3: run dynamic VOTT without any charges (but with FE)

		my_norm_epsilon = 1e-12 * (1.0 / (1.0 - theta_dict["delta"]))

		probs_route0_free, probs_hd_before, probs_dt_free,
				ss, _, _, _, _ =
			dynamic_vott_experiment(
				ss=nothing,
				theta_dict=theta_dict,
				BIGMAT=BIGMAT,
				CHARGEMAT=CHARGEMAT,
				DATAMOMS=DATAMOMS,
				hastar=hastars[1:n_resp_a, :, :],
				initial_conditions=initial_conditions,
				response_factor=0.0, # <---------------- charges switched OFF
				norm_epsilon=my_norm_epsilon,
				mydebug=(debug - 1)
			)

		# DEBUG
		# @assert all(ss.stat0 .>= 0.0)
		# @assert all(probs_route0_free .>= 0.0)

	### step 4: DT choices without charges during the experiment
		# detour sample: take from above
		probs_dt_free_full = zeros(n_resp,n_hd)
		probs_dt_free_full[sample_dynamic_vot_full, :] =
			probs_dt_free[sample_dynamic_vot, :]

		# non-detour sample: take from solve_dt() and integrate over hastar
		probs_dt_free_full[ .~sample_dynamic_vot_full, :] = sum(
			dt_choice_free[ .~sample_dynamic_vot_full,:,:] .*
			hastars[.~sample_dynamic_vot_full,:,:], dims=3)[:,:,1]

	### Step 4a: DT choices BEFORE the experiment
		probs_dt_before_full = zeros(n_resp,n_hd)
		probs_dt_before_full[sample_dynamic_vot_full, :] =
			probs_hd_before[sample_dynamic_vot, :]

		# non-detour sample: take from solve_dt() and integrate over hastar
		probs_dt_before_full[ .~sample_dynamic_vot_full, :] = sum(
			dt_choice_free[ .~sample_dynamic_vot_full,:,:] .*
			hastars[.~sample_dynamic_vot_full,:,:], dims=3)[:,:,1]

	### step 5: run dynamic VOTT WITH charges (and with FE)
		# # DEBUG
		# @assert sum(isnan.(initial_conditions)) == 0

		# no need to recompute steady state (ss)
		probs_route0_wcha, porbs_notused, probs_dt_wcha,
			_, w1_soln, w2_soln, w3_soln, w4_soln =
			dynamic_vott_experiment(
				ss=ss,
				theta_dict=theta_dict,
				BIGMAT=BIGMAT,
				CHARGEMAT=CHARGEMAT,
				DATAMOMS=DATAMOMS,
				hastar=hastars[1:n_resp_a, :, :],
				initial_conditions=initial_conditions,
				response_factor=1.0, # <---------------- charges switched ON
				norm_epsilon=my_norm_epsilon,
				mydebug=(debug - 1)
			)

		# # DEBUG
		# @assert sum(isnan.(probs_route0_wcha)) == 0
		# @assert sum(isnan.(ss.stat0)) == 0
		# @assert sum(isnan.(probs_dt_wcha)) == 0

	### step 6: DT choices WITH charges
		probs_dt_wcha_full = zeros(n_resp,n_hd)
		probs_dt_wcha_full[sample_dynamic_vot_full, :] =
			 probs_dt_wcha[sample_dynamic_vot, :]

		probs_dt_wcha_full[ .~sample_dynamic_vot_full, :] = sum(
			dt_choice_wcha[ .~sample_dynamic_vot_full,:,:] .*
			hastars[.~sample_dynamic_vot_full,:,:], dims=3)[:,:,1]

	### return
	return Dict{String, Any}(
		"dt_before" => probs_dt_before_full, # average over the relevant weeks (1-3 or 2-4)
		"dt_free" => probs_dt_free_full, # average over the relevant weeks (1-3 or 2-4)
		"dt_wcha" => probs_dt_wcha_full, # average over the relevant weeks (1-3 or 2-4)
		"probs_route0_free" => probs_route0_free,
		"probs_route0_wcha" => probs_route0_wcha,
		"w1_soln" => w1_soln,
		"w2_soln" => w2_soln,
		"w3_soln" => w3_soln,
		"w4_soln" => w4_soln
	)
end

section_dynamicvot_nodt = 1

# Dynamic VOTT without DT model: switching cost + nested logit


function nestedlogit_nodt(;
		theta_dict::Dict{String, Float64},
		e0::Vector{Float64}, e1::Vector{Float64},
		delta_u::Vector{Float64},
		sigma_adjust::Vector{Float64},
		return_probs=false,
		MAXDIFF=300.0
	)

## unpack
	mu = theta_dict["mu"]

## route level
	@. delta_u = (e1 - e0) / sigma_adjust / mu

	# bound from below relative to max
	@. delta_u = min(max(delta_u, -MAXDIFF), MAXDIFF)

	# exponentiate
	@. delta_u = exp(delta_u)

	# probabilities!
	if return_probs

		# probability to choose route 0
		prob0 = @. 1.0 / (1.0 + delta_u)

		# pi(h), pi(route=0)
		return prob0
	else
		# return just the expected utility

		# expected utility
		eu = @. e0 + log(1.0 + delta_u) * mu * sigma_adjust

		return eu
	end

end



function onestep_probs_nodt(;
		theta_dict::Dict{String, Float64},
		v0s_assumed::Vector{Float64},
		v1s_assumed::Vector{Float64},
		e0::Vector{Float64}, e1::Vector{Float64},
		delta_u::Vector{Float64},
		sigma_adjust::Vector{Float64},
		exp_utility_0::Array{Float64},
		exp_utility_1::Array{Float64},
		MAXDIFF=300.0
		)

## unpack
	gamma01 = theta_dict["gamma"] * theta_dict["gamma_01_factor"] # switch cost r_{t-1}=0
	gamma10 = theta_dict["gamma"] # switch cost r_{t-1}=1
	delta 	 = theta_dict["delta"]

## Starting with r_{t-1} = 0

	# net present value of each h x r combination
	@. e0 = exp_utility_0 + delta * v0s_assumed
	@. e1 = exp_utility_1 + delta * v1s_assumed - gamma01

	prob_route0_from0 =
	nestedlogit_nodt(
		theta_dict=theta_dict,
		e0=e0, e1=e1,
		delta_u=delta_u,
		sigma_adjust=sigma_adjust,
		return_probs=true
	)


## Starting with r_{t-1} = 1

	# net present value of each h x r combination
	@. e0 = exp_utility_0 + delta * v0s_assumed - gamma10
	@. e1 = exp_utility_1 + delta * v1s_assumed

	prob_route0_from1 =
	nestedlogit_nodt(
		theta_dict=theta_dict,
		e0=e0, e1=e1,
		delta_u=delta_u,
		sigma_adjust=sigma_adjust,
		return_probs=true
	)

	return (
		p01=1 .- prob_route0_from0,
		p10=prob_route0_from1
		)
end





"""
	Given assumed value functions for next period (two numbers):
		- compute utility from each route
		- (which is the same as) expected utility = v0 and v1 at the start of this period

	exp_utility_0 includes time costs

	exp_utility_0 = n_resp_a vector
"""
function onestep_eu_nodt(;
		theta_dict::Dict{String, Float64},
		v0s_assumed::Vector{Float64},
		v1s_assumed::Vector{Float64},
		e0::Vector{Float64}, e1::Vector{Float64},
		delta_u::Vector{Float64},
		sigma_adjust::Vector{Float64},
		exp_utility_0::Vector{Float64},
		exp_utility_1::Vector{Float64},
		MAXDIFF=300.0
		)

## unpack
	gamma01 = theta_dict["gamma"] * theta_dict["gamma_01_factor"] # switch cost r_{t-1}=0
	gamma10 = theta_dict["gamma"] # switch cost r_{t-1}=1
	delta 	 = theta_dict["delta"]

## Starting with r_{t-1} = 0

	# net present value of each h x r combination
	@. e0 = exp_utility_0 + delta * v0s_assumed
	@. e1 = exp_utility_1 + delta * v1s_assumed - gamma01

	v0_integrated =
	nestedlogit_nodt(
		theta_dict=theta_dict,
		e0=e0, e1=e1,
		delta_u=delta_u,
		sigma_adjust=sigma_adjust,
		return_probs=false
	)


## Starting with r_{t-1} = 1

	# net present value of each h x r combination
	@. e0 = exp_utility_0 + delta * v0s_assumed - gamma10
	@. e1 = exp_utility_1 + delta * v1s_assumed

	v1_integrated =
	nestedlogit_nodt(
		theta_dict=theta_dict,
		e0=e0, e1=e1,
		delta_u=delta_u,
		sigma_adjust=sigma_adjust,
		return_probs=false
	)

	return (v0_integrated, v1_integrated)
end

"""
	Compute one-step nested logit for all respondents
	Return current values
"""
function onestep_all_nodt(;
		theta_dict::Dict{String, Float64},
		v0s_assumed::Vector{Float64},
		v1s_assumed::Vector{Float64},
		e0::Vector{Float64}, e1::Vector{Float64},
		delta_u::Vector{Float64},
		sigma_adjust::Vector{Float64},
		exp_utility_0::Vector{Float64},
		exp_utility_1::Vector{Float64}
		)

	# n_resp_a = length(v0s_assumed)

	# values
	v0s, v1s = onestep_eu_nodt(
			theta_dict=theta_dict,
			v0s_assumed=v0s_assumed, v1s_assumed=v1s_assumed,
			e0=e0, e1=e1, delta_u=delta_u,
			sigma_adjust=sigma_adjust,
			exp_utility_0=exp_utility_0,
			exp_utility_1=exp_utility_1
			)

	# probabilities
	# phd0 = Pr(h| r_{it-1} = 0)
	# phd1 = Pr(h| r_{it-1} = 1)
	# p01 = Pr(r=1|r_{it-1} = 0)
	# p10 = Pr(r=0|r_{it-1} = 1)
	p01, p10 = onestep_probs_nodt(
			theta_dict=theta_dict,
			v0s_assumed=v0s_assumed, v1s_assumed=v1s_assumed,
			e0=e0, e1=e1, delta_u=delta_u,
			sigma_adjust=sigma_adjust,
			exp_utility_0=exp_utility_0,
			exp_utility_1=exp_utility_1
			)

	return (v0s=v0s, v1s=v1s, p01=p01, p10=p10, stat0=nothing)
end


function dynamic_vott_steady_state_nodt(;
		theta_dict::Dict{String, Float64},
		exp_utility_0::Vector{Float64},
		exp_utility_1::Vector{Float64},
		initial_conditions::Matrix{Float64},
		sigma_adjust::Vector{Float64},
		norm_epsilon=1e-12, lazy_param=0.1,
		max_steps=-1, mydebug=0
		)

## computation matrices
	n_resp_a = length(exp_utility_1)

	e0	= zeros(n_resp_a)
	e1	= zeros(n_resp_a)
	delta_u = zeros(n_resp_a)

	# initial conditions for value function
	v0_current = initial_conditions[:, 1]
	v1_current = initial_conditions[:, 2]

	norm_change = 1.0

	# iterate value function until we converge
	i = 0
	while (norm_change > norm_epsilon) && ((max_steps < 0) || (i < max_steps))
		# (mydebug > 0) && println("iteration ", i, " norm ", norm_change)

		v0_new, v1_new = onestep_eu_nodt(
				theta_dict=theta_dict,
				v0s_assumed=v0_current, v1s_assumed=v1_current,
				e0=e0, e1=e1,
				delta_u=delta_u,
				sigma_adjust=sigma_adjust,
				exp_utility_0=exp_utility_0,
				exp_utility_1=exp_utility_1
				)

		norm_change = maximum(
				abs.(v0_new .- v0_current) .+
				abs.(v1_new .- v1_current)
				 )

		v0_current = lazy_param * v0_new + (1.0 - lazy_param) * v0_current
		v1_current = lazy_param * v1_new + (1.0 - lazy_param) * v1_current

		i += 1
	end

	(mydebug > 0) && println("DONE after iterations: ", i)

	# error if went over the max number of steps
	if (max_steps > 0) && (i >= max_steps)
		# display(plot(1:max_steps, allv0s))
		# println("error for vi=", myvi)
		println("v0_current=", v0_current)
		println("v1_current=", v1_current)
		error("above max number of steps")
	end

	## compute probabilities
	# phd0 = Pr(h| r_{it-1} = 0)
	# phd1 = Pr(h| r_{it-1} = 1)
	# p01 = Pr(r=1|r_{it-1} = 0)
	# p10 = Pr(r=0|r_{it-1} = 1)
	p01, p10 =
	onestep_probs_nodt(
			theta_dict=theta_dict,
			v0s_assumed=v0_current, v1s_assumed=v1_current,
			e0=e0, e1=e1,
			delta_u=delta_u,
			sigma_adjust=sigma_adjust,
			exp_utility_0=exp_utility_0,
			exp_utility_1=exp_utility_1
			)

	# stationary distribution
	stat_probs_0 = p10 ./ (p01 .+ p10)

	# v0s 			= n_resp_a vector
	# probs_01 		= n_resp_a vector
	# stat_probs_0 	= n_resp_a vector

	# phd0 = Pr(h| r_{it-1} = 0)
	# phd1 = Pr(h| r_{it-1} = 1)
	# p01 = Pr(r=1|r_{it-1} = 0)
	# p10 = Pr(r=0|r_{it-1} = 1)
	# stat0 = stationary distribution of route0
	return (
		v0s=v0_current, v1s=v1_current,
		p01=p01, p10=p10, stat0=stat_probs_0
		)
end


"""
Simulate model: steady state + all weeks during experiment
"""
function dynamic_vott_experiment_nodt(;
	ss=nothing,
	theta_dict::Dict{String, Float64},
	BIGMAT::Dict{String, Array{Float64, 3}},
	CHARGEMAT::Dict{String, Array{Float64, 3}},
	DATAMOMS::Dict{String,Any},
	initial_conditions::Array{Float64},
	response_factor=1.0, ### SET THIS TO ZERO FOR NO CHARGES EVER
	my_lazy_param=1.0, norm_epsilon=1e-12,
	max_steps=50000, # maximum number of steps to find solution
	mydebug=0
	)

## unpack

	myalpha = theta_dict["alpha"]

	# charges positive, so should enter utility of route 0 negatively
	area_charges = CHARGEMAT["area_charge_mean_detour"][:,1,1]

	sigma_adjust_detour = BIGMAT["sigma_adjust_detour"][:,1,1]

	sample_early = DATAMOMS["sample_early"]
	sample_late  = DATAMOMS["sample_late"]

	# detours in HOURS
	adetourtime = convert.(Float64, DATAMOMS["data_a"].adetourtime) ./ 60.0

	n_resp_a = size(adetourtime,1)

## Make all parameters proportional to wages
	# If parameters are proportional to wages
	# equivalent to charges INVERSELY proportional to wages
	if haskey(theta_dict, "prop_to_wage") && theta_dict["prop_to_wage"] .> 0.0
		@assert theta_dict["prop_to_wage"] == 1.0
		response_factor_all = response_factor .* DATAMOMS["data_all"].inverse_income_with_pred
		response_factor_a   = response_factor_all[1:n_resp_a]
	else
		response_factor_all = response_factor
		response_factor_a   = response_factor
	end

## Step 0. Expected utilities
	exp_utility_0_free = @. adetourtime * 0.0

	exp_utility_1_free = @. (-1) * myalpha * adetourtime +
				   theta_dict["fe_late"] * sample_late +
				   theta_dict["fe_early"] * sample_early


## Step 1: solve for the steady state
	# return format: named tuple:
	# (v0s, v1s, p01, p10, stat0, phd0, phd1)

	# phd0 = Pr(h| r_{it-1} = 0)
	# phd1 = Pr(h| r_{it-1} = 1)
	# p01 = Pr(r=1|r_{it-1} = 0)
	# p10 = Pr(r=0|r_{it-1} = 1)
	# stat0 = stationary distribution of route0

	if isnothing(ss)
		ss = dynamic_vott_steady_state_nodt(
			theta_dict=theta_dict,
			exp_utility_0=exp_utility_0_free,
			exp_utility_1=exp_utility_1_free,
			initial_conditions=initial_conditions,
			sigma_adjust=sigma_adjust_detour,
			norm_epsilon=norm_epsilon, lazy_param=my_lazy_param,
			max_steps=max_steps, mydebug=max(0, mydebug-1)
			)

		# update initial conditions (value fn of route 0 and of route 1)
		initial_conditions[:, 1] = ss.v0s
		initial_conditions[:, 2] = ss.v1s
	end

## Step 2: Solve for value fn in weeks 1-4, working backward, starting from week 4
	# separately for early/late
	e0	= zeros(n_resp_a)
	e1	= zeros(n_resp_a)
	delta_u = zeros(n_resp_a)

	u0 = zeros(n_resp_a)
	u1 = zeros(n_resp_a)

## Week 4
	u0 .=  exp_utility_0_free .- theta_dict["fe_w4"] .-
		   area_charges .* response_factor_a .* sample_late
	u1 .= exp_utility_1_free

	w4_soln = onestep_all_nodt(
			theta_dict=theta_dict,
			v0s_assumed=ss.v0s,   # <---- start from steady state in week 5
			v1s_assumed=ss.v1s,
			e0=e0, e1=e1, delta_u=delta_u,
			sigma_adjust=sigma_adjust_detour,
			exp_utility_0=u0,
			exp_utility_1=u1
			)

## Week 3
	u0 .= exp_utility_0_free .- theta_dict["fe_w3"]
	u1 .= exp_utility_1_free

	w3_soln = onestep_all_nodt(
			theta_dict=theta_dict,
			v0s_assumed=w4_soln.v0s,   # <---- start from week 4
			v1s_assumed=w4_soln.v1s,
			e0=e0, e1=e1, delta_u=delta_u,
			sigma_adjust=sigma_adjust_detour,
			exp_utility_0=u0,
			exp_utility_1=u1)

## Week 2
	u0 .= exp_utility_0_free .- theta_dict["fe_w2"]
	u1 .= exp_utility_1_free

	w2_soln = onestep_all_nodt(
			theta_dict=theta_dict,
			v0s_assumed=w3_soln.v0s,   # <---- start from week 3
			v1s_assumed=w3_soln.v1s,
			e0=e0, e1=e1, delta_u=delta_u,
			sigma_adjust=sigma_adjust_detour,
			exp_utility_0=u0,
			exp_utility_1=u1)

## Week 1
	u0 .= exp_utility_0_free .- theta_dict["fe_w1"] .-
		  area_charges .* response_factor_a .* sample_early
	u1 .= exp_utility_1_free

	w1_soln = onestep_all_nodt(
			theta_dict=theta_dict,
			v0s_assumed=w2_soln.v0s,   # <---- start from week 2
			v1s_assumed=w2_soln.v1s,
			e0=e0, e1=e1, delta_u=delta_u,
			sigma_adjust=sigma_adjust_detour,
			exp_utility_0=u0,
			exp_utility_1=u1)

## Step 3: Compute fraction in each state each week using transition probabilities

	## Late: stationary probability to use route 0 (main route)
		late_route0_w0 = ss.stat0[sample_late]

		# push forward using transition probabilities
		late_route0_w1 = route0_transition(p0=late_route0_w0, prob_01=w1_soln.p01[sample_late], prob_10=w1_soln.p10[sample_late])
		late_route0_w2 = route0_transition(p0=late_route0_w1, prob_01=w2_soln.p01[sample_late], prob_10=w2_soln.p10[sample_late])
		late_route0_w3 = route0_transition(p0=late_route0_w2, prob_01=w3_soln.p01[sample_late], prob_10=w3_soln.p10[sample_late])
		late_route0_w4 = route0_transition(p0=late_route0_w3, prob_01=w4_soln.p01[sample_late], prob_10=w4_soln.p10[sample_late])

	## Early: stationary probability to use route 0 (main route)
		early_route0_w0 = ss.stat0[sample_early]

		# push forward using transition probabilities
		early_route0_w1 = route0_transition(p0=early_route0_w0, prob_01=w1_soln.p01[sample_early], prob_10=w1_soln.p10[sample_early])
		early_route0_w2 = route0_transition(p0=early_route0_w1, prob_01=w2_soln.p01[sample_early], prob_10=w2_soln.p10[sample_early])
		early_route0_w3 = route0_transition(p0=early_route0_w2, prob_01=w3_soln.p01[sample_early], prob_10=w3_soln.p10[sample_early])
		early_route0_w4 = route0_transition(p0=early_route0_w3, prob_01=w4_soln.p01[sample_early], prob_10=w4_soln.p10[sample_early])

	### Combine back route probabilities -- by week
		probs_route0 = zeros(n_resp_a, 5)
		probs_route0[sample_late, 1] = late_route0_w0
		probs_route0[sample_late, 2] = late_route0_w1
		probs_route0[sample_late, 3] = late_route0_w2
		probs_route0[sample_late, 4] = late_route0_w3
		probs_route0[sample_late, 5] = late_route0_w4

		probs_route0[sample_early, 1] = early_route0_w0
		probs_route0[sample_early, 2] = early_route0_w1
		probs_route0[sample_early, 3] = early_route0_w2
		probs_route0[sample_early, 4] = early_route0_w3
		probs_route0[sample_early, 5] = early_route0_w4

	return probs_route0, ss, w1_soln, w2_soln, w3_soln, w4_soln
end


"""
Main model: dynamic route choice model withOUT departure time choice
- uses nested logit
- performs inversion using pre-experiment route choice frequencies in data
-
"""
function solve_dynamicvot_nodt(;
		theta::Vector{Float64},
		param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}},
		BIGMAT::Dict{String, Array{Float64, 3}},
		CHARGEMAT::Dict{String, Array{Float64, 3}},
		COMPMAT::Dict{String, Array{Float64, N} where N},
		DATAMOMS::Dict{String, Any},
		max_steps::Int64=50000,
		debug::Int64=0
	)

## unpack parameters
	theta_dict = get_params_fn(theta, param_factors)

## upack matrices
	n_resp, n_resp_a, n_resp_a_ct23, n_resp_a_tr = DATAMOMS["n_resp_vec"]

	# sigma_adjust_detour = BIGMAT["sigma_adjust_detour"]

	sample_dynamic_vot = DATAMOMS["sample_dynamic_vot"]
	sample_dynamic_vot_full = vcat(sample_dynamic_vot, zeros(Bool, n_resp-n_resp_a))
	sample_early = DATAMOMS["sample_early"]
	sample_late  = DATAMOMS["sample_late"]

	initial_conditions = COMPMAT["initial_conditions_2d"]

## steps
	### step 1: dynamic VOTT without any charges (but with FE)

		my_norm_epsilon = 1e-12 * (1.0 / (1.0 - theta_dict["delta"]))

		probs_route0_free, ss, _, _, _, _ =
			dynamic_vott_experiment_nodt(
				ss=nothing,
				theta_dict=theta_dict,
				BIGMAT=BIGMAT,
				CHARGEMAT=CHARGEMAT,
				DATAMOMS=DATAMOMS,
				initial_conditions=initial_conditions,
				response_factor=0.0, # <---------------- charges switched OFF
				norm_epsilon=my_norm_epsilon,
                max_steps=max_steps,
				mydebug=(debug - 1)
			)

	### step 5: run dynamic VOTT WITH charges (and with FE)
		# no need to recompute steady state (ss)
		probs_route0_wcha, _, w1_soln, w2_soln, w3_soln, w4_soln =
			dynamic_vott_experiment_nodt(
				ss=ss,
				theta_dict=theta_dict,
				BIGMAT=BIGMAT,
				CHARGEMAT=CHARGEMAT,
				DATAMOMS=DATAMOMS,
				initial_conditions=initial_conditions,
				response_factor=1.0, # <---------------- charges switched ON
				norm_epsilon=my_norm_epsilon,
				max_steps=max_steps,
				mydebug=(debug - 1)
			)

	### return
	return Dict{String, Any}(
		"probs_route0_free" => probs_route0_free,
		"probs_route0_wcha" => probs_route0_wcha,
		"w1_soln" => w1_soln,
		"w2_soln" => w2_soln,
		"w3_soln" => w3_soln,
		"w4_soln" => w4_soln
	)
end
