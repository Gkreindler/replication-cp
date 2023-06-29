
function get_density_old(
	choice_prbs::Matrix{Float64},
	km_left::Matrix{Float64},
	i::Int64, j::Int64)

	mysum = 0.0
	myweight = 0.0

	for k=j:(i-1)
		for orig_idx=1:304
			myweight = km_left[orig_idx, k]
			if myweight .> 0.0
				# ease in "density" up to 2.5 km (=5 min at 2min/km)
				myweight = min(myweight / 2.5, 1.0)
				mysum += myweight *
					sum(
					choice_prbs[((orig_idx-1) * 10 + 1):(orig_idx * 10), k]
					)
			end
		end
	end
	return mysum
end

"""
	linear (3 parameters)
	1 = constant
	2 = slope
	3 = density conversion factor (to reach similar units to those used to
		estimate the RT relationship)
	4 = exponent (=1 for linear)
"""
function instant_delay(
		my_density::Float64,
		rt_params::Vector{Float64}
	)

	if length(rt_params) == 3
		# linear
		return rt_params[1] +
			   rt_params[2] *
			   (rt_params[3] * my_density)
	else
		# power
		@assert length(rt_params) == 4
		return rt_params[1] +
			   rt_params[2] *
			   (rt_params[3] * my_density) ^ rt_params[4]
	end
end

function make_travel_times_density(;
		choice_prbs::Matrix{Float64},  # N_agents x n_hd
		rt_params::Vector{Float64},
		extra_idx1::Int64=-1,
		extra_idx2::Int64=-1,
		km_mean::Vector{Float64},
		km_left::Matrix{Float64},
		ttimes_small::Matrix{Float64},
		ttimes::Matrix{Float64},
		delta_t::Float64,   ## time between two hdgrid pnts (MINUTES)
		extra_traffic::Matrix{Float64}=Array{Float64}(undef, 0, 0), # by default zero
		debug::Bool=false
	)
	#=
	Plan:
	1. loop over departure times
		2. get new departures
		3. get older trips that are still ongoing  <- assuming these affect equally, so approx within
		4. compute total density
		5. get instantaneous delay
		6. update travel times for all previous trips
		7. update KM left in all trips (for all j < i)
	=#

	fill!(km_left, 0.0)
	fill!(ttimes_small, 0.0)
	n_hd = size(choice_prbs, 2)
	n_agents = size(choice_prbs, 1)

	densities = zeros(n_hd)
	instant_delays = zeros(n_hd)

	# this is the smallest index where trips are still ongoing
	j=1

	for i=1:n_hd
		debug && println("i=", i)
		density_new = sum(choice_prbs[:, i])

		density_old = 0.0
		if i > 1
			density_old = get_density_old(choice_prbs, km_left, i, j)
		end

		# apply (density => instant delay) relationship to total density
		density_total = density_new + density_old + (extra_idx1 <= i <= extra_idx2) * 0.001
		instant_delay_here = instant_delay(density_total, rt_params)

		instant_dist = (1/instant_delay_here) * delta_t
		instant_dist_inv = instant_delay_here / delta_t

		# initialize current departures
		km_left[:, i] .= km_mean

		# step 6. update travel times
		@. ttimes_small[:, j:i] += instant_delay_here *
			min(instant_dist, max(0.0, km_left[:, j:i]))

		# step 7. update all leftover KM
		@. km_left[:, j:i] -= instant_dist

		densities[i] = density_new + density_old
		instant_delays[i] = instant_delay_here

		# update j
		while (sum( km_left[:, j] .> 0.0 ) == 0) && (j < n_hd)
			j += 1
			debug && println("...j=", j)
		end
	end

	## Wrap up travel times at the end -> some trips not over!
	@. ttimes_small[:, j:end] += max(0.0, km_left[:, j:end]) * instant_delays[end]

	## copy over to the large matrix
	indices = kron(1:304, ones(Int64, 10))
	ttimes .= ttimes_small[indices, :]

	return densities, instant_delays
end


function make_travel_times_density_rtcheck(;
		t_start::Vector{Float64},  # N_agents
		hdgrid::Vector{Float64},  # n_hd
		rt_params::Vector{Float64},
		km_mean::Vector{Float64},
		debug::Bool=false
	)
	#=
	Plan:
	1. loop over departure times
		2. get new departures
		3. get older trips that are still ongoing  <- assuming these affect equally, so approx within
		4. compute total density
		5. get instantaneous delay
		6. update travel times for all previous trips
		7. update KM left in all trips (for all j < i)
	=#

	km_left = zero(km_mean)
	ttimes = zero(t_start)
	n_hd = length(hdgrid)
	n_agents = length(t_start)

	delta_t_hrs = hdgrid[2] - hdgrid[1]
	delta_t_minutes = delta_t_hrs * 60

	densities = zeros(n_hd)
	instant_delays = zeros(n_hd)

	# this is the smallest index where trips are still ongoing
	j=1

	for i=1:n_hd
		debug && println("i=", i)

		new_trips = hdgrid[i] .<= t_start .< (hdgrid[i] + delta_t_hrs)

		density_new = sum(new_trips)

		# ongoing trips
		density_old = min.(max.(0.0, km_left) ./ 2.5, 1.0) |> sum

		# "activate" new trips
		km_left[new_trips] = km_mean[new_trips]

		# apply (density => instant delay) relationship to total density
		density_total = density_new + density_old # + (extra_idx1 <= i <= extra_idx2)
		instant_delay_here = instant_delay(density_total, rt_params)

		instant_dist = (1/instant_delay_here) * delta_t_minutes
		instant_dist_inv = instant_delay_here / delta_t_minutes

		# step 6. update travel times
		@. ttimes += instant_delay_here * min(instant_dist, max(0.0, km_left))

		# step 7. update all leftover KM
		@. km_left -= instant_dist

		debug && println("instant delay ", instant_delay_here)
		debug && println("instant distance ", instant_dist)

		densities[i] = density_new + density_old
		instant_delays[i] = instant_delay_here
	end

	## Wrap up travel times at the end -> some trips not over!
	@. ttimes += max(0.0, km_left) * instant_delays[end]

	return densities, instant_delays, ttimes
end
