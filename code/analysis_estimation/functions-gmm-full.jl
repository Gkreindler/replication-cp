

"""
	Put vector of parameters in a dictionary
		- combine fixed parametrs with parametrs to estimate
		- fixed parameters are included in the 3rd entry of param_factors

"""
function get_params_fn(myparams_vec, param_factors)

	myparams_dict = Dict{String, Float64}()
	idx = 1

	for i=1:length(param_factors)
		param_name = param_factors[i][1]

		if isnothing(param_factors[i][3])
			# parameter not fixed -- take from vector
			myparams_dict[param_name] = myparams_vec[idx] * param_factors[i][2]
			idx += 1
		else
			# parameter fixed -- take from param_factors
			myparams_dict[param_name] = param_factors[i][3] * param_factors[i][2]
		end
	end

	# check that we used all the parameters in the vector
	@assert idx == length(myparams_vec) + 1

	return myparams_dict
end


"""
	For debugging: print solution (for a single element)
"""
function get_1el_soln(mysoln)

	# (vgrid, e1s, e2s, v1s, v2s, probs_12, probs_21, stationary_prob1)
	vgrid, e1s, e2s, v1s, v2s, probs_12, probs_21, stationary_prob1 = mysoln

	if ~isnothing(stationary_prob1)
		display((vgrid[1], e1s[1], e2s[1], v1s[1], v2s[1], probs_12[1], probs_21[1], stationary_prob1[1]))
	else
		display((vgrid[1], e1s[1], e2s[1], v1s[1], v2s[1], probs_12[1], probs_21[1]))
	end
end


"""
	Histogram for area het moments
"""
function hist_subgrid(vals, hd_grid::Vector{Float64}, hd_bins::Vector{Float64})
	# we assume that the endpoints of hd_grid and WITHIN those of hd_bins
	n_bins = length(hd_bins) - 1

	# histogram model
    hist_mat = zeros(size(vals,1), n_bins)
    for i=0:(n_bins-1)
        # which hd are in this bin?
        h_min = hd_bins[i+1]
        h_max = hd_bins[i+2]

        if i < n_bins - 1
            hdgrid_inbin = (h_min .<= hd_grid) .& (hd_grid .< h_max)
        else
            hdgrid_inbin = (h_min .<= hd_grid) .& (hd_grid .<= h_max)
        end

        hist_mat[:,i+1] = sum(vals[:, hdgrid_inbin], dims=2)
    end

	return hist_mat
end

function histw_exp(vals::Matrix{Float64}, weights::Matrix{Float64},
					nbins::Int64)
	#  histogram function
	# %   compute nbins columns like vals,
	# %   column i = the weight for that value if val in [i/nbins, (i+1)/nbins)
	# %     (except <= at the end)
	# %   then take the *mean* of every consecutive dim1 elements

	# dim2 = div(length(vals), dim1)
	histw_out = zeros(size(vals)[1], nbins)

	for i=0:(nbins-1)
		hmin = i / nbins
		hmax = (i+1) / nbins
		if i < nbins - 1
			my_hist = ((hmin .<= vals) .& (vals .< hmax)) .* weights
		else
			my_hist = ((hmin .<= vals) .& (vals .<= hmax)) .* weights
		end

		histw_out[:,i+1] = sum(my_hist,dims=2)
	end

	return histw_out
end



"""
Moments for the full GMM:
	(1) DT shares (from old GMM)
	(2) DT heterogeneity (from old GMM)
	(3) 5+5 weekly area shares (from VOT GMM)
	(4) 2 moments: VOT response in week 1 heterogeneity by areadetourtime (from VOT GMM)
	dt_factor: 	multiply the DT diff-in-diffs by this amount. Default = 1.0. Used
				to see how estimated beta_E and beta_L change if this set of moments
				had higher or lower amplitude
"""
function moments_fullgmm(;
		theta::Vector{Float64},
		param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}},
	    SIMMODEL::Dict{String, Any},
	    DATAMOMS::Dict{String, Any},
		include_dt_main::Bool=true,
		include_dynamicarea_main::Bool=true,
		transition_moment_type::Int64=3,
		dt_factor::Float64=1.0
	)

	# always must include at least ONE moment group :-)
	@assert include_dt_main || include_dynamicarea_main
	@assert 0 <= transition_moment_type <= 3
	moms_dt_main   = nothing
	moms_dynamicarea_main = nothing

## Unpack

	cha_ct23 = DATAMOMS["cha_ct23"]
	cha_tr   = DATAMOMS["cha_tr"]

	sample_ct23   = DATAMOMS["sample_ct23"] .* DATAMOMS["sample_dynamic_vot"]
	sample_tr     = DATAMOMS["sample_tr"] .* DATAMOMS["sample_dynamic_vot"]

	nct23 = DATAMOMS["nct23"]
	ntr   = DATAMOMS["ntr"]

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])
	n_hd_bins = DATAMOMS["n_hd_bins"]
	bin_h, bin_l, n_hd_mom_l, n_hd_mom_h = DATAMOMS["hd_bins_hl"]

	n_hd_mom_l = Int64(n_hd_mom_l)
	n_hd_mom_h = Int64(n_hd_mom_h)
	n_hd_moms = (n_hd_mom_h-n_hd_mom_l+1)

	hdgrid = DATAMOMS["hdgrid"]

	sample_dt_pre = DATAMOMS["sample_dt_pre"]
	sample_dt_pos = DATAMOMS["sample_dt_pos"]
	sample_dt_treat = DATAMOMS["sample_dt_treat"]
	sample_dt_contr = DATAMOMS["sample_dt_contr"]

	mean_dt_pre   = DATAMOMS["mean_dt_pre"]
	mean_dt_pos   = DATAMOMS["mean_dt_pos"]

	shadow_rates = DATAMOMS["shadow_rates"]

	#
	prob_respond_alpha = get_params_fn(theta, param_factors)["prob_respond_alpha"]
	data_a = DATAMOMS["data_a"]


## DEPARTURE TIME MOMS

if include_dt_main
	    # for each bin:
	    bin_delta = (bin_h-bin_l)/(n_hd_bins - 1)
	    hd_bins = range(bin_l - bin_delta/2, stop=bin_h + bin_delta/2, length=n_hd_bins + 1)
		hd_bins = Array(hd_bins)

	# EMPIRICAL - POST SAMPLE
		density_empir_pos = mean_dt_pos #./ 100.0
	# EMPIRICAL - PRE SAMPLE
		density_empir_pre = mean_dt_pre #./ 100.0

	# replace zeros
	    density_empir_pos[sample_dt_pos .== 0, :] .= 0.0
	    @assert(all(.~isnan.(density_empir_pos)))
	    density_empir_pre[sample_dt_pre .== 0, :] .= 0.0
	    @assert(all(.~isnan.(density_empir_pre)))

	# check this sums up to 1.0
		@assert norm(sum(density_empir_pre[sample_dt_pre .== 1, :], dims=2) .- 1.0) .<= 1e-6
		@assert norm(sum(density_empir_pos[sample_dt_pos .== 1, :], dims=2) .- 1.0) .<= 1e-6

	# MODEL-predicted histogram - n_resp x n_bins
	    ### compute model histograms (sum within hdgrid)
	        density_model_wcha = hist_subgrid(SIMMODEL["dt_wcha"], hdgrid, hd_bins)
	        density_model_free = hist_subgrid(SIMMODEL["dt_free"], hdgrid, hd_bins)
			density_model_before = hist_subgrid(SIMMODEL["dt_before"], hdgrid, hd_bins)

		### no zeros
			c1 = all(sum(density_model_wcha, dims=2) .> 0.0)
			c2 = all(sum(density_model_free, dims=2) .> 0.0)
			c3 = all(sum(density_model_before, dims=2) .> 0.0)
			if ~c1 || ~c2 || ~c3
				println(theta)
				println(theta)
				println(theta)
				println(theta)
				println(c1)
				println(c2)
				println(c3)
				error("shares sum up to zero")
			end
			@assert c1
			@assert c2
			@assert c3

		### when using larger DT grid, need to normalize to this subsample!
			density_model_wcha = density_model_wcha ./ sum(density_model_wcha, dims=2)
			density_model_free = density_model_free ./ sum(density_model_free, dims=2)
			density_model_before = density_model_before ./ sum(density_model_before, dims=2)

		### ensure no NANs
			@assert norm(sum(density_model_wcha, dims=2) .- 1.0) <= 1e-6
			@assert norm(sum(density_model_free, dims=2) .- 1.0) <= 1e-6
			@assert norm(sum(density_model_before, dims=2) .- 1.0) <= 1e-6

		### Model during experiment attenuated
			density_model_wcha_att =
				   prob_respond_alpha  .* density_model_wcha .+
			(1.0 - prob_respond_alpha) .* density_model_free

	## DEPARTURE TIME MOMENTS
	    # four cells: treat/control pre/post
	    n_treat_pre = sum(sample_dt_treat .* sample_dt_pre)
	    n_treat_pos = sum(sample_dt_treat .* sample_dt_pos)
	    n_contr_pre = sum(sample_dt_contr .* sample_dt_pre)
	    n_contr_pos = sum(sample_dt_contr .* sample_dt_pos)

	    # treat pre
	    treat_pre_data  = density_empir_pre    .* sample_dt_treat .* sample_dt_pre * n_resp_all / n_treat_pre
		treat_pre_model = density_model_before .* sample_dt_treat .* sample_dt_pre * n_resp_all / n_treat_pre

	    # control pre
	    contr_pre_data  = density_empir_pre    .* sample_dt_contr .* sample_dt_pre * n_resp_all / n_contr_pre
	    contr_pre_model = density_model_before .* sample_dt_contr .* sample_dt_pre * n_resp_all / n_contr_pre

	    # treat post
	    treat_post_data  = density_empir_pos      .* sample_dt_treat .* sample_dt_pos * n_resp_all / n_treat_pos
	    treat_post_model = density_model_wcha_att .* sample_dt_treat .* sample_dt_pos * n_resp_all / n_treat_pos

	    # control post
	    contr_post_data  = density_empir_pos      .* sample_dt_contr .* sample_dt_pos * n_resp_all / n_contr_pos
		contr_post_model = density_model_wcha_att .* sample_dt_contr .* sample_dt_pos * n_resp_all / n_contr_pos

		# dt_factor = 1.0 by default
		# I use this to find the impact on estimated parameters from higher/lower amplitude DD DT results
		treat_data_post_minus_pre = (treat_post_data .- treat_pre_data) .* dt_factor
		contr_data_post_minus_pre = (contr_post_data .- contr_pre_data) .* dt_factor

		treat_model_post_minus_pre = treat_post_model .- treat_pre_model
		contr_model_post_minus_pre = contr_post_model .- contr_pre_model

		dd = @. (treat_model_post_minus_pre - treat_data_post_minus_pre) -
				(contr_model_post_minus_pre - contr_data_post_minus_pre)

		# dd = ((treat_post_data .- treat_pre_data) .-
		# 	  (prob_respond_alpha        .* treat_post_model .+
		# 	   (1 .- prob_respond_alpha) .* treat_post_model_noresponse .- treat_pre_model )) .-
		# 	 ((contr_post_data .- contr_pre_data) .-
		# 	   (prob_respond_alpha        .* contr_post_model .+
		# 	    (1 .- prob_respond_alpha) .* contr_post_model_noresponse .- contr_pre_model))

	    # moms[:,1:n_hd_moms] = dd[:, n_hd_mom_l:n_hd_mom_h]
		moms_dt_main = dd[:, n_hd_mom_l:n_hd_mom_h]
end


## Area moments (12 moments from VOT GMM)
if include_dynamicarea_main
	# n_resp_a x 5 matrix
	probs_route1_free = 1 .- SIMMODEL["probs_route0_free"]
	probs_route1_wcha = 1 .- SIMMODEL["probs_route0_wcha"]

	probs_route1_attenuated = probs_route1_wcha .* prob_respond_alpha .+
							  probs_route1_free .* (1 .- prob_respond_alpha)

	# sample dynamic VOT already included in .a_early and .a_late
	@assert all(data_a.a_early .== (data_a.a_early .* data_a.sample_dynamic_vot))

	nmoms = 10 + (transition_moment_type > 0)
	moms_vott	= zeros(n_resp_all, nmoms)
	moms_sample = zeros(n_resp_all, nmoms)

	# moments 1-5 are "late" means in weeks 0-4
	moms_vott[1:n_resp_a, 1] = @. (probs_route1_attenuated[:, 1] - data_a.mean_long_0) * data_a.a_late * (data_a.n_obs_0 > 0)  # late week 0
	moms_vott[1:n_resp_a, 2] = @. (probs_route1_attenuated[:, 2] - data_a.mean_long_1) * data_a.a_late * (data_a.n_obs_1 > 0)  # late week 1
	moms_vott[1:n_resp_a, 3] = @. (probs_route1_attenuated[:, 3] - data_a.mean_long_2) * data_a.a_late * (data_a.n_obs_2 > 0)  # late week 2
	moms_vott[1:n_resp_a, 4] = @. (probs_route1_attenuated[:, 4] - data_a.mean_long_3) * data_a.a_late * (data_a.n_obs_3 > 0)  # late week 3
	moms_vott[1:n_resp_a, 5] = @. (probs_route1_attenuated[:, 5] - data_a.mean_long_4) * data_a.a_late * (data_a.n_obs_4 > 0)  # late week 4

	# moments 6-10 are "early" means in weeks 0-4
	moms_vott[1:n_resp_a, 5 + 1] = @. (probs_route1_attenuated[:, 1] - data_a.mean_long_0) * data_a.a_early * (data_a.n_obs_0 > 0)  # early week 0
	moms_vott[1:n_resp_a, 5 + 2] = @. (probs_route1_attenuated[:, 2] - data_a.mean_long_1) * data_a.a_early * (data_a.n_obs_1 > 0)  # early week 1
	moms_vott[1:n_resp_a, 5 + 3] = @. (probs_route1_attenuated[:, 3] - data_a.mean_long_2) * data_a.a_early * (data_a.n_obs_2 > 0)  # early week 2
	moms_vott[1:n_resp_a, 5 + 4] = @. (probs_route1_attenuated[:, 4] - data_a.mean_long_3) * data_a.a_early * (data_a.n_obs_3 > 0)  # early week 3
	moms_vott[1:n_resp_a, 5 + 5] = @. (probs_route1_attenuated[:, 5] - data_a.mean_long_4) * data_a.a_early * (data_a.n_obs_4 > 0)  # early week 4


	# if dyanmicarea_main_12
	# 	# moment 11: early group, week 0 de-mean and correlate with detour time
	# 	mysample = data_a.a_early .* (data_a.n_obs_0 .> 0) .* data_a.sample_dynamic_vot
	# 	n_early_w0 = sum(mysample)
	# 	early_w0_mean = sum(mysample .* probs_route1_attenuated[:, 1]) / n_early_w0     ### Week 0
	# 	detour_w0_mean = sum(mysample .* data_a.adetourtime_normalized) / n_early_w0
	#
	# 	@assert sum((probs_route1_attenuated[:, 1] .- early_w0_mean) .* mysample) < 1e-8
	# 	@assert sum((data_a.adetourtime_normalized .- detour_w0_mean) .* mysample) < 1e-8
	#
	# 	moms_dynamicarea_main[1:n_resp_a, 11] =
	# 			(probs_route1_attenuated[:, 1] .- early_w0_mean) .*
	# 			(data_a.adetourtime_normalized .- detour_w0_mean) .* mysample .+
	# 			(-1) .* DATAMOMS["moment_het_w0"]
	#
	# 	# moment 12: early group, week 1 de-mean and correlate with detour time
	# 	mysample = data_a.a_early .* (data_a.n_obs_1 .> 0) .* data_a.sample_dynamic_vot
	# 	n_early_w1 = sum(mysample)
	# 	early_w1_mean = sum(mysample .* probs_route1_attenuated[:, 2]) / n_early_w1     ### Week 1
	# 	detour_w1_mean = sum(mysample .* data_a.adetourtime_normalized) / n_early_w1
	#
	# 	@assert sum((probs_route1_attenuated[:, 2] .- early_w1_mean) .* mysample) < 1e-8
	# 	@assert sum((data_a.adetourtime_normalized .- detour_w1_mean) .* mysample) < 1e-8
	#
	# 	moms_dynamicarea_main[1:n_resp_a, 12] =
	# 			(probs_route1_attenuated[:, 2] .- early_w1_mean) .*
	# 			(data_a.adetourtime_normalized .- detour_w1_mean) .* mysample .+
	# 			(-1) .* DATAMOMS["moment_het_w1"]
	# end

	# if dyanmicarea_main_p10

	# DATA:
	# switching sample: data in weeks 1 and 2 AND conditional on detour in week 1
	# sample_p10_week12 = (data_a.sample_dynamic_vot) .& (data_a.a_early .== 1) .&
	# 					(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0) .&
	# 					(data_a.mean_long_1 .>= 0.75)
	#
	# # dummy for switching to route 0 in week 2
	# switched_p10_week12 = (data_a.mean_long_1 .>= 0.5) .* (data_a.mean_long_2 .< 0.5)
	#
	# # MODEL: vector of transition probabilities from route 1 to route 0 at start of period 2
	# # same as transition from t=1 to t=2
	# moms[1:n_resp_a, 11] =
	# 	sample_p10_week12 .* (switched_p10_week12 - SIMMODEL["w2_soln"].p10)

	threshold_low = 0.5
	threshold_high = 1.0 - threshold_low

	if transition_moment_type == 1
		# DATA:
		# switching sample: data in weeks 1 and 2 AND conditional on detour in week 1
		sample_p10_week12 = (data_a.sample_dynamic_vot) .& (data_a.a_early .== 1) .&
							(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0) .&
							(data_a.mean_long_1 .>= 0.75)

		# dummy for switching to route 0 in week 2
		# switched_p10_week12 = (data_a.mean_long_1 .>= threshold_high) .* (data_a.mean_long_2 .< threshold_low)
		switched_p10_week12 = (data_a.mean_long_1 .>= 0.5) .* (data_a.mean_long_2 .< 0.5)

		# MODEL: vector of transition probabilities from route 1 to route 0 at start of period 2
		# same as transition from t=1 to t=2
		moms_sample[1:n_resp_a, 11] = sample_p10_week12
		moms_vott[1:n_resp_a, 11] = sample_p10_week12 .* (switched_p10_week12 .- SIMMODEL["w2_soln"].p10)

	elseif transition_moment_type == 2
		sample_p10_week12 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_1 .>= threshold_high) .&
							(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0)

		sample_p10_week23 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_2 .>= threshold_high) .&
							(data_a.n_obs_2 .> 0) .& (data_a.n_obs_3 .> 0)

		sample_p10_week34 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_3 .>= threshold_high) .&
							(data_a.n_obs_3 .> 0) .& (data_a.n_obs_4 .> 0)

		switched_p10_week12 = (data_a.mean_long_1 .>= threshold_high) .* (data_a.mean_long_2 .< threshold_low)
		switched_p10_week23 = (data_a.mean_long_2 .>= threshold_high) .* (data_a.mean_long_3 .< threshold_low)
		switched_p10_week34 = (data_a.mean_long_3 .>= threshold_high) .* (data_a.mean_long_4 .< threshold_low)

		moms_sample[1:n_resp_a, 11] = sample_p10_week12 .| sample_p10_week23 .| sample_p10_week34
		moms_vott[1:n_resp_a, 11] = (sample_p10_week12 .* (switched_p10_week12 .- SIMMODEL["w2_soln"].p10) .+
								sample_p10_week23 .* (switched_p10_week23 .- SIMMODEL["w3_soln"].p10) .+
								sample_p10_week34 .* (switched_p10_week34 .- SIMMODEL["w4_soln"].p10)) ./
								max.(1.0, sample_p10_week12 .+ sample_p10_week23 .+ sample_p10_week34)
	elseif transition_moment_type == 3
		sample_p01_week12 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_1 .< threshold_low) .&
							(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0)

		sample_p01_week23 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_2 .< threshold_low) .&
							(data_a.n_obs_2 .> 0) .& (data_a.n_obs_3 .> 0)

		sample_p01_week34 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_3 .< threshold_low) .&
							(data_a.n_obs_3 .> 0) .& (data_a.n_obs_4 .> 0)

		switched_p01_week12 = (data_a.mean_long_1 .< threshold_low) .* (data_a.mean_long_2 .>= threshold_high)
		switched_p01_week23 = (data_a.mean_long_2 .< threshold_low) .* (data_a.mean_long_3 .>= threshold_high)
		switched_p01_week34 = (data_a.mean_long_3 .< threshold_low) .* (data_a.mean_long_4 .>= threshold_high)

		moms_sample[1:n_resp_a, 11] = sample_p01_week12 .| sample_p01_week23 .| sample_p01_week34
		moms_vott[1:n_resp_a, 11] = (sample_p01_week12 .* (switched_p01_week12 .- SIMMODEL["w2_soln"].p01) .+
								sample_p01_week23 .* (switched_p01_week23 .- SIMMODEL["w3_soln"].p01) .+
								sample_p01_week34 .* (switched_p01_week34 .- SIMMODEL["w4_soln"].p01)) ./
								max.(1.0, sample_p01_week12 .+ sample_p01_week23 .+ sample_p01_week34)

	end

end

## combine moments
	if include_dt_main
		moms = moms_dt_main
		include_dynamicarea_main && (moms = hcat(moms, moms_vott))
		return moms
	else

		# we must be in the case with ONLY dynamic VOT moments
		@assert include_dynamicarea_main
		return moms_vott
	end

end


"""
Model with nested logit (static) VOT
DT and VOT moments:
49	(1) DT shares (from old GMM)
2	(2) AREA main
"""
function moments_nestedlogitmodel_all(;
		theta::Vector{Float64},
		param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}},
	    SIMMODEL::Dict{String, Any},
	    DATAMOMS::Dict{String, Any},
		include_dt_main::Bool=true,
		include_dt_het::Bool=true,
		include_area_main::Bool=true,
		include_area_het::Bool=true
	)

	# keep things simple: always include dt moms
	@assert include_dt_main
	moms_dt_main   = nothing
	moms_area_main = nothing

## Unpack

	cha_ct23 = DATAMOMS["cha_ct23"]
	cha_tr   = DATAMOMS["cha_tr"]

	sample_ct23   = DATAMOMS["sample_ct23"] .* DATAMOMS["sample_dynamic_vot"]
	sample_tr     = DATAMOMS["sample_tr"] .* DATAMOMS["sample_dynamic_vot"]

	nct23 = DATAMOMS["nct23"]
	ntr   = DATAMOMS["ntr"]

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])
	n_hd_bins = DATAMOMS["n_hd_bins"]
	bin_h, bin_l, n_hd_mom_l, n_hd_mom_h = DATAMOMS["hd_bins_hl"]

	n_hd_mom_l = Int64(n_hd_mom_l)
	n_hd_mom_h = Int64(n_hd_mom_h)
	n_hd_moms = (n_hd_mom_h-n_hd_mom_l+1)

	hdgrid = DATAMOMS["hdgrid"]

	sample_dt_pre = DATAMOMS["sample_dt_pre"]
	sample_dt_pos = DATAMOMS["sample_dt_pos"]
	sample_dt_treat = DATAMOMS["sample_dt_treat"]
	sample_dt_contr = DATAMOMS["sample_dt_contr"]

	mean_dt_pre   = DATAMOMS["mean_dt_pre"]
	mean_dt_pos   = DATAMOMS["mean_dt_pos"]

	# idfxnpre = DATAMOMS["idfxnpre"]
	# idfxnpos = DATAMOMS["idfxnpos"]
	# idfxbe   = DATAMOMS["idfxbe"]
	# idfx     = DATAMOMS["idfx"]

	shadow_rates = DATAMOMS["shadow_rates"]

	#
	prob_respond_alpha = get_params_fn(theta, param_factors)["prob_respond_alpha"]

	data_a = DATAMOMS["data_a"]

	# moms = zeros(n_resp_all, n_moms)


## DEPARTURE TIME MOMS

# if include_dt_main
	    # for each bin:
	    bin_delta = (bin_h-bin_l)/(n_hd_bins - 1)
	    hd_bins = range(bin_l - bin_delta/2, stop=bin_h + bin_delta/2, length=n_hd_bins + 1)
		hd_bins = Array(hd_bins)

	# EMPIRICAL - POST SAMPLE
		density_empir_pos = mean_dt_pos #./ 100.0
	# EMPIRICAL - PRE SAMPLE
		density_empir_pre = mean_dt_pre #./ 100.0

	    density_empir_pos[sample_dt_pos .== 0, :] .= 0.0
	    @assert(all(.~isnan.(density_empir_pos)))
	    density_empir_pre[sample_dt_pre .== 0, :] .= 0.0
	    @assert(all(.~isnan.(density_empir_pre)))

		@assert norm(sum(density_empir_pre[sample_dt_pre .== 1, :], dims=2) .- 1.0) .<= 1e-6
		@assert norm(sum(density_empir_pos[sample_dt_pos .== 1, :], dims=2) .- 1.0) .<= 1e-6


	# MODEL-predicted histogram - n_resp x n_bins
	    ### compute model histograms (sum within hdgrid)
	        density_model_wcha = hist_subgrid(SIMMODEL["dt_wcha"], hdgrid, hd_bins)
	        density_model_free = hist_subgrid(SIMMODEL["dt_free"], hdgrid, hd_bins)

		### no zeros
			check1 = all(sum(density_model_wcha, dims=2) .> 0.0)
			check2 = all(sum(density_model_free, dims=2) .> 0.0)
			if ~check1 || ~check2
				println(theta)
				println(theta)
				println(theta)
				error("zero departure time choice probs")
			end

		### when using larger DT grid, need to normalize to this subsample!
			density_model_wcha = density_model_wcha ./ sum(density_model_wcha, dims=2)
			density_model_free = density_model_free ./ sum(density_model_free, dims=2)

		### ensure no NANs
			@assert norm(sum(density_model_wcha, dims=2) .- 1.0) <= 1e-6
			@assert norm(sum(density_model_free, dims=2) .- 1.0) <= 1e-6

	## DEPARTURE TIME MOMENTS
	    # four cells: treat/control pre/post
	    n_treat_pre = sum(sample_dt_treat .* sample_dt_pre)
	    n_treat_pos = sum(sample_dt_treat .* sample_dt_pos)
	    n_contr_pre = sum(sample_dt_contr .* sample_dt_pre)
	    n_contr_pos = sum(sample_dt_contr .* sample_dt_pos)

	    # treat pre
	    treat_pre_data  = density_empir_pre  .* sample_dt_treat .* sample_dt_pre * n_resp_all / n_treat_pre
		treat_pre_model = density_model_free .* sample_dt_treat .* sample_dt_pre * n_resp_all / n_treat_pre

	    # control pre
	    contr_pre_data  = density_empir_pre  .* sample_dt_contr .* sample_dt_pre * n_resp_all / n_contr_pre
	    contr_pre_model = density_model_free .* sample_dt_contr .* sample_dt_pre * n_resp_all / n_contr_pre

	    # treat post
	    treat_post_data  = density_empir_pos  .* sample_dt_treat .* sample_dt_pos * n_resp_all / n_treat_pos
	    treat_post_model = density_model_wcha .* sample_dt_treat .* sample_dt_pos * n_resp_all / n_treat_pos

	    # control post
	    contr_post_data  = density_empir_pos  .* sample_dt_contr .* sample_dt_pos * n_resp_all / n_contr_pos
		contr_post_model = density_model_wcha .* sample_dt_contr .* sample_dt_pos * n_resp_all / n_contr_pos

	    # assign moments - implicit DD
	    # @assert(all(moms[:,1:n_hd_moms] .== 0))

		dd = ((treat_post_data .- treat_pre_data) .- prob_respond_alpha .* (treat_post_model .- treat_pre_model)) .-
			 ((contr_post_data .- contr_pre_data) .- prob_respond_alpha .* (contr_post_model .- contr_pre_model))

		moms_dt_main = dd[:, n_hd_mom_l:n_hd_mom_h]
	    # moms[:,1:n_hd_moms] = dd[:, n_hd_mom_l:n_hd_mom_h]
# end



## construct area moments
if include_area_main
    # control
    temp1 = (SIMMODEL["probs_route0_free"] - cha_ct23) .* sample_ct23 * n_resp_all / n_resp_a_ct23
    # moms[1:n_resp_a, end-1] = temp

    # treat
    temp2 = (SIMMODEL["probs_route0_wcha"] * prob_respond_alpha +
            SIMMODEL["probs_route0_free"] * (1 - prob_respond_alpha) - cha_tr) .*
            sample_tr * n_resp_all / n_resp_a_tr
    # moms[1:n_resp_a, end] = temp

	moms_area_main = hcat(temp1, temp2)
	moms_area_main = vcat(moms_area_main, zeros(n_resp_all - n_resp_a, 2))
end

## combine moments
	moms = moms_dt_main
	include_area_main && (moms = hcat(moms, moms_area_main))

	return moms
end


"""
Model with main DT moments only
49	(1) DT shares
"""
function moments_logitmodel_all(;
		theta::Vector{Float64},
		param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}},
	    SIMMODEL::Dict{String, Any},
	    DATAMOMS::Dict{String, Any}
		)

## Unpack

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])
	n_hd_bins = DATAMOMS["n_hd_bins"]
	bin_h, bin_l, n_hd_mom_l, n_hd_mom_h = DATAMOMS["hd_bins_hl"]

	n_hd_mom_l = Int64(n_hd_mom_l)
	n_hd_mom_h = Int64(n_hd_mom_h)
	n_hd_moms = (n_hd_mom_h-n_hd_mom_l+1)

	hdgrid = DATAMOMS["hdgrid"]

	sample_dt_pre = DATAMOMS["sample_dt_pre"]
	sample_dt_pos = DATAMOMS["sample_dt_pos"]
	sample_dt_treat = DATAMOMS["sample_dt_treat"]
	sample_dt_contr = DATAMOMS["sample_dt_contr"]

	mean_dt_pre   = DATAMOMS["mean_dt_pre"]
	mean_dt_pos   = DATAMOMS["mean_dt_pos"]

	# idfxnpre = DATAMOMS["idfxnpre"]
	# idfxnpos = DATAMOMS["idfxnpos"]
	# idfxbe   = DATAMOMS["idfxbe"]
	# idfx     = DATAMOMS["idfx"]

	shadow_rates = DATAMOMS["shadow_rates"]

	#
	prob_respond_alpha = get_params_fn(theta, param_factors)["prob_respond_alpha"]


## DEPARTURE TIME MOMS

    # for each bin:
	    bin_delta = (bin_h-bin_l)/(n_hd_bins - 1)
	    hd_bins = range(bin_l - bin_delta/2, stop=bin_h + bin_delta/2, length=n_hd_bins + 1)
		hd_bins = Array(hd_bins)

	# EMPIRICAL - POST SAMPLE
		density_empir_pos = mean_dt_pos #./ 100.0
	# EMPIRICAL - PRE SAMPLE
		density_empir_pre = mean_dt_pre #./ 100.0

	    density_empir_pos[sample_dt_pos .== 0, :] .= 0.0
	    @assert(all(.~isnan.(density_empir_pos)))
	    density_empir_pre[sample_dt_pre .== 0, :] .= 0.0
	    @assert(all(.~isnan.(density_empir_pre)))

		@assert norm(sum(density_empir_pre[sample_dt_pre .== 1, :], dims=2) .- 1.0) .<= 1e-6
		@assert norm(sum(density_empir_pos[sample_dt_pos .== 1, :], dims=2) .- 1.0) .<= 1e-6


	# MODEL-predicted histogram - n_resp x n_bins
	    ### compute model histograms (sum within hdgrid)
	        density_model_wcha = hist_subgrid(SIMMODEL["dt_wcha_1route"], hdgrid, hd_bins)
	        density_model_free = hist_subgrid(SIMMODEL["dt_free_1route"], hdgrid, hd_bins)

		### no zeros
			@assert all(sum(density_model_wcha, dims=2) .> 0.0)
			@assert all(sum(density_model_free, dims=2) .> 0.0)

		### when using larger DT grid, need to normalize to this subsample!
			density_model_wcha = density_model_wcha ./ sum(density_model_wcha, dims=2)
			density_model_free = density_model_free ./ sum(density_model_free, dims=2)

		### ensure no NANs
			@assert norm(sum(density_model_wcha, dims=2) .- 1.0) <= 1e-6
			@assert norm(sum(density_model_free, dims=2) .- 1.0) <= 1e-6

	## DEPARTURE TIME MOMENTS
	    # four cells: treat/control pre/post
	    n_treat_pre = sum(sample_dt_treat .* sample_dt_pre)
	    n_treat_pos = sum(sample_dt_treat .* sample_dt_pos)
	    n_contr_pre = sum(sample_dt_contr .* sample_dt_pre)
	    n_contr_pos = sum(sample_dt_contr .* sample_dt_pos)

	    # treat pre
	    treat_pre_data  = density_empir_pre  .* sample_dt_treat .* sample_dt_pre * n_resp_all / n_treat_pre
		treat_pre_model = density_model_free .* sample_dt_treat .* sample_dt_pre * n_resp_all / n_treat_pre

	    # control pre
	    contr_pre_data  = density_empir_pre  .* sample_dt_contr .* sample_dt_pre * n_resp_all / n_contr_pre
	    contr_pre_model = density_model_free .* sample_dt_contr .* sample_dt_pre * n_resp_all / n_contr_pre

	    # treat post
	    treat_post_data  = density_empir_pos  .* sample_dt_treat .* sample_dt_pos * n_resp_all / n_treat_pos
	    treat_post_model = density_model_wcha .* sample_dt_treat .* sample_dt_pos * n_resp_all / n_treat_pos

	    # control post
	    contr_post_data  = density_empir_pos  .* sample_dt_contr .* sample_dt_pos * n_resp_all / n_contr_pos
		contr_post_model = density_model_wcha .* sample_dt_contr .* sample_dt_pos * n_resp_all / n_contr_pos

	    # assign moments - implicit DD
	    # @assert(all(moms[:,1:n_hd_moms] .== 0))

		dd = ((treat_post_data .- treat_pre_data) .- prob_respond_alpha .* (treat_post_model .- treat_pre_model)) .-
			 ((contr_post_data .- contr_pre_data) .- prob_respond_alpha .* (contr_post_model .- contr_pre_model))

		moms = dd[:, n_hd_mom_l:n_hd_mom_h]

	return moms
end


## Data moments only
section_data_moms = 1

"""
dynamic vott data moments
"""
function moments_data_dynamic_vott(;
	    DATAMOMS::Dict{String, Any},
		transition_moment_type=1
	)

## Unpack

	# cha_ct23 = DATAMOMS["cha_ct23"]
	# cha_tr   = DATAMOMS["cha_tr"]
	#
	# sample_ct23   = DATAMOMS["sample_ct23"] .* DATAMOMS["sample_dynamic_vot"]
	# sample_tr     = DATAMOMS["sample_tr"] .* DATAMOMS["sample_dynamic_vot"]
	#
	# nct23 = DATAMOMS["nct23"]
	# ntr   = DATAMOMS["ntr"]

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])
	n_hd_bins = DATAMOMS["n_hd_bins"]
	bin_h, bin_l, n_hd_mom_l, n_hd_mom_h = DATAMOMS["hd_bins_hl"]

	data_a = DATAMOMS["data_a"]

## Area moments (11 moments from VOT GMM)

	# sample dynamic VOT already included in .a_early and .a_late
	@assert all(data_a.a_early .== (data_a.a_early .* data_a.sample_dynamic_vot))

	nmoms = 11
	moms = zeros(n_resp_all, nmoms)
	moms_sample = zeros(n_resp_all, nmoms)

	# moments 1-5 are "late" means in weeks 0-4
	moms_sample[1:n_resp_a, 1] = data_a.a_late .* (data_a.n_obs_0 .> 0)  # late week 0
	moms_sample[1:n_resp_a, 2] = data_a.a_late .* (data_a.n_obs_1 .> 0)  # late week 1
	moms_sample[1:n_resp_a, 3] = data_a.a_late .* (data_a.n_obs_2 .> 0)  # late week 2
	moms_sample[1:n_resp_a, 4] = data_a.a_late .* (data_a.n_obs_3 .> 0)  # late week 3
	moms_sample[1:n_resp_a, 5] = data_a.a_late .* (data_a.n_obs_4 .> 0)  # late week 4

	moms[1:n_resp_a, 1] = (data_a.mean_long_0) .* data_a.a_late .* (data_a.n_obs_0 .> 0)  # late week 0
	moms[1:n_resp_a, 2] = (data_a.mean_long_1) .* data_a.a_late .* (data_a.n_obs_1 .> 0)  # late week 1
	moms[1:n_resp_a, 3] = (data_a.mean_long_2) .* data_a.a_late .* (data_a.n_obs_2 .> 0)  # late week 2
	moms[1:n_resp_a, 4] = (data_a.mean_long_3) .* data_a.a_late .* (data_a.n_obs_3 .> 0)  # late week 3
	moms[1:n_resp_a, 5] = (data_a.mean_long_4) .* data_a.a_late .* (data_a.n_obs_4 .> 0)  # late week 4

	# moments 6-10 are "early" means in weeks 0-4
	moms_sample[1:n_resp_a, 5 + 1] = data_a.a_early .* (data_a.n_obs_0 .> 0)  # early week 0
	moms_sample[1:n_resp_a, 5 + 2] = data_a.a_early .* (data_a.n_obs_1 .> 0)  # early week 1
	moms_sample[1:n_resp_a, 5 + 3] = data_a.a_early .* (data_a.n_obs_2 .> 0)  # early week 2
	moms_sample[1:n_resp_a, 5 + 4] = data_a.a_early .* (data_a.n_obs_3 .> 0)  # early week 3
	moms_sample[1:n_resp_a, 5 + 5] = data_a.a_early .* (data_a.n_obs_4 .> 0)  # early week 4

	moms[1:n_resp_a, 5 + 1] = (data_a.mean_long_0) .* data_a.a_early .* (data_a.n_obs_0 .> 0)  # early week 0
	moms[1:n_resp_a, 5 + 2] = (data_a.mean_long_1) .* data_a.a_early .* (data_a.n_obs_1 .> 0)  # early week 1
	moms[1:n_resp_a, 5 + 3] = (data_a.mean_long_2) .* data_a.a_early .* (data_a.n_obs_2 .> 0)  # early week 2
	moms[1:n_resp_a, 5 + 4] = (data_a.mean_long_3) .* data_a.a_early .* (data_a.n_obs_3 .> 0)  # early week 3
	moms[1:n_resp_a, 5 + 5] = (data_a.mean_long_4) .* data_a.a_early .* (data_a.n_obs_4 .> 0)  # early week 4


if transition_moment_type == 1
	# DATA:
	# switching sample: data in weeks 1 and 2 AND conditional on detour in week 1
	sample_p10_week12 = (data_a.sample_dynamic_vot) .& (data_a.a_early .== 1) .&
						(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0) .&
						(data_a.mean_long_1 .>= 0.5)

	# dummy for switching to route 0 in week 2
	switched_p10_week12 = (data_a.mean_long_1 .>= 0.5) .* (data_a.mean_long_2 .< 0.5)

	# MODEL: vector of transition probabilities from route 1 to route 0 at start of period 2
	# same as transition from t=1 to t=2
	moms_sample[1:n_resp_a, 11] = sample_p10_week12
	moms[1:n_resp_a, 11] = sample_p10_week12 .* switched_p10_week12

elseif transition_moment_type == 2
	sample_p10_week12 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_1 .>= 0.5) .&
						(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0)

	sample_p10_week23 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_2 .>= 0.5) .&
						(data_a.n_obs_2 .> 0) .& (data_a.n_obs_3 .> 0)

	sample_p10_week34 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_3 .>= 0.5) .&
						(data_a.n_obs_3 .> 0) .& (data_a.n_obs_4 .> 0)

	switched_p10_week12 = (data_a.mean_long_1 .>= 0.5) .* (data_a.mean_long_2 .< 0.5)
	switched_p10_week23 = (data_a.mean_long_2 .>= 0.5) .* (data_a.mean_long_3 .< 0.5)
	switched_p10_week34 = (data_a.mean_long_3 .>= 0.5) .* (data_a.mean_long_4 .< 0.5)

	moms_sample[1:n_resp_a, 11] = sample_p10_week12 .| sample_p10_week23 .| sample_p10_week34
	moms[1:n_resp_a, 11] = (sample_p10_week12 .* switched_p10_week12 .+
							sample_p10_week23 .* switched_p10_week23 .+
							sample_p10_week34 .* switched_p10_week34) ./
							max.(1.0, sample_p10_week12 .+ sample_p10_week23 .+ sample_p10_week34)
elseif transition_moment_type == 3
	sample_p01_week12 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_1 .< 0.5) .&
						(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0)

	sample_p01_week23 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_2 .< 0.5) .&
						(data_a.n_obs_2 .> 0) .& (data_a.n_obs_3 .> 0)

	sample_p01_week34 = (data_a.sample_dynamic_vot) .& (data_a.mean_long_3 .< 0.5) .&
						(data_a.n_obs_3 .> 0) .& (data_a.n_obs_4 .> 0)

	switched_p01_week12 = (data_a.mean_long_1 .< 0.5) .* (data_a.mean_long_2 .>= 0.5)
	switched_p01_week23 = (data_a.mean_long_2 .< 0.5) .* (data_a.mean_long_3 .>= 0.5)
	switched_p01_week34 = (data_a.mean_long_3 .< 0.5) .* (data_a.mean_long_4 .>= 0.5)

	moms_sample[1:n_resp_a, 11] = sample_p01_week12 .| sample_p01_week23 .| sample_p01_week34
	moms[1:n_resp_a, 11] = (sample_p01_week12 .* switched_p01_week12 .+
							sample_p01_week23 .* switched_p01_week23 .+
							sample_p01_week34 .* switched_p01_week34) ./
							max.(1.0, sample_p01_week12 .+ sample_p01_week23 .+ sample_p01_week34)

end

	return moms_sample, moms
end


"""
Moments for the full GMM:
	(1) DT shares (from old GMM)
	(2) DT heterogeneity (from old GMM)
	(3) 5+5 weekly area shares (from VOT GMM)
	(4) 2 moments: VOT response in week 1 heterogeneity by areadetourtime (from VOT GMM)
"""
function moments_model_dynamic_vott(;
		theta::Vector{Float64},
		param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}},
	    SIMMODEL::Dict{String, Any},
	    DATAMOMS::Dict{String, Any}
	)

## Unpack
	prob_respond_alpha = get_params_fn(theta, param_factors)["prob_respond_alpha"]
	data_a = DATAMOMS["data_a"]

## Area moments (12 moments from VOT GMM)

	# n_resp_a x 5 matrix
	probs_route1_free = 1 .- SIMMODEL["probs_route0_free"]
	probs_route1_wcha = 1 .- SIMMODEL["probs_route0_wcha"]

	probs_route1_attenuated = probs_route1_wcha .* prob_respond_alpha .+
							  probs_route1_free .* (1 .- prob_respond_alpha)

	# sample dynamic VOT already included in .a_early and .a_late
	@assert all(data_a.a_early .== (data_a.a_early .* data_a.sample_dynamic_vot))

	nmoms = 11
	moms = zeros(n_resp_all, nmoms)
	moms_sample = zeros(n_resp_all, nmoms)

	# moments 1-5 are "late" means in weeks 0-4
	moms_sample[1:n_resp_a, 1] = data_a.a_late .* (data_a.n_obs_0 .> 0)  # late week 0
	moms_sample[1:n_resp_a, 2] = data_a.a_late .* (data_a.n_obs_1 .> 0)  # late week 1
	moms_sample[1:n_resp_a, 3] = data_a.a_late .* (data_a.n_obs_2 .> 0)  # late week 2
	moms_sample[1:n_resp_a, 4] = data_a.a_late .* (data_a.n_obs_3 .> 0)  # late week 3
	moms_sample[1:n_resp_a, 5] = data_a.a_late .* (data_a.n_obs_4 .> 0)  # late week 4

	# moments 1-5 are "late" means in weeks 0-4
	moms[1:n_resp_a, 1] = probs_route1_attenuated[:, 1] .* data_a.a_late .* (data_a.n_obs_0 .> 0)  # late week 0
	moms[1:n_resp_a, 2] = probs_route1_attenuated[:, 2] .* data_a.a_late .* (data_a.n_obs_1 .> 0)  # late week 1
	moms[1:n_resp_a, 3] = probs_route1_attenuated[:, 3] .* data_a.a_late .* (data_a.n_obs_2 .> 0)  # late week 2
	moms[1:n_resp_a, 4] = probs_route1_attenuated[:, 4] .* data_a.a_late .* (data_a.n_obs_3 .> 0)  # late week 3
	moms[1:n_resp_a, 5] = probs_route1_attenuated[:, 5] .* data_a.a_late .* (data_a.n_obs_4 .> 0)  # late week 4

	# moments 6-10 are "early" means in weeks 0-4
	moms_sample[1:n_resp_a, 5 + 1] = data_a.a_early .* (data_a.n_obs_0 .> 0)  # early week 0
	moms_sample[1:n_resp_a, 5 + 2] = data_a.a_early .* (data_a.n_obs_1 .> 0)  # early week 1
	moms_sample[1:n_resp_a, 5 + 3] = data_a.a_early .* (data_a.n_obs_2 .> 0)  # early week 2
	moms_sample[1:n_resp_a, 5 + 4] = data_a.a_early .* (data_a.n_obs_3 .> 0)  # early week 3
	moms_sample[1:n_resp_a, 5 + 5] = data_a.a_early .* (data_a.n_obs_4 .> 0)  # early week 4

	# moments 6-10 are "early" means in weeks 0-4
	moms[1:n_resp_a, 5 + 1] = probs_route1_attenuated[:, 1] .* data_a.a_early .* (data_a.n_obs_0 .> 0)  # early week 0
	moms[1:n_resp_a, 5 + 2] = probs_route1_attenuated[:, 2] .* data_a.a_early .* (data_a.n_obs_1 .> 0)  # early week 1
	moms[1:n_resp_a, 5 + 3] = probs_route1_attenuated[:, 3] .* data_a.a_early .* (data_a.n_obs_2 .> 0)  # early week 2
	moms[1:n_resp_a, 5 + 4] = probs_route1_attenuated[:, 4] .* data_a.a_early .* (data_a.n_obs_3 .> 0)  # early week 3
	moms[1:n_resp_a, 5 + 5] = probs_route1_attenuated[:, 5] .* data_a.a_early .* (data_a.n_obs_4 .> 0)  # early week 4


	# switching sample: data in weeks 1 and 2 AND conditional on detour in week 1
	sample_p10_week12 = (data_a.sample_dynamic_vot) .& (data_a.a_early .== 1) .&
						(data_a.n_obs_1 .> 0) .& (data_a.n_obs_2 .> 0) .&
						(data_a.mean_long_1 .>= 0.75)

	# dummy for switching to route 0 in week 2
	# switched_p10_week12 = (data_a.mean_long_1 .>= 0.5) .* (data_a.mean_long_2 .< 0.5)

	# MODEL: vector of transition probabilities from route 1 to route 0 at start of period 2
	# same as transition from t=1 to t=2
	moms_sample[1:n_resp_a, 11] = sample_p10_week12
	moms[1:n_resp_a, 11] = sample_p10_week12 .* SIMMODEL["w2_soln"].p10

	return moms_sample, moms
end
