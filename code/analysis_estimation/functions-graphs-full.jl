
function draw_moms_dt(;
	theta, param_factors, DATAMOMS, SIMMODEL,
	mymodelcolor=:red,
	mymodellabel="Model",
	mymodelstyle=:solid,
	drawdata=true,
	display_plot=true,
	save_path="",
	control_graph=true,
	display_plot2=true,
	save_path2="")

## Unpack

	cha_ct23 = DATAMOMS["cha_ct23"]
	cha_tr   = DATAMOMS["cha_tr"]

	sample_ct23   = DATAMOMS["sample_ct23"] .* DATAMOMS["sample_dynamic_vot"]
	sample_tr     = DATAMOMS["sample_tr"] .* DATAMOMS["sample_dynamic_vot"]

	nct23 = DATAMOMS["nct23"]
	ntr   = DATAMOMS["ntr"]

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])
	n_hd_bins = DATAMOMS["n_hd_bins"]
	bin_h, bin_l = DATAMOMS["hd_bins_hl"]
	hdgrid = DATAMOMS["hdgrid"]
	hdgrid_2525 = DATAMOMS["hdgrid_2525"]

	sample_dt_pre = DATAMOMS["sample_dt_pre"]
	sample_dt_pos = DATAMOMS["sample_dt_pos"]
	sample_dt_treat = DATAMOMS["sample_dt_treat"]
	sample_dt_contr = DATAMOMS["sample_dt_contr"]

	mean_dt_pre   = DATAMOMS["mean_dt_pre"]
	mean_dt_pos   = DATAMOMS["mean_dt_pos"]

	idfxnpre = DATAMOMS["idfxnpre"]
	idfxnpos = DATAMOMS["idfxnpos"]
	idfxbe   = DATAMOMS["idfxbe"]
	idfx     = DATAMOMS["idfx"]

	shadow_rates = DATAMOMS["shadow_rates"]

	#
	theta_dict = get_params_fn(theta, param_factors)
	prob_respond_alpha = theta_dict["prob_respond_alpha"]

	data_a = DATAMOMS["data_a"]

	# moms = zeros(n_resp_all, n_moms)


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
		density_model_wcha = hist_subgrid(SIMMODEL["dt_wcha"], hdgrid, hd_bins)
		density_model_free = hist_subgrid(SIMMODEL["dt_free"], hdgrid, hd_bins)
		density_model_before = hist_subgrid(SIMMODEL["dt_before"], hdgrid, hd_bins)

	### no zeros
		@assert all(sum(density_model_wcha, dims=2) .> 0.0)
		@assert all(sum(density_model_free, dims=2) .> 0.0)
		@assert all(sum(density_model_before, dims=2) .> 0.0)

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

## plot DATA
if drawdata
	dd_data = (treat_post_data .- treat_pre_data) .- (contr_post_data .- contr_pre_data)
	dd_data = vec(mean(dd_data, dims=1))

	Plots.scatter(hdgrid_2525, dd_data,
		label="Data",
		markercolor=:blue, markerstrokecolor=:blue)
	Plots.plot!(hdgrid_2525, dd_data,
		label="",
		linestyle=:dash,
		linewidth=1,
		color=:blue)

	model = loess(hdgrid_2525, dd_data, span=0.3)
	dd_data_smooth = Loess.predict(model, hdgrid_2525)

	Plots.plot!(hdgrid_2525, dd_data_smooth,
			linewidth=3, color=:blue, label="Data (smoothed)") #|> display
end
## plot MODEL
	dd_model = (treat_post_model .- treat_pre_model) .- (contr_post_model .- contr_pre_model)
	dd_model = vec(mean(dd_model, dims=1))

	myplotobj = Plots.plot!(
			hdgrid_2525, dd_model,
			color=mymodelcolor,
			label=mymodellabel,
			linestyle=mymodelstyle,
			linewidth=4.5,
			linealpha=0.5,
			thickness_scaling = 1.1,
			# gridlinewidth=1.0, gridalpha=0.33,
			xticks=[-2.5, -1.5, -0.5, 0.5, 1.5, 2.5],
			xlabel="Trip departure time",
			ylabel="Trip probability")

	if save_path != ""
		savefig(save_path)
	end

	if display_plot
		myplotobj |> display
	end

	# model = loess(hdgrid, dd_model, span=0.3)
	# dd_model_smooth = Loess.predict(model, hdgrid)
	#
	# Plots.plot!(hdgrid, dd_model_smooth, linewidth=1.5, color="red") |> display

## Control Post
if control_graph
	## plot DATA
	if drawdata
		dd_data = vec(mean(contr_post_data, dims=1))

		Plots.scatter(hdgrid_2525, dd_data,
			label="Data",
			markercolor=:blue, markerstrokecolor=:blue)
		Plots.plot!(hdgrid_2525, dd_data,
			label="",
			linestyle=:dash,
			linewidth=1,
			color=:blue)

		model = loess(hdgrid_2525, dd_data, span=0.3)
		dd_data_smooth = Loess.predict(model, hdgrid_2525)

		Plots.plot!(hdgrid_2525, dd_data_smooth,
				linewidth=3, color=:blue, label="Data (smoothed)") #|> display
	end
## plot MODEL
	dd_model = vec(mean(contr_post_model, dims=1))

	myplotobj2 = Plots.plot!(
			hdgrid_2525, dd_model,
			color=:red, label=mymodellabel,
			linewidth=4.5,
			linealpha=0.5,
			thickness_scaling = 1.1,
			# gridlinewidth=1.0, gridalpha=0.33,
			xticks=[-2.5, -1.5, -0.5, 0.5, 1.5, 2.5],
			xlabel="Trip departure time",
			ylabel="Trip probability")

	if save_path2 != ""
		savefig(save_path2)
	end

	if display_plot2
		myplotobj2 |> display
	end
end
end


function draw_moms_dt_het(;theta, DATAMOMS, SIMMODEL, drawdata=true)

## Unpack

	cha_ct23 = DATAMOMS["cha_ct23"]
	cha_tr   = DATAMOMS["cha_tr"]

	sample_ct23   = DATAMOMS["sample_ct23"] .* DATAMOMS["sample_dynamic_vot"]
	sample_tr     = DATAMOMS["sample_tr"] .* DATAMOMS["sample_dynamic_vot"]

	nct23 = DATAMOMS["nct23"]
	ntr   = DATAMOMS["ntr"]

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])
	n_hd_bins = DATAMOMS["n_hd_bins"]
	bin_h, bin_l = DATAMOMS["hd_bins_hl"]
	hdgrid = DATAMOMS["hdgrid"]

	sample_dt_pre = DATAMOMS["sample_dt_pre"]
	sample_dt_pos = DATAMOMS["sample_dt_pos"]
	sample_dt_treat = DATAMOMS["sample_dt_treat"]
	sample_dt_contr = DATAMOMS["sample_dt_contr"]

	mean_dt_pre   = DATAMOMS["mean_dt_pre"]
	mean_dt_pos   = DATAMOMS["mean_dt_pos"]

	idfxnpre = DATAMOMS["idfxnpre"]
	idfxnpos = DATAMOMS["idfxnpos"]
	idfxbe   = DATAMOMS["idfxbe"]
	idfx     = DATAMOMS["idfx"]

	idfx = convert.(Float64, idfx)


	shadow_rates = DATAMOMS["shadow_rates"]

	#
	prob_respond_alpha = get_params_fn(theta)["prob_respond_alpha"]

	data_a = DATAMOMS["data_a"]

	moms = zeros(n_resp_all, n_moms)

##
	# (168) respondents with pre and post data (treat/control)
	treat_idxs = (idfxnpre .* idfxnpos .* sample_dt_treat) .> 0
	contr_idxs = (idfxnpre .* idfxnpos .* sample_dt_contr) .> 0
	n_treat = sum(treat_idxs)
	n_contr = sum(contr_idxs)

	### Model predicted charges (post-pre), mean and SD
		# this is the set of departure times with rising slope
		hdgrid_early = (-2.0 .<= hdgrid .<= -0.0)
		n_early = sum(hdgrid_early)

	# probabilities to choose these departure times (w and w/o charges)
	# conditional, so normalize
		# size n_resp_all*nrc
		dt_nrc_free_early = SIMMODEL["dt_free"][:, hdgrid_early]
		# normalize and replace NAN with zero
		dt_nrc_free_early = dt_nrc_free_early ./ sum(dt_nrc_free_early, dims=2)
		dt_nrc_free_early[isnan.(dt_nrc_free_early)] .= 0.0


		dt_nrc_wcha_early = SIMMODEL["dt_wcha"][:, hdgrid_early]
		# normalize and replace NAN with zero
		dt_nrc_wcha_early = dt_nrc_wcha_early ./ sum(dt_nrc_wcha_early, dims=2)
		dt_nrc_wcha_early[isnan.(dt_nrc_wcha_early)] .= 0.0

	# charges long format
		charge_mat = repeat(shadow_rates[hdgrid_early], 1, n_resp_all) |> transpose

	### one-day charges: expectation and variance
	# expectation = E charges = \sum prob * charge
	mean_free = sum(dt_nrc_free_early .* charge_mat, dims=2)
	mean_wcha = sum(dt_nrc_wcha_early .* charge_mat, dims=2)

	# variance = E(residual) = E (charge - expected(charge))^2
	# residual
	res_free = charge_mat .- mean_free
	res_wcha = charge_mat .- mean_wcha

	# = \sum prob * residual^2
	var_free = sum(dt_nrc_free_early .* (res_free.^2), dims=2)
	var_wcha = sum(dt_nrc_wcha_early .* (res_wcha.^2), dims=2)


	### post-pre individual effect
	# model the sum of N independent random variables as normal, use E and var from above
	mean_idfx_wcha = mean_wcha .- mean_free
	mean_idfx_free = mean_wcha .* 0.0

	# with charges in post, free in pre
	var_idfx_wcha = var_wcha ./ idfxnpos .+ var_free ./ idfxnpre
	var_idfx_free = var_free .* (1 ./ idfxnpos .+ 1 ./ idfxnpre)
	var_idfx_wcha[treat_idxs .+ contr_idxs .== 0] .= 0.0
	var_idfx_free[treat_idxs .+ contr_idxs .== 0] .= 0.0
	@assert(all(.~isnan.(var_idfx_wcha)))
	@assert(all(.~isnan.(var_idfx_free)))

	### MOMENT 1 - mean - includes partial compliance in treatment
	# average over nrc
	model_mean_idfx = (mean_idfx_wcha .* prob_respond_alpha .+
					   mean_idfx_free .* (1.0 .- prob_respond_alpha)) .* treat_idxs .+
					   mean_idfx_free .* contr_idxs

	moment_mean_idfx = (model_mean_idfx .* treat_idxs ./ n_treat .-
								  (idfx .* treat_idxs ./ n_treat .- idfx .* contr_idxs ./ n_contr))
	moment_mean_idfx = n_resp_all * moment_mean_idfx


	### MOMENT 2 - variance

	# model in control (and treatment not responding)
	# variance only from sampling variance (on avearge = 0)
	model_var_idfx_contr = var_idfx_free

	# model in treat (responding)
	# variance
	mean_mean_idfx_treat = mean(mean_idfx_wcha[treat_idxs])
	model_var_idfx_treat = var_idfx_wcha .+ (mean_idfx_wcha .- mean_mean_idfx_treat) .^ 2

	# data
	idfx_treat_resid = (idfx .- mean(idfx[treat_idxs])) .^2
	idfx_contr_resid = (idfx .- mean(idfx[contr_idxs])) .^2

	# build moment
	moment_var_idfx_treat =
		((model_var_idfx_treat .* prob_respond_alpha .+
		  model_var_idfx_contr .* (1.0 .- prob_respond_alpha) ) .* treat_idxs ./ n_treat .-
			  idfx_treat_resid .* treat_idxs ./ n_treat)
	moment_var_idfx_treat = n_resp_all .* moment_var_idfx_treat

	moment_var_idfx_contr =
		(model_var_idfx_contr .* contr_idxs ./ n_contr .-
			 idfx_contr_resid .* contr_idxs ./ n_contr)
	moment_var_idfx_contr = n_resp_all * moment_var_idfx_contr

	## MOMENTS
	# mean tr-control
	# @assert(all(moms[:, n_hd_bins + 1] .== 0))

	# variance treatment
	# @assert(all(moms[:, n_hd_moms + 1] .== 0))
	# moms[:, n_hd_moms + 1] = (moment_var_idfx_treat) ./ 10000.0

	# variance control
	# @assert(all(moms[:, n_hd_moms + 2] .== 0))
	# moms[:, n_hd_moms + 2] = moment_var_idfx_contr ./ 10000.0

## Draw

	# simulate model idfx draws using mean and SD
	mean_contr = mean_idfx_free[contr_idxs]
	sd_contr   = sqrt.(var_idfx_free[contr_idxs])

	mean_treat_resp = mean_idfx_wcha[treat_idxs]
	  sd_treat_resp = sqrt.(var_idfx_wcha[treat_idxs])
	mean_treat_nors = mean_idfx_free[treat_idxs]
	  sd_treat_nors = sqrt.(var_idfx_free[treat_idxs])

	  # Draw random ns
	draws_treat_resp = randn(size(mean_treat_resp)) .* sd_treat_resp .+ mean_treat_resp
	draws_treat_nors = randn(size(mean_treat_nors)) .* sd_treat_nors .+ mean_treat_nors
	draws_contr = randn(size(mean_contr)) .* sd_contr .+ mean_contr

	# % model
	# [fm_tr_resp, xim_tr] = ksdensity(draws_treat_resp,         'Kernel', 'epanechnikov', 'bandwidth', 6);
	# [fm_tr_nors, xim_tr] = ksdensity(draws_treat_nors, xim_tr, 'Kernel', 'epanechnikov', 'bandwidth', 6);
	# [fm_ct, xim_ct] = ksdensity(draws_contr, 'Kernel', 'epanechnikov', 'bandwidth', 6);

	# % data
	# [f_tr, xi_tr] = ksdensity(idfx[treat_idxs], 'Kernel', 'epanechnikov', 'bandwidth', 6);
	# [f_ct, xi_ct] = ksdensity(idfx[contr_idxs], 'Kernel', 'epanechnikov', 'bandwidth', 6);

	U_model_treat_resp = kde(draws_treat_resp, bandwidth=5, boundary=(-100.0,100.0), npoints=2048)
	U_model_treat_nors = kde(draws_treat_nors, bandwidth=5, boundary=(-100.0,100.0), npoints=2048)
	U_model_contr = kde(draws_contr, bandwidth=5)

	# % data
	U_data_treat = kde(idfx[treat_idxs], bandwidth=5)
	U_data_contr = kde(idfx[contr_idxs], bandwidth=5)

	# plot control
	Plots.plot(U_data_contr.x, U_data_contr.density, label="data control")
	Plots.plot!(U_model_contr.x, U_model_contr.density, label="model control") |> display

	# plot treatment
	Plots.plot(U_data_treat.x, U_data_treat.density, label="data treatment")

	umodelx = U_model_treat_resp.x
	umodeld = U_model_treat_resp.density .* prob_respond_alpha .+
			  U_model_treat_nors.density .* (1.0 .- prob_respond_alpha)

	Plots.plot!(umodelx, umodeld, label="model treatment") |> display

	# display(U_model_treat_nors.x)
	# display(U_model_treat_resp.x)

	# fdsfd
	# figure
	# hold on
	# grid on
	# plot(xi_tr, f_tr, 'b--', 'linewidth', 1.5);
	# plot(xim_tr, prob_respond_alpha * fm_tr_resp + (1-prob_respond_alpha) * fm_tr_nors, 'r', 'linewidth', 1.5);
	# legend('Data', 'Model', 'Location','northeast')
	# xlabel('Individual-level Change in Shadow Charges (Rs.)')
	# set(gca,'fontsize', 14)
	#
	# # if debug==2
	# # 	saveas(gca, [figsavepath 'dep_time_idfx_treat.eps'],'epsc');
	# # end
	#
	# figure
	# hold on
	# grid on
	# plot(xi_ct, f_ct, 'b--', 'linewidth', 1.5);
	# plot(xim_ct, fm_ct, 'r', 'linewidth', 1.5);
	# legend('Data', 'Model', 'Location','northeast')
	# xlabel('Individual-level Change in Shadow Charges (Rs.)')
	# set(gca,'fontsize', 14)

	# if debug==2
	# 	saveas(gca, [figsavepath 'dep_time_idfx_contr.eps'],'epsc');
	# end

end

function draw_moms_vot_het(;theta, DATAMOMS, SIMMODEL)

## UNPACK

	# SIMMODEL = sim_model
	# theta=theta_optimal

	# probability to choose SHORT route
	cha_ct23 = DATAMOMS["cha_ct23"]
	cha_tr   = DATAMOMS["cha_tr"]

	sample_ct23   = DATAMOMS["sample_ct23"] .* DATAMOMS["sample_dynamic_vot"]
	sample_tr     = DATAMOMS["sample_tr"] .* DATAMOMS["sample_dynamic_vot"]

	nct23 = DATAMOMS["nct23"]
	ntr   = DATAMOMS["ntr"]

	n_resp_all, n_resp_a, n_resp_a_ct23, n_resp_a_tr = Int64.(DATAMOMS["n_resp_vec"])
	n_hd_bins = DATAMOMS["n_hd_bins"]
	bin_h, bin_l = DATAMOMS["hd_bins_hl"]
	hdgrid = DATAMOMS["hdgrid"]

	sample_dt_pre = DATAMOMS["sample_dt_pre"]
	sample_dt_pos = DATAMOMS["sample_dt_pos"]
	sample_dt_treat = DATAMOMS["sample_dt_treat"]
	sample_dt_contr = DATAMOMS["sample_dt_contr"]

	mean_dt_pre   = DATAMOMS["mean_dt_pre"]
	mean_dt_pos   = DATAMOMS["mean_dt_pos"]

	idfxnpre = DATAMOMS["idfxnpre"]
	idfxnpos = DATAMOMS["idfxnpos"]
	idfxbe   = DATAMOMS["idfxbe"]
	idfx     = DATAMOMS["idfx"]

	probs_route1_free = 1 .- SIMMODEL["probs_route0_free"]
	probs_route1_wcha = 1 .- SIMMODEL["probs_route0_wcha"]

	data_a = DATAMOMS["data_a"]

	prob_respond_alpha = get_params_fn(theta)["prob_respond_alpha"]

##
	## Week 0

	## weeks 23 (w/o charges)
		temp = (probs_route1_free[:,3] .+ probs_route1_free[:,4]) ./ 2
		area_ct23_free = 1 .- temp[sample_ct23] |> vec

	## weeks 23 (with charges)
		temp = (probs_route1_wcha[:,3] .+ probs_route1_wcha[:,4]) ./ 2
		area_ct23_wcha = 1 .- temp[sample_ct23] |> vec

	## Weeks 1 and 4 (treated)
		temp = probs_route1_free[:,2] .* data_a.a_early .+
			   probs_route1_free[:,5] .* data_a.a_late
		area_tr_free   = 1 .- temp[sample_tr]  |> vec

	## treatment but no response (week 0)
		# area_wcha_trnr = 1 .- probs_route1[sample_tr,1]   |> vec
		temp = probs_route1_wcha[:,2] .* data_a.a_early .+
			   probs_route1_wcha[:,5] .* data_a.a_late
		area_tr_wcha   = 1 .- temp[sample_tr]  |> vec

	# checks
	# if ~all(-1e-6 .<= area_free_ct23 .<= 1.0 .+ 1e-6) |
	#    ~all(-1e-6 .<= area_wcha_tr .<= 1.0 .+ 1e-6) |
	#    ~all(-1e-6 .<= area_wcha_trnr .<= 1.0 .+ 1e-6)
	#    println("ERROR ASSERT area bounds theta=", theta)
	#    println("BRR ", theta)
	#    println("BRR ", theta)
	#    println("BRR ", theta)
	#    println("BRR ", theta)
	#    println("BRR ", theta)
	#    println("BRR ", theta)
	#    println("BRR ", theta)
	#    println(maximum(area_free_ct23))
	#    println(maximum(area_wcha_tr))
	#    println(maximum(area_wcha_trnr))
	#    display(theta)
	#    # display(SIMMODEL["r0_free"])
	#    @assert all(0.0 .<= area_free_ct23 .<= 1.0)
	#    @assert all(0.0 .<= area_wcha_tr .<= 1.0)
	#    @assert all(0.0 .<= area_wcha_trnr .<= 1.0)
   # end

   # make sure EXACTLY between 0.0 and 1.0
   area_ct23_free = max.(0.0, min.(1.0, area_ct23_free))
   area_ct23_wcha = max.(0.0, min.(1.0, area_ct23_wcha))
   area_tr_free   = max.(0.0, min.(1.0, area_tr_free))
   area_tr_wcha   = max.(0.0, min.(1.0, area_tr_wcha))

	# probability to choose area - data
		data_tr   = cha_tr[sample_tr]
		data_ct23 = cha_ct23[sample_ct23]

	# the corresponding number of observations - control
		nct23_nrc = nct23[sample_ct23]
		# nct23_nrc = kron(nct23_nrc, ones(nrc, 1))

		# the corresponding number of observations - treat
		ntr_nrc = ntr[sample_tr]
		# ntr_nrc = kron(ntr_nrc, ones(nrc,1));

	### compute discrete probabilities (based on n days)
	# maximum number of days in control
		maxn_ct23 = maximum(nct23_nrc)
		maxn1_ct23 = maxn_ct23 + 1
	# maximum number of days in treatment
		maxn_tr = maximum(ntr_nrc)
		maxn1_tr = maxn_tr + 1

	# prob to observe each frequency - based on n days and binomial dist'n
	# # Binomial(n, p)

		# control grid, and freqs

		# ct23 free
		grid_ct23_free_nrc = repeat((0:maxn_ct23), 1, length(nct23_nrc)) |> transpose |> Matrix
		nreps_ct23 = repeat(nct23_nrc, 1, maxn1_ct23)
		probs_ct23_free = repeat(area_ct23_free, 1, maxn1_ct23)
		bindists = Binomial.(nreps_ct23, probs_ct23_free)
		model_ct23_free_nrc = pdf.(bindists, grid_ct23_free_nrc)
		grid_ct23_free_nrc = grid_ct23_free_nrc ./ nreps_ct23

		# ct23 wcha
		grid_ct23_wcha_nrc = repeat((0:maxn_ct23), 1, length(nct23_nrc)) |> transpose |> Matrix
		# nreps_ct23 = repeat(nct23_nrc, 1, maxn1_ct23)
		probs_ct23_wcha = repeat(area_ct23_wcha, 1, maxn1_ct23)
		bindists = Binomial.(nreps_ct23, probs_ct23_wcha)
		model_ct23_wcha_nrc = pdf.(bindists, grid_ct23_wcha_nrc)
		grid_ct23_wcha_nrc = grid_ct23_wcha_nrc ./ nreps_ct23

		# tr free
		grid_tr_free_nrc = repeat((0:maxn_tr), 1, length(ntr_nrc)) |> transpose |> Matrix
		nreps_tr = repeat(ntr_nrc, 1, maxn1_tr)
		probs_tr_free = repeat(area_tr_free, 1, maxn1_tr)
		bindists = Binomial.(nreps_tr, probs_tr_free)
		model_tr_free_nrc = pdf.(bindists, grid_tr_free_nrc)
		grid_tr_free_nrc = grid_tr_free_nrc ./ nreps_tr

		# treatment grid and freqs -- no response
		grid_tr_wcha_nrc = repeat((0:maxn_tr), 1, length(ntr_nrc)) |> transpose |> Matrix
		# nreps_tr = repeat(ntr_nrc, 1, maxn1_tr)
		probs_tr_wcha = repeat(area_tr_wcha, 1, maxn1_tr)
		bindists = Binomial.(nreps_tr, probs_tr_wcha)
		model_tr_wcha_nrc = pdf.(bindists, grid_tr_wcha_nrc)
		grid_tr_wcha_nrc = grid_tr_wcha_nrc ./ nreps_tr

	### compute histograms
		# syntax: histw_exp(vals, weights, nbins, dim1)

		# histogram model
		n_area_bins = 3
		# n_area_bins = DATAMOMS["n_area_bins"]
		hm_ct23_free = histw_exp(grid_ct23_free_nrc, model_ct23_free_nrc, n_area_bins)
		hm_ct23_wcha = histw_exp(grid_ct23_wcha_nrc, model_ct23_wcha_nrc, n_area_bins)
		hm_tr_free   = histw_exp(grid_tr_free_nrc,   model_tr_free_nrc,   n_area_bins)
		hm_tr_wcha   = histw_exp(grid_tr_wcha_nrc,   model_tr_wcha_nrc,   n_area_bins)
		@assert norm(sum(hm_ct23_free,dims=2) .- 1.0) <= 1.0e-5
		@assert norm(sum(hm_ct23_wcha,dims=2) .- 1.0) <= 1.0e-5
		@assert norm(sum(hm_tr_free  ,dims=2) .- 1.0) <= 1.0e-5
		@assert norm(sum(hm_tr_wcha  ,dims=2) .- 1.0) <= 1.0e-5

		# histogram for data
		hd_ct23 = histw_exp(reshape(data_ct23, length(data_ct23), 1), ones(length(data_ct23),1), n_area_bins)
		hd_tr   = histw_exp(reshape(data_tr  , length(data_tr  ), 1), ones(length(data_tr),  1), n_area_bins)
		@assert norm(sum(hd_ct23,dims=2) .- 1) <= 1.0e-5
		@assert norm(sum(hd_tr  ,dims=2) .- 1) <= 1.0e-5

	### MOMENTS
	# control going backwards
	ahet_ct_moms = zeros(n_resp_a, n_area_bins)
	ahet_ct_moms[sample_ct23, :] =
		(prob_respond_alpha .* hm_ct23_wcha .+
		(1 .- prob_respond_alpha) .* hm_ct23_free .- hd_ct23) .* n_resp_all ./ n_resp_a_ct23

	# ix1 = size(moms,2) + 1 - (n_area_bins-1)
	# ix2 = size(moms,2)
	# @assert all(moms[:,ix1:ix2] .== 0)
	# moms[1:n_resp_a,ix1:ix2] = ahet_ct_moms[:,2:n_area_bins]

	# treatment going backwards
	ahet_tr_moms = zeros(n_resp_a,n_area_bins)
	# %     ahet_tr_moms(sample_tr  , :) = (prob_respond * hm_tr   - hd_tr  ) * n_resp_all / n_resp_a_tr;
	# %     ahet_tr_moms(sample_ct23  , :) = ahet_tr_moms(sample_ct23  , :) + ...
	# %                 ((1-prob_respond ) * hm_ct23) * n_resp_all / n_resp_a_ct23;
	ahet_tr_moms[sample_tr  , :] =
		(prob_respond_alpha .* hm_tr_wcha .+
		 (1 .- prob_respond_alpha) .* hm_tr_free .- hd_tr  ) .* n_resp_all ./ n_resp_a_tr
	# ahet_tr_moms = cat(1,zeros(n_resp_na,n_area_bins), ahet_tr_moms);

	# ix1 = size(moms,2) + 1 - (n_area_bins-1)*2
	# ix2 = size(moms,2)     - (n_area_bins-1)
	# @assert all(moms[:,ix1:ix2] .== 0)
	# moms[1:n_resp_a,ix1:ix2] = ahet_tr_moms[:, 2:n_area_bins]

## Draw

	binwidth = 1/n_area_bins
	ahetgrid = vec(Array((binwidth/2):binwidth:(1-binwidth/2)))

	# barplot(ahetgrid, hd_ct23)

	ahetgrid_datamodel = vcat(ahetgrid, ahetgrid)
	ct23_datamodel = vcat(
			vec(mean(hd_ct23, dims=1)),
			vec(mean(      prob_respond_alpha  .* hm_ct23_wcha .+
					 (1 .- prob_respond_alpha) .* hm_ct23_free .- hd_ct23, dims=1)))

	hm_tr_final  = prob_respond_alpha .* hm_tr_wcha .+ (1 .- prob_respond_alpha) .* hm_tr_free
	tr_datamodel = vcat(vec(mean(hd_tr, dims=1)), vec(mean(hm_tr_final, dims=1)))

	### Data and model -- control
	barplot(ahetgrid_datamodel, ct23_datamodel,
        dodge = [1, 1, 1, 2, 2, 2],
        color = [1, 1, 1, 2, 2, 2],
        axis = (xticks = (ahetgrid, ["[0, 0.33)", "[0.33, 0.66)", "[0.66, 1]"]),
                title = "Control Group"),
        ) |> display

	### Data and model -- treatment
	barplot(ahetgrid_datamodel, tr_datamodel,
        dodge = [1, 1, 1, 2, 2, 2],
        color = [1, 1, 1, 2, 2, 2],
        axis = (xticks = (ahetgrid, ["[0, 0.33)", "[0.33, 0.66)", "[0.66, 1]"]),
                title = "Treatment Group"),
        ) |> display

	### Only data
	# barplot(ahetgrid, vec(mean(hd_ct23, dims=1)),
    #     dodge = [1, 1, 1],
    #     color = [1, 1, 1],
    #     axis = (xticks = (ahetgrid, ["[0, 0.33)", "[0.33, 0.66)", "[0.66, 1]"]),
    #             title = "Control Group"),
    #     )


end

function draw_vot_dynamic_model(moms;
		mylinealpha=1.0, mylinewidth=1.5, mylabel="",
		mylinecolor=:black, hold=false, offset=0, mytitle="")

	if hold
		myplotfn = Plots.plot!
	else
		myplotfn = Plots.plot
	end

   myplotfn([0,1,2,3,4], moms[6:10],
		   label=mylabel,
		   linewidth=mylinewidth, linestyle=:solid,
		   color=mylinecolor, linealpha=mylinealpha)

   ## pretty
   xticklabels = ["before exp", "week 1", "week 2", "week 3", "week 4", "1->0\ntransition\nprobability"]

	# mom11_format = @sprintf "%4.2f" moms[11]
	offsetwidth = 0.25
	Plots.bar!([5 + offset * offsetwidth], [moms[11]],
				bar_width=0.25, color=mylinecolor,
				fillalpha=0.5, label="")

   myplotobj = Plots.plot!([0,1,2,3,4], moms[1:5],
	   linestyle=:dash,
	   linealpha=mylinealpha,
	   linewidth=mylinewidth,
	   ylims=(0,1.0),
	   yticks=0:0.1:1.0,
	   xticks = (0:5, xticklabels),
	   ylabel="Fraction detour route",
	   title=mytitle,
	   thickness_scaling = 1.1,
	   legend=(0.55, 0.95),
	   label="",
	   color=mylinecolor)

	myplotobj |> display
end

function draw_moms_vot_dynamic(;
		theta, param_factors, DATAMOMS, SIMMODEL,
		display_plot=false, save_path="")

	data_a = DATAMOMS["data_a"]
	probs_route1_free = 1 .- SIMMODEL["probs_route0_free"]
	probs_route1_wcha = 1 .- SIMMODEL["probs_route0_wcha"]

	theta_dict = get_params_fn(theta, param_factors)
	prob_respond_alpha = theta_dict["prob_respond_alpha"]
	probs_route1_attenuated = probs_route1_wcha .* prob_respond_alpha .+
							  probs_route1_free .* (1 .- prob_respond_alpha)


## Weekly average samples
	mysample_early_0 = @. data_a.a_early * (data_a.n_obs_0 > 0)
	mysample_early_1 = @. data_a.a_early * (data_a.n_obs_1 > 0)
	mysample_early_2 = @. data_a.a_early * (data_a.n_obs_2 > 0)
	mysample_early_3 = @. data_a.a_early * (data_a.n_obs_3 > 0)
	mysample_early_4 = @. data_a.a_early * (data_a.n_obs_4 > 0)

	mysample_late_0 = @. data_a.a_late * (data_a.n_obs_0 > 0)
	mysample_late_1 = @. data_a.a_late * (data_a.n_obs_1 > 0)
	mysample_late_2 = @. data_a.a_late * (data_a.n_obs_2 > 0)
	mysample_late_3 = @. data_a.a_late * (data_a.n_obs_3 > 0)
	mysample_late_4 = @. data_a.a_late * (data_a.n_obs_4 > 0)

## Data Moments
	data_early_w0_mean = mean(data_a.mean_long_0[mysample_early_0])
	data_early_w1_mean = mean(data_a.mean_long_1[mysample_early_1])
	data_early_w2_mean = mean(data_a.mean_long_2[mysample_early_2])
	data_early_w3_mean = mean(data_a.mean_long_3[mysample_early_3])
	data_early_w4_mean = mean(data_a.mean_long_4[mysample_early_4])

	data_late_w0_mean = mean(data_a.mean_long_0[mysample_late_0])
	data_late_w1_mean = mean(data_a.mean_long_1[mysample_late_1])
	data_late_w2_mean = mean(data_a.mean_long_2[mysample_late_2])
	data_late_w3_mean = mean(data_a.mean_long_3[mysample_late_3])
	data_late_w4_mean = mean(data_a.mean_long_4[mysample_late_4])

## Model Moments
	model_early_w0_mean = mean(probs_route1_attenuated[mysample_early_0, 1])
	model_early_w1_mean = mean(probs_route1_attenuated[mysample_early_1, 2])
	model_early_w2_mean = mean(probs_route1_attenuated[mysample_early_2, 3])
	model_early_w3_mean = mean(probs_route1_attenuated[mysample_early_3, 4])
	model_early_w4_mean = mean(probs_route1_attenuated[mysample_early_4, 5])

	model_late_w0_mean = mean(probs_route1_attenuated[mysample_late_0, 1])
	model_late_w1_mean = mean(probs_route1_attenuated[mysample_late_1, 2])
	model_late_w2_mean = mean(probs_route1_attenuated[mysample_late_2, 3])
	model_late_w3_mean = mean(probs_route1_attenuated[mysample_late_3, 4])
	model_late_w4_mean = mean(probs_route1_attenuated[mysample_late_4, 5])

## graph 1 with weekly take-up rates

	 Plots.plot([0,1,2,3,4], [data_early_w0_mean,
	 					data_early_w1_mean,
						data_early_w2_mean,
						data_early_w3_mean,
						data_early_w4_mean],
			label="Data: charges week 1", linewidth=3, color=:blue)

	Plots.plot!([0,1,2,3,4], [data_late_w0_mean,
	 					data_late_w1_mean,
						data_late_w2_mean,
						data_late_w3_mean,
						data_late_w4_mean],
			label="Data: charges week 4", linewidth=3, linestyle=:dash, color=:blue)

	Plots.plot!([0,1,2,3,4], [model_early_w0_mean,
						model_early_w1_mean,
						model_early_w2_mean,
						model_early_w3_mean,
						model_early_w4_mean],
			label="Model: charges week 1",
			linewidth=4.5,
			color=:red, linealpha=0.5)


	## pretty
	xticklabels = ["before exp", "week 1", "week 2", "week 3", "week 4"]

	myplotobj = Plots.plot!([0,1,2,3,4],
					   [model_late_w0_mean,
						model_late_w1_mean,
						model_late_w2_mean,
						model_late_w3_mean,
						model_late_w4_mean],
		linestyle=:dash,
		linealpha=0.5,
		linewidth=4.5,
		ylims=(0,0.5),
		yticks=0:0.1:0.5,
		xticks = (0:4, xticklabels),
		ylabel="Fraction detour route",
		thickness_scaling = 1.1,
		legend=(0.55, 0.95),
		label="Model: charges week 4", color=:red)

	if save_path != ""
		savefig(save_path)
	end

	if display_plot
		myplotobj |> display
		return
	end

# legend=(2.5, 0.475),
end


function draw_moms_vot_dynamic_het(;
		theta, param_factors, DATAMOMS, SIMMODEL,
		display_plot=true,
		save_path="")

## unpack
	data_a = DATAMOMS["data_a"]
	probs_route1_free = 1 .- SIMMODEL["probs_route0_free"]
	probs_route1_wcha = 1 .- SIMMODEL["probs_route0_wcha"]

	theta_dict = get_params_fn(theta, param_factors)
	prob_respond_alpha = theta_dict["prob_respond_alpha"]
	probs_route1_attenuated =
			probs_route1_wcha .* prob_respond_alpha .+
			probs_route1_free .* (1 .- prob_respond_alpha)

## Weekly average samples
	mysample_early_0 = @. data_a.a_early * (data_a.n_obs_0 > 0)
	mysample_early_1 = @. data_a.a_early * (data_a.n_obs_1 > 0)
	mysample_early_2 = @. data_a.a_early * (data_a.n_obs_2 > 0)
	mysample_early_3 = @. data_a.a_early * (data_a.n_obs_3 > 0)
	mysample_early_4 = @. data_a.a_early * (data_a.n_obs_4 > 0)

	mysample_late_0 = @. data_a.a_late * (data_a.n_obs_0 > 0)
	mysample_late_1 = @. data_a.a_late * (data_a.n_obs_1 > 0)
	mysample_late_2 = @. data_a.a_late * (data_a.n_obs_2 > 0)
	mysample_late_3 = @. data_a.a_late * (data_a.n_obs_3 > 0)
	mysample_late_4 = @. data_a.a_late * (data_a.n_obs_4 > 0)

## graph 2 with take-up as fn of detourtime in "treat" group
	# 	 before (week 0), charges (week 1) and persistence (weeks 2 and 3 -- not in GMM)

	### DATA ###
		### Treat as a function of the detour -- ONLY TREATMENT
		# mysample = (atdf.study_cycle .== 0) .& (atdf.a_early .== 1) .& (atdf.adetourtime .< 15)
		# vgrid_control = atdf[mysample, :adetourtime]
		# rout2_control = atdf[mysample, :is_long_route]
		vgrid_control = data_a[mysample_early_0 .& (data_a.adetourtime .< 15), :adetourtime]
		rout1_control = data_a[mysample_early_0 .& (data_a.adetourtime .< 15), :mean_long_0]

		mysample = mysample_early_1 .& (data_a.adetourtime .< 15)
		vgrid_treat = data_a[mysample, :adetourtime]
		rout1_treat = data_a[mysample, :mean_long_1]

		mysample2 = mysample_early_2 .& (data_a.adetourtime .< 15)
		mysample3 = mysample_early_3 .& (data_a.adetourtime .< 15)
		vgrid_persi = vcat(	data_a[mysample2, :adetourtime],
							data_a[mysample3, :adetourtime])
		rout1_persi = vcat( data_a[mysample2, :mean_long_2],
							data_a[mysample3, :mean_long_3])

		vgrid_control = convert.(Float64, vgrid_control)
		vgrid_treat = convert.(Float64, vgrid_treat)
		vgrid_persi = convert.(Float64, vgrid_persi)

		rout1_control = convert.(Float64, rout1_control)
		# rout1_control = convert.(Float64, rout1_control)
		# rout1_control = convert.(Float64, rout1_control)


		# plot!(dtgrid, pi2_after, color="green", linestyle=:dashdot, label="model after")
		data_control = loess(vgrid_control, rout1_control; span=0.8, degree=1)
		us = range(extrema(vgrid_control)...; step = 0.1)
		vs = Loess.predict(data_control, us)
		Plots.plot(us, vs,
				color=:blue,
				linestyle=:solid,
				linewidth=3,
				label="Data before experiment")

		data_treat = loess(vgrid_treat, rout1_treat; span=0.8, degree=1)
		us = range(extrema(vgrid_treat)...; step = 0.1)
		vs = Loess.predict(data_treat, us)
		Plots.plot!(us, vs,
				color=:blue,
				linestyle=:dash,
				linewidth=3,
				label="Data week 1 (charges)")

		data_persi = loess(vgrid_persi, rout1_persi; span=0.8, degree=1)
		us = range(extrema(vgrid_persi)...; step = 0.1)
		vs = Loess.predict(data_persi, us)
		Plots.plot!(us, vs,
				color=:blue,
				linestyle=:dot,
				linewidth=3,
				label="Data week 2 (persistence)")

	### MODEL ###
		detours_vec = data_a.adetourtime
		# detours_vec = convert.(Float64, detours_vec)

		mysample0 = (detours_vec .< 15) .& mysample_early_0
		mysample1 = (detours_vec .< 15) .& mysample_early_1
		mysample2 = (detours_vec .< 15) .& mysample_early_2
		mysample3 = (detours_vec .< 15) .& mysample_early_3
		mysample4 = (detours_vec .< 15) .& mysample_early_4

		### Treat and control as a function of the detour
		vgrid_control = detours_vec[mysample0]
		rout1_control = probs_route1_attenuated[mysample0, 1]

		vgrid_treat = detours_vec[mysample1]
		rout1_treat = probs_route1_attenuated[mysample1, 2]

		vgrid_persi = vcat(detours_vec[mysample2], detours_vec[mysample3])
		rout1_persi = vcat(
					probs_route1_attenuated[mysample2, 3],
					probs_route1_attenuated[mysample3, 4]
						)

		vgrid_control = convert.(Float64, vgrid_control)
		vgrid_treat = convert.(Float64, vgrid_treat)
		vgrid_persi = convert.(Float64, vgrid_persi)

		# plot!(dtgrid, pi2_after, color="green", linestyle=:dashdot, label="model after")
		data_control = loess(vgrid_control, rout1_control; span=0.8, degree=1)
		us = range(extrema(vgrid_control)...; step = 0.1)
		vs = Loess.predict(data_control, us)
		# scatter!(vgrid_control, rout2_control)
		Plots.plot!(us, vs,
				color=:red, linestyle=:solid,
				linewidth=4.5,
				linealpha=0.5,
				label="Model before experiment")

		data_treat = loess(vgrid_treat, rout1_treat; span=0.8, degree=1)
		us = range(extrema(vgrid_treat)...; step = 0.1)
		vs = Loess.predict(data_treat, us)
		Plots.plot!(us, vs,
				color=:red,
				linestyle=:dash,
				linewidth=4.5,
				linealpha=0.5,
				label="Model week 1 (charges)")

		data_persi = loess(vgrid_persi, rout1_persi; span=0.8, degree=1)
		us = range(extrema(vgrid_persi)...; step = 0.1)
		vs = Loess.predict(data_persi, us)
		myplotobj = Plots.plot!(us, vs,
				color=:red,
				linestyle=:dot,
				linewidth=4.5,
				linealpha=0.5,
				label="Model week 2 (persistence)",
				xlabel="Detour route (r=1) extra detour (minutes)",
				ylabel="Fraction detour route",
				thickness_scaling = 1.1,
				yticks=0:0.1:0.7,
				ylims=(0,0.7)
				)

	if save_path != ""
		savefig(save_path)
	end

	if display_plot
		myplotobj |> display
	end

end


# linestyle=:dash,
# linewidth=3,
# ylims=(0,0.5),
# yticks=0:0.1:0.5,
# xticks = (0:4, xticklabels),
# ylabel="Fraction detour route",
# thickness_scaling = 1.1,
# legend=(0.55, 0.95),
# label="Model: charges week 4", color=:red)
