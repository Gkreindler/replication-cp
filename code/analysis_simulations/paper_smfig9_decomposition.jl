
include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using Plots

using CSV
using DelimitedFiles
using Glob
using BSON

using SpecialFunctions

using DataFrames
using StatFiles

using GLM
using FixedEffectModels

using FreqTables
using Statistics
using StatsBase
using Printf

using InlineStrings

using Loess
using Random

## All functions to read and organize GMM results
    include("functions-table.jl")
	include("sims-functions-model.jl")
	include("sims-functions-loadprep.jl")


function get_ci(mymat)
    ncol = size(mymat)[2]

    cis = zeros(2, ncol)

    for icol=1:ncol
        cis[1, icol] = percentile(mymat[:,icol] |> vec, 2.5)
        cis[2, icol] = percentile(mymat[:,icol] |> vec, 97.5)
    end
    return cis
end

function loess_bootstrap(xs, ys;
            nboot=200, step_x=0.1,
            fresh_graph=true,
			mycolor=:black, mylinewidth=2, mylinestyle=:solid,
			mylabel="y var")

    nx = length(xs)

    # min and max of the x range
    idx_xmin = argmin(xs)
    idx_xmax = argmax(xs)

    # main model
    model = loess(xs, ys, span=0.5)
    us = range(extrema(xs)...; step = step_x)
    vs_main = Loess.predict(model, us)

    n_us = length(us)

    # Bootstrap
    vs_boot = zeros(nboot, n_us) # will hold all bootstrap predictions here
    for iboot=1:nboot
        print(".")

        # sample
        idx_boot = sample( (1:nx) |> Array |> vec, nx-2, replace=true)

        # always include x range endpoints (otherwise can't predict)
        # this point is probably sketchy
        push!(idx_boot, idx_xmin)
        push!(idx_boot, idx_xmax)

        # run model on bootstrap sample
        model = loess(xs[idx_boot], ys[idx_boot], span=0.5)

        vs_boot[iboot, :] .= Loess.predict(model, us)
    end

    # get CI percentiles for each value on the x axis
    vs_ci_low  = zeros(n_us)
    vs_ci_high = zeros(n_us)
    for icol=1:n_us
        vs_ci_low[icol]  = percentile(vs_boot[:, icol], 2.5)
        vs_ci_high[icol] = percentile(vs_boot[:, icol], 97.5)
    end

    # plot

	# main line
    if fresh_graph
        plot(us, vs_main, color=mycolor, label=mylabel,
			linewidth=mylinewidth, linestyle=mylinestyle)
    else
        plot!(us, vs_main, color=mycolor, label=mylabel,
			linewidth=mylinewidth, linestyle=mylinestyle)
    end

	# CI
    plot!(us, vs_ci_low , color=mycolor, alpha=0.5, linestyle=:dash, label="95% CI")
    myplot = plot!(us, vs_ci_high, color=mycolor, alpha=0.5,
			linestyle=:dash, label="")

	return myplot
end

### PATHS
	sm9b_path = rootpath_base * "paper/figures/smfigure9/"
	isdir(sm9b_path) || mkdir(sm9b_path)

    # simulation results
    rootpath_sims = rootpath_base * "model/simulation_results/main/"
    rootpath_boot = rootpath_base * "model/simulation_results/main_boot/"

# plot settings
    Plots.scalefontsizes()
    # Plots.scalefontsizes(1.1)

## separate VOTT and schedule cost changes from Nash to SO
#=
1. load
2. get average departure time under Nash at the person level
3. run make_choice_prbs() to get schedule costs
4. run loess -> save
=#

path = rootpath_sims * "so_eqm.bson"

### 1. Load and definitions
	main_so_eqm = BSON.load(path)

	hdgrid = (5:(1/12):14) |> transpose |> Matrix
	agents = main_so_eqm[:agents]

	# bigmats that we reuse
	t_exact_mat 	  = agents.ha .- hdgrid
	log_t_exact_mat = log.(max.(0, t_exact_mat))

	BIGMAT = Dict{String, Array}()
	BIGMAT["t_exact_mat"] = t_exact_mat
	BIGMAT["log_t_exact_mat"] = log_t_exact_mat
	BIGMAT["t_exact_mat_flip"] = Matrix(transpose(BIGMAT["t_exact_mat"]))
	BIGMAT["log_t_exact_mat_flip"] = Matrix(transpose(BIGMAT["log_t_exact_mat"]))
	BIGMAT["ttimes"] = copy(main_so_eqm[:so_ttimes])

### 2. avearge departure time
	avg_deptime_nash = sum(main_so_eqm[:nash_choice_probs] .* hdgrid, dims=2) |> vec
	avg_deptime_so   = sum(main_so_eqm[:so_choice_probs]   .* hdgrid, dims=2) |> vec

	change_avg_deptime = (avg_deptime_so .- avg_deptime_nash) .* 60.0

	# describe(change_avg_deptime)
	percentile(change_avg_deptime, [10, 90]) |> display
	# percentile(change_avg_deptime, [ 5, 95]) |> display
	# percentile(change_avg_deptime, [ 1, 99]) |> display

	# change in average departure time (SO - NASH) by original average departure time
	loess_bootstrap(avg_deptime_nash, change_avg_deptime)

	### arrival time before/after ideal arrival time in Nash
		arrival_relative = sum(main_so_eqm[:nash_choice_probs] .* (hdgrid .+ main_so_eqm[:nash_ttimes] ./ 60 .- agents.ha ), dims=2) |> vec
		arrival_relative .*= 60.0
		loess_bootstrap(avg_deptime_nash, arrival_relative)


### 3. get expected schedule cost
	nash_eu = make_choice_prbs(
				agents_mean_km=agents.mean_km,
				agents_a=agents.a,
				agents_be=agents.be,
				agents_bl=agents.bl,
				agents_sig=agents.sig,
				ttimes=main_so_eqm[:nash_ttimes],
				dt_charges=main_so_eqm[:so_charges] .* 0.0,
				t_exact_mat=BIGMAT["t_exact_mat_flip"],
				log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
				compute_logsum=1,
				return_utility=true)

	### define approximate average experienced schedule costs
	# add back time
	nash_eu_sch = nash_eu .+ agents.a .* main_so_eqm[:nash_ttimes] ./ 60.0
	nash_eu_sch_avg = sum(nash_eu_sch .* main_so_eqm[:nash_choice_probs], dims=2) |> vec

	nash_eu_time_avg = (-1) .* sum(agents.a .* main_so_eqm[:nash_ttimes] ./ 60.0 .*
							main_so_eqm[:nash_choice_probs], dims=2) |> vec

	so_eu = make_choice_prbs(
				agents_mean_km=agents.mean_km,
				agents_a=agents.a,
				agents_be=agents.be,
				agents_bl=agents.bl,
				agents_sig=agents.sig,
				ttimes=main_so_eqm[:so_ttimes],
				dt_charges=main_so_eqm[:so_charges] .* 0.0,
				t_exact_mat=BIGMAT["t_exact_mat_flip"],
				log_t_exact_mat=BIGMAT["log_t_exact_mat_flip"],
				compute_logsum=1,
				return_utility=true)

	### define approximate average experienced schedule costs
	so_eu_sch = so_eu .+ agents.a .* main_so_eqm[:so_ttimes] ./ 60.0
	so_eu_sch_avg = sum(so_eu_sch .* main_so_eqm[:so_choice_probs], dims=2) |> vec

	so_eu_time_avg = (-1) .* sum(agents.a .* main_so_eqm[:so_ttimes] ./ 60.0 .*
							main_so_eqm[:so_choice_probs], dims=2) |> vec


	### Take difference
	change_eu_sch = so_eu_sch_avg .- nash_eu_sch_avg

	change_eu_time = so_eu_time_avg .- nash_eu_time_avg

	## How about if charges are redistributed propto KM

		# utility difference (without charges)
		change_logsum = main_so_eqm[:so_logsum_vec] .- main_so_eqm[:nash_logsum_vec]

		# paid charges
	    charges_pay = sum(main_so_eqm[:so_charges] .* main_so_eqm[:so_choice_probs], dims=2) |> vec
	    total_charges = sum(charges_pay)

		# rebate proportional to KM
	    mean_km = agents.mean_km
	    charges_get_km = -(mean_km .* total_charges) ./ sum(mean_km)
	    @assert abs(sum(charges_get_km) + total_charges) < 1e-6

		ΔEU_km = change_logsum .+ charges_pay .+ charges_get_km
		ΔEU_km = ΔEU_km |> vec

		# flat rebate
		nagents = size(agents, 1)
		charges_get_flat = zeros(nagents) .- total_charges / nagents
	    @assert abs(sum(charges_get_flat) + total_charges) < 1e-6

	    ΔEU_flat = change_logsum .+ charges_pay .+ charges_get_flat
		ΔEU_flat = ΔEU_flat |> vec

## Distributional consequences with flat/per-km rebate
	# read original income
	path = rootpath_base * "data/coded_cto/recruitment coded.dta"
	df_coded = DataFrame(load(path))
	df_coded = select(df_coded, [:uidp_original, :log_income_with_pred])
	rename!(df_coded, [:uidp, :log_income_with_pred])

	agents_with_income = leftjoin(agents, df_coded, on=:uidp, indicator="_m", validate=(false, true))
	@assert all(agents_with_income._m .== "both")
	@assert sum(ismissing.(agents_with_income.log_income_with_pred)) == 0

	agents_with_income.income_with_pred = exp.(agents_with_income.log_income_with_pred) ./ 20.0 ./ 8.0

	describe(agents_with_income.income_with_pred)

	low_income = agents_with_income.income_with_pred .< median(agents_with_income.income_with_pred)

	mydf = DataFrame("uidp" => agents_with_income.uidp, "ΔEU_km" => ΔEU_km, "ΔEU_flat" => ΔEU_flat, "low_income" => low_income)
	low_income_model = reg(mydf, @formula(ΔEU_km ~ low_income), Vcov.cluster(:uidp))

	describe(ΔEU_km[low_income])
	describe(ΔEU_km[.~low_income])


##
	nboot = 200

### 4. run loess -> save
	minu, maxu = extrema(avg_deptime_nash)
	xs = avg_deptime_nash
	ys = change_eu_time
	loess_bootstrap(xs, ys, fresh_graph=true, mycolor=:black,
			mylabel="Δ VOTT EU", nboot=nboot)

	ys = change_eu_sch
	myplot = loess_bootstrap(xs, ys, nboot=nboot, fresh_graph=false,
			mycolor=:green, mylinestyle=:dash,
			mylabel="Δ Schedule EU")

	# ys = change_logsum_wcharges .≥ 0.0
	# myplot = loess_bootstrap(xs, ys, nboot=20, fresh_graph=false, mycolor=:black,
	# 		mylabel="ΔEU (Flat Rebate)")

	# histogram!(twinx(myplot), avg_deptime_nash, linealpha=0.5, fillcolor=:white,
	# 			fillalpha=0.0, label="Avg Dep Time (Nash)")
	hline!([0.0], color=:gray, linestyle=:dash, label="")
	plot!(myplot, xlabel="Average commuter departure time in Nash equilibrium")
	plot!(myplot, xticks=6:12)
	plot!(myplot, xlims=(5, 12))  # , ylims=(-50, 50)
	plot!(myplot, ylabel="Change in Expected Utility (INR)")
	plot!(myplot, legend=:topleft,
		legendtitlefontsize=10,
		legendfontsize=10, xtickfontsize=10, ytickfontsize=10,
		xguidefontsize=10, yguidefontsize=10)

	# save to file
	savefig(myplot, sm9b_path * "cost_change_decomposition.pdf")

	myplot |> display

## Distributional
	xs = agents_with_income.income_with_pred
	ys = ΔEU_km #.> 0.0
	myplot = loess_bootstrap(xs, ys, nboot=nboot, fresh_graph=true,
			mycolor=:black,
			mylabel="ΔEU (Rebate ∝ KM)")
	hline!([0.0], color=:gray, linestyle=:dash, label="")



## Including congestion charges
	xs = avg_deptime_nash
	ys = ΔEU_km #.> 0.0
	# scatter(avg_deptime_nash, ys)
	myplot = loess_bootstrap(xs, ys, nboot=nboot, fresh_graph=true,
			mycolor=:black,
			mylabel="ΔEU (Rebate ∝ KM)")
	hline!([0.0], color=:gray, linestyle=:dash, label="") 

	ys = ΔEU_flat #.> 0.0
	# scatter(avg_deptime_nash, ys)
	myplot = loess_bootstrap(xs, ys, nboot=nboot, fresh_graph=false,
			mycolor=:green, mylinestyle=:dash,
			mylabel="ΔEU (Flat Rebate)")

	plot!(myplot, xlabel="Average commuter departure time in Nash equilibrium")
	plot!(myplot, xticks=6:12)
	plot!(myplot, xlims=(5, 12))  # , ylims=(-50, 50)
	plot!(myplot, ylabel="Change in Expected Utility (INR)")
	plot!(myplot, legend=:topright,
		yticks=-50:50:350,
		legendtitlefontsize=10,
		legendfontsize=10, xtickfontsize=10, ytickfontsize=10,
		xguidefontsize=10, yguidefontsize=10)

	# save to file
	savefig(myplot, sm9b_path * "cost_change_decomposition_wcharges.pdf")

	myplot |> display

## Compute stat in footnote on bus users:
	# 2.4 = 2.4 minutes gained on average per trip in social optimum vs Nash
	# 611 = value of time
	# /2  = assume bus users have VOTT/2
	# 8.9 INR = SO vs Nash gains
	2.4 / 60 * 611 / 2
	2.4 / 60 * 611 / 2 / 8.9
