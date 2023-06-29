
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

using DataFrames
using FreqTables
using Statistics
using StatsBase
using Printf

using InlineStrings


## All functions to read and organize GMM results
    include("functions-table.jl")


## Read all extensive margin sims
	rootpath_sims = rootpath_base * "model/simulation_results/main_ext/"

## 
	table_path = rootpath_base * "paper/tables/table5/"
	isdir(table_path) || mkdir(table_path)

function f0(mynumber)
	return string(Int64(floor(mynumber)))
end
function f1(mynumber)
	return string(Int64(floor(mynumber * 10)) / 10  )
end
function f2(mynumber)
	return string(Int64(floor(mynumber * 100)) / 100)
end


tex_contents =
"\\begin{tabular}{cccccc} \n" *
"\\toprule \n" *
"\\addlinespace \n" *
"\\thead{Trip value \\\\ \$ \\omega \$} & " *
	"\\thead{Logit trip \\\\ parameter \\\\ \$ \\eta \$} & " *
	"\\thead{ \\% Traffic Reduction \\\\ from flat 100 INR fee }  & " *
	"\\thead{Implied elasticity \\\\ at flat 100 INR fee} & " *
	"\\thead{ Nash \\\\ Welfare (INR)} & " *
	"\\thead{Improvement \\\\ at Social Optimum \\\\ (\\% of Nash)} \\\\  \n" *
"\\addlinespace\\addlinespace \n"

all_so_eqm = []
all_elastic = []
all_avg_charges = []
for (idx_eta, eta) in enumerate([20.0, 100.0])
	# for (idx_delta, delta) in enumerate(600.0:100.0:1100.0)
	for (idx_delta, delta) in enumerate(600.0:200.0:1000.0)

		theta_add_dict = Dict{String, Float64}(
			"delta_ext" => -delta, # penalty for no trip
			"eta" => eta
		)

		# Paths
		rootpath_output = rootpath_sims *
						 "eta_" * string(Int64(eta)) * "_" *
						"delta_" * string(Int64(delta)) * "/"

		# println(rootpath_output)

		path = rootpath_output * "so_eqm.bson"
		so_eqm = BSON.load(path)
		push!(all_so_eqm, so_eqm)

		# display(so_eqm)
		# fdsfd

		# average charges
		avg_charges = sum(so_eqm[:so_charges] .* so_eqm[:so_choice_probs], dims=2) |> vec
		push!(all_avg_charges, mean(avg_charges))

		path = rootpath_output * "elasticity_results.bson"
		elastics = BSON.load(path)
		push!(all_elastic, elastics)

		## get level reduction and elasticity from 100 Rs
		level_resp100 = elastics[:level_ch_100]
		elasticity100 = elastics[:elastic_100]

		# no congestion benchmark
		nocon_welfare = so_eqm[:nocon_logsum]
		nocon_trip_prob = so_eqm[:nocon_trip_prob]

		# nash welfare
		nash_welfare = so_eqm[:nash_logsum]
		nash_trip_prob = so_eqm[:nash_trip_prob] |> mean

		# nash welfare
		so_welfare = so_eqm[:so_logsum]
		so_trip_prob = so_eqm[:so_trip_prob]

		deadweight_loss = (so_welfare - nash_welfare) / abs(nash_welfare)

		## print key points
		println(
			"delta=",f0(-delta), " eta=", f1(eta),
			"  level=", f2(level_resp100),
			"  elast=", f2(elasticity100),
			"  nash=", f0(nash_welfare),
			"  dwl=", f1(deadweight_loss * 100),
			"  nash π=", f2(nash_trip_prob))

		this_row1 = f0(delta) * " INR & " *
					f0(eta) * " & " *
					f2(level_resp100) * " & " *
					f2(elasticity100) * " & " *
					f1(nash_welfare) * " & " *
					f1(deadweight_loss * 100) * "\\%"

		global tex_contents = tex_contents * this_row1 * " \\\\ \n\\addlinespace\n"
	end
end

tex_contents *=
"   \\bottomrule \n " *
"\\end{tabular}  \n"

open(table_path * "table_sims_ext.tex","w") do io
   println(io,tex_contents)
end

## Average charges
describe(convert.(Float64, all_avg_charges))
