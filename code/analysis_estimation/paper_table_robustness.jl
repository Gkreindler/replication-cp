include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

# ] add CSV, Glob, DataFrames, FreqTables, Statistics, StatsBase, Printf
using CSV
using Glob
using DataFrames
using FreqTables
using Statistics
using StatsBase
using Printf

## All functions to read and organize GMM results
    include("functions-table.jl")

##

rootpath = string(rootpath_base, "model/")


all_results = Vector{Any}(undef,0)
all_specs = Vector{Any}(undef,0)

ci_levels = [2.5, 97.5]


## Static Route Choice
    rootpath_results      =  string(rootpath,"estimation_results/results_dt_static/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_dt_static/")]
    param_names = ["alpha", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)

    println("Spec diagnostics for DT static")
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    spec_dict_df["delta"] = (90.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)

## Asymmetric switching costs
    rootpath_results      =  string(rootpath,"estimation_results/results_asym_switchcost/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_asym_switchcost/")]
    param_names = ["alpha", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    
    println("Spec diagnostics for asymetric switching cost")
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    spec_dict_df["delta"] = (90.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)


## With FEs
    rootpath_results      =  string(rootpath,"estimation_results/results_withFE/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_withFE/")]
    param_names = ["alpha", "gamma", "fe_early", "fe_w1", "fe_w2", "fe_w3", "fe_w4", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    
    println("Spec diagnostics for with FE")
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    spec_dict_df["delta"] = (90.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)


## Assuming half attention (p=0.5)
    rootpath_results      =  string(rootpath,"estimation_results/results_halfattention/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_halfattention/")]
    param_names = ["alpha", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    
    
    println("Spec diagnostics for half attention")
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    # spec_dict_df["delta"] = (90.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)

## All parameters proportional to wages
    rootpath_results      =  string(rootpath,"estimation_results/results_propwage/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_propwage/")]
    param_names = ["alpha", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    
    
    println("Spec diagnostics for prop to wage")
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    spec_dict_df["delta"] = (90.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)

## Single h_A*
    rootpath_results      =  string(rootpath,"estimation_results/results_singlehastar/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_singlehastar/")]
    param_names = ["alpha", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
            
    println("Spec diagnostics for single ideal arrival time")
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    spec_dict_df["delta"] = (90.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)


## Build the table

row_names = Dict(
    "beta_early"=> "\$\\beta_E\$: Schedule cost early (INR/hour)",
    "beta_late" => "\$\\beta_L\$: Schedule cost late (INR/hour)",
    "alpha"     => "\$\\alpha\$: Value of travel time (INR/hour)",
    "gamma"     => "\$\\gamma\$: Route switching cost (INR)",
    "mu"        => "\$\\sigma^R\$: Logit route nest",
    "sigma_dt"  => "\$\\sigma^{DT}\$: Logit departure time",
)

row_order = ["beta_early", "beta_late", "alpha", "gamma", "sigma_dt", "mu"]
row_factor = [6.0, 6.0, 6.0, 1.0, 1.0, 1.0, 0.01]
row_d      = [0,   0,   0,   1,   1,   1,   2]

tex_contents =
"\\resizebox*{\\textwidth}{!}{" *
"\\begin{tabular}{lcccccc}\n" *
"   \\toprule\n" *
"   & (1) & (2) & (3) & (4) & (5) & (6) \\\\\n" *
"   & \\thead{Static \\\\ Route Choice} & \\thead{ Asymmetric \\\\ switching cost} & Time FE & \\thead{ Half \\\\ attention} & \\thead{ Parameters \\\\ prop. to \\\\ wage} & \\thead{ Single ideal \\\\ arrival time} \\\\\n" *
"   \\addlinespace\\addlinespace\n"


for irow=1:length(row_order)

    this_row1 = row_names[row_order[irow]]
    this_row2 = "   "
    this_row3 = "   "

    for jcol=1:length(all_specs)
        if haskey(all_specs[jcol], row_order[irow])

            temp = all_specs[jcol][row_order[irow]]

            myfactor = row_factor[irow]

            this_row1 *= " & " * pn1(temp[1] * myfactor, d=row_d[irow])
            if isnothing(temp[2])
                this_row2 *= " & "
                this_row3 *= " & "
            else
                this_row2 *= " & (" * pn1(temp[4] * myfactor, d=row_d[irow]) * ")"  # standard deviation boot
                this_row3 *= " & [" * pn1(temp[2] * myfactor, d=row_d[irow]) * ", " * pn1(temp[3] * myfactor, d=row_d[irow]) * "]"
            end
        else
            this_row1 *= " & "
            this_row2 *= " & "
            this_row3 *= " & "
        end
    end

    this_row1 *= " \\\\ \n"
    this_row2 *= " \\\\ \n"
    this_row3 *= " \\\\ \n"

    # global tex_contents = tex_contents * this_row1 * this_row2 * this_row3 * "   \\addlinespace\n"
    global tex_contents = tex_contents * this_row1 * this_row3 * "   \\addlinespace\n"
end

tex_contents *=
"\\addlinespace \n" *
"\\multicolumn{7}{l}{\\textit{Model Components:}}  \\\\ \n \\addlinespace \n " *
" Route choice model           & Static & Dynamic & Dynamic  & Dynamic & Dynamic  & Dynamic \\\\ \n \\addlinespace \n " *
" Fixed discount factor (\$\\delta\$)  & - & 0.90 & 0.90 & 0.90 & 0.90 & 0.90 \\\\ \n \\addlinespace \n " *
" Asymmetric switch cost (\$\\gamma_{01}=2\\gamma_{10}\$)" *
                                    "  & - & Yes &  -  &  -  &  -  &  -  \\\\ \n \\addlinespace \n " *
" Route Choice Time FE                 & - & -   & Yes &  -  &  -  &  -  \\\\ \n \\addlinespace \n " *
"\\addlinespace \n" *
"\\multicolumn{7}{l}{\\textit{Moments:}}  \\\\ \n \\addlinespace \n " *
" Departure Time (49)           & Yes & Yes & Yes & Yes & Yes & Yes  \\\\ \n \\addlinespace \n " *
" Dynamic route choice (10)     & -   & Yes & Yes & Yes & Yes & Yes\\\\ \n \\addlinespace \n " *
" Static route choice (2)       & Yes & -   & -   & -   & -   & -  \\\\ \n \\addlinespace \n " *
"   \\bottomrule \n " *
"\\end{tabular} \n } \n"


table_path = rootpath_base * "paper/tables/smtable9/"
isdir(table_path) || mkdir(table_path)

open(table_path * "table_gmm_appendix.tex","w") do io
   println(io,tex_contents)
end
