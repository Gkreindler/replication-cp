include("../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

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

### column 1 -> (full model) dynamic VOTT and DT momements
    rootpath_results      =  rootpath * "estimation_results/results_delta_90/"
    rootpath_boot_results = [rootpath * "estimation_results/results_boot_delta_90/"]
    param_names = ["alpha", "gamma", "fe_early",
            "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)




### column 2 -> static VOTT model with fixed alpha = wage /2
    rootpath_results       = string(rootpath,"estimation_results/results_dt_halfwage/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_dt_halfwage/")]
    param_names = ["beta_early", "beta_late", "sigma_dt"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120
            )
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df["alpha"] = (165.0 / 6.0 / 2.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)



### column 3 -> static VOTT model with fixed alpha = wage
    rootpath_results       = string(rootpath,"estimation_results/results_dt_wage/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_dt_wage/")]
    param_names = ["beta_early", "beta_late", "sigma_dt"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120
            )
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df["alpha"] = (165.0 / 6.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)




### column 4 -> dynamic VOTT model without DT
    rootpath_results =       string(rootpath,"estimation_results/results_nodt/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_nodt/")]
    param_names = ["alpha", "gamma", "fe_early", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120
            )
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(
            main2_df, boot2_df, param_names;
            ci_levels=ci_levels, keep_only_converged=true)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)

## Build the table
function pn1(mynumber; d=1)
    (d==0) && return @sprintf("%2.0f", mynumber)
    (d==1) && return @sprintf("%2.1f", mynumber)
    (d==2) && return @sprintf("%2.2f", mynumber)
    (d==3) && return @sprintf("%2.3f", mynumber)
end

row_names = Dict(
    "beta_early"=> "\$\\beta_E\$: Schedule cost early (INR/hour)",
    "beta_late" => "\$\\beta_L\$: Schedule cost late (INR/hour)",
    "alpha"     => "\$\\alpha\$: Value of travel time (INR/hour)",
    "gamma"     => "\$\\gamma\$: Route switching cost (INR)",
    "delta"     => "\$\\delta\$: Discount factor",
    "mu"        => "\$\\sigma^R\$: Logit route (upper nest)",
    "sigma_dt"  => "\$\\sigma^{DT}\$: Logit departure time",
)

row_order = ["beta_early", "beta_late", "alpha", "gamma", "sigma_dt", "mu"]
row_factor = [6.0, 6.0, 6.0, 1.0, 1.0, 1.0]
row_d      = [0,   0,   0,   1,   1,   1]


tex_contents =
"\\resizebox*{\\textwidth}{!}{" *
"\\begin{tabular}{lcccc}\n" *
"   \\toprule\n" *
"   & (1) & (2) & (3) & (4) \\\\\n" *
" \\textit{Model:}  & Benchmark & \\multicolumn{2}{c}{\\thead{Only Departure Time \\\\ Calibrated VOTT}} & \\thead{Only \\\\ Route Choice}  \\\\\n" *
"   \\addlinespace\\addlinespace\n"


for irow=1:length(row_order)

    this_row1 = row_names[row_order[irow]]
    this_row2 = "   "
    this_row3 = "   "

    for jcol=1:4
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
"\\multicolumn{5}{l}{\\textit{Model Components:}}  \\\\ \n \\addlinespace \n " *
"Route choice model         & Dynamic & -  & - & Dynamic  \\\\ \n \\addlinespace \n " *
"Calibrated VOTT \$ \\alpha \$ &  -   & \$50\\%\$ wage    & \$100\\%\$ wage    &    -     \\\\ \n \\addlinespace \n " *
"Departure time model       & Yes     & Yes    & Yes    &     -     \\\\ \n \\addlinespace\\addlinespace \n " *
"\\multicolumn{5}{l}{\\textit{Moments:}}  \\\\ \n \\addlinespace \n " *
" Departure time (49)        & Yes & Yes & Yes & -   \\\\ \n \\addlinespace \n " *
" Dynamic route choice (10)  & Yes &  -  &  -  & Yes \\\\ \n \\addlinespace \n " *
"   \\bottomrule \n " *
"\\end{tabular} \n } \n"

table_path = rootpath_base * "paper/tables/table2/"
isdir(table_path) || mkdir(table_path)

open(table_path * "table_main_gmm_4col.tex","w") do io
   println(io,tex_contents)
end
