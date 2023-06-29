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

table_path = rootpath_base * "paper/tables/smtable10/"
isdir(table_path) || mkdir(table_path)


all_results = Vector{Any}(undef,0)
all_specs = Vector{Any}(undef,0)

ci_levels = [2.5, 97.5]

### delta's
deltas     = [ 0.0, 50.0, 80.0, 90.0, 95.0, 99.0]
deltas_str = ["00", "50", "80", "90", "95", "99"]

deltas     = [ 0.0, 50.0, 90.0, 99.0]
deltas_str = ["00", "50", "90", "99"]


for i=1:length(deltas)
    rootpath_results      =  string(rootpath,"estimation_results/results_delta_", deltas_str[i], "/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_delta_", deltas_str[i], "/")]
    param_names = ["alpha", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    spec_dict_df["delta"] = (deltas[i], nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)
end

## Estimate Delta
    rootpath_results      =  string(rootpath,"estimation_results/results_with_delta/")
    rootpath_boot_results = [string(rootpath,"estimation_results/results_boot_with_delta/")]
    param_names = ["alpha", "delta", "gamma", "fe_early", "beta_early", "beta_late", "sigma_dt", "mu"]

    main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df =
        get_all_results(
            rootpath_results=rootpath_results,
            rootpath_boot_results=rootpath_boot_results,
            pm_range=0.1, boot_runs=120, debug=false)
    spec_diagnostics_df |> display

    spec_dict_df = main_table_column(main2_df, boot2_df, param_names; ci_levels=ci_levels)
    spec_dict_df |> display
    # spec_dict_df["delta"] = (90.0, nothing, nothing, nothing)

    push!(all_results, [main1_df, boot1_df, main2_df, boot2_df, spec_diagnostics_df, full_df])
    push!(all_specs, spec_dict_df)




## Build the table

row_names = Dict(
    "beta_early"=> "\$\\beta_E\$: Schedule cost early (INR/hour)",
    "beta_late" => "\$\\beta_L\$: Schedule cost late (INR/hour)",
    "alpha"     => "\$\\alpha\$: Value of travel time (INR/hour)",
    "gamma"     => "\$\\gamma\$: Route switching cost (INR)",
    "mu"        => "\$\\sigma^R\$: Logit route (upper nest)",
    "sigma_dt"  => "\$\\sigma^{DT}\$: Logit departure time",
    "delta"     => "\$\\delta\$: discount factor"
)

row_order = ["beta_early", "beta_late", "alpha", "gamma", "sigma_dt", "mu", "delta"]
row_factor = [6.0, 6.0, 6.0, 1.0, 1.0, 1.0, 0.01]
row_d      = [0,   0,   0,   1,   1,   1,   2]

tex_contents =
"\\resizebox*{\\textwidth}{!}{" *
"\\begin{tabular}{lccccc}\n" *
"   \\toprule\n" *
"   & (1) & (2) & (3) & (4) & (5) \\\\\n \\addlinespace \n" *
"   & \\multicolumn{4}{c}{Varying Discount Factor \$\\delta\$} & Estimate \$\\delta\$  \\\\\n" *
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
"\\multicolumn{6}{l}{\\textit{Model:}}  \\\\ \n \\addlinespace \n " *
" Dynamic route choice model           & Dynamic & Dynamic  & Dynamic  & Dynamic  & Dynamic \\\\ \n \\addlinespace \n " *
" Fixed discount factor (\$\\delta\$)  & 0.0 & 0.50 & 0.90 & 0.99 &  -  \\\\ \n \\addlinespace \n " *
"\\addlinespace  \n" *
"\\multicolumn{6}{l}{\\textit{Moments:}}  \\\\ \n \\addlinespace \n " *
" Dynamic route choice (10)     & Yes & Yes & Yes & Yes & Yes \\\\ \n \\addlinespace \n " *
" Route choice transition (1)   & -   &  -  &  -  &  -  & Yes \\\\ \n \\addlinespace \n " *
"   \\bottomrule \n " *
"\\end{tabular} \n } \n"


open(table_path * "table_gmm_appendix_deltas.tex","w") do io
   println(io,tex_contents)
end
