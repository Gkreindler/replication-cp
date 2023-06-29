
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

    # Rootpath
    rootpath_sims = rootpath_base * "model/simulation_results/main/"
    rootpath_boot = rootpath_base * "model/simulation_results/main_boot/"

    table_path = rootpath_base * "paper/tables/table4/"
    isdir(table_path) || mkdir(table_path)
    table_path = table_path * "table_policy_main.tex"

## BOOT main simulation - read BSON objects

    boot_tt         = zeros(120,4)
    boot_tt_rel     = zeros(120,4)
    boot_logsum     = zeros(120,4)
    boot_logsum_rel = zeros(120,4)

    for i=1:120
        # Paths
        path = rootpath_boot * "boot_" * string(i) * "/so_eqm.bson"
        main_so_eqm = BSON.load(path)

        # Average travel time
        nash_tt = sum(main_so_eqm[:nash_ttimes] .* main_so_eqm[:nash_choice_probs], dims=2) |> mean
        so_tt   = sum(main_so_eqm[:so_ttimes] .* main_so_eqm[:so_choice_probs], dims=2) |> mean
        level_tt  = nash_tt - so_tt
        # improve_perc = (nash_tt - so_tt) / nash_tt
        perc_tt = (nash_tt - so_tt) / nash_tt * 100.0

        # Average travel time relative to no congestion benchmark
        ttime_nocon = minimum(main_so_eqm[:nash_ttimes], dims=2)
        nash_tt_rel = sum((main_so_eqm[:nash_ttimes] .- ttime_nocon) .* main_so_eqm[:nash_choice_probs], dims=2) |> mean
        so_tt_rel   = sum((main_so_eqm[:so_ttimes]   .- ttime_nocon) .* main_so_eqm[:so_choice_probs], dims=2) |> mean
        level_tt_rel = nash_tt_rel - so_tt_rel
        perc_tt_rel  = (nash_tt_rel - so_tt_rel) / nash_tt_rel * 100.0

        # Welfare
        nash_logsum = main_so_eqm[:nash_logsum]
        so_logsum = main_so_eqm[:so_logsum]
        level_logsum = so_logsum - nash_logsum
        perc_logsum  = (nash_logsum - so_logsum) / nash_logsum * 100.0

        # Welfare relative to no congestion benchmark
        nash_logsum_rel = mean(main_so_eqm[:nash_logsum]) - mean(main_so_eqm[:logsum_nocon])
        so_logsum_rel = mean(main_so_eqm[:so_logsum]) - mean(main_so_eqm[:logsum_nocon])
        level_logsum_rel = so_logsum_rel - nash_logsum_rel
        perc_logsum_rel  = (nash_logsum_rel - so_logsum_rel) / nash_logsum_rel * 100.0

        boot_tt[i,:]         .= [nash_tt        , so_tt        , level_tt        , perc_tt        ]
        boot_tt_rel[i,:]     .= [nash_tt_rel    , so_tt_rel    , level_tt_rel    , perc_tt_rel    ]
        boot_logsum[i,:]     .= [nash_logsum    , so_logsum    , level_logsum    , perc_logsum    ]
        boot_logsum_rel[i,:] .= [nash_logsum_rel, so_logsum_rel, level_logsum_rel, perc_logsum_rel]
    end

    function get_ci(mymat)
        ncol = size(mymat)[2]

        cis = zeros(2, ncol)

        for icol=1:ncol
            cis[1, icol] = percentile(mymat[:,icol] |> vec, 2.5)
            cis[2, icol] = percentile(mymat[:,icol] |> vec, 97.5)
        end
        return cis
    end

    boot_tt_cis         = get_ci(boot_tt)
    boot_tt_rel_cis     = get_ci(boot_tt_rel)
    boot_logsum_cis     = get_ci(boot_logsum)
    boot_logsum_rel_cis = get_ci(boot_logsum_rel)

    path = rootpath_sims * "so_eqm.bson"
    main_so_eqm = BSON.load(path)

## Load Panel A (benchmark specification)
    # columns: Nash, SO, Improvement, Improvement %

    # Average travel time
    nash_tt = sum(main_so_eqm[:nash_ttimes] .* main_so_eqm[:nash_choice_probs], dims=2) |> mean
    so_tt   = sum(main_so_eqm[:so_ttimes] .* main_so_eqm[:so_choice_probs], dims=2) |> mean
    level_tt  = nash_tt - so_tt
    # improve_perc = (nash_tt - so_tt) / nash_tt
    perc_tt = (nash_tt - so_tt) / nash_tt * 100.0

    # Average travel time relative to no congestion benchmark
    ttime_nocon = minimum(main_so_eqm[:nash_ttimes], dims=2)
    nash_tt_rel = sum((main_so_eqm[:nash_ttimes] .- ttime_nocon) .* main_so_eqm[:nash_choice_probs], dims=2) |> mean
    so_tt_rel   = sum((main_so_eqm[:so_ttimes]   .- ttime_nocon) .* main_so_eqm[:so_choice_probs], dims=2) |> mean
    level_tt_rel = nash_tt_rel - so_tt_rel
    perc_tt_rel  = (nash_tt_rel - so_tt_rel) / nash_tt_rel * 100.0

    # Welfare
    nash_logsum = main_so_eqm[:nash_logsum]
    so_logsum = main_so_eqm[:so_logsum]
    level_logsum = so_logsum - nash_logsum
    perc_logsum  = (nash_logsum - so_logsum) / nash_logsum * 100.0

    # Welfare relative to no congestion benchmark
    nash_logsum_rel = mean(main_so_eqm[:nash_logsum]) - mean(main_so_eqm[:logsum_nocon])
    so_logsum_rel = mean(main_so_eqm[:so_logsum]) - mean(main_so_eqm[:logsum_nocon])
    level_logsum_rel = so_logsum_rel - nash_logsum_rel
    perc_logsum_rel  = (nash_logsum_rel - so_logsum_rel) / nash_logsum_rel * 100.0

## Load Heterogeneity (Panel B)
    rootpath_sims = rootpath * "model/simulation_results/main_het/"
    rootpath_boot = rootpath * "model/simulation_results/propwage_boot/" # main_het_boot/"

    ### boot
    het_boot_logsum = zeros(120,4)
    for i=1:120
        # Paths
        path = rootpath_boot * "boot_" * string(i) * "/so_eqm.bson"
        het_main_so_eqm = BSON.load(path)

        # Welfare
        het_nash_logsum = het_main_so_eqm[:nash_logsum]
        het_so_logsum = het_main_so_eqm[:so_logsum]
        het_level_logsum = het_so_logsum - het_nash_logsum
        het_perc_logsum  = (het_nash_logsum - het_so_logsum) / het_nash_logsum * 100.0

        het_boot_logsum[i,:] .=
            [het_nash_logsum, het_so_logsum, het_level_logsum, het_perc_logsum]
    end

    het_boot_logsum_cis     = get_ci(het_boot_logsum)

    # Read main het sim
    path = rootpath_sims * "so_eqm.bson"
    het_main_so_eqm = BSON.load(path)

    # Welfare
    het_nash_logsum = het_main_so_eqm[:nash_logsum]
    het_so_logsum = het_main_so_eqm[:so_logsum]
    het_level_logsum = het_so_logsum - het_nash_logsum
    het_perc_logsum  = (het_nash_logsum - het_so_logsum) / het_nash_logsum * 100.0


## Load Preferences (Panel C)
    rootpath_sims = rootpath_base * "model/simulation_results/varyparams/"

    myfolders = ["sim_vary_linear_alpha_0_25/",
                "sim_vary_linear_alpha_4_00/",
                "sim_vary_linear_beta_e_0_25/",
                "sim_vary_linear_beta_e_4_00/"]

    v_logsum = zeros(4,4)
    for i=1:4
        # Paths
        path = rootpath_sims * myfolders[i] * "so_eqm.bson"
        v_so_eqm = BSON.load(path)

        # Welfare
        v_nash_logsum  = v_so_eqm[:nash_logsum]
        v_so_logsum    = v_so_eqm[:so_logsum]
        v_level_logsum = v_so_logsum - v_nash_logsum
        v_perc_logsum  = (v_nash_logsum - v_so_logsum) / v_nash_logsum * 100.0

        v_logsum[i,:] .= [v_nash_logsum, v_so_logsum, v_level_logsum, v_perc_logsum]
    end

## Load Power Road Technology (panel D)
    rootpath_sims = rootpath_base * "model/simulation_results/power/"

    myfolders = ["power_5/", "power_15/"]

    pow_logsum = zeros(2,4)
    for i=1:2
        # Paths
        path = rootpath_sims * myfolders[i] * "so_eqm.bson"
        pow_so_eqm = BSON.load(path)

        # Welfare
        pow_nash_logsum  = pow_so_eqm[:nash_logsum]
        pow_so_logsum    = pow_so_eqm[:so_logsum]
        pow_level_logsum = pow_so_logsum - pow_nash_logsum
        pow_perc_logsum  = (pow_nash_logsum - pow_so_logsum) / pow_nash_logsum * 100.0

        pow_logsum[i,:] .= [pow_nash_logsum, pow_so_logsum, pow_level_logsum, pow_perc_logsum]
    end

## Load 2 routes simulation (Panel E)
    rootpath_sims = rootpath * "model/simulation_results/main_2routes/"

    r2_logsum = zeros(2,4)

    # Paths
    path = rootpath_sims * "v15_low_mu/so_eqm.bson"
    r2_so_eqm = BSON.load(path)

    # Welfare
    r2_nash_logsum  = r2_so_eqm[:nash_logsum]
    r2_so_logsum    = r2_so_eqm[:so_logsum]
    r2_level_logsum = r2_so_logsum - r2_nash_logsum
    r2_perc_logsum  = (r2_nash_logsum - r2_so_logsum) / r2_nash_logsum * 100.0
    r2_logsum[1, :] = [r2_nash_logsum, r2_so_logsum, r2_level_logsum, r2_perc_logsum]

    # Paths
    path = rootpath_sims * "v30/so_eqm.bson"
    r2_so_eqm = BSON.load(path)

    # Welfare
    r2_nash_logsum  = r2_so_eqm[:nash_logsum]
    r2_so_logsum    = r2_so_eqm[:so_logsum]
    r2_level_logsum = r2_so_logsum - r2_nash_logsum
    r2_perc_logsum  = (r2_nash_logsum - r2_so_logsum) / r2_nash_logsum * 100.0
    r2_logsum[2, :] = [r2_nash_logsum, r2_so_logsum, r2_level_logsum, r2_perc_logsum]

## Main main table
    ### HEADING
    tex_contents =
        "\\begin{tabular}{lcccc}\n" *
        "\\toprule\n" *
        "\\addlinespace\n" *
        "& (1) & (2) & (3) & (4) \\\\ \n" *
        "& Nash & \\thead{Social \\\\ Optimum} & " *
        " Improvement &" *
        " \\thead{Improvement \\\\ (\\% of Nash)} \\\\ \n" *
        "\\addlinespace\\addlinespace \n"


    ### PANEL A
        tex_contents *= "\\multicolumn{5}{l}{\\emph{Panel A. Benchmark Model:" *
                        " Travel Time and Commuter Welfare}} \\\\ \n" *
                        "\\addlinespace \n"
        # Travel time
        tex_contents *= "Travel Time (minutes) & " *
                        f1(nash_tt) * " & " *
                        f1(so_tt) * " & " *
                        f1(level_tt) * " & " *
                        f1(perc_tt) * "\\%  \\\\ \n"

        tex_contents *= " & (" * f1(boot_tt_cis[1,1]) * "," * f1(boot_tt_cis[2,1]) * ") " *
                        " & (" * f1(boot_tt_cis[1,2]) * "," * f1(boot_tt_cis[2,2]) * ") " *
                        " & (" * f1(boot_tt_cis[1,3]) * "," * f1(boot_tt_cis[2,3]) * ") " *
                        " & (" * f1(boot_tt_cis[1,4]) * "," * f1(boot_tt_cis[2,4]) * ") \\\\ \n" * "\\addlinespace \n"


        # Travel time relative
        tex_contents *= "Travel Time above Free Flow & " *
                        f1(nash_tt_rel) * " & " *
                        f1(so_tt_rel) * " & " *
                        f1(level_tt_rel) * " & " *
                        f1(perc_tt_rel) * "\\%  \\\\ \n"

        tex_contents *= " & (" * f1(boot_tt_rel_cis[1,1]) * "," * f1(boot_tt_rel_cis[2,1]) * ") " *
                        " & (" * f1(boot_tt_rel_cis[1,2]) * "," * f1(boot_tt_rel_cis[2,2]) * ") " *
                        " & (" * f1(boot_tt_rel_cis[1,3]) * "," * f1(boot_tt_rel_cis[2,3]) * ") " *
                        " & (" * f1(boot_tt_rel_cis[1,4]) * "," * f1(boot_tt_rel_cis[2,4]) * ") \\\\ \n" * "\\addlinespace \n"

        # Welfare (INR)
        tex_contents *= "Welfare (INR) & " *
                        f0(nash_logsum) * " & " *
                        f0(so_logsum) * " & " *
                        f1(level_logsum) * " & " *
                        f1(perc_logsum) * "\\%  \\\\ \n"

        tex_contents *= " & (" * f0(boot_logsum_cis[1,1]) * "," * f0(boot_logsum_cis[2,1]) * ") " *
                        " & (" * f0(boot_logsum_cis[1,2]) * "," * f0(boot_logsum_cis[2,2]) * ") " *
                        " & (" * f1(boot_logsum_cis[1,3]) * "," * f1(boot_logsum_cis[2,3]) * ") " *
                        " & (" * f1(boot_logsum_cis[1,4]) * "," * f1(boot_logsum_cis[2,4]) * ") \\\\ \n" * "\\addlinespace \n"

        # Welfare relative (INR)
        tex_contents *= "Welfare above Free Flow (INR) & " *
                        f0(nash_logsum_rel) * " & " *
                        f0(so_logsum_rel) * " & " *
                        f1(level_logsum_rel) * " & " *
                        f1(perc_logsum_rel) * "\\%  \\\\ \n"

        tex_contents *= " & (" * f0(boot_logsum_rel_cis[1,1]) * "," * f0(boot_logsum_rel_cis[2,1]) * ") " *
                        " & (" * f0(boot_logsum_rel_cis[1,2]) * "," * f0(boot_logsum_rel_cis[2,2]) * ") " *
                        " & (" * f1(boot_logsum_rel_cis[1,3]) * "," * f1(boot_logsum_rel_cis[2,3]) * ") " *
                        " & (" * f1(boot_logsum_rel_cis[1,4]) * "," * f1(boot_logsum_rel_cis[2,4]) * ") \\\\ \n" * "\\addlinespace \n"

        tex_contents *= "\\addlinespace \n"

## Panel B - Heterogeneity

    tex_contents *= "\\multicolumn{5}{l}{\\emph{Panel B. " *
                " Preference Heterogeneity (Commuter Welfare, INR) }} \\\\ \n " *
                "\\addlinespace \n"

    # Welfare (INR)
    tex_contents *= "Preferences \$\\propto\$ wage & " *
                    f0(het_nash_logsum) * " & " *
                    f0(het_so_logsum) * " & " *
                    f1(het_level_logsum) * " & " *
                    f1(het_perc_logsum) * "\\%  \\\\ \n"

    tex_contents *= " & (" * f0(het_boot_logsum_cis[1,1]) * "," * f0(het_boot_logsum_cis[2,1]) * ") " *
                    " & (" * f0(het_boot_logsum_cis[1,2]) * "," * f0(het_boot_logsum_cis[2,2]) * ") " *
                    " & (" * f1(het_boot_logsum_cis[1,3]) * "," * f1(het_boot_logsum_cis[2,3]) * ") " *
                    " & (" * f1(het_boot_logsum_cis[1,4]) * "," * f1(het_boot_logsum_cis[2,4]) * ") \\\\ \n" * "\\addlinespace \n"

    tex_contents *= "\\addlinespace \n"

## Panel C - varying preferences and technology
    tex_contents *= "\\multicolumn{5}{l}{\\emph{Panel C. " *
                    "Varying Preferences (Commuter Welfare, INR) }} \\\\ \n " *
                    "\\addlinespace \n"

    # VOTT
    tex_contents *= "VOTT\$\\approx 100\\%\$ wage (\$4\\times\$ smaller \$\\alpha\$) & " *
                    f0(v_logsum[1,1]) * " & " *
                    f0(v_logsum[1,2]) * " & " *
                    f1(v_logsum[1,3]) * " & " *
                    f1(v_logsum[1,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"
    tex_contents *= "VOTT\$\\approx 1,600\\%\$ wage (\$4\\times\$ larger \$\\alpha\$) & " *
                    f0(v_logsum[2,1]) * " & " *
                    f0(v_logsum[2,2]) * " & " *
                    f1(v_logsum[2,3]) * " & " *
                    f1(v_logsum[2,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"
    # Schedule cost early
    tex_contents *= "High schedule cost (\$4\\times\$ larger \$\\beta_E\$) & " *
                    f0(v_logsum[4,1]) * " & " *
                    f0(v_logsum[4,2]) * " & " *
                    f1(v_logsum[4,3]) * " & " *
                    f1(v_logsum[4,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"
    tex_contents *= "Low schedule cost (\$4\\times\$ smaller \$\\beta_E\$) & " *
                    f0(v_logsum[3,1]) * " & " *
                    f0(v_logsum[3,2]) * " & " *
                    f1(v_logsum[3,3]) * " & " *
                    f1(v_logsum[3,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"

    tex_contents *= "\\addlinespace \n"

## Panel D - varying road technology
    tex_contents *= "\\multicolumn{5}{l}{\\emph{Panel D. " *
                    "Varying Road Technology (Commuter Welfare, INR) }} \\\\ \n " *
                    "\\addlinespace \n"

    # POWER
    tex_contents *= "Concave Road Technology (\$\\nu=0.5\$) & " *
                    f0(pow_logsum[1,1]) * " & " *
                    f0(pow_logsum[1,2]) * " & " *
                    f1(pow_logsum[1,3]) * " & " *
                    f1(pow_logsum[1,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"
    tex_contents *= "Convex Road Technology (\$\\nu=1.5\$) & " *
                    f0(pow_logsum[2,1]) * " & " *
                    f0(pow_logsum[2,2]) * " & " *
                    f1(pow_logsum[2,3]) * " & " *
                    f1(pow_logsum[2,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"

    tex_contents *= "\\addlinespace \n"

## Panel E - two routes equilibrium
    tex_contents *= "\\multicolumn{5}{l}{\\emph{Panel E." *
                    "Equilibrium with Two Routes (Commuter Welfare, INR) }} \\\\ \n " *
                    "\\addlinespace \n"

    # POWER
    tex_contents *= "One route 15\\% steeper & " *
                    f0(r2_logsum[1,1]) * " & " *
                    f0(r2_logsum[1,2]) * " & " *
                    f1(r2_logsum[1,3]) * " & " *
                    f1(r2_logsum[1,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"

    tex_contents *= "One route 30\\% steeper & " *
                    f0(r2_logsum[2,1]) * " & " *
                    f0(r2_logsum[2,2]) * " & " *
                    f1(r2_logsum[2,3]) * " & " *
                    f1(r2_logsum[2,4]) * "\\%  \\\\ \n" * "\\addlinespace \n"

    tex_contents *= "\\addlinespace \n"


## Wrap up
    tex_contents *=
    "   \\bottomrule \n " *
    "\\end{tabular} \n"

    open(table_path,"w") do io
       println(io, tex_contents)
    end
