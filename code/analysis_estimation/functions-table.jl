

using CSV
# using Glob
using DataFrames
using FreqTables
using Statistics
using StatsBase
using Printf

function pn1(mynumber; d=1)
    (d==0) && return @sprintf("%2.0f", mynumber)
    (d==1) && return @sprintf("%2.1f", mynumber)
    (d==2) && return @sprintf("%2.2f", mynumber)
    (d==3) && return @sprintf("%2.3f", mynumber)
end


function f0(mynumber)
	return string(Int64(floor(mynumber)))
end
function f1(mynumber)
	return string(Int64(floor(mynumber * 10)) / 10  )
end
function f2(mynumber)
	return string(Int64(floor(mynumber * 100)) / 100)
end

function compute_reliability(gmm_df; pm_range=0.1)
    gmm_optimum = gmm_df[gmm_df.is_best_vec .== 1, :]
    @assert size(gmm_optimum, 1) == 1

    for mycol in names(gmm_df)
        if startswith(mycol, "param_")
            theta_star = gmm_optimum[1, mycol]

            gmm_df[!, string("inrange_", mycol)] =
                (abs.(gmm_df[:, mycol] .- theta_star) ./ abs.(theta_star)) .< pm_range

            gmm_df[!, string("inrange_",mycol)] .=
                mean(gmm_df[:, string("inrange_",mycol)])
        end
    end
end

function read_results_file(; path, boot_idx=-1, pm_range=0.1)
    # read stage 1
    if isfile(string(path, "gmm-stage1-all.csv"))
        filepath = string(path, "gmm-stage1-all.csv")
        stage1 = CSV.read(filepath, DataFrame)
        stage1[!, "boot_idx"] .= boot_idx
        stage1[!, "stage"] .= 1
        stage1[!, "file_exists"] .= 1
        compute_reliability(stage1, pm_range=pm_range)
    else
        stage1 = DataFrame("file_exists" => [0])
        stage1[!, "boot_idx"] .= boot_idx
        stage1[!, "is_best_vec"] .= 1
        stage1[!, "opt_converged"] .= 0
        stage1[!, "stage"] .= 1
    end

    # read stage 1
    if isfile(string(path, "gmm-stage2-all.csv"))
        filepath = string(path, "gmm-stage2-all.csv")
        stage2 = CSV.read(filepath, DataFrame)
        stage2[!, "boot_idx"] .= boot_idx
        stage2[!, "stage"] .= 2
        stage2[!, "file_exists"] .= 1
        compute_reliability(stage2, pm_range=pm_range)
    else
        stage2 = DataFrame("file_exists" => [0])
        stage2[!, "boot_idx"] .= boot_idx
        stage2[!, "is_best_vec"] .= 1
        stage2[!, "opt_converged"] .= 0
        stage2[!, "stage"] .= 2
    end

return vcat(stage1, stage2, cols = :union)
end


function get_all_results(;
            rootpath_results="",
            rootpath_boot_results=[],
            boot_runs=120, pm_range=0.1, debug=false)

    ### Read main results
        main_df = read_results_file(path=rootpath_results)
        main_df[!, :dir_exists] .= 1

## Loop over bootstrap directories
boot_dfs = []
for mybootfolder=rootpath_boot_results
    ### parse folder
    dir = readdir(mybootfolder)
    files = [elm for elm in dir if !isdir(string(mybootfolder,elm))]
    dirs  = [elm for elm in dir if isdir(string(mybootfolder,elm))]

    ### Read bootstrap results
    for idx_bootdir = 1:boot_runs

        mydir = string(mybootfolder, "boot_run_", idx_bootdir, "/")

        if ~isdir(mydir)
            debug && println(mydir, " does not exist")
            # boot_stage1_df[idx_bootdir,:dir_exists] = 0
            # boot_stage2_df[idx_bootdir,:dir_exists] = 0
            tempdf = DataFrame(
                "boot_idx" => idx_bootdir,
                "dir_exists" => 0,
                "stage" => 1,
                "file_exists" => 0,
                "opt_converged" => 0,
                "is_best_vec" => 1
                )
            push!(boot_dfs, copy(tempdf))
            tempdf[!, "stage"] .= 2
            push!(boot_dfs, copy(tempdf))
        else

            # println(mydir[10:end])
            debug && println("reading folder ", mydir)

            # folder_idx = parse(Int64, mydir[10:end])
            # string(rootpath_results, mydir) |> display

            tempdf = read_results_file(
                path=mydir, #string(mybootfolder, mydir, "/"),
                boot_idx=idx_bootdir,
                pm_range=pm_range)
            tempdf[!, :dir_exists] .= 1
            push!(boot_dfs, tempdf)

        end
    end

end
    # display(boot_dfs)

    full_df = vcat(main_df, boot_dfs..., cols=:union)
    debug && freqtable(full_df.file_exists) |> display
    # myregex = r"boot_run_"
    # match(myregex, dirs[1])

    ## parse (1) main (2nd stage), and (2) boot (2nd stage)
    main_df1 = full_df[(full_df.is_best_vec .== 1) .& (full_df.stage .== 1) .& (full_df.boot_idx .== -1), :] |> copy
    @assert size(main_df1, 1) == 1
    boot_df1 = full_df[(full_df.is_best_vec .== 1) .& (full_df.stage .== 1) .& (full_df.boot_idx .> 0), :] |> copy
    @assert size(boot_df1, 1) == boot_runs

    main_df2 = full_df[(full_df.is_best_vec .== 1) .& (full_df.stage .== 2) .& (full_df.boot_idx .== -1), :] |> copy
    @assert size(main_df2, 1) == 1
    boot_df2 = full_df[(full_df.is_best_vec .== 1) .& (full_df.stage .== 2) .& (full_df.boot_idx .> 0), :] |> copy
    @assert size(boot_df2, 1) == boot_runs

    # display(main_df1)
    # display(main_df2)
    # display(names(boot_df1))
    # display(boot_df1[:, :inrange_param_1])
    # display(boot_df2[:, :inrange_param_1])
    # fdsfd

    ## df with results
    diag_df = DataFrame(
            "statistic" => [],
            "main stage 1" => [],
            "main stage 2" => [],
            "boot stage 1" => [],
            "boot stage 2" => [])

    push!(diag_df, ["frac_dir_exists",
                mean(main_df1.dir_exists |> skipmissing),
                mean(main_df2.dir_exists |> skipmissing),
                mean(boot_df1.dir_exists |> skipmissing),
                mean(boot_df2.dir_exists |> skipmissing)])

    # diag_df[!, "spec"] = ["main stage 1", "main stage 2", "boot stage 1", "boot stage 2"]
    push!(diag_df, ["frac_file_exists",
                mean(main_df1.file_exists |> skipmissing),
                mean(main_df2.file_exists |> skipmissing),
                mean(boot_df1.file_exists |> skipmissing),
                mean(boot_df2.file_exists |> skipmissing)])

    push!(diag_df, ["frac_converged",
                mean(main_df1.opt_converged |> skipmissing),
                mean(main_df2.opt_converged |> skipmissing),
                mean(boot_df1.opt_converged |> skipmissing),
                mean(boot_df2.opt_converged |> skipmissing)])

    for mycol in names(full_df)
        if startswith(mycol, "inrange_")

            # debug && display(mycol)
            # debug && display(boot_df2[:, mycol])

            push!(diag_df, [string("min_", mycol),
                minimum(main_df1[:, mycol] |> skipmissing),
                minimum(main_df2[:, mycol] |> skipmissing),
                minimum(boot_df1[:, mycol] |> skipmissing),
                minimum(boot_df2[:, mycol] |> skipmissing)])

            push!(diag_df, [string("med_", mycol),
                median(main_df1[:, mycol] |> skipmissing),
                median(main_df2[:, mycol] |> skipmissing),
                median(boot_df1[:, mycol] |> skipmissing),
                median(boot_df2[:, mycol] |> skipmissing)])
        end
    end

    return main_df1, boot_df1, main_df2, boot_df2, diag_df, full_df
end


function main_table_column(main_df, boot_df, param_names;
        ci_levels=(5, 95), keep_only_converged=true)


    # check that we have exactly the right number of parameters
    param_cols = [mycol for mycol in names(main_df) if startswith(mycol, "param_")]
    @assert length(param_cols) == length(param_names)

    results = Dict()
    for i=1:length(param_names)
        param_name = param_names[i]

        if keep_only_converged
            boot_vec = boot_df[boot_df.opt_converged .== 1, string("param_", i)]
        else
            boot_vec = boot_df[:, string("param_", i)]
        end
        boot_vec = skipmissing(boot_vec) |> collect
        boot_vec = sort(boot_vec)

        # nboot = length(boot_vec)
        # bootidx_low  = floor(nboot * (100.0 - ci_level)/100.0) |> Int64
        # bootidx_high =  ceil(nboot * ci_level/100.0) |> Int64
        boot_cilow, boot_cihigh = percentile(boot_vec, ci_levels)

        results[param_name] = (
            main_df[1, string("param_", i)],
            # boot_vec[bootidx_low],
            # boot_vec[bootidx_high],
            boot_cilow, boot_cihigh,
            std(boot_vec)
            )
    end
    return results
end
