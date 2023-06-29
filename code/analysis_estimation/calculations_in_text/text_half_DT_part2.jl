"""

    Using estimated parameters, analyze the departure time response in minutes to
    a congestion charge that is equal to one hour's wage, and other values.

    1) load estimated parameters
    2) load data
    3) invert ideal arrival time distributions
    4) define individual specific charge profile
	5) compute and summarize response!

"""

include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

using CSV
using DataFrames


## Read baseline results
	rootpath_params = rootpath_base * "model/estimation_results/results_delta_90/"

	main_df = read_results_file(path=rootpath_params)
	theta_optimum = select(main_df[(main_df.stage .== 2) .& (main_df.is_best_vec .== 1),:], r"^param_") |> Array |> vec


## Read new results
    rootpath_params = rootpath_base * "model/estimation_results/results_DD_05x/"

	main_df = read_results_file(path=rootpath_params)
	theta_05x = select(main_df[(main_df.stage .== 2) .& (main_df.is_best_vec .== 1),:], r"^param_") |> Array |> vec

## Compare
	beta_e, beta_l = (theta_optimum[4], theta_optimum[5]) .* 6.0
	beta_e_x05, beta_l_x05 = (theta_05x[4], theta_05x[5]) .* 6.0

	println("With x0.5 DT response we get beta_e, beta_l: ", beta_e_x05, ", ", beta_l_x05)
	println("With x0.5 DT response we get ratios: ", beta_e_x05/beta_e , ", ", beta_l_x05/beta_l)

## Save output

    mytext = ""
    mytext *= "With x0.5 DT response we get beta_e, beta_l: " * string(beta_e_x05) * ", " * string(beta_l_x05) * "\n"
    mytext *= "With x0.5 DT response we get ratios: " * string(beta_e_x05/beta_e) * ", " * string(beta_l_x05/beta_l)

    text_results_path = rootpath_base * "paper/text_results/model_DT_05.txt"
    file = open(text_results_path, "w")
    write(file, mytext)
    close(file)