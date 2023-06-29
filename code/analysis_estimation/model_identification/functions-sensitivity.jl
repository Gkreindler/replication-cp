using Distributed
using DelimitedFiles

function compute_sensitivity(;
                DATA,
				ESTIMATION_PARAMS,
                mymomfunction::Function,
                rootpath_results,
                rootpath_boot_results,
                optimal_W_path,
				parameter_names_list="",
				moms_names_list="",
				myfactor=1e15)

## read all bootstrap runs
	println("... reading estimates parameters (including bootstrap runs)")
	main_df1, boot_df1, main_df2, boot_df2, diag_df, full_df = get_all_results(
			rootpath_results=rootpath_results,
			rootpath_boot_results=rootpath_boot_results,
			boot_runs=120, pm_range=0.1, debug=false)

	# get matrix or vector of estimates ("^" anchors to start of the word)
	theta_optimum = select(main_df2, r"^param_") |> Array |> vec
	theta_boot    = select(boot_df2, r"^param_") |> Array

## Read optimal weighting matrix
	println("... reading optimal weighting matrix")
	# optimal_W_path = string(rootpath,"gmm37/results/gmm-stage1-optW.csv")
	optW = readdlm(optimal_W_path, ',', Float64)

## prep function
	println("... run moments at the optimum")
	mymomfunction_main = theta -> mymomfunction(
		theta, DATA, ESTIMATION_PARAMS
	)

	# adding minus so it's data minus model like in Andrews et al
	function mymomfunction_main_avg(theta)
		return mean( mymomfunction_main(theta), dims=1) |> vec
	end

	moms = mymomfunction_main(theta_optimum)
	temp = mymomfunction_main_avg(theta_optimum)
	n_moms = length(temp)

## S matrix = E(f*transpose(f)) where f are moms:
# https://ocw.mit.edu/courses/sloan-school-of-management/15-450-analytics-of-finance-fall-2010/lecture-notes/MIT15_450F10_lec08.pdf

	n_obs = size(moms,1)

	S = transpose(moms) * moms / n_obs
	S = Hermitian(S)

## bootstrapped standard error of moments
	println("... bootstrap moment values")

	n_boot = size(theta_boot, 1)
	boot_moms = zeros(n_boot, n_moms)
	for i=1:n_boot
		print(".")
		one_theta_boot = theta_boot[i, :]

		boot_folder = rootpath_boot_results[1] * "boot_run_" * string(i)* "/"

		## Load bootstrap samples
			bootstrap_samples = Dict{String, Vector{Int64}}()

			## Load boot samples to file
			path = string(boot_folder,"boot_sample_a.csv")
			bootstrap_samples["a"] = readdlm(path, ',', Float64) |> vec
			# display(bootstrap_samples["a"])

			path = string(boot_folder,"boot_sample_na.csv")
			bootstrap_samples["na"] = readdlm(path, ',', Float64) |> vec

		## Load boot data
		BIGMAT_BOOT, CHARGEMAT_BOOT, COMPMAT_BOOT, DATAMOMS_BOOT =
			sample_bigmats(
			DATA["BIGMAT"], DATA["CHARGEMAT"], DATA["COMPMAT"], DATA["DATAMOMS"],
			bootstrap_samples["a"], bootstrap_samples["na"])

		DATA_BOOT=Dict(
			"BIGMAT" => BIGMAT_BOOT,
			"CHARGEMAT" => CHARGEMAT_BOOT,
			"COMPMAT" => COMPMAT_BOOT,
			"DATAMOMS" => DATAMOMS_BOOT
		)

		boot_moms[i, :] = mean(
			mymomfunction(
				one_theta_boot,
				DATA_BOOT,
				ESTIMATION_PARAMS), dims=1) |> vec
	end
	println(".")

## write moments to file
	path = string(rootpath_results, "boot_moms.csv")
	writedlm(path, boot_moms, ",")

	boot_moms_std = std(boot_moms, dims=1) |> vec

## Run Jacobian
	println("... computing jacobian")

	# higher factor = larger changes in parameters
	myjac = jacobian(central_fdm(5, 1, factor=myfactor), mymomfunction_main_avg, theta_optimum)
	myjac_df = DataFrame(myjac[1], :auto)
	myjac_df[!, "moment"] = moms_names_list

	select!(myjac_df, vcat(["moment"], "x" .* string.(1:length(theta_optimum))) )

	if parameter_names_list != ""
		rename!(myjac_df, vcat(["moment"], parameter_names_list))
	end

	path = string(rootpath_results, "jac.csv")
	CSV.write(path, myjac_df)

## Calculate and save the Gentzkow-Shapiro sensitivity matrix
	println("... computing AGS sensitivity measure")

	G = myjac[1]
	W = Hermitian(optW)

# Andrews, Gentzkow, and Shapiro sensitivity parameter - NOT SCALED
	Lambda = - Hermitian(transpose(G) * W * G) \ transpose(G) * W

	path = string(rootpath_results, "ags_lambda.csv")
	writedlm(path, Lambda, ",")

# Gentzkow Shapiro sensitivity parameter - SCALED
# Lambdatilde = diag(boot_ests_std)\Lambda*diag(boot_moms_std);
    Lambdatilde = Lambda * diagm(boot_moms_std)

	path = string(rootpath_results, "ags_lambda_scaled.csv")
	writedlm(path, Lambdatilde, ",")

## Analytic standard errors
	vcov = inv(transpose(G) * inv(S) * G) / n_obs

	ses = sqrt.(diag(vcov))

	display(theta_optimum)
	display(ses)

	return theta_optimum, theta_boot, boot_moms, G, W, Lambda, Lambdatilde, ses
end
