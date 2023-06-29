"""
model definitions for GMM
"""

##
model_1_dt_static = 1
"""
Model 1: DT + static nested logit
"""
function model_dt_static_params(;
		theta_fixed_dict, gmm_options=nothing, main_n_start_pts=1, boot_n_start_pts=1)

		if ~isnothing(gmm_options)
			main_n_start_pts = gmm_options["main_n_start_pts"]
			boot_n_start_pts = gmm_options["boot_n_start_pts"]
		end

        # full parameter list with expansion factors
	# this is the order that will be used in optimization
	# will next drop from optimization any parameters held fixed
        param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}} = [
		("alpha", 	6.0, nothing),   ## x 6
		("beta_e", 	6.0, nothing),  ## x 6
		("beta_l", 	6.0, nothing),  ## x 6
		("sigma_dt", 	1.0, nothing),
		("mu", 		1.0, nothing),
		("prob_respond_alpha", 0.01, nothing) ## /100
	]
	all_param_names = [el[1] for el in param_factors]

	# update the given dictionary (overwrite)
	param_included_vec = ones(Bool, length(all_param_names))
	for mykey in keys(theta_fixed_dict)
		if ~in(mykey, all_param_names)
		    error("Fixed parameter name is not in the current model:", mykey)
		end

		idx = findfirst(all_param_names .== mykey)

		# do not include the parameter in estimation
		param_included_vec[idx] = false

		# set paramter value
		param_factors[idx] = (param_factors[idx][1], param_factors[idx][2], theta_fixed_dict[mykey])
	end

	# initial conditions and bounds
    	theta_matrix = [
            100.0 0.0 500.0;  # alpha in INR per 10 MINUTES
            50.0 1.0 Inf;   # beta_early per 10 MINUTES
            50.0 1.0 Inf;   # beta_late  per 10 MINUTES
            25.0 1.0 Inf;   # sigma_dt = logit parameter for departure time
            25.0 1.0 Inf;   # mu = logit parameter for route choice
            50.0 0.0 100.0  # probability to respond x 100
    		]
	theta_matrix = theta_matrix[param_included_vec, :]

    	theta_initial = theta_matrix[:, 1] |> vec
    	theta_lower   = theta_matrix[:, 2] |> vec
    	theta_upper   = theta_matrix[:, 3] |> vec

    	theta_initials = zeros(main_n_start_pts, length(theta_initial))
    	for i=1:main_n_start_pts
    		theta_initials[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials[i, :] = min.(vec(theta_matrix[:,3]), theta_initials[i, :])
    		theta_initials[i, :] = max.(vec(theta_matrix[:,2]), theta_initials[i, :])
    	end

    	theta_initials_boot = zeros(boot_n_start_pts, length(theta_initial))
    	for i=1:boot_n_start_pts
    		theta_initials_boot[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials_boot[i, :] = min.(vec(theta_matrix[:,3]), theta_initials_boot[i, :])
    		theta_initials_boot[i, :] = max.(vec(theta_matrix[:,2]), theta_initials_boot[i, :])
    	end

    	# see inside this function for parameter names (order)
    	theta_dict = get_params_fn(theta_initial, param_factors)

	ESTIMATION_PARAMS = Dict{String, Any}(
		"param_factors" => param_factors,
		"theta_lower" => theta_lower,
		"theta_upper" => theta_upper,
		"theta_initials" => theta_initials,
		"theta_initials_boot" => theta_initials_boot
	)
    println("Model parameters and initial conditions loaded")
    return ESTIMATION_PARAMS
end

### moment function (only takes a vector of parameters)
function model_dt_static(theta, DATA, PARAMS)
	# BIGMAT, CHARGEMAT, COMPMAT, DATAMOMS
	sim_model = solve_nestedlogit(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		BIGMAT=DATA["BIGMAT"],
		CHARGEMAT=DATA["CHARGEMAT"],
		COMPMAT=DATA["COMPMAT"],
		DATAMOMS=DATA["DATAMOMS"]
	)

	return moments_nestedlogitmodel_all(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		SIMMODEL=sim_model,
		DATAMOMS=DATA["DATAMOMS"],
		include_dt_het=false,
		include_area_main=true,
		include_area_het=false
	)
end






##
model_2_dt_dynamic = 1
"""
Model 2: DT + dynamic nested logit
 - hastar_inversion_type="distribution": how to invert
 - transition_moment_type=0: include 11th moment on transition probabilities?
 	0 = no
	1 = week 1 early 1->0
	2 = all weeks 1->0
	3 = all weeks 0->1
"""
function model_dt_dynamic_params(;
	theta_fixed_dict, gmm_options=nothing, main_n_start_pts=1, boot_n_start_pts=1)

	if ~isnothing(gmm_options)
		main_n_start_pts = gmm_options["main_n_start_pts"]
		boot_n_start_pts = gmm_options["boot_n_start_pts"]
	end

        # full parameter list with expansion factors
	# this is the order that will be used in optimization
	# will next drop from optimization any parameters held fixed
        param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}} = [
		("alpha", 				6.0, nothing),   ## x 6
		("delta",      			0.01, nothing),  ## /100
		("gamma", 				1.0, nothing),
		("gamma_01_factor", 	1.0, nothing),
		("fe_early", 			1.0, nothing),
		("fe_late", 			1.0, nothing),
		("fe_w1", 				1.0, nothing),
		("fe_w2", 				1.0, nothing),
		("fe_w3", 				1.0, nothing),
		("fe_w4", 				1.0, nothing),
		("beta_e", 				6.0, nothing),  ## x 6
		("beta_l", 				6.0, nothing),  ## x 6
		("sigma_dt", 			1.0, nothing),
		("mu", 					1.0, nothing),
		("prob_respond_alpha", 0.01, nothing), ## /100
		("prop_to_wage", 		1.0, nothing)
	]
	all_param_names = [el[1] for el in param_factors]

	# update the given dictionary (overwrite)
	param_included_vec = ones(Bool, length(all_param_names))
	for mykey in keys(theta_fixed_dict)
		if ~in(mykey, all_param_names)
		    error("Fixed parameter name is not in the current model:", mykey)
		end

		idx = findfirst(all_param_names .== mykey)

		# do not include the parameter in estimation
		param_included_vec[idx] = false

		# set paramter value
		param_factors[idx] = (param_factors[idx][1], param_factors[idx][2], theta_fixed_dict[mykey])
	end

	# initial conditions and bounds
		theta_matrix = [
			100.0 0.0 500.0;  # alpha in INR per 10 MINUTES
			50.0 0.0 99.0; # delta x 100   <------------------------ 95 limit!
			100.0 0.0 Inf;  # gamma = switching cost
			1.0 0.0 Inf;  # gamma = switching cost fact  <- never really estimated
			0.0 -Inf Inf; # fe_early = FE of early treatment
			0.0 -Inf Inf; # fe_late = FE of late treatment
			0.0 -Inf Inf;  # time week 1
			0.0 -Inf Inf;  # time week 2
			0.0 -Inf Inf;  # time week 3
			0.0 -Inf Inf;  # time week 4
			50.0 1.0 Inf;   # beta_early per 10 MINUTES
			50.0 1.0 Inf;   # beta_late  per 10 MINUTES
			25.0 1.0 Inf;   # sigma_dt = logit parameter for departure time
			25.0 1.0 Inf;   # mu = logit parameter for route choice
			50.0 0.0 100.0; # probability to respond x 100
			0.0 0.0 1.0 # prop_to_wage (never estimated)
		]

	theta_matrix = theta_matrix[param_included_vec, :]

    	theta_initial = theta_matrix[:, 1] |> vec
    	theta_lower   = theta_matrix[:, 2] |> vec
    	theta_upper   = theta_matrix[:, 3] |> vec

    	theta_initials = zeros(main_n_start_pts, length(theta_initial))
    	for i=1:main_n_start_pts
    		theta_initials[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials[i, :] = min.(vec(theta_matrix[:,3]), theta_initials[i, :])
    		theta_initials[i, :] = max.(vec(theta_matrix[:,2]), theta_initials[i, :])
    	end

    	theta_initials_boot = zeros(boot_n_start_pts, length(theta_initial))
    	for i=1:boot_n_start_pts
    		theta_initials_boot[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials_boot[i, :] = min.(vec(theta_matrix[:,3]), theta_initials_boot[i, :])
    		theta_initials_boot[i, :] = max.(vec(theta_matrix[:,2]), theta_initials_boot[i, :])
    	end

    	# see inside this function for parameter names (order)
    	theta_dict = get_params_fn(theta_initial, param_factors)

	ESTIMATION_PARAMS = Dict{String, Any}(
		"param_factors" => param_factors,
		"theta_lower" => theta_lower,
		"theta_upper" => theta_upper,
		"theta_initials" => theta_initials,
		"theta_initials_boot" => theta_initials_boot
	)
	println("Model parameters and initial conditions loaded")
    return ESTIMATION_PARAMS
end

### moment function (only takes a vector of parameters)
function model_dt_dynamic(theta, DATA, PARAMS; dt_factor=1.0)
	sim_model = solve_dynamicvot(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		BIGMAT=DATA["BIGMAT"],
		CHARGEMAT=DATA["CHARGEMAT"],
		COMPMAT=DATA["COMPMAT"],
		DATAMOMS=DATA["DATAMOMS"],
		hastar_inversion_type=PARAMS["hastar_inversion_type"]
	)

	return moments_fullgmm(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		SIMMODEL=sim_model,
		DATAMOMS=DATA["DATAMOMS"],
		include_dynamicarea_main=true,
		transition_moment_type=PARAMS["transition_moment_type"],
        dt_factor=dt_factor
	)
end


##
model_3_nodt_dynamic = 1
"""
Model 3: dynamic nested logit without departure time
 - transition_moment_type=0: include 11th moment on transition probabilities?
 	0 = no
	1 = week 1 early 1->0
	2 = all weeks 1->0
	3 = all weeks 0->1
"""
function model_nodt_dynamic_params(;
	theta_fixed_dict, gmm_options=nothing, main_n_start_pts=1, boot_n_start_pts=1)

	if ~isnothing(gmm_options)
		main_n_start_pts = gmm_options["main_n_start_pts"]
		boot_n_start_pts = gmm_options["boot_n_start_pts"]
	end

        # full parameter list with expansion factors
	# this is the order that will be used in optimization
	# will next drop from optimization any parameters held fixed
        param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}} = [
		("alpha", 		6.0, nothing),   ## x 6
		("delta",      0.01, nothing),  ## /100
		("gamma", 		1.0, nothing),
		("gamma_01_factor", 1.0, nothing),
		("fe_early", 	1.0, nothing),
		("fe_late", 	1.0, nothing),
		("fe_w1", 		1.0, nothing),
		("fe_w2", 		1.0, nothing),
		("fe_w3", 		1.0, nothing),
		("fe_w4", 		1.0, nothing),
		# ("beta_e", 	6.0, nothing),  ## x 6
		# ("beta_l", 	6.0, nothing),  ## x 6
		# ("sigma_dt", 	1.0, nothing),
		("mu", 			1.0, nothing),
		("prob_respond_alpha", 0.01, nothing), ## /100
		("prop_to_wage", 1.0, nothing)
	]
	all_param_names = [el[1] for el in param_factors]

	# update fixed parameter names
	param_included_vec = ones(Bool, length(all_param_names))
	for mykey in keys(theta_fixed_dict)
		# println("processing ", mykey)
		if ~in(mykey, all_param_names)
		    error("Fixed parameter name is not in the current model:", mykey)
		end

		idx = findfirst(all_param_names .== mykey)
		@assert ~isnothing(idx)

		# do not include the parameter in estimation
		param_included_vec[idx] = false

		# set paramter value
		param_factors[idx] = (param_factors[idx][1], param_factors[idx][2], theta_fixed_dict[mykey])
	end

	# initial conditions and bounds
		theta_matrix = [
			100.0 0.0 500.0;  # alpha in INR per 10 MINUTES
			50.0 0.0 99.0; # delta x 100   <------------------------ 95 limit!
			100.0 0.0 Inf;  # gamma = switching cost
			1.0 0.0 Inf;  # gamma = switching cost fact  <- never really estimated
			0.0 -Inf Inf; # fe_early = FE of early treatment
			0.0 -Inf Inf; # fe_late = FE of late treatment
			0.0 -Inf Inf;  # time week 1
			0.0 -Inf Inf;  # time week 2
			0.0 -Inf Inf;  # time week 3
			0.0 -Inf Inf;  # time week 4
			# 50.0 1.0 Inf;   # beta_early per 10 MINUTES
			# 50.0 1.0 Inf;   # beta_late  per 10 MINUTES
			# 25.0 1.0 Inf;   # sigma_dt = logit parameter for departure time
			25.0 1.0 Inf;   # mu = logit parameter for route choice
			50.0 0.0 100.0  # probability to respond x 100
			0.0 0.0 1.0 # prop_to_wage
		]

	theta_matrix = theta_matrix[param_included_vec, :]

    	theta_initial = theta_matrix[:, 1] |> vec
    	theta_lower   = theta_matrix[:, 2] |> vec
    	theta_upper   = theta_matrix[:, 3] |> vec

    	theta_initials = zeros(main_n_start_pts, length(theta_initial))
    	for i=1:main_n_start_pts
    		theta_initials[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials[i, :] = min.(vec(theta_matrix[:,3]), theta_initials[i, :])
    		theta_initials[i, :] = max.(vec(theta_matrix[:,2]), theta_initials[i, :])
    	end

    	theta_initials_boot = zeros(boot_n_start_pts, length(theta_initial))
    	for i=1:boot_n_start_pts
    		theta_initials_boot[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials_boot[i, :] = min.(vec(theta_matrix[:,3]), theta_initials_boot[i, :])
    		theta_initials_boot[i, :] = max.(vec(theta_matrix[:,2]), theta_initials_boot[i, :])
    	end

    	# see inside this function for parameter names (order)
    	theta_dict = get_params_fn(theta_initial, param_factors)

	ESTIMATION_PARAMS = Dict{String, Any}(
		"param_factors" => param_factors,
		"theta_lower" => theta_lower,
		"theta_upper" => theta_upper,
		"theta_initials" => theta_initials,
		"theta_initials_boot" => theta_initials_boot
	)
	println("Model parameters and initial conditions loaded")
    return ESTIMATION_PARAMS
end

### moment function (only takes a vector of parameters)
function model_nodt_dynamic(theta, DATA, PARAMS)
	sim_model = solve_dynamicvot_nodt(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		BIGMAT=DATA["BIGMAT"],
		CHARGEMAT=DATA["CHARGEMAT"],
		COMPMAT=DATA["COMPMAT"],
		DATAMOMS=DATA["DATAMOMS"]
	)

	return moments_fullgmm(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		SIMMODEL=sim_model,
		DATAMOMS=DATA["DATAMOMS"],
		include_dt_main=false,
		include_dynamicarea_main=true,
		transition_moment_type=PARAMS["transition_moment_type"]
	)
end




##
model_4_dt = 1
"""
Model 4: departure time without
for convenience we use the nested logit from the static route choicemodel
but in the moments we only use the direct route
hence, the value of mu doesn't matter
 - dt_free_1route and dt_wcha_1route
"""
function model_dt_params(;
	theta_fixed_dict, gmm_options=nothing, main_n_start_pts=1, boot_n_start_pts=1)

	if ~isnothing(gmm_options)
		main_n_start_pts = gmm_options["main_n_start_pts"]
		boot_n_start_pts = gmm_options["boot_n_start_pts"]
	end
        # full parameter list with expansion factors
	# this is the order that will be used in optimization
	# will next drop from optimization any parameters held fixed
        param_factors::Vector{Tuple{String,Float64,Union{Nothing, Float64}}} = [
		("alpha",	 6.0, nothing),  ## x 6
		("beta_e", 	 6.0, nothing),  ## x 6
		("beta_l", 	 6.0, nothing),  ## x 6
		("sigma_dt", 1.0, nothing),
		("mu", 		 1.0, nothing),   # SET THIS SO WE CAN COMPUTE ROUTE BUT DOESN'T AFFECT RESULTS
		("prob_respond_alpha", 0.01, nothing) ## /100
	]
	all_param_names = [el[1] for el in param_factors]

	# update the given dictionary (overwrite)
	param_included_vec = ones(Bool, length(all_param_names))
	for mykey in keys(theta_fixed_dict)
		if ~in(mykey, all_param_names)
		    error("Fixed parameter name is not in the current model:", mykey)
		end

		idx = findfirst(all_param_names .== mykey)

		# do not include the parameter in estimation
		param_included_vec[idx] = false

		# set paramter value
		param_factors[idx] = (param_factors[idx][1], param_factors[idx][2], theta_fixed_dict[mykey])
	end

	# initial conditions and bounds
    	theta_matrix = [
		   100.0 0.0 500.0;  # alpha in INR per 10 MINUTES
            50.0 1.0 Inf;   # beta_early per 10 MINUTES
            50.0 1.0 Inf;   # beta_late  per 10 MINUTES
            25.0 1.0 Inf;   # sigma_dt = logit parameter for departure time
			25.0 1.0 Inf;   # mu
            50.0 0.0 100.0  # probability to respond x 100
    		]

	theta_matrix = theta_matrix[param_included_vec, :]

    	theta_initial = theta_matrix[:, 1] |> vec
    	theta_lower   = theta_matrix[:, 2] |> vec
    	theta_upper   = theta_matrix[:, 3] |> vec

    	theta_initials = zeros(main_n_start_pts, length(theta_initial))
    	for i=1:main_n_start_pts
    		theta_initials[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials[i, :] = min.(vec(theta_matrix[:,3]), theta_initials[i, :])
    		theta_initials[i, :] = max.(vec(theta_matrix[:,2]), theta_initials[i, :])
    	end

    	theta_initials_boot = zeros(boot_n_start_pts, length(theta_initial))
    	for i=1:boot_n_start_pts
    		theta_initials_boot[i, :] = theta_initial .* (0.5 .+ rand(length(theta_initial)))
    		theta_initials_boot[i, :] = min.(vec(theta_matrix[:,3]), theta_initials_boot[i, :])
    		theta_initials_boot[i, :] = max.(vec(theta_matrix[:,2]), theta_initials_boot[i, :])
    	end

    	# see inside this function for parameter names (order)
    	theta_dict = get_params_fn(theta_initial, param_factors)

	ESTIMATION_PARAMS = Dict{String, Any}(
		"param_factors" => param_factors,
		"theta_lower" => theta_lower,
		"theta_upper" => theta_upper,
		"theta_initials" => theta_initials,
		"theta_initials_boot" => theta_initials_boot
	)
	println("Model loaded")
    return ESTIMATION_PARAMS
end

### moment function (only takes a vector of parameters)
function model_dt(theta, DATA, PARAMS)
	sim_model = solve_nestedlogit(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		BIGMAT=DATA["BIGMAT"],
		CHARGEMAT=DATA["CHARGEMAT"],
		COMPMAT=DATA["COMPMAT"],
		DATAMOMS=DATA["DATAMOMS"]
	)

	return moments_logitmodel_all(
		theta=theta,
		param_factors=PARAMS["param_factors"],
		SIMMODEL=sim_model,
		DATAMOMS=DATA["DATAMOMS"],
	)
end
