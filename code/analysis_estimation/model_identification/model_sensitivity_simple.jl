
## Prep
include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"

using Pkg
Pkg.activate(env_path) # the dollar sign takes the variable from the local (main) worker
Pkg.instantiate()

using Optim
using FiniteDifferences
using ForwardDiff
using LinearAlgebra
using DataFrames
using CSV

"""
    1. compute transition probabilities
    2. compute route 0 choice probabilities
    3. compare with data and compute objective value (loss)

    constant detour ΔT = 6.4 minutes = 0.106667 hours
    constant price  p = 144 INR (avg of 96 and 192)
"""
function vott_3periods(θ; ΔT=6.4/60.0, p=144)

    α,γ,μ = θ

    # exponential factors
    A = exp(-α*ΔT/μ)
    G = exp(-γ/μ)
    M = exp(-p/μ)

    # transition probabilities
    π00 = 1/(1+G*A)
    π01 = G/(G+A)
    π00p = M/(M+G*A)
    π01p = G*M/(G*M+A)

    # route 0 choice probabilities
    s0 = π01 / (1 - π00 + π01)
    s1 = s0 * π00p + (1-s0)*π01p
    s2 = s1 * π00  + (1-s1)*π01

    # use prob(route=1)
    return 1 .- [s0, s1, s2]
end

function obj_value(θ;s_vec_data=[0.1, 0.4, 0.25])

    s_vec_model = vott_3periods(θ)

    obj_val = s_vec_model .- s_vec_data
    obj_val = dot(obj_val, obj_val)

    return obj_val
end

# Optimize
    θ0 = [100.0, 100.0, 50.0]

    optim_results = optimize(obj_value, θ0; autodiff = :forward)
    optim_results |> display
    θ = Optim.minimizer(optim_results)

# jacobian
    # myjac = jacobian(central_fdm(5, 1, factor=myfactor), mymomfunction_main_avg, theta_optimum)

    myjac = ForwardDiff.jacobian(vott_3periods, θ)

    # normalize by θ
    myjac_norm = myjac .* reshape(θ, (1,3))

    myjac_norm |> display

# save
    ## Paths
	rootpath_gmm = rootpath_base * "model/"

	savepath = rootpath_base * "paper/figures/smfigure3/"
	isdir(savepath) || mkdir(savepath)

	myjacobian_df = DataFrame(myjac_norm, :auto)
	push!(myjacobian_df, θ)
	myjacobian_df[!, :moment] = ["w0", "w1", "w2", "optimal"]
	select!(myjacobian_df, [:moment, :x1, :x2, :x3])
	rename!(myjacobian_df, [:moment, :alpha, :gamma, :mu])
	CSV.write(savepath * "jacobian_simple.csv", myjacobian_df)
