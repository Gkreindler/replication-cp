include("paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

include("analysis_estimation/estimation_main_benchmark.jl")
include("analysis_estimation/estimation_main_fixed_vott.jl")
include("analysis_estimation/estimation_main_no_dt.jl")

include("analysis_estimation/estimation_robust1_static.jl")
include("analysis_estimation/estimation_robust2_asym_switch.jl")
include("analysis_estimation/estimation_robust3_FE.jl")
include("analysis_estimation/estimation_robust4_half_attention.jl")
include("analysis_estimation/estimation_robust5_propto_wage.jl")
include("analysis_estimation/estimation_robust6_single_arrival_time.jl")

include("analysis_estimation/estimation_delta_other.jl")
include("analysis_estimation/estimation_delta_d99.jl")
include("analysis_estimation/estimation_delta_estimation.jl")

include("analysis_estimation/numerical_identification/simulate_data_numerical_check.jl")
include("analysis_estimation/numerical_identification/estimation_simulated_data_asy.jl")
include("analysis_estimation/numerical_identification/estimation_simulated_data_fin.jl")