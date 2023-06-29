include("paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

include("analysis_estimation/paper_table_II.jl")
include("analysis_estimation/paper_table_robustness.jl")
include("analysis_estimation/paper_table_delta.jl")
include("analysis_estimation/paper_figure_model_fit.jl")

include("analysis_estimation/model_identification/model_sensitivity_full.jl")
include("analysis_estimation/model_identification/model_sensitivity_vott.jl")
include("analysis_estimation/model_identification/model_sensitivity_simple.jl")
include("analysis_estimation/model_identification/table_sensitivity_vott.jl")

include("analysis_estimation/numerical_identification/paper_figure_sm_4.jl")
include("analysis_estimation/numerical_identification/paper_table_sm_12_part1.jl")