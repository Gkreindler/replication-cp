include("paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

include("analysis_simulations/paper_table_IV.jl")
include("analysis_simulations/paper_table_V.jl")

include("analysis_simulations/paper_figure5.jl")

include("analysis_simulations/paper_smfig8A_charges.jl")

include("analysis_simulations/paper_smfig8D_factor.jl")
include("analysis_simulations/paper_smfig9_decomposition.jl")

include("analysis_simulations/paper_smfig10_2routes.jl")