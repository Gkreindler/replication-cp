include("paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"
using Pkg
Pkg.activate(env_path)
Pkg.instantiate()

include("analysis_estimation/runsims-main.jl")
include("analysis_estimation/runsims-main-boot.jl")

include("analysis_estimation/runsims-het.jl")
include("analysis_estimation/runsims-het-boot.jl")

include("analysis_estimation/runsims-varyparams.jl")

include("analysis_estimation/runsims-power.jl")

include("analysis_estimation/runsims-2routes.jl")

include("analysis_estimation/runsims-ext.jl")

include("analysis_estimation/runsims-factor.jl")




