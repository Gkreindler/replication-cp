## Prep
include("../../paths.jl") # defines the var `rootpath_base` that points to the root of the replication package
env_path = rootpath_base * "code/analysis_estimation/cp_env/"

using Pkg
Pkg.activate(env_path) # the dollar sign takes the variable from the local (main) worker
Pkg.instantiate()

using CSV
using DataFrames

include("../functions-table.jl")

## Paths
rootpath_gmm = rootpath_base * "model/"

savepath = rootpath_base * "paper/figures/smfigure3/"
isdir(savepath) || mkdir(savepath)

## Read full and vott-only (nodt) model jacobians and parameters

    G_full = CSV.read(savepath * "jacobian_full.csv", DataFrame)
    G1 = Matrix(G_full)[1:3,2:4]
    θfull = Matrix(G_full)[4,2:4]

    Gnodt = CSV.read(savepath * "jacobian_nodt.csv", DataFrame)
    G2 = Matrix(Gnodt)[1:3,2:4]
    θnodt = Matrix(Gnodt)[4,2:4]

    G_simple = CSV.read(savepath * "jacobian_simple.csv", DataFrame)
    G3 = Matrix(G_simple)[1:3,2:4]
    θsimp = Matrix(G_simple)[4,2:4]

## make table
tex_contents =
"\\resizebox*{\\textwidth}{!}{" *
"\\begin{tabular}{lccc@{\\hskip 1cm}ccc@{\\hskip 1cm}ccc}\n" *
"   \\toprule\n" *
"   & \\multicolumn{3}{c}{(1)} & \\multicolumn{3}{c}{(2)} & \\multicolumn{3}{c}{(3)}  \\\\ \n" *
"   \\addlinespace\n" *
"   & \\multicolumn{3}{c}{Full Model} & \\multicolumn{3}{c}{No departure time} & \\multicolumn{3}{c}{Simple Model (\$\\delta=0\$)} \\\\ \n" *
"   \\addlinespace\n" *
"   & \$\\alpha\$ & \$\\sigma^R\$ & \$\\sigma^R\$   & \$\\alpha\$ & \$\\gamma\$ & \$\\sigma^R\$   & \$\\alpha\$ & \$\\gamma\$ & \$\\sigma^R\$ \\\\ \n" *
"   \\addlinespace\n"

# estimated values
    model_full = f1(θfull[1] * 6.0) * " & " *
                 f1(θfull[2]) * " & " *
                 f1(θfull[3]) * " & "

    model_nodt = f1(θnodt[1] * 6.0) * " & " *
                 f1(θnodt[2]) * " & " *
                 f1(θnodt[3]) * " & "

    model_simp = f1(θsimp[1]) * " & " *
                 f1(θsimp[2]) * " & " *
                 f1(θsimp[3])

    tex_contents *=
    " Estimated values  & " * model_full * model_nodt * model_simp * "  \\\\ \n" *
    "   \\addlinespace\\addlinespace\n"

# Jacobian
    tex_contents *=
    " \\multicolumn{10}{l}{\\textit{Jacobian: Change in Route 1 take-up Due to Change in Parameter}}  \\\\ \n" *
    "   \\addlinespace\n"

# Jacobian
    tex_contents *=
    " Before Experiment & " *
                f2(G1[1,1]) * " & " *
                f2(G1[1,2]) * " & " *
                f2(G1[1,3]) * " & " *
                f2(G2[1,1]) * " & " *
                f2(G2[1,2]) * " & " *
                f2(G2[1,3]) * " & " *
                f2(G3[1,1]) * " & " *
                f2(G3[1,2]) * " & " *
                f2(G3[1,3]) * " \\\\ \n \\addlinespace\n" *
    " Week 1 (Charges) & " *
                f2(G1[2,1]) * " & " *
                f2(G1[2,2]) * " & " *
                f2(G1[2,3]) * " & " *
                f2(G2[2,1]) * " & " *
                f2(G2[2,2]) * " & " *
                f2(G2[2,3]) * " & " *
                f2(G3[2,1]) * " & " *
                f2(G3[2,2]) * " & " *
                f2(G3[2,3]) * " \\\\ \n \\addlinespace\n" *
    " Week 2 (After Charges) & " *
                f2(G1[3,1]) * " & " *
                f2(G1[3,2]) * " & " *
                f2(G1[3,3]) * " & " *
                f2(G2[3,1]) * " & " *
                f2(G2[3,2]) * " & " *
                f2(G2[3,3]) * " & " *
                f2(G3[3,1]) * " & " *
                f2(G3[3,2]) * " & " *
                f2(G3[3,3]) * " \\\\ \n \\addlinespace\n"

# Close
    tex_contents *=
        "   \\bottomrule \n " *
        "\\end{tabular} \n } \n"

table_path = string(savepath, "table_jacobian_vott.tex")
open(table_path,"w") do io
   println(io,tex_contents)
end

println("done!")
