module OPF

using Distributed
using JuMP
using Gurobi
using Distributions
using Random
using Statistics
import Base: @kwdef

const GRB_ENV = Ref{Gurobi.Env}()
function __init__()
    GRB_ENV[] = Gurobi.Env()
    return
end

include("utils.jl")
include("data_structures.jl")
include("tree_structure.jl")
include("deterministic_model.jl")
include("nodal_stochastic_model.jl")
include("master_problem.jl")
include("subproblem.jl")
include("CG.jl")

export Tree,
    TestNetwork,
    TestSystemData,
    ModelDataNodalStochastic,
    SolutionInfo,
    CGSolutionInfo,
    solve_stoch_model_FS,
    solve_CG
end
