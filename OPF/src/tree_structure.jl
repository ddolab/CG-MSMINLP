#code adapted from: https://github.com/kirui93/ScenTrees.jl

"""
MIT License

Copyright (c) 2019 Kipngeno Kirui <kipngenokirui1993@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

rng = MersenneTwister(03102023);

#This defines the tree structure.
mutable struct Tree
    name::String                        # name of the tree
    parent::Vector{Int64}               # parents of nodes in the tree
    children::Vector{Vector{Int64}}     # successor nodes of each parent
    probability::Matrix{Float64}        # probability to go from one node to another

    function child(parent::Vector{Int64})
        allchildren = Vector{Vector{Int64}}([])
        for node in unique(parent)
            push!(
                allchildren,
                [i for i in 1:length(parent) if parent[i] == node],
            )
        end
        return allchildren
    end

    function Tree(bstructure::Vector{Int64})
        self = new()
        leaves = bstructure
        for stage in 1:length(bstructure)
            if stage == 1
                leaves = 1
                self.name = "Tree $(bstructure[1])"
                self.parent = zeros(Int64, bstructure[1])
                self.children = child(self.parent)
                self.probability = ones(bstructure[1], 1)
            else
                leaves = leaves * bstructure[stage-1]
                newleaves = vec(
                    ones(Int64, bstructure[stage]) .*
                    transpose(length(self.parent) .+ (1 .- leaves:0)),
                )
                self.parent = vcat(self.parent, newleaves)
                self.children = child(self.parent)
                self.name = "$(self.name)x$(bstructure[stage])"
                
                # Filler probabilities. Modified later when data instance is created
                tmp = rand(rng, Uniform(0.3, 1.0), bstructure[stage], leaves)
                tmp =
                    tmp ./ (
                        transpose(ones(1, bstructure[stage])) .*
                        sum(tmp; dims = 1)
                    )
                self.probability = vcat(self.probability, vec(tmp))
            end
        end
        return self
    end
end

function stage(trr::Tree, node = Int64[])
    if isempty(node)
        node = 1:length(trr.parent)
    elseif isa(node, Int64)
        node = Int64[node]
    end
    stage = zero(node)
    for i in 1:length(node)
        pred = node[i]
        while pred > 0 && trr.parent[pred] > 0
            pred = trr.parent[pred]
            stage[i] += 1
        end
    end
    return stage
end

function leaves(trr::Tree, node = Int64[])
    nodes = 1:length(trr.parent)
    leaves = setdiff(nodes, trr.parent)
    omegas = 1:length(leaves)
    if !isempty(node) && isa(node, Int64)
        node = Int64[node]
    end
    if !isempty(node) && (0 ∉ node)
        omegas = Set{Int64}()
        nodes = leaves
        while any(j != 0 for j in nodes)
            omegas = union(
                omegas,
                (ind for (ind, j) in enumerate(nodes) if j in node),
            )
            nodes = Int64[trr.parent[max(0, j)] for j in nodes]
        end
        omegas = collect(omegas)
    end
    leaves = Int64[leaves[j] for j in omegas]
    prob = ones(Float64, length(leaves))
    nodes = leaves
    while any(j != 0 for j in nodes)
        prob = prob .* trr.probability[[j for j in nodes]]
        nodes = Int64[trr.parent[max(0, j)] for j in nodes]
    end
    return leaves, omegas, prob
end

function nodes(trr::Tree, t = Int64[])
    nodes = 1:length(trr.parent)
    if isempty(t) #if stage t is not given, return all nodes of the tree
        return nodes
    else # else return nodes at the given stage t
        stg = stage(trr)
        return Int64[i for i in nodes if stg[i] == t]
    end
end

function root(trr::Tree, nodes = Int64[])
    if isempty(nodes)
        nodes = trr.children[1]
    elseif isa(nodes, Int64)
        nodes = Int64[nodes]
    end
    root = 1:length(trr.parent)
    for i in nodes
        iroot = Vector{Int64}([])
        tmp = i
        while tmp > 0
            push!(iroot, tmp)
            tmp = trr.parent[tmp]
        end
        root = Int64[i for i in reverse(iroot) if i in root]
    end
    return root
end

function total_node_probability(tree::Tree, n::Int64)
    path = root(tree, n)

    p = 1.0
    for node in path
        p *= tree.probability[node]
    end
    return p
end

function sibling_nodes(tree::Tree)
    N = maximum(nodes(tree))

    siblings = Dict{Int64, Vector{Int64}}()
    for n in 1:(minimum(leaves(tree)[1])-1)
        siblings[n] = Int64[]
        populating = true
        encountered_siblings = false
        for n_ in n+1:N 
            if tree.parent[n_] == n
                push!(siblings[n], n_)
                encountered_siblings = true
                continue
            elseif encountered_siblings && !populating
                break
            end
            populating = false
        end
    end

    return siblings
end

function pairs_of_sibling_nodes(tree::Tree)
    sibling_dictionary = sibling_nodes(tree)
    sibling_pairs = []
    for n in 1:(minimum(leaves(tree)[1])-1)
        push!(sibling_pairs, create_all_possible_pairs(sibling_dictionary[n]))
    end
    return vcat(sibling_pairs...)
end

function pairs_of_sibling_nodes(tree::Tree, sibling_pairs, sp_solutions::Vector{SPSolutionInfo}, columns)
    pairs_with_column_info = []
    
    for (source, target) in sibling_pairs
        col_index = 1
        for col in sp_solutions[source].columns_info
            if col.add_column && !(col.column in columns[target])
                push!(pairs_with_column_info, (source, target, col_index))
            end
            col_index += 1
        end
    end

    return pairs_with_column_info
end

function nodes_in_stage(tree::Tree, stage_number::Int64)
    member_nodes = Int64[]
    for node in nodes(tree)
        if stage(tree, node)[1] + 1 == stage_number
           push!(member_nodes, node) 
        end
    end
    return member_nodes
end