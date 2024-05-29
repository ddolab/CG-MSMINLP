# Create all possible pairs of elements in array vec (order matters)
# For example: vec = [3,5,4]
# returns [[3,5],[5,3],[3,4],[4,3],[5,4],[4,5]]
function create_all_possible_pairs(vec)
    all_pairs = []
    for i in 1:length(vec)
        for j in i:length(vec)
            if i!=j
                push!(all_pairs, [vec[i], vec[j]])
                push!(all_pairs, [vec[j], vec[i]])
            end
        end
    end
    return all_pairs
end

function is_solution_binary(x::Real; kwargs...)
    return isapprox(x, 0; kwargs...) || isapprox(x, 1; kwargs...)
end