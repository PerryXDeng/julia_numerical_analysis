#=
Workshop #8 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
# Question 1
function generate_x_basis(n::Int)::Vector{Float64}
    # returns a vector of size n for basis, equally spaced between 0 and 1
    # including 0 and 1
    return Vector{Float64}(0:n-1) ./ Float64(n-1)
end
function generate_x_graph(M::Int)::Vector{Float64}
    # returns a vector of size M for graphing
    g = Vector{Float64}(1:M)
    return (g .- 0.5) ./ Float64(M)
end

function incorporate_new_data_into_triangle!(triangle::Vector{Vector{Float64}},
                                          x::Vector{Float64},
                                          y::Vector{Float64},
                                          new_data_index::Int)
    num_existing = length(triangle)
    push!(triangle, Vector{Float64}(undef, new_data_index))
    new_row = triangle[num_existing + 1]
    new_row[num_existing + 1] = y[new_data_index]

    if num_existing > 0
        new_x_i = x[new_data_index]
        previous_row = triangle[num_existing]
        for i in num_existing:-1:1
            # has to be calculated sequentially
            numerator = new_row[i+1] - previous_row[i]
            denominator = new_x_i - x[i]
            new_row[i] = numerator/denominator
        end
    end
end
function construct_triangle(x::Vector{Float64}, y::Vector{Float64})::Vector{Vector{Float64}}
    # x, y should have the same size
    triangle::Vector{Vector{Float64}} = []
    for i in 1:length(x)
        # has to be inserted one at a time
        incorporate_new_data_into_triangle!(triangle, x, y, i)
    end
    return triangle
end
function get_triangle_poly_coefs(triangle::Vector{Vector{Float64}})::Vector{Float64}
    return getindex.(triangle, 1) # gets the first element of all rows
end


# test triangle, using one from the slides, should print [1, 0.5, 0,5]. passes
println(get_triangle_poly_coefs(construct_triangle([0.0, 2.0, 3.0], [1.0, 2.0, 4.0])))
