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
    return (Vector{Float64}(0:n-1) ./ Float64(n-1))
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

# x: x values of interpolating points
# c: array of coefficients taken from top edge of triangle
# x_in: x value to evaluate polynomial at
function horners_method(x::Vector{Float64}, c::Vector{Float64}, x_in::Float64)::Float64
    n = length(x)
    # loop through x_in to evaluate
    p = c[n]
    for i in length(c)-1:-1:1
        p = p*(x_in - x[i]) + c[i]
    end
    return p
end
function horners_method_vector(x::Vector{Float64}, c::Vector{Float64},
                               x_in::Vector{Float64})::Vector{Float64}
    # if confused about why we have (x, )
    # see https://discourse.julialang.org/t/how-to-broadcast-over-only-certain-function-arguments/19274
    return horners_method.((x,), (c,), x_in)
end

# test triangle, using one from the slides, should print [1, 0.5, 0,5]. passes
println(get_triangle_poly_coefs(construct_triangle([0.0, 2.0, 3.0], [1.0, 2.0, 4.0])))
c = get_triangle_poly_coefs(construct_triangle([0.0, 2.0, 3.0], [1.0, 2.0, 4.0]))
println(c)
println(horners_method_vector([0.0, 2.0, 3.0], c, [0.0, 2.0])) # returns 1,2 as expected

# Question 2
using LinearAlgebra, Statistics, Compat
using Plots
using LaTeXStrings
A = zeros(6,6)
D = Diagonal([1,1,1,1,1,1])
M = A + D
println(M)
kx = [1.0,2.0,3.0,4.0,5.0,6.0]
#println(k1)

function x_g(M::Int)::Vector{Float64}
    # returns a vector of size M for graphing
    g = Vector{Float64}(1:M)
    return (g .- 0.5) ./ Float64(M)
end
function x_i(n::Int)::Vector{Float64}
    # returns a vector of size n for basis, equally spaced between 0 and 1
    # including 0 and 1
    return Vector{Float64}(0:n-1) ./ Float64(n-1)
end

function generate_lagrange_basis_vector_func(x::Vector{Float64})
    n = size(x)[1]
    function basis_k(x_interp::Float64, x_data::Vector{Float64}, k::Int)::Float64
        non_k_indices = cat(1:(k-1), (k+1):n, dims=[1]) # dims = [n-1]
        x_data_not_k = x_data[non_k_indices] # dims = [n-1]
        numerators = x_interp .- x_data_not_k # dims = [n-1]
        denominators = x_data[k] .- x_data_not_k # dims = [n-1]
        fractions = numerators ./ denominators # dims = [n-1]
        product = prod(fractions) # dims = [1]
        return product
    end
    # dims = [n]
    # if confused about why we have (x, )
    # see https://discourse.julialang.org/t/how-to-broadcast-over-only-certain-function-arguments/19274
    basis_vector_func = (x_interp::Float64 -> basis_k.(x_interp, (x,), 1:n))
    return basis_vector_func
end

function validate(M, kx)
    p = plot(title="Question 2")
    for i in 1:length(kx)
        ky = M[i,:]
        println(ky)
        coeffs = get_triangle_poly_coefs(construct_triangle(kx/6.0, ky))
        x = x_i(6)
        x_g_vec = x_g(100) # 100 equally spaced points between 0 and 1 for graphing
        x_basis_vec = x_i(6) # 6 points between 0 and 1 for basis
        #println(x_g_vec)
        println(i, " ", x_basis_vec)
        #println(coeffs)
        # function that returns 6 values for a scalar input
        #basis_vec_f = generate_lagrange_basis_vector_func(x_basis_vec)
        # should see a [100, 6] output, which is 6 output for each input value
        out = horners_method_vector(x_basis_vec, coeffs, x_g_vec)
        plot!(x_g_vec, out, legend=:bottomright, xaxis=L"x", yaxis=L"f(x)")
    end
    display(p)
end

validate(M, kx)
