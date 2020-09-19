#=
Workshop #7 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#

# Question 1
# data generation function as instructed by Q1
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

# Perry: I dislike complex for-loops so here's a vectorized implementation
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

# test run for basis function, passes all tests
test_x = Vector{Float64}([1.618, 2.71828, 3.14159, 420.69])
basis_vec_f = generate_lagrange_basis_vector_func(test_x)
println(basis_vec_f.(test_x)) # should return identity matrix
println(basis_vec_f(20.20)) # should return a vector of size 4

# Question 2
using Plots
using LaTeXStrings
x_g_vec = x_g(100) # 100 equally spaced points between 0 and 1 for graphing
x_basis_vec = x_i(6) # 6 points between 0 and 1 for basis
# function that returns 6 values for a scalar input
basis_vec_f = generate_lagrange_basis_vector_func(x_basis_vec)
# should see a [100, 6] output, which is 6 output for each input value
out = basis_vec_f.(x_g_vec)
# turns out to be a vector of 100 vectors of size 6
# we need to turn it into mult of [100, 6] for plotting
# https://stackoverflow.com/a/52257481/11652747
out = transpose(hcat(out...))
p = plot(x_g_vec, out, title="Question 2",
     label=[L"f_{k=1}" L"f_{k=2}" L"f_{k=3}" L"f_{k=4}" L"f_{k=5}" L"f_{k=6}"],
     legend=:bottomright, xaxis=L"x", yaxis=L"f(x)")
display(p)

# Question 3
function do_question_3()
    abs_outs::Array{Vector{Float64}, 1} = []
    labels::Vector{String} = []
    maximums::Vector{Float64} = []
    ks = [5, 7, 9, 11]
    ns = [9, 13, 17, 21]
    for i in 1:4
        k = ks[i]
        n = ns[i]
        x_basis_vec = x_i(n)
        # function that returns n values for 1 input
        basis_vec_f_n = generate_lagrange_basis_vector_func(x_basis_vec)
        # function that returns 1 value for 1 input
        basis_f_k_n = (x->basis_vec_f_n(x)[k])
        abs_out = abs.(basis_f_k_n.(x_g_vec))
        m = maximum(abs_out)
        lab = string("k=", k, ", n=", n, ", max=", m)
        push!(abs_outs, abs_out)
        push!(labels, lab)
        push!(maximums, m)
    end
    p = plot(x_g_vec, abs_outs, title="Question 3",
             xaxis=L"x", yaxis=L"|f_{k;n}(x)|",
             label=permutedims(labels))
    display(p)
    p = plot(ns, maximums, title="Question 3",
             xaxis=L"n", yaxis=L"max(|f_{k;n}(x)|)",
             lw=2, marker = ([:hex :d], 3, 0.8), legend=false)
    display(p)
end
do_question_3()

# Question 4
x_basis_vec = x_i(21)
# function that returns n values for 1 input
basis_vec_f_n21 = generate_lagrange_basis_vector_func(x_basis_vec)
# function that returns 1 value for 1 input
basis_f_k11_n21 = (x->basis_vec_f_n21(x)[11])
x_g_vec_5000 = Vector(range(-2, 3, length=5000)) # 5000 datapoints to graph
maximums = abs.(basis_f_k11_n21.(x_g_vec_5000)) # dims = [5000, 1]
p = plot(x_g_vec_5000, maximums, title="Question 4",
         xaxis=L"x_{extrapolation}", yaxis=L"|f_{k;n}(x)|", legend=false)
display(p)
