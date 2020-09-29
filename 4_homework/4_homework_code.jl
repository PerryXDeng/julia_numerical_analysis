#=
Workshop #8 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
# Newton's Divided Differences and Horner's Method implementation
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
        p += p*(x_in - x[i]) + c[i]
    end
    return p
end
function horners_method_vector(x::Vector{Float64}, c::Vector{Float64},
                               x_in::Vector{Float64})::Vector{Float64}
    # if confused about why we have (x, )
    # see https://discourse.julialang.org/t/how-to-broadcast-over-only-certain-function-arguments/19274
    return horners_method.((x,), (c,), x_in)
end

# test triangle, using one from the slides
c = get_triangle_poly_coefs(construct_triangle([0.0, 2.0, 3.0], [1.0, 2.0, 4.0]))
out = horners_method_vector([0.0, 2.0, 3.0], c, [0.0, 2.0])
println(c) #should print [1, 0.5, 0,5]. passes
println(out) #should print [2,4]. passes


# Question 3 Numerical Part (Need someone to analytically estimate error bound)
x = Vector{Float64}(1:10)
y = log.(x) # we are approximating the f(x)=ln(x) function

# 3a. what are the divided differences coefficients?
c = get_triangle_poly_coefs(construct_triangle(x, y))
println(c)
#=
[0.0, 0.6931471805599453, -0.1438410362258904, 0.028316506132566182, -0.004860605094857686, 0.0007260544531257324, -9.536660491298434e-5, 1.1144297217017655e-5, -1.1710707694472781e-6, 1.1169422368964362e-7]
=#

# 3b. numerical interpolation
x_new = Vector{Float64}(1:9) .+ 0.5
y_estimated = horners_method_vector(x, c, x_new)
println(y_estimated)
#=
[0.9165188922906382, 1.2527066051592177, 1.5041015312414938, 1.7047316940606938, 1.871819317912857, 2.0148752083875547, 2.140140249284577, 2.250914263723907, 2.3578923136608574]
=#

# 3c. theoretical error bound
nth_derivative_abs(x::Float64)::Float64 = abs(-362880/(x^10))
nth_d_max = maximum(nth_derivative_abs.(x_new))
error_bound(in::Float64) = abs((nth_d_max*prod(in .- x)/factorial(10)))
println(error_bound.(x_new))

# 3d. empirical errors
y_expected = log.(x_new)
errors = abs.(y_expected - y_estimated)
println(errors)
using Plots
using LaTeXStrings
p1 = plot(x_new, [y_expected y_estimated], label=[L"ln(x)" L"L(x)"],
          xaxis=L"x", marker=([:hex :d], 3, 0.8))
p2 = plot(x_new, errors, xaxis=L"x", yaxis=L"|ln(x) - L(x)|")
display(p1)
display(p2)

# Question 4
function brightness(t::Float64)::Float64
    return 1/(1+t^2)
end

using Statistics
function calculate_lagrange_rms_error(n::Int, f)::Float64
    # f should take in a scalar and output a scalar
    start = -5
    stop = 5
    x = range(-5, stop=5, length=n)
    step_size = Float64(x.step)
    x = Vector{Float64}(x)
    y = f.(x)

    c = get_triangle_poly_coefs(construct_triangle(x, y))

    x_new = (x .+ step_size)[1:n-1] # halfway betwee ngridpoints
    y_actual = f.(x_new)
    y_estim = horners_method_vector(x, c, x_new)
    return sqrt(mean((y_actual .- y_estim).^2))
end

n_options = [5, 10, 15, 20, 25]
rms_for_each_n = calculate_lagrange_rms_error.(n_options, brightness)
# as expected, errors blow up as n increases
p = plot(n_options, rms_for_each_n, xaxis=L"n", yaxis=L"RMS_{error}",
         marker=([:hex :d], 3, 0.8))
display(p)
# this one plots partial data (excluding n = 25)
p = plot(n_options[1:4], rms_for_each_n[1:4], xaxis=L"n", yaxis=L"RMS_{error}",
         marker=([:hex :d], 3, 0.8))
display(p)

# Question 5
t = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
magnitude = [0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561, 0.468, 0.302]
c = get_triangle_poly_coefs(construct_triangle(t, magnitude))

y_estimates = horners_method_vector(t, c, t)
plot([t, t], [magnitude, y_estimates], xaxis=L"t", yaxis="magnitude",
     label=["recorded" "estimated"])
println(y_estimates)
#=
Something very wrong with these estimates past x_1
Output:
[0.3020000000000009, -259.6910000000002, -999.2598750000009, -3020.0496000000003, -7799.168250000004, -17980.1004, -38000.16881250001, -74938.89700000004, -248110.627]
=#
y_estimate_at_point_one = horners_method(t, c, 0.1)
println(y_estimate_at_point_one)
#=
Estimate is way off the given data distribution
Output:
-42.32244999999994
=#
