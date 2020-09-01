#====
Numerical Analysis
Workshop #3

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#

using Plots
using Statistics

function make_plots(forward_abs_errors::Array{Float64},
                    backward_abs_errors::Array{Float64},
                    meta_info::String, question::String)
    title = string(question, "\n", meta_info)
    num_iterations = size(forward_abs_errors, 1)
    p = plot(1:num_iterations, forward_abs_errors,
         title=title, legend=false)
    ylabel!(p, "|f(x)|")
    xlabel!(p, "n")
    display(p)
    p = plot(1:num_iterations, backward_abs_errors,
         title=title, legend=false)
    ylabel!(p, "|x-7|")
    xlabel!(p, "n")
    display(p)
    log_forward_errors = log.(forward_abs_errors)
    log_errors_n = log_forward_errors[2:num_iterations]
    log_errors_n_minus_one = log_forward_errors[1:(num_iterations-1)]
    delta_y = log_errors_n[2:end] .- log_errors_n[1:end-1]
    delta_x = log_errors_n_minus_one[2:end] .- log_errors_n_minus_one[1:end-1]
    dy_dx = delta_y ./ delta_x
    dy_dx_avg = mean(dy_dx)
    p = plot(log_errors_n_minus_one, log_errors_n,
         title=string(title, ", avg_slope=", dy_dx_avg), legend=false)
    ylabel!(p, "ln|f(x_n)|")
    xlabel!(p, "ln|f(x_n_minus_one)|")
    display(p)
    p = plot(2:num_iterations-1, dy_dx,
         title=title, legend=false)
    ylabel!(p, "dy/dx")
    xlabel!(p, "n")
    display(p)
end

function secant_method(f, x0::Float64, x1::Float64, max_iter::Int,
                       threshold::Float64, question::String,
                       visualize::Bool=true)
    s = string("Secant Method x0=", x0, " x1=", x1)
    println(s)
    x2::Float64 = 0.0
    forward_abs_errors::Array{Float64} = Float64[] # plot 1, plot 3
    backward_abs_errors::Array{Float64} = Float64[] # plot 2
    for n in 0:max_iter
        fx1::Float64 = f(x1)
        denominator::Float64 = fx1 - f(x0)
        if denominator == 0.0
            info = string("Singularity, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)
            end
            return 1
        end
        x2 = x1 - fx1 * (x1 - x0) / denominator
        x0, x1 = x1, x2
        forward_abs_error = abs(f(x2))
        push!(forward_abs_errors, forward_abs_error)
        backward_abs_error = abs(x2 - 7.0)
        push!(backward_abs_errors, backward_abs_error)
        if forward_abs_error < threshold
            info = string(":) Converged, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)

            end
            return 0
        end
    end
    info = string("Max iter reached, Root: ", x2)
    println(info)
    s = string(s, "\n", info)
    if visualize
        make_plots(forward_abs_errors, backward_abs_errors, s, question)
    end
    return 1
end

function false_position(f, x0::Float64, x1::Float64, max_iter::Int,
                        threshold::Float64, question::String,
                        visualize::Bool=true)
    s = string("False Position Method x0=", x0, " x1=", x1)
    println(s)
    forward_abs_errors::Array{Float64} = Float64[] # plot 1, plot 3
    backward_abs_errors::Array{Float64} = Float64[] # plot 2
    if f(x0) * f(x1) >= 0.0
        info = string("Bad x0 & x1: ", x0, " & ", x1)
        println(info)
        s = string(s, "\n", info)
        if visualize
            make_plots(forward_abs_errors, backward_abs_errors, s, question)
        end
        return 1
    end
    x2 = x0
    for n in 0:max_iter
        fx0::Float64 = f(x0)
        fx1::Float64 = f(x1)
        denominator::Float64 = fx1 - fx0
        if denominator == 0.0
            info = string("Singularity, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)
            end
            return 1
        end
        x2::Float64 = x1 - fx1 * (x1 - x0) / denominator
        fx2::Float64 = f(x2)
        forward_abs_error = abs(f(x2))
        push!(forward_abs_errors, forward_abs_error)
        backward_abs_error = abs(x2 - 7.0)
        push!(backward_abs_errors, backward_abs_error)
        if forward_abs_error < threshold
            info = string(":) Converged, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)
            end
            return 0
        elseif f(x2) * f(x0) < 0
            x1 = x2
        else
            x0 = x2
        end
    end
    info = string("Max iter reached, Root: ", x2)
    println(info)
    s = string(s, "\n", info)
    if visualize
        make_plots(forward_abs_errors, backward_abs_errors, s, question)
    end
    return 1
end

thresh = .00001
num_iter = 9999
function f(x::Float64) :: Float64
    return (exp(x-7) - 1)
end

# Question 2
# 2.a
println()
secant_method(f, 1.0, 10.0, num_iter, thresh, "Question 2a")
false_position(f, 1.0, 10.0, num_iter, thresh, "Question 2a")
#= Output 2a
Secant Method x0=1.0 x1=10.0
Singularity, Root: 1.8712868878076563
False Position Method x0=1.0 x1=10.0
:) Converged, Root: 6.9999910450421465
=#
# 2.b
println()
secant_method(f, 6.0, 9.0, num_iter, thresh, "Question 2b")
false_position(f, 6.0, 9.0, num_iter, thresh, "Question 2b")
# 2.c
#= Output 2b
Secant Method x0=6.0 x1=9.0
:) Converged, Root: 6.999999068985341
False Position Method x0=6.0 x1=9.0
:) Converged, Root: 6.999992411394308
=#
println()
secant_method(f, 10.0, 1.0, num_iter, thresh, "Question 2c")
false_position(f, 10.0, 1.0, num_iter, thresh, "Question 2c")
#= Output 2c
Secant Method x0=10.0 x1=1.0
Singularity, Root: 1.4470280947235779
False Position Method x0=10.0 x1=1.0
:) Converged, Root: 6.999991045042148
=#
#=
For all cases in the log-log plot, we see a linear relationship that indicates convergence.
For the secant method the slope should be 1.618 around convergence and for false position it
should be 1. This mostly agrees with our results in the plots. Since we took the average
which isn't robust against outliers, we plotted a dy/dx vs n. plot to show that it converges
to the power or alpha.
=#

# Question 3
println()
secant_method(f, 4.0, 5.0, num_iter, thresh, "Question 3")
false_position(f, 9.0, 10.0, num_iter, thresh, "Question 3")
#= Output 3
Secant Method x0=4.0 x1=5.0
:) Converged, Root: 6.999999945795481
False Position Method x0=9.0 x1=10.0
Bad x0 & x1: 9.0 & 10.0
Explanation:
Although both cases involve an x0 & x1 that are close together, only the first case
converges, while the second case not only fails to converge but actually breaks during
its calculations, resulting in an empty plot for all types. When you look at the graph
of f(x), you can see that the slope if pretty flat from 0 and starts to increase rapidly
around x=8 & x=9. Since the algorithm requires taking a slope of the secant line, it is able to converge at
x0,x1 = 4,5 because there is enough difference in the secant when the function is more flat.
For x0,x1 = 9,10 the secant line is incredibly steep and following the algorithm becomes steeper
until the slope is basically dividing by 0, which causes it to break.
=#

# Question 4
evals = []
for i in -50.0:1:50.0
    for j in -50.0:1:50.0
        if secant_method(f, i, j, num_iter, thresh, "", false) == 0
            push!(evals, (i, j, 1.0::Float64))
        else
            push!(evals, (i, j, 0.0::Float64))
        end
    end
end
#println(evals)

x = [x[1] for x in evals]
y = [y[2] for y in evals]
z = [z[3] for z in evals]

using Plots
gr()

scatter(x,y,z)
