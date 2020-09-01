#====
Numerical Analysis
Workshop #3

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#

using Plots

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
    p = plot(log_errors_n_minus_one, log_errors_n,
         title=title, legend=false)
    ylabel!(p, "ln|f(x_n)|")
    xlabel!(p, "ln|f(x_n_minus_one)|")
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
# 2.b
println()
secant_method(f, 6.0, 9.0, num_iter, thresh, "Question 2b")
false_position(f, 6.0, 9.0, num_iter, thresh, "Question 2b")
# 2.c
println()
secant_method(f, 10.0, 1.0, num_iter, thresh, "Question 2c")
false_position(f, 10.0, 1.0, num_iter, thresh, "Question 2c")

# Question 3
println()
secant_method(f, 4.0, 5.0, num_iter, thresh, "Question 3")
false_position(f, 9.0, 10.0, num_iter, thresh, "Question 3")

# Question 4
evals = []
for i in -50.0:1:50.0
    for j in -50.0:1:50.0
        if secant_method(f, i, j, max_iter, thresh, "", false) == 0
            push!(evals, (i, j, 1.0::Float64))
        else
            push!(evals, (i, j, 0.0::Float64))
        end
    end
end
println(evals)

x = [x[1] for x in evals]
y = [y[2] for y in evals]
z = [z[3] for z in evals]

using Plots
gr()

scatter(x,y,z)
