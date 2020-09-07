#=
Homework #2 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
using Plots
using Statistics


# ========== Functions ===========
function bisection(f, a, b, tolerance)
    if f(a) == 0 return a end
    if f(b) == 0 return b end
    if f(a) * f(b) > 0
        return "no root in range [a,b]:[" * a * "," * b * "]"
    end

    f_t = tolerance
    n = log2((b-a)/f_t)
    counter = 0
    while true
        dx = b - a
        x_m = a + (dx/2.0)
        if abs(f(x_m)) < f_t
            return string("x=", x_m, ", f(x)=", f(x_m), ", steps=", counter)
        end
        f(a) * f(x_m) <= 0 ? b = x_m : a = x_m
        counter += 1
        if counter >= (n + 2)
            return string("x=", x_m, ", f(x)=", f(x_m), ", n=", n, " iterations reached.")
        end
    end
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


FACC = 10.0^-10
MAXACC = 10.0^20
max_iter = 50

function fpi(func, xold)
    diff = 10*FACC
    step = 0

    while step < max_iter && abs(diff) > FACC && abs(diff) < MAXACC
        xnew = func(xold)
        diff = xnew-xold
        xold = xnew

        println("step=", step, "   xold,new= ", xold, ", ", xnew, "  diff=", diff)
        step += 1
    end
end



# ============= Question 1 =============
function f(x)
    return x^3 - 2
end

println(bisection(f, 0.0, 2.0, 10.0^-8))























# =========== Plotting Function =============

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
    ylabel!(p, "d ln|f(x_n)| / d ln|f(x_n_minus_one)|")
    xlabel!(p, "n")
    display(p)
end
