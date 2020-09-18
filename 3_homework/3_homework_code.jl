#=
Homework #3 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#

# Question 1
using ForwardDiff # for efficient differentiation with automatic chainrule
# TODO: modified_newton implementation, how it might work for other functions, last priority
function newton_method(f, x::Vector, num_iter::Int64;
                       forward_err_threshold::Float64=1e-8)
    f_prime = (x_in -> ForwardDiff.gradient(f, x_in)[1])
    #f_2prime = (x_in -> ForwardDiff.gradient(f_prime, x_in)[1])
    f_x = f(x)
    err = abs(f_x)
    println("step:", 0, ", ", x, " error:", err)

    f_prime_x = f_prime(x)
    if iszero(f_prime_x)
        println("SINGULARITY")
        break
    end
    if isinf(x[1])
        println("SINGULARITY")
        break
    end
    if isinf(err) || isnan(err)
        println("DIVERGENCE/ERROR")
        break
    end
    if err < forward_err_threshold
        println("CONVERGENCE")
        break
    end
    for i in 1:num_iter
        x = [x[1] - f_x/f_prime_x] # conversion between scalar and vector gets complicated
        f_x = f(x)
        err = abs(f_x)
        f_prime_x = f_prime(x)
        println("step:",i, ", ", x, " error:", err)
        if iszero(f_prime_x)
            println("SINGULARITY")
            break
        end
        if isinf(x[1])
            println("SINGULARITY")
            break
        end
        if isinf(err) || isnan(err)
            println("DIVERGENCE/ERROR")
            break
        end
        if err < forward_err_threshold
            println("CONVERGENCE")
            break
        end
    end
    return x
end

function modified_newton()
end


# Question 2
using Plots
using Statistics

function make_loglog_error_plots(backward_abs_errors::Array{Float64},
                                 meta_info::String, question::String)
    title = string(question, "\n", meta_info)
    num_iterations = size(backward_abs_errors, 1)

    # log log plot slope correspond to the order of the error convergence
    eps = 1e-15 # prevent breaking when error is 0
    log_backward_errors = log.(backward_abs_errors .+ eps)
    log_errors_n = log_backward_errors[2:num_iterations]
    log_errors_n_minus_one = log_backward_errors[1:(num_iterations-1)]

    # slope for the log log plot over iterations
    delta_y = log_errors_n[2:end] .- log_errors_n[1:end-1]
    delta_x = log_errors_n_minus_one[2:end] .- log_errors_n_minus_one[1:end-1]
    dy_dx = delta_y ./ delta_x
    mean_slope = mean(dy_dx)

    p = plot(log_errors_n_minus_one, log_errors_n,
             title=string(title, ", mean dy/dx=", mean_slope), legend=false,
             lw=2, marker = ([:hex :d], 3, 0.8))
    # epsilon_n = |x_n - true_root|
    ylabel!(p, "ln|epsilon_n|")
    xlabel!(p, "ln|epsilon_n_minus_one|")
    display(p)

    if size(dy_dx, 1) > 0
        p = plot(2:num_iterations-1, dy_dx,
             title=title, legend=false, lw=2, marker = ([:hex :d], 3, 0.8))
        ylabel!(p, "d ln|epsilon_n| / d ln|epsilon_n_minus_one|")
        xlabel!(p, "n")
        display(p)
    end
end

function halley(f, x::Vector;
                r::Float64=0,
                num_iter::Int64=50,
                forward_err_threshold::Float64=1e-8,
                question::String)
    s = string("Halley's Method x0=", x[1])
    println(s)
    f_prime = (x_in -> ForwardDiff.gradient(f, x_in)[1])
    f_primeprime = (x_in -> ForwardDiff.gradient(f_prime, x_in)[1])

    backward_errs::Array{Float64} = Float64[]

    f_x = f(x)
    forward_err = abs(f_x)
    # backward_err = abs(x[1] - r)
    # push!(backward_errs, backward_err)
    println("step:",0, ", ", x, " forward_error:", forward_err)

    f_prime_x = f_prime(x)
    f_primeprime_x = f_primeprime(x)
    numerator = 2*f_x*f_prime_x
    denominator = 2*(f_prime_x^2) - f_x * f_primeprime_x
    if iszero(denominator)
        println("SINGULARITY")
        return
    end
    if isinf(forward_err) || isnan(forward_err)
        println("DIVERGENCE/ERROR")
        return
    end
    if forward_err < forward_err_threshold
        println("CONVERGENCE")
        return
    end
    for i in 1:num_iter
        x = [x[1] - numerator/denominator] # conversion between scalar and vector gets complicated
        f_x = f(x)
        forward_err = abs(f_x)
        backward_err = abs(x[1] - r)
        push!(backward_errs, backward_err)
        println("step:",i, ", ", x, " forward_error:", forward_err, " backward_error:", backward_err)

        f_prime_x = f_prime(x)
        f_primeprime_x = f_primeprime(x)
        numerator = 2*f_x*f_prime_x
        denominator = 2*(f_prime_x^2) - f_x * f_primeprime_x
        if iszero(denominator)
            println("SINGULARITY")
            break
        end
        if isinf(x[1])
            println("SINGULARITY")
            break
        end
        if isinf(forward_err) || isnan(forward_err)
            println("DIVERGENCE/ERROR")
            break
        end
        if forward_err < forward_err_threshold
            println("CONVERGENCE")
            break
        end
    end
    println(backward_errs)
    make_loglog_error_plots(backward_errs, s, question)
end

# cubic convergence right away
function f(x::Vector)
    return x[1]^3 - 2
end
halley(f, [2.0], r=2^(1/3), question="f(x)=x^3 - 2")

# cubic convergence right away
function f(x::Vector)
    return x[1]^2 - x[1] - 1
end
halley(f, [3.0], r=((1+sqrt(5))/2), question="f(x)=x^2 - x - 1")


# Question 3
# working on this one rn -clarissa
# x vector is in form: [V, W] for 2-D Newton
function generalized_newton_method(f, x::Vector, num_iter::Int64=50;
                       forward_err_threshold::Float64=1e-8)
    # f should be a vector valued function (returns multiple values)
    # this allows efficient generalization to finding roots of many scalar functions
    jacobian = (x_in -> ForwardDiff.jacobian(f, x_in))
    f_x = f(x)
    err_vec = abs.(f_x)
    # for the other error question
    #err_vec = x_in -> sqrt(((e^(x[1]+x[2])-e^7)^2)+((ln(3*x[1]-2*x[2])-ln(11))^2))
    println("step:", 0, ", ", x, " error:", err_vec)
    for i in 1:num_iter
        j = jacobian(x)
        #inv_j = inv(j)
        inv_j = C_NULL
        try
            inv_j = inv(j)
        catch e
            println("SINGULARITY")
            break
        end
        x = x - inv_j * f_x
        f_x = f(x)
        err_vec = abs.(f_x)
        println("step:",i, ", ", x, " error:", err_vec)
        if count(x_in->isinf(x_in), x) > 0
            println("SINGULARITY")
            break
        end
        if count(err->(isnan(err)||isinf(err)), err_vec) > 0
            println("DIVERGENCE")
            break
        end
        if count(err->(err>forward_err_threshold), err_vec) == 0
            println("CONVERGENCE")
            break
        end
    end
    return x
end

function eq1(V::Float64, W::Float64) ::Float64
    return (V * (0.04 + W)^2) - 0.0011111 * ((0.04 + 2W)/(W^2)) - 0.083333
end

function eq2(V::Float64, W::Float64) ::Float64
    return ((-0.04/2)(1 + V)) + ((0.0011111/2)*(W^(-2))) - W + (1/2)*((1-V)W - 1.29074*(sqrt(1-V))) + 1.774598
end

function V(W::Float64) ::Float64
    return (0.083333 + (0.0011111)*(((0.04) + 2W)/(W^2)))/((0.04 + W)^2)
end

function combined_eq(W::Float64) ::Float64
    V = V(W)
    return ((-0.04/2)(1 + V)) + ((0.0011111/2)*(W^(-2))) - W + (1/2)*((1-V)W - 1.29074*(sqrt(1-V))) + 1.774598
end

# 2D Newton-Raphson Solver
generalized_newton_method()

# 1D Newton-Raphson Solver
generalized_newton_method(combined_eq, [1.0])

# Question 4
