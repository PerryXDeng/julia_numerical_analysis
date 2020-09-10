#=
Workshop #5 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
using Plots
using Statistics

# ========== Question 1 ===========
using ForwardDiff # for efficient differentiation with automatic chainrule
#= implement Newton code (univariate) =#
function newton_method(f, x::Vector, num_iter::Int64;
                       forward_err_threshold::Float64=10^-8)
    f_prime = (x_in -> ForwardDiff.gradient(f, x_in)[1])
    println(0, ", ", x)
    for i in 1:num_iter
        f_x = f(x)
        x = [x[1] - f_x/f_prime(x)] # conversion between scalar and vector gets complicated
        println(i, ", ", x)
        if abs(f_x) < forward_err_threshold
            println("CONVERGENCE")
            break
        end
    end
    return x
end

function linear_func(x::Vector)
    return sum(x*5)
end
function affine_func(x::Vector)
    return sum(x*5 .+ 1)
end
function poly_func(x::Vector)
    return sum(x^5 .+ 1)
end
println("f(x) = x*5")
root = newton_method(linear_func, [420.69], 50)
println("f(x) = x*5 + 1")
root = newton_method(affine_func, [420.69], 50)
println("f(x) = x^5 + 1")
root = newton_method(poly_func, [420.69], 50)

# ========== Question 3 ===========

# x_new function
# y_new function
