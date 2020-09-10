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
    println("step:", 0, ", ", x)
    for i in 1:num_iter
        f_x = f(x)
        x = [x[1] - f_x/f_prime(x)] # conversion between scalar and vector gets complicated
        println("step:",i, ", ", x)
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
    return sum(x[1]^5 .+ 1)
end
println("f(x) = x*5")
root = newton_method(linear_func, [420.69], 50)
println("f(x) = x*5 + 1")
root = newton_method(affine_func, [420.69], 50)
println("f(x) = x^5 + 1")
root = newton_method(poly_func, [420.69], 50)
#= OUTPUT
f(x) = x*5
step: 0, [420.69]
step: 1, [5.684341886080802e-14]
step: 2, [0.0]
CONVERGENCE

f(x) = x*5 + 1
step: 0, [420.69]
step: 1, [-0.19999999999998863]
step: 2, [-0.2]
CONVERGENCE

f(x) = x^5 + 1
step: 0, [420.69]
step: 1, [336.55199999999365]
step: 2, [269.2415999999793]
step: 3, [215.3932799999454]
step: 4, [172.31462399986341]
step: 5, [137.85169919966387]
step: 6, [110.28135935917726]
step: 7, [88.22508748598968]
step: 8, [70.58006998549061]
step: 9, [56.46405598033311]
step: 10, [45.17124476459026]
step: 11, [36.13699576363456]
step: 12, [28.909596493628225]
step: 13, [23.12767690857586]
step: 14, [18.502140827820813]
step: 15, [14.801710955616077]
step: 16, [11.841364597887972]
step: 17, [9.473081505920845]
step: 18, [7.578440369694619]
step: 19, [6.062691662534009]
step: 20, [4.850005293778907]
step: 21, [3.8796427742733397]
step: 22, [3.102831417964252]
step: 23, [2.4801074035032844]
step: 24, [1.9787996686353164]
step: 25, [1.5699953777982345]
step: 26, [1.223078070336621]
step: 27, [0.8890880276925446]
step: 28, [0.39119601698412665]
step: 29, [-8.226933613726082]
step: 30, [-6.581590550407136]
step: 31, [-5.265379027893607]
step: 32, [-4.2125634247990424]
step: 33, [-3.3706858428693875]
step: 34, [-2.6980980505215157]
step: 35, [-2.1622524159745624]
step: 36, [-1.7389515817290095]
step: 37, [-1.4130328392552782]
step: 38, [-1.180593599978475]
step: 39, [-1.0474253413848853]
step: 40, [-1.0041045598332308]
step: 41, [-1.0000334201918948]
step: 42, [-1.0000000022336692]
step: 43, [-1.0]
step: 44, [-1.0]
CONVERGENCE
# ========== Question 3 ===========

# x_new function
# y_new function
