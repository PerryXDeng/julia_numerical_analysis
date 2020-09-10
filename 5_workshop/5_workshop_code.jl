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
=#
# ========== Question 2 ===========
#= PART 1
golden ratio = 1.61803 =#

function q2_f1(x::Vector)
    return x[1]^2 - x[1] - 1
end

println("Initial Guess: 1.5")
newton_method(q2_f1, [1.5], 50)
println("Initial Guess: 0.0")
newton_method(q2_f1, [0.0], 50)
println("Initial Guess: 1.0")
newton_method(q2_f1, [1.0], 50)
println("Initial Guess: 0.5")
newton_method(q2_f1, [0.5], 50)
println("Initial Guess: 0.49")
newton_method(q2_f1, [0.49], 50)
println("Initial Guess: 0.51")
newton_method(q2_f1, [0.51], 50)
#= OUTPUT
Initial Guess: 1.5
step: 0, [1.5]
step: 1, [1.625]
step: 2, [1.6180555555555556]
step: 3, [1.6180339889579018]
step: 4, [1.618033988749895]
CONVERGENCE

Initial Guess: 0.0
step: 0, [0.0]
step: 1, [-1.0]
step: 2, [-0.6666666666666667]
step: 3, [-0.6190476190476191]
step: 4, [-0.6180344478216819]
step: 5, [-0.6180339887499892]
step: 6, [-0.6180339887498948]
CONVERGENCE

Initial Guess: 1.0
step: 0, [1.0]
step: 1, [2.0]
step: 2, [1.6666666666666667]
step: 3, [1.619047619047619]
step: 4, [1.618034447821682]
step: 5, [1.618033988749989]
step: 6, [1.618033988749895]
CONVERGENCE

Initial Guess: 0.5
step: 0, [0.5]
step: 1, [Inf]
step: 2, [NaN]
step: 3, [NaN]
step: 4, [NaN]
step: 5, [NaN]
step: 6, [NaN]
step: 7, [NaN]
step: 8, [NaN]
step: 9, [NaN]
step: 10, [NaN]
step: 11, [NaN]
step: 12, [NaN]
step: 13, [NaN]
step: 14, [NaN]
step: 15, [NaN]
step: 16, [NaN]
step: 17, [NaN]
step: 18, [NaN]
step: 19, [NaN]
step: 20, [NaN]
step: 21, [NaN]
step: 22, [NaN]
step: 23, [NaN]
step: 24, [NaN]
step: 25, [NaN]
step: 26, [NaN]
step: 27, [NaN]
step: 28, [NaN]
step: 29, [NaN]
step: 30, [NaN]
step: 31, [NaN]
step: 32, [NaN]
step: 33, [NaN]
step: 34, [NaN]
step: 35, [NaN]
step: 36, [NaN]
step: 37, [NaN]
step: 38, [NaN]
step: 39, [NaN]
step: 40, [NaN]
step: 41, [NaN]
step: 42, [NaN]
step: 43, [NaN]
step: 44, [NaN]
step: 45, [NaN]
step: 46, [NaN]
step: 47, [NaN]
step: 48, [NaN]
step: 49, [NaN]
step: 50, [NaN]
Initial Guess: 0.49
step: 0, [0.49]
step: 1, [-62.004999999999946]
step: 2, [-30.762499200063967]
step: 3, [-15.151241603742253]
step: 4, [-7.365553736000171]
step: 5, [-3.512237262141112]
step: 6, [-1.6618920709659004]
step: 7, [-0.8700446488659763]
step: 8, [-0.6412118365903691]
step: 9, [-0.618269358123613]
step: 10, [-0.618034013519751]
step: 11, [-0.6180339887498951]
step: 12, [-0.6180339887498948]
CONVERGENCE

Initial Guess: 0.51
step: 0, [0.51]
step: 1, [63.004999999999946]
step: 2, [31.76249920006397]
step: 3, [16.151241603742257]
step: 4, [8.365553736000173]
step: 5, [4.512237262141113]
step: 6, [2.6618920709659006]
step: 7, [1.870044648865976]
step: 8, [1.6412118365903692]
step: 9, [1.6182693581236132]
step: 10, [1.6180340135197508]
step: 11, [1.6180339887498951]
step: 12, [1.618033988749895]
CONVERGENCE

EXPLANATION: It did converge at the golden ratio! We also found another root
at -0.61803. After some other tries of initial guess between 0.0 and 1.5, we found
that the dividing line is 0.5. Less than 0.5 yields the negative root, and more
than 0.5 yields the positive root. =#

#= PART 2 =#
function
# ========== Question 3 ===========

# x_new function
# y_new function
