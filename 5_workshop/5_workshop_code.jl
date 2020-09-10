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
                       forward_err_threshold::Float64=1e-8)
    f_prime = (x_in -> ForwardDiff.gradient(f, x_in)[1])
    println("step:", 0, ", ", x)
    for i in 1:num_iter
        f_x = f(x)
        x = [x[1] - f_x/f_prime(x)] # conversion between scalar and vector gets complicated

        println("step:",i, ", ", x, " error:", abs(f_x))
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
step:0, [420.69]
step:1, [5.684341886080802e-14] error:2103.45
step:2, [0.0] error:2.8421709430404007e-13
CONVERGENCE

f(x) = x*5 + 1
step:0, [420.69]
1, [-0.19999999999998863] error:2104.45
step:2, [-0.2] error:5.684341886080802e-14
CONVERGENCE

f(x) = x^5 + 1
step:0, [420.69]
step:1, [336.55199999999365] error:1.3176830024935041e13
step:2, [269.2415999999793] error:4.31778366257098e12
step:3, [215.3932799999454] error:1.4148513505515208e12
step:4, [172.31462399986341] error:4.636184905489852e11
step:5, [137.85169919966387] error:1.5191850698335425e11
step:6, [110.28135935917726] error:4.9780656368568214e10
step:7, [88.22508748598968] error:1.6312125479115152e10
step:8, [70.58006998549061] error:5.345157277259175e9
step:9, [56.46405598033311] error:1.7515011368750057e9
step:10, [45.17124476459026] error:5.739318927939218e8
step:11, [36.13699576363456] error:1.880660028934323e8
step:12, [28.909596493628225] error:6.162546809083989e7
step:13, [23.12767690857586] error:2.019343364672642e7
step:14, [18.502140827820813] error:6.616984600079321e6
step:15, [14.801710955616077] error:2.1682537764740223e6
step:16, [11.841364597887972] error:710493.6601951023
step:17, [9.473081505920845] error:232814.82529301933
step:18, [7.578440369694619] error:76289.02467289627
step:19, [6.062691662534009] error:24998.650327499203
step:20, [4.850005293778907] error:8191.820467507624
step:21, [3.8796427742733397] error:2684.558475795735
step:22, [3.102831417964252] error:879.938917658215
step:23, [2.4801074035032844] error:288.6013374802062
step:24, [1.9787996686353164] error:94.83231774361398
step:25, [1.5699953777982345] error:31.339550687638173
step:26, [1.223078070336621] error:10.538758840230674
step:27, [0.8890880276925446] error:3.736975418614904
step:28, [0.39119601698412665] error:1.5555508427614377
step:29, [-8.226933613726082] error:1.009161616775151
step:30, [-6.581590550407136] error:37685.86085236034
step:31, [-5.265379027893607] error:12348.640169535725
step:32, [-4.2125634247990424] error:4046.1397073372827
step:33, [-3.3706858428693875] error:1325.5763899070469
step:34, [-2.6980980505215157] error:434.10230585620286
step:35, [-2.1622524159745624] error:141.98439454738553
step:36, [-1.7389515817290095] error:46.2641612359542
step:37, [-1.4130328392552782] error:14.901476426894032
step:38, [-1.180593599978475] error:4.633279185531521
step:39, [-1.0474253413848853] error:1.293517841299023
step:40, [-1.0041045598332308] error:0.26071054379612724
step:41, [-1.0000334201918948] error:0.02069196621282332
step:42, [-1.0000000022336692] error:0.000167112128939495
step:43, [-1.0] error:1.1168346247814043e-8
step:44, [-1.0] error:0.0
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
step:0, [1.5]
step:1, [1.625] error:0.25
step:2, [1.6180555555555556] error:0.015625
step:3, [1.6180339889579018] error:4.822530864223573e-5
step:4, [1.618033988749895] error:4.651177221148828e-10
CONVERGENCE

Initial Guess: 0.0
step:0, [0.0]
step:1, [-1.0] error:1.0
step:2, [-0.6666666666666667] error:1.0
step:3, [-0.6190476190476191] error:0.11111111111111116
step:4, [-0.6180344478216819] error:0.0022675736961450532
step:5, [-0.6180339887499892] error:1.0265159331446227e-6
step:6, [-0.6180339887498948] error:2.1094237467877974e-13
CONVERGENCE

Initial Guess: 1.0
step:0, [1.0]
step:1, [2.0] error:1.0
step:2, [1.6666666666666667] error:1.0
step:3, [1.619047619047619] error:0.11111111111111138
step:4, [1.618034447821682] error:0.0022675736961450532
step:5, [1.618033988749989] error:1.0265159333666674e-6
step:6, [1.618033988749895] error:2.1049828546892968e-13
CONVERGENCE

Initial Guess: 0.5
step:0, [0.5]
step:1, [Inf] error:1.25
step:2, [NaN] error:NaN
step:3, [NaN] error:NaN
step:4, [NaN] error:NaN
step:5, [NaN] error:NaN
step:6, [NaN] error:NaN
step:7, [NaN] error:NaN
step:8, [NaN] error:NaN
step:9, [NaN] error:NaN
step:10, [NaN] error:NaN
step:11, [NaN] error:NaN
step:12, [NaN] error:NaN
step:13, [NaN] error:NaN
step:14, [NaN] error:NaN
step:15, [NaN] error:NaN
step:16, [NaN] error:NaN
step:17, [NaN] error:NaN
step:18, [NaN] error:NaN
step:19, [NaN] error:NaN
step:20, [NaN] error:NaN
step:21, [NaN] error:NaN
step:22, [NaN] error:NaN
step:23, [NaN] error:NaN
step:24, [NaN] error:NaN
step:25, [NaN] error:NaN
step:26, [NaN] error:NaN
step:27, [NaN] error:NaN
step:28, [NaN] error:NaN
step:29, [NaN] error:NaN
step:30, [NaN] error:NaN
step:31, [NaN] error:NaN
step:32, [NaN] error:NaN
step:33, [NaN] error:NaN
step:34, [NaN] error:NaN
step:35, [NaN] error:NaN
step:36, [NaN] error:NaN
step:37, [NaN] error:NaN
step:38, [NaN] error:NaN
step:39, [NaN] error:NaN
step:40, [NaN] error:NaN
step:41, [NaN] error:NaN
step:42, [NaN] error:NaN
step:43, [NaN] error:NaN
step:44, [NaN] error:NaN
step:45, [NaN] error:NaN
step:46, [NaN] error:NaN
step:47, [NaN] error:NaN
step:48, [NaN] error:NaN
step:49, [NaN] error:NaN
step:50, [NaN] error:NaN

Initial Guess: 0.49
step:0, [0.49]
step:1, [-62.004999999999946] error:1.2499
step:2, [-30.762499200063967] error:3905.6250249999935
step:3, [-15.151241603742253] error:976.0938562340002
step:4, [-7.365553736000171] error:243.71136373871238
step:5, [-3.512237262141112] error:60.616935573906254
step:6, [-1.6618920709659004] error:14.848047847713609
step:7, [-0.8700446488659763] error:3.4237773265052294
step:8, [-0.6412118365903691] error:0.6270223398862964
step:9, [-0.618269358123613] error:0.05236445597396333
step:10, [-0.618034013519751] error:0.0005263573181975012
step:11, [-0.6180339887498951] error:5.538708269803294e-8
step:12, [-0.6180339887498948] error:6.661338147750939e-16
CONVERGENCE

Initial Guess: 0.51
step:0, [0.51]
step:1, [63.004999999999946] error:1.2499
step:2, [31.76249920006397] error:3905.625024999993
step:3, [16.151241603742257] error:976.0938562340004
step:4, [8.365553736000173] error:243.7113637387125
step:5, [4.512237262141113] error:60.61693557390628
step:6, [2.6618920709659006] error:14.848047847713616
step:7, [1.870044648865976] error:3.423777326505231
step:8, [1.6412118365903692] error:0.6270223398862957
step:9, [1.6182693581236132] error:0.05236445597396333
step:10, [1.6180340135197508] error:0.0005263573181979453
step:11, [1.6180339887498951] error:5.5387082475988336e-8
step:12, [1.618033988749895] error:6.661338147750939e-16
CONVERGENCE

EXPLANATION: It did converge at the golden ratio! We also found another root
at -0.61803. After some other tries of initial guess between 0.0 and 1.5, we found
that the dividing line is 0.5. Less than 0.5 yields the negative root, and more
than 0.5 yields the positive root. =#

#= PART 2 =#
n = 7 # value between 3 & 8
function q2_f2(x::Vector)
    return (x[1]-1)^n
end

println("n = 7")
newton_method(q2_f2, [1.5], 50)

#= OUTPUT
n = 5
step:0, [1.5]
step:1, [1.4] error:0.03125
step:2, [1.3199999999999998] error:0.010239999999999989
step:3, [1.2559999999999998] error:0.0033554431999999914
step:4, [1.2047999999999999] error:0.0010995116277759953
step:5, [1.16384] error:0.0003602879701896385
step:6, [1.131072] error:0.00011805916207174108
step:7, [1.1048576] error:3.8685626227668245e-5
step:8, [1.08388608] error:1.2676506002282358e-5
step:9, [1.0671088640000002] error:4.153837486827883e-6
step:10, [1.0536870912] error:1.3611294676837697e-6
step:11, [1.04294967296] error:4.46014903970614e-7
step:12, [1.034359738368] error:1.4615016373309154e-7
step:13, [1.0274877906944] error:4.789048565205975e-8
step:14, [1.02199023255552] error:1.5692754338466684e-8
step:15, [1.017592186044416] error:5.14220174162866e-9

n = 7
step:0, [1.5]
step:1, [1.4] error:0.03125
step:2, [1.3199999999999998] error:0.010239999999999989
step:3, [1.2559999999999998] error:0.0033554431999999914
step:4, [1.2047999999999999] error:0.0010995116277759953
step:5, [1.16384] error:0.0003602879701896385
step:6, [1.131072] error:0.00011805916207174108
step:7, [1.1048576] error:3.8685626227668245e-5
step:8, [1.08388608] error:1.2676506002282358e-5
step:9, [1.0671088640000002] error:4.153837486827883e-6
step:10, [1.0536870912] error:1.3611294676837697e-6
step:11, [1.04294967296] error:4.46014903970614e-7
step:12, [1.034359738368] error:1.4615016373309154e-7
step:13, [1.0274877906944] error:4.789048565205975e-8
step:14, [1.02199023255552] error:1.5692754338466684e-8
step:15, [1.017592186044416] error:5.14220174162866e-9

EXPLANATION
For both initial guesses of 5 and 7, we can see that Newton's Method ends up bouncing
back and forth and fails to quadratically converge.
=#


# ========== Question 3 ===========

# x_new function
# y_new function
