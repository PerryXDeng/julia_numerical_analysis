#=
Workshop #6 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#

# Question 1
# TODO: Perry - not tested on anything yet, maybe broken
function steffensen(g, x0::Float64;
                    max_iter::Int64=50,
                    x_convergence_threshold::Float64=1e-8)
    loop_iters = max_iter÷3
    for i in 0:loop_iters
        x1 = g(x0)
        println("i=", 3i + 1, " x1:", x1)
        if abs(x1 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 1)
            x0=x1
            break
        end
        if 3i+1 == max_iter
            break
        end

        x2 = g(x1)
        println("i=", 3i + 2, " x2:", x2)
        if abs(x2 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 2)
            x0=x2
            break
        end
        if 3i+2 == max_iter
            break
        end

        x3 = x0 - ((x1 - x0)^2) / (x2 - 2x1 + x0)
        println("i=", 3i + 3, " x0:", x0," x1:", x1, " x2:", x2, " x3:", x3)
        if abs(x3 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 3)
            x0=x3
            break
        end
        if 3i+3 == max_iter
            break
        end
        x0 = x3
    end
    println("x=",x0)
end

function plain_fpi(g, x0::Float64;
                   max_iter::Int64=50,
                   x_convergence_threshold::Float64=1e-8)
    for i in 1:max_iter
        x1 = g(x0)
        if abs(x1 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", i)
            x0=x1
            break
        end
        x0=x1
    end
    println("x=",x0)
end



# Question 3
println("\ng(x)=cos(x); x0 = 1; expected root approx 0.7")
g(x) = cos(x)
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
Improves convergence
Output:
g(x)=cos(x); x0 = 1; expected root approx 0.7
FPI
CONVERGENCE at i=46
x=0.7390851366465718
STEFFENSEN
CONVERGENCE at i=10
x=0.7390851332482249
=#
println("\ng(x)=x/2 + 1/x^2; x0 = 1; expected root approx 1.26")
g(x) = x/2 + 1/x^2
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
Improves convergence
Output:
g(x)=x/2 + 1/x^2; x0 = 1; expected root approx 1.26
FPI
CONVERGENCE at i=27
x=1.2599210520301471
STEFFENSEN
CONVERGENCE at i=13
x=1.2599210498948539
=#
println("\ng(x)=x-0.5((x^2)-x-1); x0 = 1; expected root 1.618")
g(x) = x-0.5((x^2)-x-1)
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
No improvement, possibly due to quadraticly converging FPI on g(x)
Output:
g(x)=x-0.5((x^2)-x-1); x0 = 1; expected root 1.618
FPI
CONVERGENCE at i=10
x=1.6180339890192184
STEFFENSEN
CONVERGENCE at i=10
x=1.618033988885942
=#
println("\ng(x)=(-2.5x+((x^2)-1))/-1.5; x0 = 1; expected root 1.618")
g(x) = (-2.5x+((x^2)-1))/(-1.5)
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
Improvement
Output:
g(x)=(-2.5x+((x^2)-1))/-1.5; x0 = 1; expected root 1.618
FPI
CONVERGENCE at i=25
x=1.6180339906726358
STEFFENSEN
CONVERGENCE at i=10
x=1.6180339892762732
=#

# Question 4
# part 1
function false_position(f, x0::Float64, x1::Float64, max_iter::Int64=50, threshold::Float64=1e-8)
    if f(x0) * f(x1) >= 0.0
        info = string("Bad x0 & x1: ", x0, " & ", x1)
        println(info)
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
            return 1
        end
        x2::Float64 = x1 - fx1 * (x1 - x0) / denominator
        fx2::Float64 = f(x2)
        if abs(f(x2)) < threshold
            info = string(":) Converged, Root: ", x2, " n:", n)
            println(info)
            return 0
        elseif f(x2) * f(x0) < 0
            x1 = x2
        else
            x0 = x2
        end
    end
    info = string("Max iter reached, Root: ", x2)
    println(info)
    return 1
end

function false_aitkens(g, x0::Float64, x1::Float64, max_iter::Int64=50, threshold::Float64=1e-8)
    if g(x0) * g(x1) >= 0.0
        info = string("Bad x0 & x1: ", x0, " & ", x1)
        println(info)
        s = string(s, "\n", info)
        return 1
    end
    loop_iters = max_iter÷3
    x2 = x0
    for i in 0:loop_iters
        gx0::Float64 = g(x0)
        gx1::Float64 = g(x1)
        denominator::Float64 = gx1 - gx0
        if denominator == 0.0
            info = string("Singularity, Root: ", x2)
            println(info)
            return 1
        end

        x2::Float64 = x1 - gx1 * (x1 - x0) / denominator
        gx2::Float64 = g(x2)
        println("i=", 3i + 2, " x0:", x0," x1:", x1, " x2:", x2, " error:", abs(g(x2)))
        if abs(g(x2)) < threshold
            info = string(":) Converged, Root: ", x2, " error:", abs(g(x2)))
            println(info)
            x0=x2
            break
        end
        if 3i+2 == max_iter
            break
        end

        x3 = x0 - ((x1 - x0)^2) / (x2 - 2x1 + x0)
        println("i=", 3i + 3, " x0:", x0," x1:", x1, " x2:", x2, " x3:", x3, " error:", abs(x3 - x0))
        if abs(x3 - x0) < threshold
            println("CONVERGENCE at i=", 3i + 3, " error:", abs(x3 - x0))
            x0=x3
            break
        elseif g(x3) * g(x1) < 0
            x0 = x3
        else
            x1 = x3
            x0 = x2
        end
        if 3i+3 == max_iter
            break
        end
    end
    println("x=",x0)
end
println("\ng(x)=x^3-2; initial guess x0 = 0, x1 = 2; expected root approx 1.26")
g(x) = x^3 - 2
x0 = 0.0
x1 = 2.0
false_aitkens(g, x0, x1)
false_position(g, x0, x1)
#=
Output:
False Position with Aitkens
i=2 x0:0.0 x1:2.0 x2:0.5 error:1.875
i=3 x0:0.0 x1:2.0 x2:0.5 x3:1.1428571428571428 error:1.1428571428571428
i=5 x0:1.1428571428571428 x1:2.0 x2:1.2096774193548385 error:0.22985549326978028
i=6 x0:1.1428571428571428 x1:2.0 x2:1.2096774193548385 x3:1.5888111888111887 error:0.44595404595404586
i=8 x0:1.2096774193548385 x1:1.5888111888111887 x2:1.2485727596125022 error:0.053557553416991865
i=9 x0:1.2096774193548385 x1:1.5888111888111887 x2:1.2485727596125022 x3:1.4094938918442212 error:0.19981647248938272
i=11 x0:1.2485727596125022 x1:1.4094938918442212 x2:1.2586675537891296 error:0.005963466117995253
i=12 x0:1.2485727596125022 x1:1.4094938918442212 x2:1.2586675537891296 x3:1.3316387449592773 error:0.083065985346775
i=14 x0:1.2586675537891296 x1:1.3316387449592773 x2:1.2598522864999055 error:0.00032744738460155354
i=15 x0:1.2586675537891296 x1:1.3316387449592773 x2:1.2598522864999055 x3:1.2954517565888992 error:0.03678420279976957
i=17 x0:1.2598522864999055 x1:1.2954517565888992 x2:1.2599191465968447 error:9.063878185777696e-6
i=18 x0:1.2598522864999055 x1:1.2954517565888992 x2:1.2599191465968447 x3:1.2776687522798063 error:0.017816465779900703
i=20 x0:1.2599191465968447 x1:1.2776687522798063 x2:1.2599210233343145 error:1.2648677349957893e-7
i=21 x0:1.2599191465968447 x1:1.2776687522798063 x2:1.2599210233343145 x3:1.2687944186474986 error:0.00887527205065397
i=23 x0:1.2599210233343145 x1:1.2687944186474986 x2:1.2599210497086877 error:8.866529732642903e-10
:) Converged, Root: 1.2599210497086877 error:8.866529732642903e-10
x=1.2599210497086877

False Position
:) Converged, Root: 1.2599210480979497 n:23
=#

# part 2
println("\ng(x)=2x/3 + 2/3x^2; initial guess x0 = 1; expected root approx 1.26")
g(x) = (2x)/3 + 2/(3x^2)
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
Not an improvement.
Output:
g(x)=2x/3 + 2/3x^2; initial guess x0 = 1; expected root approx 1.26
FPI
CONVERGENCE at i=5
x=1.2599210498948732
STEFFENSEN
CONVERGENCE at i=10
x=1.259921049894873
=#
