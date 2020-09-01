#====
Numerical Analysis
Workshop #3

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#

# Question 1
FACC = .00001
function secant_method(f, x0::Float64, x1::Float64, iter)
    println("Secant Method x0=", x0, " x1=", x1)
    x2::Float64 = 0.0
    for i in 0:iter
        fx1::Float64 = f(x1)
        denominator::Float64 = fx1 - f(x0)
        if denominator == 0.0
            println("Singularity, Root: ", x2)
            return 1
        end
        x2 = x1 - fx1 * (x1 - x0) / denominator
        x0, x1 = x1, x2
        if abs(f(x2)) < FACC
            println("Converged, Root: ", x2)
            return 0
        end
    end
    println("Max Iteration, Root: ", x2)
    return 1
end

function false_position(f, x0::Float64, x1::Float64, iter)
    println("False Position Method x0=", x0, " x1=", x1)
    if f(x0) * f(x1) >= 0.0
        println("Bad x0 & x1: ", x0, " & ", x1)
        return 1
    end
    x2 = x0
    for i in 0:iter
        fx0::Float64 = f(x0)
        fx1::Float64 = f(x1)
        denominator::Float64 = fx1 - fx0
        if denominator == 0.0
            println("Singularity, Root: ", x2)
            return 1
        end
        x2::Float64 = x1 - fx1 * (x1 - x0) / denominator
        fx2::Float64 = f(x2)
        if abs(fx2) < FACC
            println("Converged, Root: ", x2)
            return 0
        elseif f(x2) * f(x0) < 0
            x1 = x2
        else
            x0 = x2
        end
    end
    println("Max Iteration, Root: ", x2)
    return 1
end


# Question 2
function f(x::Float64) :: Float64
    return (exp(x-7) - 1)
end
max_iter = 9999
# 2.a
secant_method(f, 1.0, 10.0, max_iter)
false_position(f, 1.0, 10.0, max_iter)
# 2.b
secant_method(f, 6.0, 9.0, max_iter)
false_position(f, 6.0, 9.0, max_iter)
# 2.c
secant_method(f, 10.0, 1.0, max_iter)
false_position(f, 10.0, 1.0, max_iter)

# Question 3
secant_method(f, 4.0, 5.0, max_iter)
false_position(f, 9.0, 10.0, max_iter)

# Question 4
evals = []
for i in -50.0:1:50.0
    for j in -50.0:1:50.0
        if secant_method(f, i, j, max_iter) == 0
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
