#====
Numerical Analysis
Workshop #3

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#

# Question 1

function secant_method(f, x0, x1, iter)
    x2 = 0
    for i in 0:iter
        x2 = x1 - f(x1) * (x1 - x0) / (float(f(x1) - f(x0)))
        x0, x1 = x1, x2
    end
    return x2
end

function false_postition(f, x0, x1, iter)
    if f(x0) * f(x1) >= 0
        println("Bad x0 & x1: ", x0, " & ", x1)
        return -1
    end
    x2 = x0
    for i in 0:iter
        x2 = x1 - f(x1) * (x1 - x0) / (float(f(x1) - f(x0)))
        if f(x2) == 0
            break
        elseif f(x2) * f(x0) < 0
            x1 = x2
        else
            x0 = x2
        end
    end
    return x2
end

function f(x)
    return x^2 - 612
end

root_sec = secant_method(f, 10, 30, 5)
root_fal = false_postition(f, 10, 30, 5)

println("Secant Method: ", root_sec)
println("False Position: ", root_fal)
