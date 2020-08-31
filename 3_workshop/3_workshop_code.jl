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
    x2 = 0
    for i in 0:iter
        x2 = x1 - f(x1) * (x1 - x0) / (float(f(x1) - f(x0)))
        x0, x1 = x1, x2
        if abs(f(x2)) < FACC
            break
        end
    end
    return x2
end

function false_postition(f, x0::Float64, x1::Float64, iter)
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
        if abs(f(x2)) < FACC
            break
        end
    end
    return x2
end


# Question 2
function f(x::Float64) :: Float64
    return (exp(x-7) - 1)
end
max_iter = 9999
# 2.a
secant_method(f, 1, 10, max_iter)
false_position(f, 1, 10, max_iter)
# 2.b
secant_method(f, 6, 9, max_iter)
false_position(f, 6, 9, max_iter)
# 2.c
secant_method(f, 10, 1, max_iter)
false_position(f, 10, 1, max_iter)

# Question 3
secant_method(f, 4, 5, max_iter)
false_position(f, 9, 10, max_iter)

# Question 4 ??? pen and paper?
