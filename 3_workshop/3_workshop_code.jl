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

function f(x)
    return x^2 - 612
end

root = secant_method(f, 10, 30, 5)

println("Root: ", root)
