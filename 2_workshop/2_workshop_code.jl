#====
Numerical Analysis
Workshop #2

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#

# Question 1
function bisection(f, a, b)
    if f(a) == 0 return a end
    if f(b) == 0 return b end
    if f(a) * f(b) > 0
        return "no root in range [a,b]:[" * a * "," * b * "]"
    end

    f_t = .001
    n = log2((b-a)/f_t)
    counter = 0
    while true
        dx = b - a
        x_m = a + (dx/2)
        if abs(f(x_m)) < f_t
            return x_m
        end
        f(a) * f(x_m) <= 0 ? b = x_m : a = x_m
        counter += 1
        if counter >= (n + 4)
            return string("n=", n, " iterations reached.")
        end
    end
end

function f(x)
    return x^2 - 10
end

println(bisection(f, 0, 10))
