#====
Numerical Analysis
Workshop #2

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#

# Question 1
function bisection(f, a, b, tolerance)
    if f(a) == 0 return a end
    if f(b) == 0 return b end
    if f(a) * f(b) > 0
        return "no root in range [a,b]:[" * a * "," * b * "]"
    end

    f_t = tolerance
    n = log2((b-a)/f_t)
    counter = 0
    while true
        dx = b - a
        x_m = a + (dx/2.0)
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
    return -128 + 448*x - 672*x^2 + 560*x^3 - 280*x^4 + 84*x^5 - 14*x^6 + x^7
end

println("[1.5, 2.5]")
for i in 8:16
    println(string("tol = 10^-", i , ": val = ", bisection(f, 1.5, 2.5, 10.0^-i)))
end
println("[1.5, 3]")
for i in 8:16
    println(string("tol = 10^-", i , ": val = ", bisection(f, 1.5, 3, 10.0^-i)))
end

#=
Output:
[1.5, 2.5]
tol = 10^-8: val = 2.0
tol = 10^-9: val = 2.0
tol = 10^-10: val = 2.0
tol = 10^-11: val = 2.0
tol = 10^-12: val = 2.0
tol = 10^-13: val = 2.0
tol = 10^-14: val = 2.0
tol = 10^-15: val = 2.0
tol = 10^-16: val = 2.0
[1.5, 3]
tol = 10^-8: val = 2.0625
tol = 10^-9: val = 1.96875
tol = 10^-10: val = 1.96875
tol = 10^-11: val = 2.015625
tol = 10^-12: val = 2.015625
tol = 10^-13: val = 1.9921875
tol = 10^-14: val = 1.9921875
tol = 10^-15: val = 1.9921875
tol = 10^-16: val = 1.9921875
=#
