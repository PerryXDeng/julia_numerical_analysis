#====
Numerical Analysis
Workshop #2

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#

# Question 1 & 2
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
            return string("x=", x_m, ", f(x)=", f(x_m))
        end
        f(a) * f(x_m) <= 0 ? b = x_m : a = x_m
        counter += 1
        if counter >= (n + 2)
            return string("x=", x_m, ", f(x)=", f(x_m), ", n=", n, " iterations reached.")
        end
    end
end

function f(x)
    return -128 + 448*x - 672*x^2 + 560*x^3 - 280*x^4 + 84*x^5 - 14*x^6 + x^7
end

println("Question 2")
println("[1.5, 2.5]")
for i in 8:16
    println(string("tol = 10^-", i , ":", bisection(f, 1.5, 2.5, 10.0^-i)))
end
println("[1.5, 3]")
for i in 8:16
    println(string("tol = 10^-", i , ":", bisection(f, 1.5, 3, 10.0^-i)))
end

#=
Output:
[1.5, 2.5]
tol = 10^-8:x=2.0, f(x)=0.0
tol = 10^-9:x=2.0, f(x)=0.0
tol = 10^-10:x=2.0, f(x)=0.0
tol = 10^-11:x=2.0, f(x)=0.0
tol = 10^-12:x=2.0, f(x)=0.0
tol = 10^-13:x=2.0, f(x)=0.0
tol = 10^-14:x=2.0, f(x)=0.0
tol = 10^-15:x=2.0, f(x)=0.0
tol = 10^-16:x=2.0, f(x)=0.0
[1.5, 3]
tol = 10^-8:x=2.0625, f(x)=3.725290298461914e-9
tol = 10^-9:x=1.96875, f(x)=-2.9103830456733704e-11
tol = 10^-10:x=1.96875, f(x)=-2.9103830456733704e-11
tol = 10^-11:x=2.015625, f(x)=2.2737367544323206e-13
tol = 10^-12:x=2.015625, f(x)=2.2737367544323206e-13
tol = 10^-13:x=1.9921875, f(x)=0.0
tol = 10^-14:x=1.9921875, f(x)=0.0
tol = 10^-15:x=1.9921875, f(x)=0.0
tol = 10^-16:x=1.9921875, f(x)=0.0
Question 2 comments: smaller tolerance leads to increasingly accurate solutions
=#


# Question 3
# compulsive type hinting
using LinearAlgebra
function wilkinson_polynomials_expanded_eval(x::Float64) :: Float64
    poly_constants::Array{Float64} = [1.0, -210.0, 20615.0,-1256850.0,
          53327946.0,-1672280820.0, 40171771630.0, -756111184500.0,
          11310276995381.0, -135585182899530.0,
          1307535010540395.0,     -10142299865511450.0,
          63030812099294896.0,     -311333643161390640.0,
          1206647803780373360.0,     -3599979517947607200.0,
          8037811822645051776.0,      -12870931245150988800.0,
          13803759753640704000.0,      -8752948036761600000.0,
          2432902008176640000.0]
    poly_vars::Array{Float64} = x .^ reverse(0:20)
    return dot(poly_constants, poly_vars)
end
fn = wilkinson_polynomials_expanded_eval
i = 8
println("Question 3")
println("[1.9, 2.1]")
println(string("tol = 10^-", i , ": ", bisection(fn, 1.9, 2.1, 10.0^-i)))
println("[9.9, 10.1]")
println(string("tol = 10^-", i , ": ", bisection(fn, 9.9, 10.1, 10.0^-i)))
println("[15.9, 16.1]")
println(string("tol = 10^-", i , ": ", bisection(fn, 15.9, 16.1, 10.0^-i)))
#=
[1.9, 2.1]
tol = 10^-8: x=1.999999998509884, f(x)=-9.527296e6, n=24.25349666421154 iterations reached.
[9.9, 10.1]
tol = 10^-8: x=9.999821965396404, f(x)=-9.55604992e9, n=24.253496664211532 iterations reached.
[15.9, 16.1]
tol = 10^-8: x=15.990879048407082, f(x)=-7.70704719872e11, n=24.253496664211543 iterations reached.
Question 3 comments: It appears that both forward and backward errors become larger
on the later two trials compared to the first trial.
=#
using Winston
xs = range(1, stop=20, step=0.01)
ys = map(wilkinson_polynomials_expanded_eval, xs)
plot(xs, ys, "b")
#=
Question 3 comments: It appears that the function has steep tangents near
the first two and last two roots, but much smoother tangents otherwise.
This should cause much larger condition number
=#
