#=
Project #1
Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
using Plots

@enum SOLVE_TYPE begin
    bisect_type
    fpi_type
    newton_type
end



# =================== Part 1 ===================
#=
INPUTS:
M = number of different binding molecules
ξ = total concentration of the ligand
ks = equilibrium constants
ns = total concentrations fo the various binding molecules
UB: upperbound for plot
LB: lowerbound for plot, default 0
OUTPUT:
x = concentration of unbound lignand in solution
=#
function concentrate_unbound_lignand(M, ξ, ks, ns)
    return function f(x)
        sum = 0
        for j in 1:M
            if isa(x, Vector)
                sum += (ks[j] * ns[j]) / (1 + (ks[j] * x[1]))
            else
                sum += (ks[j] * ns[j]) / (1 + (ks[j] * x))
            end
        end
        if isa(x, Vector)
            return x[1] * (1 + sum) - ξ
        else
            return x * (1 + sum) - ξ
        end
    end
end

function concentrate_unbounc_lignand_fpi(M, ξ, ks, ns)
    return function f(x)
        sum = 0
        for j in 1:M
            sum += (ks[j] * ns[j]) / (1 + (ks[j] * x))
        end
        return ξ - (x * sum)
    end
end

function concentrate_unbound_lignand_solver(M, ξ, ks, ns, type, guess)
    func = concentrate_unbound_lignand(M, ξ, ks, ns)
    func_fpi = concentrate_unbounc_lignand_fpi(M, ξ, ks, ns)
    if type == bisect_type
        println(bisection(func, guess, 10, 10.0^-6))
    elseif type == fpi_type
        fpi(func_fpi, guess)
    elseif type == newton_type
        newton_method(func, [guess], 50)
    end

    return "done solving"
end

function plot_concentrate_unbound_lignand(M, ξ, ks, ns, UB, LB=0)
    f = concentrate_unbound_lignand(M, ξ, ks, ns)
    plot(f, LB, UB, title = "Ligands and Binding Molecules", label = "f(x)")
    xlabel!("x")
    ylabel!("f(x)")
end

function derivative_concentrate_unbound_lignand(M, ξ, ks, ns)
    return function f_prime(x)
        gsum = 0
        gsum_prime = 0
        for j in 1:M
            gsum_prime += -((ks[j]^2) * (ns[j])) / ((1 + (ks[j] * x))^2)
            gsum += (ks[j] * ns[j]) / (1 + (ks[j] * x[1]))
        end
        return gsum + (x[1] * gsum_prime)
    end
end

function plot_derivative_concentrate_unbound_lignand(M, ξ, ks, ns, UB, LB=0)
    f_prime = derivative_concentrate_unbound_lignand(M, ξ, ks, ns)
    plot(f_prime, LB, UB, title = "Derivative", label = "f'(x)")
    xlabel!("x")
    ylabel!("f'(x)")
end

M = 1
ξ = 3
ks = [1]
ns = [1]
UB = 50
plot_concentrate_unbound_lignand(M, ξ, ks, ns, UB)
plot_derivative_concentrate_unbound_lignand(M, ξ, ks, ns, UB)


# =================== Fixed-point Iteration ===================

M = 1
ξ = 3
ks = [1]
ns = [1]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)


M = 1
ξ = 3
ks = [1]
ns = [5]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)



# =================== Newton's Method ===================

using ForwardDiff

M = 1
ξ = 1
ks = [1]
ns = [10]
guess = 1.6
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)


guess = 1.4999
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# this takes a while to converge


guess = 1.5
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# this takes sooooo long to converge, not even sure if it will

guess = 1.0
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# 6 steps to converge, nice



# ================= Find Positive Root ==================

M = 3
ξ = 9
ks = [1, 2, 6]
ns = [2, 3, 1]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)




# Root Finding Algorithms
# ========== Functions ===========
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
            return string("x=", x_m, ", f(x)=", f(x_m), ", steps=", counter)
        end
        f(a) * f(x_m) <= 0 ? b = x_m : a = x_m
        counter += 1
        if counter >= (n + 2)
            return string("x=", x_m, ", f(x)=", f(x_m), ", n=", n, " iterations reached.")
        end
        println( "step=", counter, ", x=", x_m, ", f(x)=", f(x_m))
    end
end

function fpi(func, xold; max_iter=50, FACC=10.0^-10, DIVERGENCE_THRESHHOLD=10.0^20)
    diff = 10*FACC
    step = 0
    xnew = 420.69
    while step < max_iter && abs(diff) > FACC && abs(diff) < DIVERGENCE_THRESHHOLD
        println("step=", step, "   xold,xnew= ", xold, ", ", xnew, "  diff=", diff)
        xnew = func(xold)
        diff = xnew-xold
        xold = xnew
        step += 1
    end
end

function newton_method(f, x::Vector, num_iter::Int64;
                       forward_err_threshold::Float64=1e-8)
    f_prime = (x_in -> ForwardDiff.gradient(f, x_in)[1])
    f_x = f(x)
    err = abs(f_x)
    println("step:", 0, ", ", x, " error:", err)

    f_prime_x = f_prime(x)
    if iszero(f_prime_x)
        println("SINGULARITY")
        return
    end
    if isinf(x[1])
        println("SINGULARITY")
        return
    end
    if isinf(err) || isnan(err)
        println("DIVERGENCE/ERROR")
        return
    end
    if err < forward_err_threshold
        println("CONVERGENCE")
        return
    end
    for i in 1:num_iter
        x = [x[1] - f_x/f_prime_x] # conversion between scalar and vector gets complicated
        f_x = f(x)
        err = abs(f_x)
        f_prime_x = f_prime(x)
        println("step:",i, ", ", x, " error:", err)
        if iszero(f_prime_x)
            println("SINGULARITY")
            break
        end
        if isinf(x[1])
            println("SINGULARITY")
            break
        end
        if isinf(err) || isnan(err)
            println("DIVERGENCE/ERROR")
            break
        end
        if err < forward_err_threshold
            println("CONVERGENCE")
            break
        end
    end
    return x
end
