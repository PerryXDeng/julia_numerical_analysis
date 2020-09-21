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

# Fixed-point Iteration scheme
function concentrate_unbounc_lignand_fpi(M, ξ, ks, ns)
    return function f(x)
        sum = 0
        for j in 1:M
            sum += (ks[j] * ns[j]) / (1 + (ks[j] * x))
        end
        return ξ - (x * sum)
    end
end

function concentrate_unbound_lignand_solver(M, ξ, ks, ns, type, guess, root=nothing)
    func = concentrate_unbound_lignand(M, ξ, ks, ns)
    func_fpi = concentrate_unbounc_lignand_fpi(M, ξ, ks, ns)

    if type == bisect_type
        return bisection(func, guess, 10, 10.0^-6, root)
    elseif type == fpi_type
        return fpi(func_fpi, guess, root)
    elseif type == newton_type
        return newton_method(func, [guess], 80, root)
    end

    return "done solving"
end

function plot_concentrate_unbound_lignand(M, ξ, ks, ns, UB, LB=0)
    f = concentrate_unbound_lignand(M, ξ, ks, ns)
    plot(f, LB, UB, title = "Ligands and Binding Molecules", label = "f(x)")
    xlabel!("Concentration")
    ylabel!("")
end

M = 1
ξ = 3
ks = [1]
ns = [1]
UB = 50
plot_concentrate_unbound_lignand(M, ξ, ks, ns, UB)

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



# ================== Find Positive Root & Error Calculations ==================

M = 3
ξ = 9
ks = [1, 2, 6]
ns = [2, 3, 1]
guess = 2
root = 3.806382815465563265304167

errors = concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess, root)
plot([1:length(errors)-1], errors[2:end] ./ errors[1:end-1])


errors = concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess, root)
plot([1:length(errors)-1], errors[2:end] ./ errors[1:end-1])


errors = abs.(concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess, root))
plot(errors[1:end-1], errors[2:end], yscale=:log10, xscale=:log10)
title!("Newton's Method ")

plot([1:length(errors)-1], errors[2:end] ./ ( errors[1:end-1]))



# =================== 3 Different Cases ===================

M = 5
ξ = 4
ks = [4, 7, 9, 3, 2]
ns = [2, 3, 1, 2, 2]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)


M = 3
ξ = 9
ks = [1, 5, 7]
ns = [2, 7, 9]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)


M = 3
ξ = 2
ks = [1, 2, 6]
ns = [2, 3, 1]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)




# =================== Root Finding Algorithms ===================

function bisection(f, a, b, tolerance, root)
    if f(a) == 0 return a end
    if f(b) == 0 return b end
    if f(a) * f(b) > 0
        return string("no root in range [a,b]:[", a, ",", b, "]")
    end

    ϵ_n = []
    f_t = tolerance
    n = log2((b-a)/f_t)
    counter = 0
    while true
        dx = b - a
        x_m = a + (dx/2.0)
        if abs(f(x_m)) < f_t
            println(string(" -- DONE -- x=", x_m, ", f(x)=", f(x_m), ", steps=", counter))
            return ϵ_n
        end
        f(a) * f(x_m) <= 0 ? b = x_m : a = x_m
        counter += 1
        if counter >= (n + 2)
            println(string("x=", x_m, ", f(x)=", f(x_m), ", n=", n, " iterations reached."))
            return ϵ_n
        end
        if root != nothing
            append!(ϵ_n, x_m - root)
        end
        println( "step=", counter, ", x=", x_m, ", f(x)=", f(x_m))
    end
end

function fpi(func, xold, root;
    max_iter=50, FACC=10.0^-10, DIVERGENCE_THRESHHOLD=10.0^20)
    diff = 10*FACC
    step = 0
    xnew = 420.69
    ϵ_n = []
    while step < max_iter && abs(diff) > FACC && abs(diff) < DIVERGENCE_THRESHHOLD
        println("step=", step, "   xold,xnew= ", xold, ", ", xnew, "  diff=", diff)
        xnew = func(xold)
        diff = xnew-xold
        xold = xnew
        step += 1
        if root != nothing
            append!(ϵ_n, xold - root)
        end
    end
    return ϵ_n
end

function newton_method(f, x::Vector, num_iter::Int64, root::Float64, forward_err_threshold::Float64=1e-8)
    f_prime = (x_in -> ForwardDiff.gradient(f, x_in)[1])
    #f_2prime = (x_in -> ForwardDiff.gradient(f_prime, x_in)[1])
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
    ϵ_n = []
    for i in 1:num_iter
        x = [x[1] - f_x/f_prime_x] # conversion between scalar and vector gets complicated
        f_x = f(x)
        err = abs(f_x)
        f_prime_x = f_prime(x)
        println("step:",i, ", ", x, " error:", err)
        if root != nothing
            append!(ϵ_n, x[1] - root)
        end
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
    return ϵ_n
end
