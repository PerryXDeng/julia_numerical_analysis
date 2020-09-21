#=
Project #1
Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
ENV["GKS_ENCODING"]="utf8"
@enum SOLVE_TYPE begin
    bisect_type
    fpi_type
    newton_type
end


# =================== Part 1 ===================
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

function concentrate_unbound_lignand_solver(M, ξ, ks, ns, type, guess)
    func = concentrate_unbound_lignand(M, ξ, ks, ns)
    if type == bisect_type
        println(bisection(func, guess, 10, 10.0^-6))
    elseif type == fpi_type
        fpi(func, guess)
    elseif type == newton_type
        newton_method(func, [guess], 50)
    end

    return "done solving"
end


# =================== Fixed-point Iteration ===================
function f(x)
    return 3 - x * (1/(1+x))
end

fpi(f, 2)
#= -- Output --
step=0   xold,xnew= 2, 420.69  diff=1.0e-9
step=1   xold,xnew= 2.25, 2.25  diff=0.25
step=2   xold,xnew= 2.2941176470588234, 2.2941176470588234  diff=0.04411764705882337
step=3   xold,xnew= 2.3013698630136985, 2.3013698630136985  diff=0.00725221595487513
step=4   xold,xnew= 2.3025477707006368, 2.3025477707006368  diff=0.001177907686938262
step=5   xold,xnew= 2.302738712065137, 2.302738712065137  diff=0.0001909413645000413
step=6   xold,xnew= 2.3027696542232925, 2.3027696542232925  diff=3.0942158155689015e-5
step=7   xold,xnew= 2.3027746681592833, 2.3027746681592833  diff=5.013935990838547e-6
step=8   xold,xnew= 2.3027754806218117, 2.3027754806218117  diff=8.124625283656428e-7
step=9   xold,xnew= 2.302775612273765, 2.302775612273765  diff=1.3165195333897373e-7
step=10   xold,xnew= 2.3027756336067275, 2.3027756336067275  diff=2.1332962507614184e-8
step=11   xold,xnew= 2.302775637063534, 2.302775637063534  diff=3.456806396684442e-9
step=12   xold,xnew= 2.302775637623677, 2.302775637623677  diff=5.601430430601795e-10
Converges to 2.30277
=#

function f2(x)
    return 3 - x * (5/(1+x))
end

fpi(f2, 2)





# =================== Newton's Method ===================
using ForwardDiff
function g(x::Vector)
    return x[1] * (1 + (10 /(1 + x[1]))) - 1
end


concentrate_unbound_lignand_solver(1, 1, [1], [10], newton_type, 1.6)
#= -- Output --
step:0, [1.6] error:6.753846153846154
step:1, [-1.1241050119331741] error:88.45281806498997
step:2, [-1.2601310118670497] error:46.18203675669911
step:3, [-1.5705358283808182] error:24.95684720798347
step:4, [-2.3572989696013047] error:14.010274589330532
step:5, [-4.5368300976196] error:7.290560453803925
step:6, [-8.588461092752654] error:1.729329151905843
step:7, [-10.061914568142107] error:0.04160499439667453
step:8, [-10.099003087922368] error:1.8409638565586306e-5
step:9, [-10.099019513589589] error:3.582023566650605e-12
CONVERGENCE
=#

concentrate_unbound_lignand_solver(1, 1, [1], [10], newton_type, 1.4999)
# this takes a while to converge


concentrate_unbound_lignand_solver(1, 1, [1], [10], newton_type, 1.5)
# this takes sooooo long to converge, not even sure if it will

concentrate_unbound_lignand_solver(1, 1, [1], [10], newton_type, 1.0)
# 6 steps to converge, nice


concentrate_unbound_lignand_solver(3, 9, [1 2 6], [2 3 1], bisect_type, 2)
concentrate_unbound_lignand_solver(3, 9, [1 2 6], [2 3 1], newton_type, 2)





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
