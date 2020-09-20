#=
Project #1

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
@enum SOLVE_TYPE begin
    bisect
    fpi
    newton
end


# Part 1
function concentrate_unbound_lignand(M, ξ, ks, ns, x)
    return function f(x)
        sum = 0
        for j in 1:M
            sum += (ks[j] * ns[j]) / (1 + (k[j] * x))
        end
        return x * (1 + sum) - ξ
    end
end

function concentrate_unbound_lignand_solver(M, ξ, ks, ns, type)
    if type == bisect
        bisection(concentrate_unbound_lignand(M, ξ, ks, ns, x), 0, 10, 1e^-5)
    elseif type = fpi
        fpi(concentrate_unbound_lignand(M, ξ, ks, ns, x), 0)
    elseif type = newton
        newton(concentrate_unbound_lignand(M, ξ, ks, ns, x), [1], 10)
    end

    return "done solving"
end

























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
