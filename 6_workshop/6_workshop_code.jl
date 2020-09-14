#=
Workshop #6 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#

# Question 1
# TODO: Perry - not tested on anything yet, maybe broken
function steffensen(g, x0::Float64;
                    max_iter::Int64=50,
                    x_convergence_threshold::Float64=1e-8)
    loop_iters = max_iter÷3
    for i in 0:loop_iters:
        x1 = g(x0)
        if abs(x1 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 1)
            x0=x1
            break
        end

        x2 = g(x1)
        if abs(x2 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 2)
            x0=x2
            break
        end

        x3 = x0 - ((x1 - x0)^2) / (x2 - 2x1 + x0)
        if abs(x3 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 3)
            x0=x3
            break
        end
        x0 = x3
    end
    println("x=",x0)
end

# Question 2

# Question 3
println("g(x)=cos(x); x0 = 1; expected root approx 0.7")
#steffensen(x::Float64->cos(x), )
