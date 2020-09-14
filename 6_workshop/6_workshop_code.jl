#=
Workshop #6 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#

# TODO: Perry - not tested on anything yet, maybe broken
function steffensen(g, x0::Float64;
                    max_iter::Int64=50,
                    x_convergence_threshold::Float64=1e-8)
    for i in 1:max_iter:
        x1 = g(x0)
        x2 = g(x1)
        x3 = x0 - ((x1 - x0)^2) / (x2 - 2x1 + x0)
        if abs(x3 - x0) < x_convergence_threshold
            println("CONVERGENCE at x=",x3)
            break
        end
        x0 = x3
    end
    println("EXCEEDED max_iter at x=",x3)
end
