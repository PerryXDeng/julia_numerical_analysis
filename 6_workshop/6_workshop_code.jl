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
    for i in 0:loop_iters
        x1 = g(x0)
        if abs(x1 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 1)
            x0=x1
            break
        end
        if 3i+1 == max_iter
            break
        end

        x2 = g(x1)
        if abs(x2 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 2)
            x0=x2
            break
        end
        if 3i+2 == max_iter
            break
        end

        x3 = x0 - ((x1 - x0)^2) / (x2 - 2x1 + x0)
        if abs(x3 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", 3i + 3)
            x0=x3
            break
        end
        if 3i+3 == max_iter
            break
        end

        x0 = x3
    end
    println("x=",x0)
end

function plain_fpi(g, x0::Float64;
                   max_iter::Int64=50,
                   x_convergence_threshold::Float64=1e-8)
    for i in 1:max_iter
        x1 = g(x0)
        if abs(x1 - x0) < x_convergence_threshold
            println("CONVERGENCE at i=", i)
            x0=x1
            break
        end
        x0=x1
    end
    println("x=",x0)
end


# Question 2



# Question 3
println("\ng(x)=cos(x); x0 = 1; expected root approx 0.7")
g(x) = cos(x)
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
Improves convergence
Output:
g(x)=cos(x); x0 = 1; expected root approx 0.7
FPI
CONVERGENCE at i=46
x=0.7390851366465718
STEFFENSEN
CONVERGENCE at i=10
x=0.7390851332482249
=#
println("\ng(x)=x/2 + 1/x^2; x0 = 1; expected root approx 1.26")
g(x) = x/2 + 1/x^2
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
Improves convergence
Output:
g(x)=x/2 + 1/x^2; x0 = 1; expected root approx 1.26
FPI
CONVERGENCE at i=27
x=1.2599210520301471
STEFFENSEN
CONVERGENCE at i=13
x=1.2599210498948539
=#
println("\ng(x)=x-0.5((x^2)-x-1); x0 = 1; expected root 1.618")
g(x) = x-0.5((x^2)-x-1)
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
No improvement, possibly due to quadraticly converging FPI on g(x)
Output:
g(x)=x-0.5((x^2)-x-1); x0 = 1; expected root 1.618
FPI
CONVERGENCE at i=10
x=1.6180339890192184
STEFFENSEN
CONVERGENCE at i=10
x=1.618033988885942
=#
println("\ng(x)=(-2.5x+((x^2)-1))/-1.5; x0 = 1; expected root 1.618")
g(x) = (-2.5x+((x^2)-1))/(-1.5)
println("FPI")
plain_fpi(g, 1.0)
println("STEFFENSEN")
steffensen(g, 1.0)
#=
Improvement
Output:
g(x)=x-0.5((x^2)-x-1); x0 = 1; expected root 1.618
FPI
CONVERGENCE at i=25
x=1.6180339890192184
STEFFENSEN
CONVERGENCE at i=10
x=1.618033988885942
=#
