#====
Numerical Analysis
Workshop #1

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#
using Plots

# Question 1
println("hello world")
#=
Output:
> hello world

Source:
https://juliabyexample.helpmanual.io/
=#


# Question 2
function count_collatz(num)
    counter = 0
    while num > 1
        if num % 2 == 0
            num /= 2
        else
            num = (num * 3) + 1
        end
        counter += 1
    end
    return counter
end

function plotcollatz(min, max)
    x = collect(min:2:max)
    y = []
    for i in min:2:max
        push!(y, count_collatz(i))
    end
    gr()
    plot(x,y, title="Collatz Conjecture Steps", lw=2)
    xlabel!("Input Value")
    ylabel!("# of steps")
end

plotcollatz(3, 99)


# Question 3
function find_roots(n)
    sqrt_out::Real = sqrt((10)^(2n) + 4)
    pos = (-(10^n) + sqrt_out ) / 2
    neg = (-(10^n) - sqrt_out ) / 2
    return (pos, neg)
end

function find_roots_part2(n)
    sqrt_out::Real = sqrt((10)^(2n) + 4)
    pos = 2 / ((10^n) + sqrt_out)
    neg = 2 / ((10^n) - sqrt_out)
    return (pos, neg)
end


for i in 1:16
    println("$i - $(find_roots(i))")
    println("$i - $(find_roots_part2(i))")
end


# Question 4
function break_math()

end
