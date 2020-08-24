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
function determinant_(n)
    return (10.0)^(2n) + 4
end

function find_roots(n)
    sqrt_out = sqrt(determinant_(n))
    pos = (-(10.0^n) + sqrt_out ) / 2.0
    neg = (-(10.0^n) - sqrt_out ) / 2.0
    return (pos, neg)
end

function find_roots_part2(n)
    sqrt_out = sqrt(determinant_(n))
    pos_denominator = (10.0^n) + sqrt_out
    neg_denominator = (10.0^n) - sqrt_out
    pos = 2.0 / pos_denominator
    neg = 2.0 / neg_denominator
    return (pos, neg)
end


for i in 1:16
    println(i)
    println("first $(find_roots(i))")
    println("second $(find_roots_part2(i))")
    println()
end

#=
Answer:
The first solution is better because it is not yielding infinity for one of the
roots. Infinity is infinitely far away from the actual root, so the second
solution is infinitely wrong.

Output:
1
first (0.09901951359278449, -10.099019513592784)
second (0.09901951359278484, -10.09901951359282)

2
first (0.009999000199947261, -100.00999900019994)
second (0.009999000199950016, -100.00999900022748)

3
first (0.0009999989999869285, -1000.000999999)
second (0.000999999000002, -1000.0010000140715)

4
first (9.999999929277692e-5, -10000.0001)
second (9.999999900000002e-5, -10000.000070722308)

5
first (1.0000003385357559e-5, -100000.00001)
second (9.999999999e-6, -99999.96614643588)

6
first (1.00000761449337e-6, -1.000000000001e6)
second (9.99999999999e-7, -999992.38556461)

7
first (9.96515154838562e-8, -1.00000000000001e7)
second (9.999999999999899e-8, -1.0034970317757009e7)

8
first (7.450580596923828e-9, -1.0e8)
second (1.0e-8, -1.34217728e8)

9
first (0.0, -1.0e9)
second (1.0e-9, Inf)

10
first (0.0, -1.0e10)
second (1.0e-10, Inf)

11
first (0.0, -1.0e11)
second (1.0e-11, Inf)

12
first (0.0, -1.0e12)
second (1.0e-12, Inf)

13
first (0.0, -1.0e13)
second (1.0e-13, Inf)

14
first (0.0, -1.0e14)
second (1.0e-14, Inf)

15
first (0.0, -1.0e15)
second (1.0e-15, Inf)

16
first (0.0, -1.0e16)
second (1.0e-16, Inf)
=#



# Question 4

function break_math()
    println(10^19)
end

function fix_math()
    println(10.0^19)
end

break_math()
fix_math()

#=
Answers:
If we use a float, rather than an int, with a large exponent then we have the
greater range given to use by floating point numbers.

Output:
-8446744073709551616
1.0e19
=#
