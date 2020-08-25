#=
Homework #1 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#

# Question 1
function question_1()
    println("Hello World!")
    input = "Invalid Input"
    try
        input = parse(Int, readline())
    catch err
        if isa(err, ArgumentError)
            input = "Invalid Input"
        end
    end
    println(input)

    # calculate pi
    p = 4 * (4*atan(1/5) - atan(1/239))
    println(p)
end

question_1()

#=
Output:
> Hello World!
4
> 4
> 3.1415926535897936

Fun Fact: 'pi' is a keyword in Julia and is stored with type Irrational, which
is automatically rounded to the correct precision in arithmetic operations.
https://julialang.org/blog/2017/03/piday/
=#


# Question 2
m = rand(1:10,(4,4))
show(stdout, "text/plain", m)
show(stdout, "text/plain", inv(m))

#=
Output:
4×4 Array{Int64,2}:
 7  4  7   5
 8  5  9  10
 3  1  4   9
 7  6  2  10
Array{Float64,2}:
  1.0       -1.0         0.454545   0.0909091
 -0.884615   0.961538   -0.646853   0.0629371
 -0.269231   0.423077   -0.157343  -0.146853
 -0.115385   0.0384615   0.101399   0.027972
=#


# Question 3
# a)
