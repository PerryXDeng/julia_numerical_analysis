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




# Question 4
function question_4(x, max_n, threshold)
    println("x: ", x)
    sum = 0
    for n in 0:max_n
        term = ((3.0*x)^n)/factorial(big(n))
        sum += term
        #println(term)
        if abs(term) < threshold
            println("Reached threashold at n=", i)
            println("sum: ", sum)
            return
        end
    end
    println("Stopped at n=", max_n)
    println("sum: ", sum)
    return
end

println("First Part")
question_4(2,100, 10.0^-9)
question_4(-2,100, 10.0^-9)
question_4(-12,100, 10.0^-9)
println()
println("Second Part")
question_4(2,99999, 10.0^-15)
question_4(-2,99999, 10.0^-15)
question_4(-12,99999, 10.0^-15)
println()
println("Third Part")
question_4(20,99999, 10.0^-15)
question_4(-20,99999, 10.0^-15)
#=
Output:
First Part
x: 2
Reached threashold at n=30
sum: 403.4287934925368953042380035039178855952930700943860429870714844947427734992497
x: -2
Reached threashold at n=30
sum: 0.002478752312403788906848224047953192464984510657938127357093666775546201459982776
x: -12
Stopped at n=100
sum: 0.02693711718205233623544164316681056498927302324551674853214105955994607046619373

Second Part
x: 2
Reached threashold at n=38
sum: 403.4287934927351224799854552166752647226141726622076519934716990286399703491546
x: -2
Reached threashold at n=38
sum: 0.002478752176666453359615404993428939203360429655544956872785297639258892548241384
x: -12
Reached threashold at n=126
sum: 0.02573302635381157790151388191216509352840576198434274992371366140787095476048348

Third Part
x: 20
Stopped at n=100
sum: 1.142006397937806098761365236566413971041926087530146112196237928690682336949673e+26
x: -20
Stopped at n=100
sum: 2.614878889445793964482521864727837075302150645365233360179877459749452006797551e+19
=#
#=
Comments:
The taylor sums are compared to the function evaluations in wolfram alpha for
accuracy.

part one - the results on x=2 and x=-2 appears to mostly agree with the output
from wolfram alpha. For x=-12, the result is off - the sum should be smaller
than the sum for x=-2, but it turns out to be larger. The number of iterations
is also much higher on x=-12, reaching the maximum iterations without triggering
the convergence difference threshold. We hypothesize that this is because
the polynomial terms are larger in absolute terms with larger values of x,
and thus takes more iterations to converge to a small number. Perhaps with the
relaxing of the max-iteration stopping criterion in part two, higher accuracy
can be achieved.

part two - for this part, we set maximum iterations to a much larger number.
for n=2, -2, and -12, we gain better accuracy. However, for x=-12, the accuracy
is still much worse than x=-2 even afte convergence. Our new hypothesis is
that in x=-12 we have large terms of alternating signs, which degrades the
floating point accuracy.

part three - we have positive infinity for x=20, which is infinitely inaccurate.
According to wolfram alpha, the number for x=20 should not exceed the floating
point representation range. For x=-20, the sum gives NaN, which is wrong by
an immusurable amount. We suspect that like x=-12, this is because of large
alternating signs of floating point terms.
=#

# Question 5
# notes:
# 1 + hv/kt + ((hv/kt)^2)/2
# then replace v with λ
