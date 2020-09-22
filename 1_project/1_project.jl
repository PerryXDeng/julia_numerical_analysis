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
function concentrate_unbound_lignand_fpi(M, ξ, ks, ns, genAlpha::Bool=true, α=1)
    return function f(x)
        sum = 0
        if genAlpha
            α = 0
        end
        for j in 1:M
            sum += (ks[j] * ns[j]) / (1 + (ks[j] * x))
            if genAlpha
                α += ks[j] * ns[j]
            end
        end
        if genAlpha
            return (1 / (1 + α)) * (ξ + (x * (α - sum)))
        else
            return ξ - (x * sum)
        end
    end
end

function concentrate_unbound_lignand_solver(M, ξ, ks, ns, type, guess, genAlpha::Bool=true, root=nothing)
    func = concentrate_unbound_lignand(M, ξ, ks, ns)
    func_fpi = concentrate_unbound_lignand_fpi(M, ξ, ks, ns, genAlpha)

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
    plot(f, LB, UB, title = "Lignands and Binding Molecules", label = "f(x)")
    xlabel!("x")
    ylabel!("f(x)")
end

function derivative_concentrate_unbound_lignand(M, ξ, ks, ns)
    return function f_prime(x)
        gsum = 0
        gsum_prime = 0
        for j in 1:M
            gsum_prime += -((ks[j]^2) * (ns[j])) / ((1 + (ks[j] * x))^2)
            gsum += (ks[j] * ns[j]) / (1 + (ks[j] * x[1]))
        end
        return (1 + gsum) + (x[1] * gsum_prime)
    end
end

function plot_derivative_concentrate_unbound_lignand(M, ξ, ks, ns, UB, LB=0)
    f_prime = derivative_concentrate_unbound_lignand(M, ξ, ks, ns)
    plot(f_prime, LB, UB, title = "Derivative", label = "f'(x)")
    xlabel!("x")
    ylabel!("f'(x)")
end

M = 1
ξ = 3
ks = [1]
ns = [1]
UB = 5
plot_concentrate_unbound_lignand(M, ξ, ks, ns, UB)
plot_derivative_concentrate_unbound_lignand(M, ξ, ks, ns, UB)


# =================== Fixed-point Iteration ===================

M = 1
ξ = 3
ks = [1]
ns = [1]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess, false) # α = 0
#= OUTPUT
step=0   xold,xnew= 2, 420.69  diff=1.0e-9
step=1   xold,xnew= 2.3333333333333335, 2.3333333333333335  diff=0.3333333333333335
step=2   xold,xnew= 2.3, 2.3  diff=-0.03333333333333366
step=3   xold,xnew= 2.303030303030303, 2.303030303030303  diff=0.003030303030303383
step=4   xold,xnew= 2.302752293577982, 2.302752293577982  diff=-0.00027800945232137053
step=5   xold,xnew= 2.302777777777778, 2.302777777777778  diff=2.5484199796199647e-5
step=6   xold,xnew= 2.302775441547519, 2.302775441547519  diff=-2.33623025902574e-6
step=7   xold,xnew= 2.3027756557168324, 2.3027756557168324  diff=2.1416931339501843e-7
step=8   xold,xnew= 2.302775636083269, 2.302775636083269  diff=-1.963356321965648e-8
step=9   xold,xnew= 2.3027756378831383, 2.3027756378831383  diff=1.799869142615762e-9
step=10   xold,xnew= 2.3027756377181388, 2.3027756377181388  diff=-1.6499956956295136e-10
=#

M = 1
ξ = 3
ks = [1]
ns = [5]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess, false) # α = 0
#= OUTPUT
gives a negative solution, so this doesn't work
step=60   xold,xnew= -3.791287847538623, -3.791287847538623  diff=-1.552935557924684e-10

concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess) # α = generalized
#= OUTPUT
step=0   xold,xnew= 2, 420.69  diff=1.0e-9
step=1   xold,xnew= 1.611111111111111, 1.611111111111111  diff=-0.38888888888888906
step=2   xold,xnew= 1.3284081954294717, 1.3284081954294717  diff=-0.2827029156816392
step=3   xold,xnew= 1.131571795539733, 1.131571795539733  diff=-0.19683639988973867
step=4   xold,xnew= 1.0005909767072425, 1.0005909767072425  diff=-0.13098081883249058
step=5   xold,xnew= 0.9170360634785159, 0.9170360634785159  diff=-0.08355491322872666
step=6   xold,xnew= 0.8655622159563672, 0.8655622159563672  diff=-0.051473847522148675
step=7   xold,xnew= 0.8346613795041341, 0.8346613795041341  diff=-0.03090083645223307
step=8   xold,xnew= 0.8164342414983934, 0.8164342414983934  diff=-0.018227138005740717
step=9   xold,xnew= 0.8058028267060566, 0.8058028267060566  diff=-0.010631414792336757
step=10   xold,xnew= 0.7996442846752565, 0.7996442846752565  diff=-0.006158542030800129
step=11   xold,xnew= 0.7960913755099789, 0.7960913755099789  diff=-0.003552909165277618
step=12   xold,xnew= 0.7940466016015011, 0.7940466016015011  diff=-0.0020447739084777927
step=13   xold,xnew= 0.7928714360745759, 0.7928714360745759  diff=-0.0011751655269252037
step=14   xold,xnew= 0.7921965947538673, 0.7921965947538673  diff=-0.0006748413207086079
step=15   xold,xnew= 0.7918092461585348, 0.7918092461585348  diff=-0.0003873485953325151
step=16   xold,xnew= 0.7915869735419161, 0.7915869735419161  diff=-0.00022227261661866837
step=17   xold,xnew= 0.7914594462595641, 0.7914594462595641  diff=-0.0001275272823519913
step=18   xold,xnew= 0.79138628488139, 0.79138628488139  diff=-7.316137817414692e-5
step=19   xold,xnew= 0.7913443149141133, 0.7913443149141133  diff=-4.196996727667823e-5
step=20   xold,xnew= 0.791320239002048, 0.791320239002048  diff=-2.4075912065324445e-5
step=21   xold,xnew= 0.7913064281772002, 0.7913064281772002  diff=-1.3810824847748648e-5
step=22   xold,xnew= 0.791298505858158, 0.791298505858158  diff=-7.922319042230619e-6
step=23   xold,xnew= 0.7912939613943079, 0.7912939613943079  diff=-4.544463850120195e-6
step=24   xold,xnew= 0.7912913545709143, 0.7912913545709143  diff=-2.6068233935916396e-6
step=25   xold,xnew= 0.7912898592314975, 0.7912898592314975  diff=-1.4953394167349643e-6
step=26   xold,xnew= 0.7912890014681146, 0.7912890014681146  diff=-8.577633828776854e-7
step=27   xold,xnew= 0.7912885094342826, 0.7912885094342826  diff=-4.920338320113515e-7
step=28   xold,xnew= 0.7912882271918669, 0.7912882271918669  diff=-2.822424157722381e-7
step=29   xold,xnew= 0.7912880652908747, 0.7912880652908747  diff=-1.6190099216828457e-7
step=30   xold,xnew= 0.7912879724206082, 0.7912879724206082  diff=-9.28702664770853e-8
step=31   xold,xnew= 0.7912879191480142, 0.7912879191480142  diff=-5.32725940027845e-8
step=32   xold,xnew= 0.7912878885895882, 0.7912878885895882  diff=-3.055842601185077e-8
step=33   xold,xnew= 0.7912878710605491, 0.7912878710605491  diff=-1.752903910912096e-8
step=34   xold,xnew= 0.791287861005476, 0.791287861005476  diff=-1.0055073107473333e-8
step=35   xold,xnew= 0.7912878552376471, 0.7912878552376471  diff=-5.767828881175774e-9
step=36   xold,xnew= 0.7912878519290836, 0.7912878519290836  diff=-3.3085635342544606e-9
step=37   xold,xnew= 0.791287850031213, 0.791287850031213  diff=-1.897870527400869e-9
step=38   xold,xnew= 0.7912878489425494, 0.7912878489425494  diff=-1.0886636037099606e-9
step=39   xold,xnew= 0.7912878483180662, 0.7912878483180662  diff=-6.244832428947689e-10
step=40   xold,xnew= 0.7912878479598477, 0.7912878479598477  diff=-3.5821845489891757e-10
step=41   xold,xnew= 0.7912878477543651, 0.7912878477543651  diff=-2.054826309105806e-10
step=42   xold,xnew= 0.7912878476364954, 0.7912878476364954  diff=-1.1786971398919377e-10
=#



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
title!("Bisection Method Errors")

errors = concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess, root)
plot([1:length(errors)-1], errors[2:end] ./ errors[1:end-1])
title!("FPI Method Errors")

errors = concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess, root)
plot([1:length(errors)-1], errors[2:end] ./ ( errors[1:end-1]))
title!("Newton's Method Errors")



# =================== 3 Different Cases ===================

M = 5
ξ = 4
ks = [4, 7, 9, 3, 2]
ns = [2, 3, 1, 2, 2]
guess = 2.0
println("\nCase 1: Bisection Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
println("\nCase 1: FPI Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
println("\nCase 1: Newton Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# === OUTPUT === #
#=
Case 1: Bisection Method
no root in range [a,b]:[2.0,10]

Case 1: FPI Method
step=0   xold,xnew= 2.0, 420.69  diff=1.0e-9
step=1   xold,xnew= -4.839431913116124, -4.839431913116124  diff=-6.839431913116124
step=2   xold,xnew= -6.602089475828411, -6.602089475828411  diff=-1.7626575627122865
step=3   xold,xnew= -6.432408054553257, -6.432408054553257  diff=0.16968142127515407
step=4   xold,xnew= -6.444463938599245, -6.444463938599245  diff=-0.012055884045988563
step=5   xold,xnew= -6.443585214085237, -6.443585214085237  diff=0.0008787245140080557
step=6   xold,xnew= -6.443649144736874, -6.443649144736874  diff=-6.393065163656786e-5
step=7   xold,xnew= -6.443644492909572, -6.443644492909572  diff=4.651827302026845e-6
step=8   xold,xnew= -6.443644831390193, -6.443644831390193  diff=-3.3848062130914514e-7
step=9   xold,xnew= -6.443644806761332, -6.443644806761332  diff=2.4628860728626023e-8
step=10   xold,xnew= -6.443644808553401, -6.443644808553401  diff=-1.7920687156447457e-9
step=11   xold,xnew= -6.443644808423006, -6.443644808423006  diff=1.3039525015301479e-10

Case 1: Newton Method
step:0, [2.0] error:6.839431913116124
step:1, [-2.561207076519719] error:4.662087851270968
step:2, [-5.5591765709406875] error:0.9597522852464397
step:3, [-6.432018247892473] error:0.01247416091732756
step:4, [-6.443643298511636] error:1.6197865440048531e-6
step:5, [-6.443644808431824] error:2.842170943040401e-14
CONVERGENCE
=#

M = 3
ξ = 9
ks = [1, 5, 7]
ns = [2, 7, 9]
guess = 2
println("\nCase 2: Bisection Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
println("\nCase 2: FPI Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
println("\nCase 2: Newton Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# === OUTPUT === #
#=
Case 2: Bisection Method
no root in range [a,b]:[2,10]

Case 2: FPI Method
step=0   xold,xnew= 2, 420.69  diff=1.0e-9
step=1   xold,xnew= -7.096969696969698, -7.096969696969698  diff=-9.096969696969698
step=2   xold,xnew= -9.71590496502775, -9.71590496502775  diff=-2.6189352680580527
step=3   xold,xnew= -9.510893303631672, -9.510893303631672  diff=0.20501166139607818
step=4   xold,xnew= -9.52259927526886, -9.52259927526886  diff=-0.011705971637187673
step=5   xold,xnew= -9.521916422899729, -9.521916422899729  diff=0.0006828523691311261
step=6   xold,xnew= -9.521956207068555, -9.521956207068555  diff=-3.978416882688407e-5
step=7   xold,xnew= -9.52195388900681, -9.52195388900681  diff=2.318061746109379e-6
step=8   xold,xnew= -9.521954024070276, -9.521954024070276  diff=-1.3506346618896714e-7
step=9   xold,xnew= -9.52195401620071, -9.52195401620071  diff=7.869566331919486e-9
step=10   xold,xnew= -9.521954016659233, -9.521954016659233  diff=-4.5852388552702905e-10

Case 2: Newton Method
step:0, [2] error:9.096969696969698
step:1, [-3.07791195948578] error:7.809111139149351
step:2, [-7.461382657410847] error:2.2166287782524545
step:3, [-9.479322863532037] error:0.045126990282222224
step:4, [-9.521942720871907] error:1.1953918431117927e-5
step:5, [-9.521954016633202] error:8.313350008393172e-13
CONVERGENCE
=#

M = 3
ξ = 2
ks = [1, 2, 6]
ns = [2, 3, 1]
guess = 2
println("\nCase 3: Bisection Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
println("\nCase 3: FPI Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
println("\nCase 3: Newton Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# === OUTPUT === #
#=
Case 3: Bisection Method
no root in range [a,b]:[2,10]

Case 3: FPI Method
step=0   xold,xnew= 2, 420.69  diff=1.0e-9
step=1   xold,xnew= -2.6564102564102567, -2.6564102564102567  diff=-4.656410256410257
step=2   xold,xnew= -5.96997211381284, -5.96997211381284  diff=-3.3135618574025836
step=3   xold,xnew= -4.705360442520642, -4.705360442520642  diff=1.2646116712921982
step=4   xold,xnew= -4.9331674584822425, -4.9331674584822425  diff=-0.22780701596160036
step=5   xold,xnew= -4.881820796575109, -4.881820796575109  diff=0.05134666190713322
step=6   xold,xnew= -4.892892606665089, -4.892892606665089  diff=-0.011071810089979373
step=7   xold,xnew= -4.890481664598472, -4.890481664598472  diff=0.002410942066616961
step=8   xold,xnew= -4.891005545286776, -4.891005545286776  diff=-0.0005238806883038905
step=9   xold,xnew= -4.890891657088721, -4.890891657088721  diff=0.00011388819805446815
step=10   xold,xnew= -4.8909164131432865, -4.8909164131432865  diff=-2.4756054565422403e-5
step=11   xold,xnew= -4.890911031763738, -4.890911031763738  diff=5.381379548197174e-6
step=12   xold,xnew= -4.890912201542562, -4.890912201542562  diff=-1.1697788240283558e-6
step=13   xold,xnew= -4.890911947261315, -4.890911947261315  diff=2.5428124761361914e-7
step=14   xold,xnew= -4.890912002535814, -4.890912002535814  diff=-5.527449964404241e-8
step=15   xold,xnew= -4.890911990520493, -4.890911990520493  diff=1.201532118244586e-8
step=16   xold,xnew= -4.89091199313233, -4.89091199313233  diff=-2.611836968924308e-9
step=17   xold,xnew= -4.890911992564581, -4.890911992564581  diff=5.677494030464914e-10
step=18   xold,xnew= -4.890911992687996, -4.890911992687996  diff=-1.234150559525915e-10

Case 3: Newton Method
step:0, [2] error:4.656410256410257
step:1, [-1.1089884286491891] error:23.88155260125967
step:2, [-1.2465517929290493] error:13.028909149471417
step:3, [-1.6012186766881678] error:7.203665522324506
step:4, [-2.518769654860344] error:3.611971127865367
step:5, [-4.113290700232235] error:0.9864806189457216
step:6, [-4.853928121817551] error:0.04509676358156023
step:7, [-4.890851256040428] error:7.393947344258578e-5
step:8, [-4.890911992504588] error:1.964490792261131e-10
CONVERGENCE
=#



# =================== Root Finding Algorithms ===================

function bisection(f, a, b, tolerance, root)
    if f(a) == 0 return a end
    if f(b) == 0 return b end
    if f(a) * f(b) > 0
        println("no root in range [a,b]:[", a, ",", b, "]")
        return 0
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
    max_iter=150, FACC=10.0^-10, DIVERGENCE_THRESHHOLD=10.0^20)
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

function newton_method(f, x::Vector, num_iter::Float64, root, forward_err_threshold::Float64=1e-8)
    f_prime = (x_in -> ForwardDiff.gradient(f, x_in)[1])
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
