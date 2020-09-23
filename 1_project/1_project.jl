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
            println(string("x=", x_m, ", f(x)=", f(x_m), ", n=", counter, " iterations reached."))
            return ϵ_n
        end
        if root != nothing
            append!(ϵ_n, x_m - root)
        end
        println( "step=", counter, ", x=", x_m, ", f(x)=", f(x_m))
    end
end

function fpi(func, xold, root;
    max_iter=100, FACC=10.0^-10, DIVERGENCE_THRESHHOLD=10.0^20)
    diff = 10*FACC
    step = 0
    xnew = 30
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

function newton_method(f, x::Vector, num_iter::Int64, root, forward_err_threshold::Float64=1e-8)
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
# use genAlpha = false when need to use α = 0 or guess α
function concentrate_unbound_lignand_fpi(M, ξ, ks, ns, genAlpha::Bool=true, α=0)
    return function f(x)
        sum = 0
        if genAlpha
            α = Base.sum(ks .* ns)
        end
        for j in 1:M
            sum += (ks[j] * ns[j]) / (1 + (ks[j] * x))
        end
        return (1 / (1 + α)) * (ξ + (x * (α - sum)))
    end
end

# use genAlpha = false when need to use α = 0
function concentrate_unbound_lignand_solver(M, ξ, ks, ns, type, guess, genAlpha::Bool=true, root=nothing, α=0)
    func = concentrate_unbound_lignand(M, ξ, ks, ns)
    func_fpi = concentrate_unbound_lignand_fpi(M, ξ, ks, ns, genAlpha, α)

    if type == bisect_type
        return bisection(func, guess, 100, 10.0^-6, root)
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
UB = 50
plot_concentrate_unbound_lignand(M, ξ, ks, ns, UB)



# =================== Fixed-point Iteration ===================

M = 1
ξ = 3
ks = [1]
ns = [1]
guess = 2
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess, false) # α = 0
#= OUTPUT
step=0   xold,xnew= 2, 30  diff=1.0e-9
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
=#
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
plotly()

errors = concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess, root)
plot([1:length(errors)-1], errors[2:end] ./ errors[1:end-1])
xaxis!("n")
yaxis!("ϵ_{n+1} / ϵ_n")
plot([1:length(errors)-1], abs.(errors[2:end] .- errors[1:end-1]))
xaxis!("n")
yaxis!("|ϵ_{n+1} - ϵ_n|")
plot([1:length(errors)-1], abs.(errors[2:end] .- root) ./ (abs.(errors[1:end-1] .- root) .^ (1)))
xaxis!("n")
yaxis!("|ϵ_{n+1} - r| / |ϵ_n - r|^α")

errors = concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess, root)
plot([1:length(errors)-1], errors[2:end] ./ errors[1:end-1])
xaxis!("n")
yaxis!("ϵ_{n+1} / ϵ_n")
plot([1:length(errors)-1], abs.(errors[2:end] - errors[1:end-1]))
xaxis!("n")
yaxis!("|ϵ_{n+1} - ϵ_n|")
plot([1:length(errors)-1], abs.(errors[2:end] .- root) ./ (abs.(errors[1:end-1] .- root) .^ (1)))
xaxis!("n")
yaxis!("|ϵ_{n+1} - r| / |ϵ_n - r|^α")

errors = concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess, root)
plot([1:length(errors)-1], errors[2:end] ./ ( errors[1:end-1]))
xaxis!("n")
yaxis!("ϵ_{n+1} / ϵ_n")
plot([1:length(errors)-1], abs.(errors[2:end] - errors[1:end-1]))
xaxis!("n")
yaxis!("|ϵ_{n+1} - ϵ_n|")
plot([1:length(errors)-1], abs.(errors[2:end] .- root) ./ (abs.(errors[1:end-1] .- root) .^ (2)))
xaxis!("n")
yaxis!("|ϵ_{n+1} - r| / |ϵ_n - r|^α")


# =================== 3 Different Cases ===================

M = 5
ξ = 4
ks = [4, 7, 9, 3, 2]
ns = [2, 3, 1, 2, 2]
guess = 0
println("\nCase 1: Bisection Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
println("\nCase 1: FPI Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
println("\nCase 1: Newton Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# === OUTPUT === #
#=
Case 1: Bisection Method
step=1, x=50.0, f(x)=55.94623843448594
step=2, x=25.0, f(x)=30.89319631074727
step=3, x=12.5, f(x)=18.289204307082393
step=4, x=6.25, f(x)=11.839156597090415
step=5, x=3.125, f(x)=8.342736393616882
step=6, x=1.5625, f(x)=6.13244215252031
step=7, x=0.78125, f(x)=4.329501134532899
step=8, x=0.390625, f(x)=2.5416663006936133
step=9, x=0.1953125, f(x)=0.74328652645488
step=10, x=0.09765625, f(x)=-0.8747034905607389
step=11, x=0.146484375, f(x)=0.036621014635512594
step=12, x=0.1220703125, f(x)=-0.3873940389297781
step=13, x=0.13427734375, f(x)=-0.1683863326913575
step=14, x=0.140380859375, f(x)=-0.06423062789954548
step=15, x=0.1434326171875, f(x)=-0.013403208361766339
step=16, x=0.14495849609375, f(x)=0.011707925128036756
step=17, x=0.144195556640625, f(x)=-0.0008227160203047035
step=18, x=0.1445770263671875, f(x)=0.0054488145613520444
step=19, x=0.14438629150390625, f(x)=0.002314604438213408
step=20, x=0.14429092407226562, f(x)=0.0007463333348143308
step=21, x=0.1442432403564453, f(x)=-3.809401949306235e-5
step=22, x=0.14426708221435547, f(x)=0.00035414398324817853
step=23, x=0.1442551612854004, f(x)=0.0001580310639273108
step=24, x=0.14424920082092285, f(x)=5.9970042810775226e-5
step=25, x=0.14424622058868408, f(x)=1.0938391817205684e-5
step=26, x=0.1442447304725647, f(x)=-1.3577718797286309e-5
step=27, x=0.1442454755306244, f(x)=-1.3196397294912288e-6
step=28, x=0.14424584805965424, f(x)=4.809381983328365e-6
x=0.1442456617951393, f(x)=1.744872611730841e-6, n=29 iterations reached.

Case 1: FPI Method
step=0   xold,xnew= 0, 30  diff=1.0e-9
step=1   xold,xnew= 0.08163265306122448, 0.08163265306122448  diff=0.08163265306122448
step=2   xold,xnew= 0.10688728111470411, 0.10688728111470411  diff=0.025254628053479625
step=3   xold,xnew= 0.12080539796100603, 0.12080539796100603  diff=0.013918116846301926
step=4   xold,xnew= 0.12919224922081127, 0.12919224922081127  diff=0.008386851259805234
step=5   xold,xnew= 0.13445396087025846, 0.13445396087025846  diff=0.005261711649447193
step=6   xold,xnew= 0.13782792470248206, 0.13782792470248206  diff=0.003373963832223603
step=7   xold,xnew= 0.14001942393909242, 0.14001942393909242  diff=0.002191499236610356
step=8   xold,xnew= 0.14145421221488236, 0.14145421221488236  diff=0.0014347882757899388
step=9   xold,xnew= 0.14239831084044935, 0.14239831084044935  diff=0.0009440986255669914
step=10   xold,xnew= 0.14302154804754727, 0.14302154804754727  diff=0.0006232372070979209
step=11   xold,xnew= 0.14343384030442494, 0.14343384030442494  diff=0.0004122922568776666
step=12   xold,xnew= 0.14370696274713463, 0.14370696274713463  diff=0.000273122442709689
step=13   xold,xnew= 0.14388805713187192, 0.14388805713187192  diff=0.00018109438473729034
step=14   xold,xnew= 0.14400820434970643, 0.14400820434970643  diff=0.00012014721783451532
step=15   xold,xnew= 0.14408794784311785, 0.14408794784311785  diff=7.974349341141607e-5
step=16   xold,xnew= 0.1441408887397391, 0.1441408887397391  diff=5.294089662125079e-5
step=17   xold,xnew= 0.14417604180932236, 0.14417604180932236  diff=3.5153069583265806e-5
step=18   xold,xnew= 0.14419938636395474, 0.14419938636395474  diff=2.334455463237628e-5
step=19   xold,xnew= 0.14421489027864665, 0.14421489027864665  diff=1.5503914691905596e-5
step=20   xold,xnew= 0.14422518748337146, 0.14422518748337146  diff=1.029720472481288e-5
step=21   xold,xnew= 0.144232026789855, 0.144232026789855  diff=6.839306483541474e-6
step=22   xold,xnew= 0.14423656949526248, 0.14423656949526248  diff=4.542705407478698e-6
step=23   xold,xnew= 0.1442395868306324, 0.1442395868306324  diff=3.0173353699203265e-6
step=24   xold,xnew= 0.14424159101149894, 0.14424159101149894  diff=2.0041808665438587e-6
step=25   xold,xnew= 0.14424292224153484, 0.14424292224153484  diff=1.3312300359002993e-6
step=26   xold,xnew= 0.14424380648367677, 0.14424380648367677  diff=8.842421419241742e-7
step=27   xold,xnew= 0.14424439382502469, 0.14424439382502469  diff=5.873413479184908e-7
step=28   xold,xnew= 0.14424478395631365, 0.14424478395631365  diff=3.901312889686981e-7
step=29   xold,xnew= 0.1442450430945794, 0.1442450430945794  diff=2.591382657346486e-7
step=30   xold,xnew= 0.14424521522303968, 0.14424521522303968  diff=1.7212846029068807e-7
step=31   xold,xnew= 0.14424532955668923, 0.14424532955668923  diff=1.1433364954793745e-7
step=32   xold,xnew= 0.14424540550106424, 0.14424540550106424  diff=7.594437501090567e-8
step=33   xold,xnew= 0.14424545594596447, 0.14424545594596447  diff=5.0444900229207335e-8
step=34   xold,xnew= 0.14424548945323054, 0.14424548945323054  diff=3.350726607287413e-8
step=35   xold,xnew= 0.14424551170993039, 0.14424551170993039  diff=2.2256699844236394e-8
step=36   xold,xnew= 0.14424552649361114, 0.14424552649361114  diff=1.4783680751806472e-8
step=37   xold,xnew= 0.14424553631345094, 0.14424553631345094  diff=9.819839802416297e-9
step=38   xold,xnew= 0.14424554283613353, 0.14424554283613353  diff=6.522682588494533e-9
step=39   xold,xnew= 0.14424554716872853, 0.14424554716872853  diff=4.332595004985507e-9
step=40   xold,xnew= 0.14424555004659045, 0.14424555004659045  diff=2.877861920991265e-9
step=41   xold,xnew= 0.14424555195816754, 0.14424555195816754  diff=1.9115770910627106e-9
step=42   xold,xnew= 0.14424555322790428, 0.14424555322790428  diff=1.2697367324232545e-9
step=43   xold,xnew= 0.14424555407130804, 0.14424555407130804  diff=8.434037634952318e-10
step=44   xold,xnew= 0.14424555463152652, 0.14424555463152652  diff=5.602184827147028e-10
step=45   xold,xnew= 0.1442455550036433, 0.1442455550036433  diff=3.721167818326876e-10
step=46   xold,xnew= 0.14424555525081637, 0.14424555525081637  diff=2.471730597974897e-10
step=47   xold,xnew= 0.14424555541499742, 0.14424555541499742  diff=1.6418105763804647e-10
step=48   xold,xnew= 0.14424555552405227, 0.14424555552405227  diff=1.090548484850018e-10

Case 1: Newton Method
step:0, [0] error:4.0
step:1, [0.08163265306122448] error:1.237476774620502
step:2, [0.13347761344319795] error:0.18228837756813698
step:3, [0.14394125899657323] error:0.0050104019767109875
step:4, [0.14424531478433483] error:3.964316485305375e-6
step:5, [0.1442455557396417] error:2.4851232183209504e-12
CONVERGENCE
=#

M = 3
ξ = 9
ks = [1, 5, 7]
ns = [2, 7, 9]
guess = 0
println("\nCase 2: Bisection Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
println("\nCase 2: FPI Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
println("\nCase 2: Newton Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# === OUTPUT === #
#=
Case 2: Bisection Method
step=1, x=50.0, f(x)=58.90725484186933
step=2, x=25.0, f(x)=33.816385003885
step=3, x=12.5, f(x)=21.139920716125175
step=4, x=6.25, f(x)=14.555966349032836
step=5, x=3.125, f(x)=10.825656260621749
step=6, x=1.5625, f(x)=8.233759252416519
step=7, x=0.78125, f(x)=5.840387042273367
step=8, x=0.390625, f(x)=3.1720105414344353
step=9, x=0.1953125, f(x)=0.178627711316528
step=10, x=0.09765625, f(x)=-2.7735275370428445
step=11, x=0.146484375, f(x)=-1.082152514829704
step=12, x=0.1708984375, f(x)=-0.40963111156902166
step=13, x=0.18310546875, f(x)=-0.10606828296489468
step=14, x=0.189208984375, f(x)=0.03851733086211517
step=15, x=0.1861572265625, f(x)=-0.03320163924025188
step=16, x=0.18768310546875, f(x)=0.0027994617009916567
step=17, x=0.186920166015625, f(x)=-0.015165456894365548
step=18, x=0.1873016357421875, f(x)=-0.006174118270775963
step=19, x=0.18749237060546875, f(x)=-0.001685112020142654
step=20, x=0.18758773803710938, f(x)=0.0005577284616187228
step=21, x=0.18754005432128906, f(x)=-0.00056355331837743
step=22, x=0.18756389617919922, f(x)=-2.877820108437845e-6
step=23, x=0.1875758171081543, f(x)=0.00027743397195401087
step=24, x=0.18756985664367676, f(x)=0.00013728023883174956
step=25, x=0.187566876411438, f(x)=6.720175010244134e-5
step=26, x=0.1875653862953186, f(x)=3.2162100183086295e-5
step=27, x=0.1875646412372589, f(x)=1.4642173834289451e-5
step=28, x=0.18756426870822906, f(x)=5.882185313055288e-6
x=0.18756408244371414, f(x)=1.5021847143970035e-6, n=29 iterations reached.

Case 2: FPI Method
step=0   xold,xnew= 0, 30  diff=1.0e-9
step=1   xold,xnew= 0.0891089108910891, 0.0891089108910891  diff=0.0891089108910891
step=2   xold,xnew= 0.12012280599729182, 0.12012280599729182  diff=0.031013895106202713
step=3   xold,xnew= 0.13920931659365218, 0.13920931659365218  diff=0.019086510596360365
step=4   xold,xnew= 0.1520987790860459, 0.1520987790860459  diff=0.012889462492393727
step=5   xold,xnew= 0.16119803879061817, 0.16119803879061817  diff=0.00909925970457226
step=6   xold,xnew= 0.16778919560270844, 0.16778919560270844  diff=0.006591156812090271
step=7   xold,xnew= 0.17264279107576283, 0.17264279107576283  diff=0.0048535954730543895
step=8   xold,xnew= 0.17625685569671531, 0.17625685569671531  diff=0.003614064620952484
step=9   xold,xnew= 0.17896900389086676, 0.17896900389086676  diff=0.0027121481941514425
step=10   xold,xnew= 0.18101574123330344, 0.18101574123330344  diff=0.002046737342436683
step=11   xold,xnew= 0.18256665329442284, 0.18256665329442284  diff=0.0015509120611194027
step=12   xold,xnew= 0.18374541574989567, 0.18374541574989567  diff=0.001178762455472826
step=13   xold,xnew= 0.18464335368853557, 0.18464335368853557  diff=0.0008979379386399033
step=14   xold,xnew= 0.18532853179755365, 0.18532853179755365  diff=0.0006851781090180808
step=15   xold,xnew= 0.18585203260021377, 0.18585203260021377  diff=0.0005235008026601151
step=16   xold,xnew= 0.1862523949809608, 0.1862523949809608  diff=0.00040036238074703245
step=17   xold,xnew= 0.18655880998340632, 0.18655880998340632  diff=0.00030641500244552033
step=18   xold,xnew= 0.186793454946439, 0.186793454946439  diff=0.00023464496303268745
step=19   xold,xnew= 0.18697321740190648, 0.18697321740190648  diff=0.00017976245546746927
step=20   xold,xnew= 0.18711097935200452, 0.18711097935200452  diff=0.00013776195009804204
step=21   xold,xnew= 0.18721658048616724, 0.18721658048616724  diff=0.00010560113416271943
step=22   xold,xnew= 0.18729754436671314, 0.18729754436671314  diff=8.096388054590342e-5
step=23   xold,xnew= 0.18735962811640447, 0.18735962811640447  diff=6.208374969132735e-5
step=24   xold,xnew= 0.18740723979689522, 0.18740723979689522  diff=4.761168049075004e-5
step=25   xold,xnew= 0.1874437560823199, 0.1874437560823199  diff=3.6516285424670336e-5
step=26   xold,xnew= 0.18747176449068192, 0.18747176449068192  diff=2.8008408362034665e-5
step=27   xold,xnew= 0.1874932483506204, 0.1874932483506204  diff=2.148385993847035e-5
step=28   xold,xnew= 0.1875097281947552, 0.1875097281947552  diff=1.6479844134803523e-5
step=29   xold,xnew= 0.1875223699346845, 0.1875223699346845  diff=1.2641739929292184e-5
step=30   xold,xnew= 0.18753206767491062, 0.18753206767491062  diff=9.697740226133345e-6
step=31   xold,xnew= 0.18753950714270076, 0.18753950714270076  diff=7.4394677901379325e-6
step=32   xold,xnew= 0.1875452142894395, 0.1875452142894395  diff=5.7071467387259656e-6
step=33   xold,xnew= 0.18754959254086315, 0.18754959254086315  diff=4.3782514236589964e-6
step=34   xold,xnew= 0.18755295135365388, 0.18755295135365388  diff=3.3588127907324683e-6
step=35   xold,xnew= 0.18755552811110043, 0.18755552811110043  diff=2.5767574465540566e-6
step=36   xold,xnew= 0.18755750491371548, 0.18755750491371548  diff=1.976802615044626e-6
step=37   xold,xnew= 0.18755902145636372, 0.18755902145636372  diff=1.516542648238861e-6
step=38   xold,xnew= 0.18756018490480303, 0.18756018490480303  diff=1.1634484393119315e-6
step=39   xold,xnew= 0.18756107747127526, 0.18756107747127526  diff=8.925664722325699e-7
step=40   xold,xnew= 0.18756176222546592, 0.18756176222546592  diff=6.84754190655168e-7
step=41   xold,xnew= 0.18756228755203522, 0.18756228755203522  diff=5.253265693005993e-7
step=42   xold,xnew= 0.18756269057003874, 0.18756269057003874  diff=4.030180035252684e-7
step=43   xold,xnew= 0.18756299975605586, 0.18756299975605586  diff=3.091860171222649e-7
step=44   xold,xnew= 0.18756323695649318, 0.18756323695649318  diff=2.3720043731967344e-7
step=45   xold,xnew= 0.18756341893131903, 0.18756341893131903  diff=1.8197482584692004e-7
step=46   xold,xnew= 0.18756355853834722, 0.18756355853834722  diff=1.39607028187827e-7
step=47   xold,xnew= 0.18756366564177301, 0.18756366564177301  diff=1.0710342579489662e-7
step=48   xold,xnew= 0.18756374780916948, 0.18756374780916948  diff=8.216739647015636e-8
step=49   xold,xnew= 0.1875638108462012, 0.1875638108462012  diff=6.303703170562613e-8
step=50   xold,xnew= 0.18756385920684046, 0.18756385920684046  diff=4.836063927093903e-8
step=51   xold,xnew= 0.187563896308074, 0.187563896308074  diff=3.7101233529845956e-8
step=52   xold,xnew= 0.18756392477133757, 0.18756392477133757  diff=2.8463263579414644e-8
step=53   xold,xnew= 0.18756394660773856, 0.18756394660773856  diff=2.183640099295836e-8
step=54   xold,xnew= 0.1875639633601545, 0.1875639633601545  diff=1.6752415926690745e-8
step=55   xold,xnew= 0.18756397621224707, 0.18756397621224707  diff=1.2852092584614283e-8
step=56   xold,xnew= 0.18756398607209585, 0.18756398607209585  diff=9.859848770776836e-9
step=57   xold,xnew= 0.18756399363635973, 0.18756399363635973  diff=7.564263881931765e-9
step=58   xold,xnew= 0.18756399943950028, 0.18756399943950028  diff=5.803140551430275e-9
step=59   xold,xnew= 0.187564003891545, 0.187564003891545  diff=4.452044732872196e-9
step=60   xold,xnew= 0.18756400730705805, 0.18756400730705805  diff=3.41551303906229e-9
step=61   xold,xnew= 0.18756400992736627, 0.18756400992736627  diff=2.6203082204023787e-9
step=62   xold,xnew= 0.18756401193761044, 0.18756401193761044  diff=2.0102441655733827e-9
step=63   xold,xnew= 0.18756401347982662, 0.18756401347982662  diff=1.5422161880884744e-9
step=64   xold,xnew= 0.18756401466298173, 0.18756401466298173  diff=1.183155101669442e-9
step=65   xold,xnew= 0.18756401557067293, 0.18756401557067293  diff=9.076911999805759e-10
step=66   xold,xnew= 0.18756401626703414, 0.18756401626703414  diff=6.963612186883239e-10
step=67   xold,xnew= 0.1875640168012675, 0.1875640168012675  diff=5.34233352000868e-10
step=68   xold,xnew= 0.18756401721111984, 0.18756401721111984  diff=4.0985234650570135e-10
step=69   xold,xnew= 0.18756401752554977, 0.18756401752554977  diff=3.1442992654007185e-10
step=70   xold,xnew= 0.18756401776677367, 0.18756401776677367  diff=2.412239019644602e-10
step=71   xold,xnew= 0.18756401795183542, 0.18756401795183542  diff=1.850617437071378e-10
step=72   xold,xnew= 0.18756401809381085, 0.18756401809381085  diff=1.4197543141136748e-10
step=73   xold,xnew= 0.1875640182027314, 0.1875640182027314  diff=1.0892053925459777e-10

Case 2: Newton Method
step:0, [0] error:9.0
step:1, [0.0891089108910891] error:3.132403405726474
step:2, [0.16140051853273918] error:0.660218450329376
step:3, [0.18577150615662155] error:0.042347486474454854
step:4, [0.18755569290051247] error:0.000195781877604162
step:5, [0.18756401838216843] error:4.2202525918355605e-9
CONVERGENCE
=#

M = 3
ξ = 2
ks = [1, 2, 6]
ns = [2, 3, 1]
guess = 0
println("\nCase 3: Bisection Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, bisect_type, guess)
println("\nCase 3: FPI Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, fpi_type, guess)
println("\nCase 3: Newton Method")
concentrate_unbound_lignand_solver(M, ξ, ks, ns, newton_type, guess)
# === OUTPUT === #
#=
Case 3: Bisection Method
step=1, x=50.0, f(x)=53.92775908429225
step=2, x=25.0, f(x)=28.857630877108868
step=3, x=12.5, f(x)=16.223309341730392
step=4, x=6.25, f(x)=9.725941682838235
step=5, x=3.125, f(x)=6.1757255003108344
step=6, x=1.5625, f(x)=3.95835392568055
step=7, x=0.78125, f(x)=2.311887099314891
step=8, x=0.390625, f(x)=0.9691468059324517
step=9, x=0.1953125, f(x)=-0.09562513984186083
step=10, x=0.29296875, f(x)=0.491909107379902
step=11, x=0.244140625, f(x)=0.21515204188892323
step=12, x=0.2197265625, f(x)=0.06455045131777037
step=13, x=0.20751953125, f(x)=-0.014262242188606278
step=14, x=0.213623046875, f(x)=0.025452480017329027
step=15, x=0.2105712890625, f(x)=0.0056734720291404805
step=16, x=0.20904541015625, f(x)=-0.004274635940175742
step=17, x=0.209808349609375, f(x)=0.0007043350640718593
step=18, x=0.2094268798828125, f(x)=-0.0017839186599419854
step=19, x=0.20961761474609375, f(x)=-0.0005394841694279506
step=20, x=0.20971298217773438, f(x)=8.250231500017691e-5
step=21, x=0.20966529846191406, f(x)=-0.0002284717053659957
step=22, x=0.20968914031982422, f(x)=-7.297989033694385e-5
step=23, x=0.2097010612487793, f(x)=4.762413465808635e-6
step=24, x=0.20969510078430176, f(x)=-3.410843814255493e-5
step=25, x=0.20969808101654053, f(x)=-1.4672937266091424e-5
step=26, x=0.2096995711326599, f(x)=-4.955243132043208e-6
 -- DONE -- x=0.2097003161907196, f(x)=-9.641014120376212e-8, steps=26

Case 3: FPI Method
step=0   xold,xnew= 0, 30  diff=1.0e-9
step=1   xold,xnew= 0.13333333333333333, 0.13333333333333333  diff=0.13333333333333333
step=2   xold,xnew= 0.17035661048044948, 0.17035661048044948  diff=0.037023277147116146
step=3   xold,xnew= 0.18840093806386501, 0.18840093806386501  diff=0.018044327583415537
step=4   xold,xnew= 0.19792691801074824, 0.19792691801074824  diff=0.009525979946883228
step=5   xold,xnew= 0.20312532988015464, 0.20312532988015464  diff=0.005198411869406394
step=6   xold,xnew= 0.2060085320410449, 0.2060085320410449  diff=0.0028832021608902747
step=7   xold,xnew= 0.20762131184799185, 0.20762131184799185  diff=0.0016127798069469346
step=8   xold,xnew= 0.20852763208610176, 0.20852763208610176  diff=0.0009063202381099089
step=9   xold,xnew= 0.20903825188454384, 0.20903825188454384  diff=0.0005106197984420802
step=10   xold,xnew= 0.2093263450698887, 0.2093263450698887  diff=0.0002880931853448676
step=11   xold,xnew= 0.20948901826351907, 0.20948901826351907  diff=0.0001626731936303638
step=12   xold,xnew= 0.20958091387108949, 0.20958091387108949  diff=9.189560757041915e-5
step=13   xold,xnew= 0.20963283975743016, 0.20963283975743016  diff=5.192588634067574e-5
step=14   xold,xnew= 0.20966218484506646, 0.20966218484506646  diff=2.934508763630106e-5
step=15   xold,xnew= 0.2096787700983602, 0.2096787700983602  diff=1.6585253293732727e-5
step=16   xold,xnew= 0.20968814417927234, 0.20968814417927234  diff=9.374080912144533e-6
step=17   xold,xnew= 0.20969344260106132, 0.20969344260106132  diff=5.298421788979546e-6
step=18   xold,xnew= 0.20969643742096816, 0.20969643742096816  diff=2.994819906843027e-6
step=19   xold,xnew= 0.20969813019304007, 0.20969813019304007  diff=1.6927720719039918e-6
step=20   xold,xnew= 0.209699087008731, 0.209699087008731  diff=9.568156909278702e-7
step=21   xold,xnew= 0.20969962783681992, 0.20969962783681992  diff=5.408280889240835e-7
step=22   xold,xnew= 0.2096999335335818, 0.2096999335335818  diff=3.0569676187641726e-7
step=23   xold,xnew= 0.2097001063252526, 0.2097001063252526  diff=1.727916708049726e-7
step=24   xold,xnew= 0.20970020399385583, 0.20970020399385583  diff=9.766860323279225e-8
step=25   xold,xnew= 0.20970025919998106, 0.20970025919998106  diff=5.520612522680324e-8
step=26   xold,xnew= 0.20970029040465313, 0.20970029040465313  diff=3.1204672068518846e-8
step=27   xold,xnew= 0.20970030804276188, 0.20970030804276188  diff=1.7638108751327763e-8
step=28   xold,xnew= 0.2097003180125156, 0.2097003180125156  diff=9.969753717031793e-9
step=29   xold,xnew= 0.2097003236478132, 0.2097003236478132  diff=5.6352976174345315e-9
step=30   xold,xnew= 0.2097003268331055, 0.2097003268331055  diff=3.1852922799391337e-9
step=31   xold,xnew= 0.20970032863355825, 0.20970032863355825  diff=1.8004527591042319e-9
step=32   xold,xnew= 0.2097003296512451, 0.2097003296512451  diff=1.017686851456645e-9
step=33   xold,xnew= 0.20970033022648185, 0.20970033022648185  diff=5.752367471245634e-10
step=34   xold,xnew= 0.20970033055162832, 0.20970033055162832  diff=3.2514646530756863e-10
step=35   xold,xnew= 0.20970033073541394, 0.20970033073541394  diff=1.8378562560705802e-10
step=36   xold,xnew= 0.2097003308392968, 0.2097003308392968  diff=1.0388284676920989e-10

Case 3: Newton Method
step:0, [0] error:2.0
step:1, [0.13333333333333333] error:0.5553491572067424
step:2, [0.20148639717757902] error:0.05414522486472473
step:3, [0.20961216350918793] error:0.0005750417885419523
step:4, [0.20970032089839608] error:6.570945787487403e-8
step:5, [0.20970033097435653] error:8.881784197001252e-16
CONVERGENCE
=#
