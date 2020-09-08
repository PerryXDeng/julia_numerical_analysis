#=
Homework #2 Solutions

Contributors:
Owen Miller
Perry Deng
Clarissa Xue
=#
using Plots
using Statistics

# =========== Plotting Function =============

function make_plots(forward_abs_errors::Array{Float64},
                    backward_abs_errors::Array{Float64},
                    meta_info::String, question::String)
    title = string(question, "\n", meta_info)
    num_iterations = size(forward_abs_errors, 1)
    p = plot(1:num_iterations, forward_abs_errors,
         title=title, legend=false)
    ylabel!(p, "|f(x)|")
    xlabel!(p, "n")
    display(p)
    p = plot(1:num_iterations, backward_abs_errors,
         title=title, legend=false)
    ylabel!(p, "|x-7|")
    xlabel!(p, "n")
    display(p)
    log_forward_errors = log.(forward_abs_errors)
    log_errors_n = log_forward_errors[2:num_iterations]
    log_errors_n_minus_one = log_forward_errors[1:(num_iterations-1)]
    delta_y = log_errors_n[2:end] .- log_errors_n[1:end-1]
    delta_x = log_errors_n_minus_one[2:end] .- log_errors_n_minus_one[1:end-1]
    dy_dx = delta_y ./ delta_x
    dy_dx_avg = mean(dy_dx)
    p = plot(log_errors_n_minus_one, log_errors_n,
         title=string(title, ", avg_slope=", dy_dx_avg), legend=false)
    ylabel!(p, "ln|f(x_n)|")
    xlabel!(p, "ln|f(x_n_minus_one)|")
    display(p)
    p = plot(2:num_iterations-1, dy_dx,
         title=title, legend=false)
    ylabel!(p, "d ln|f(x_n)| / d ln|f(x_n_minus_one)|")
    xlabel!(p, "n")
    display(p)
end

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

function secant_method(f, x0::Float64, x1::Float64, max_iter::Int,
                       threshold::Float64, question::String,
                       visualize::Bool=true)
    s = string("Secant Method x0=", x0, " x1=", x1)
    println(s)
    x2::Float64 = 0.0
    forward_abs_errors::Array{Float64} = Float64[] # plot 1, plot 3
    backward_abs_errors::Array{Float64} = Float64[] # plot 2
    for n in 0:max_iter
        fx1::Float64 = f(x1)
        denominator::Float64 = fx1 - f(x0)
        if denominator == 0.0
            info = string("Singularity, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)
            end
            return 1
        end
        x2 = x1 - fx1 * (x1 - x0) / denominator
        x0, x1 = x1, x2
        forward_abs_error = abs(f(x2))
        push!(forward_abs_errors, forward_abs_error)
        backward_abs_error = abs(x2 - 7.0)
        push!(backward_abs_errors, backward_abs_error)
        if forward_abs_error < threshold
            info = string(":) Converged, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)

            end
            return 0
        end
    end
    info = string("Max iter reached, Root: ", x2)
    println(info)
    s = string(s, "\n", info)
    if visualize
        make_plots(forward_abs_errors, backward_abs_errors, s, question)
    end
    return 1
end

function false_position(f, x0::Float64, x1::Float64, max_iter::Int,
                        threshold::Float64, question::String,
                        visualize::Bool=false)
    s = string("False Position Method x0=", x0, " x1=", x1)
    println(s)
    forward_abs_errors::Array{Float64} = Float64[] # plot 1, plot 3
    backward_abs_errors::Array{Float64} = Float64[] # plot 2
    if f(x0) * f(x1) >= 0.0
        info = string("Bad x0 & x1: ", x0, " & ", x1)
        println(info)
        s = string(s, "\n", info)
        if visualize
            make_plots(forward_abs_errors, backward_abs_errors, s, question)
        end
        return 1
    end
    x2 = x0
    for n in 0:max_iter
        fx0::Float64 = f(x0)
        fx1::Float64 = f(x1)
        denominator::Float64 = fx1 - fx0
        if denominator == 0.0
            info = string("Singularity, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)
            end
            return 1
        end
        x2::Float64 = x1 - fx1 * (x1 - x0) / denominator
        fx2::Float64 = f(x2)
        forward_abs_error = abs(f(x2))
        push!(forward_abs_errors, forward_abs_error)
        backward_abs_error = abs(x2 - 7.0)
        push!(backward_abs_errors, backward_abs_error)
        if forward_abs_error < threshold
            info = string(":) Converged, Root: ", x2)
            println(info)
            s = string(s, "\n", info)
            if visualize
                make_plots(forward_abs_errors, backward_abs_errors, s, question)
            end
            return 0
        elseif f(x2) * f(x0) < 0
            x1 = x2
        else
            x0 = x2
        end
    end
    info = string("Max iter reached, Root: ", x2)
    println(info)
    s = string(s, "\n", info)
    if visualize
        make_plots(forward_abs_errors, backward_abs_errors, s, question)
    end
    return 1
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



# ============= Question 1 =============
function f(x)
    return x^3 - 2
end

println(bisection(f, 0.0, 2.0, 10.0^-8))
#=
x=1.2599210515618324, f(x)=7.938398383089407e-9, steps=27
=#






# ============= Question 2 =============
# f(x) = 2x^3 -8x - 1
# according to wolfram alpha, for f(x) = 0:
#   x_1 \approx -1.9343
#   x_2 \approx -0.12549
#   x_3 \approx 2.0598
function q2_f(x::Float64)::Float64
    return 2x^3 - 8x - 1
end

function q2_f_prime(x::Float64)::Float64
    return 6x^2 - 8
end
-1/q2_f_prime(-1.9343) # -0.06920846788803289
-1/q2_f_prime(2.0598) # -0.057284739199286644

function g_static_alpha(x::Float64, alpha::Float64)::Float64
    return x + alpha * q2_f(x)
end

function g_1(x::Float64)::Float64
    # alpha = - 1 / f'(x=-1.9343)
    # g(x) = x + alpha * f(x)
    return g_static_alpha(x, -0.06920846788803289)
end
fpi(g_1, -2.0)
#=
CONVERGENCE! x=-1.9343
Output:
step=1   xold,xnew= -2.0, -1.930791532111967  diff=0.06920846788803292
step=2   xold,xnew= -1.930791532111967, -1.9342879946383285  diff=-0.003496462526361377
step=3   xold,xnew= -1.9342879946383285, -1.9342978756539178  diff=-9.881015589296993e-6
step=4   xold,xnew= -1.9342978756539178, -1.9342978757660596  diff=-1.1214185136054766e-10
step=5   xold,xnew= -1.9342978757660596, -1.93429787576606  diff=-4.440892098500626e-16
=#
fpi(g_1, -10.0)
#=
DIVERGENCE!
Output:
step=1   xold,xnew= -10.0, 122.9494668129112  diff=132.9494668129112
step=2   xold,xnew= 122.9494668129112, -257067.0804771323  diff=-257190.0299439452
step=3   xold,xnew= -257067.0804771323, 2.351411438540143e15  diff=2.35141143879721e15
step=4   xold,xnew= 2.351411438540143e15, -1.7995963775364438e45  diff=-1.7995963775364438e45
=#
fpi(g_1, 0.0)
#=
CONVERGENCE! x=2.0598
Output:
step=1   xold,xnew= 0.0, 0.06920846788803289  diff=0.06920846788803289
step=2   xold,xnew= 0.06920846788803289, 0.17668954739682102  diff=0.10748107950878813
step=3   xold,xnew= 0.17668954739682102, 0.34296179578213193  diff=0.1662722483853109
step=4   xold,xnew= 0.34296179578213193, 0.5964733906983445  diff=0.25351159491621256
step=5   xold,xnew= 0.5964733906983445, 0.9665559777510444  diff=0.3700825870526999
step=6   xold,xnew= 0.9665559777510444, 1.4459267525635433  diff=0.4793707748124989
step=7   xold,xnew= 1.4459267525635433, 1.8972631199357948  diff=0.45133636737225147
step=8   xold,xnew= 1.8972631199357948, 2.071620042903614  diff=0.17435692296781902
step=9   xold,xnew= 2.071620042903614, 2.057210245175736  diff=-0.014409797727877649
step=10   xold,xnew= 2.057210245175736, 2.060323617980693  diff=0.0031133728049570486
step=11   xold,xnew= 2.060323617980693, 2.059681073951929  diff=-0.0006425440287642914
step=12   xold,xnew= 2.059681073951929, 2.0598150409190765  diff=0.00013396696714762513
step=13   xold,xnew= 2.0598150409190765, 2.0597871678084103  diff=-2.787311066621001e-5
step=14   xold,xnew= 2.0597871678084103, 2.059792969605372  diff=5.801796961524275e-6
step=15   xold,xnew= 2.059792969605372, 2.0597917620689943  diff=-1.2075363775565506e-6
step=16   xold,xnew= 2.0597917620689943, 2.059792013400029  diff=2.5133103465080353e-7
step=17   xold,xnew= 2.059792013400029, 2.059791961089356  diff=-5.231067312649884e-8
step=18   xold,xnew= 2.059791961089356, 2.059791971977023  diff=1.0887667212244878e-8
step=19   xold,xnew= 2.059791971977023, 2.059791969710922  diff=-2.2661010845581586e-9
step=20   xold,xnew= 2.059791969710922, 2.059791970182576  diff=4.716542711946659e-10
step=21   xold,xnew= 2.059791970182576, 2.0597919700844085  diff=-9.816769619419574e-11
=#
fpi(g_1, 4.0)
#=
CONVERGENCE! x=-1.9343
Output:
step=1   xold,xnew= 4.0, -2.574804449363125  diff=-6.574804449363125
step=2   xold,xnew= -2.574804449363125, -1.5684092361652102  diff=1.0063952131979146
step=3   xold,xnew= -1.5684092361652102, -1.8335462193328433  diff=-0.2651369831676331
step=4   xold,xnew= -1.8335462193328433, -1.9262857034085419  diff=-0.09273948407569854
step=5   xold,xnew= -1.9262857034085419, -1.9342463570660187  diff=-0.00796065365747678
step=6   xold,xnew= -1.9342463570660187, -1.9342978734583898  diff=-5.1516392371109276e-5
step=7   xold,xnew= -1.9342978734583898, -1.9342978757660523  diff=-2.3076625055296063e-9
step=8   xold,xnew= -1.9342978757660523, -1.93429787576606  diff=-7.771561172376096e-1
=#
fpi(g_1, 2000.0)
#=
DIVERGENCE!
likely due to large f(x) value at initial guess
Output:
step=1   xold,xnew= 2000.0, -1.1073323788038316e9  diff=-1.1073343788038316e9
step=2   xold,xnew= -1.1073323788038316e9, 1.8794173330572365e26  diff=1.8794173330572365e26
=#

function g_2(x::Float64)::Float64
    # alpha = - 1 / f'(x=-1.9343)
    return g_static_alpha(x, 2.0598)
end
fpi(g_2, -2.0)
#=
DIVERGENCE!
Output:
step=1   xold,xnew= -2.0, -4.0598  diff=-2.0598
step=2   xold,xnew= -4.0598, -214.8775547375002  diff=-210.8177547375002
step=3   xold,xnew= -214.8775547375002, -4.086889437292179e7  diff=-4.086867949536705e7
step=4   xold,xnew= -4.086889437292179e7, -2.8121191388446592e23  diff=-2.812119138844659e23
=#
fpi(g_2, -10.0)
#=
DIVERGENCE!
Output:
step=1   xold,xnew= -10.0, -3966.8758000000003  diff=-3956.8758000000003
step=2   xold,xnew= -3966.8758000000003, -2.5715842399546158e11  diff=-2.571584200285858e11
step=3   xold,xnew= -2.5715842399546158e11, -7.005793240620509e34  diff=-7.005793240620509e34
=#
fpi(g_2, 0.0)
#=
DIVERGENCE!
Output:
step=1   xold,xnew= 0.0, -2.0598  diff=-2.0598
step=2   xold,xnew= -2.0598, -6.179688731196172  diff=-4.119888731196172
step=3   xold,xnew= -6.179688731196172, -878.6063908870058  diff=-872.4267021558096
step=4   xold,xnew= -878.6063908870058, -2.794061813340066e9  diff=-2.794060934733675e9
step=5   xold,xnew= -2.794061813340066e9, -8.985930989577029e28  diff=-8.985930989577029e28
=#
fpi(g_2, 4.0)
#=
DIVERGENCE!
Output:
step=1   xold,xnew= 4.0, 199.681  diff=195.681
step=2   xold,xnew= 199.681, 3.2796260304849505e7  diff=3.2796060623849504e7
step=3   xold,xnew= 3.2796260304849505e7, 1.4532088155969493e23  diff=1.453208815596949e23
=#
fpi(g_2, 2000.0)
#=
DIVERGENCE!
Output:
step=1   xold,xnew= 2000.0, 3.29567690411402e10  diff=3.29567670411402e10
step=2   xold,xnew= 3.29567690411402e10, 1.4746499314525168e32  diff=1.4746499314525168e32
=#
function q2_f_square_prime(x::Float64)::Float64
    # d(f(x)^2)/dx
    # 0 in extremas of f(x)^2, including local maximas
    # also 0 in global minimas of f(x)^2, including when f(x) = 0
    return 4 * (3 * x^2 - 4) * (2 * x^3 - 8 * x - 1)
end

function g_gradient_descent(x::Float64, alpha::Float64)::Float64
    # gradient descent
    return x - alpha * q2_f_square_prime(x)
end

function g_3(x::Float64)::Float64
    return g_gradient_descent(x, 0.001)
end
fpi(g_3, -2.0)
#=
CONVERGENCE! x=-1.9342978759447438
Output:
step=0   xold,xnew= -2.0, 420.69  diff=1.0e-9
step=1   xold,xnew= -1.968, -1.968  diff=0.03200000000000003
step=2   xold,xnew= -1.9527550761230663, -1.9527550761230663  diff=0.015244923876933658
step=3   xold,xnew= -1.9447006372137023, -1.9447006372137023  diff=0.00805443890936397
step=4   xold,xnew= -1.940247210336428, -1.940247210336428  diff=0.004453426877274236
step=5   xold,xnew= -1.9377272903855418, -1.9377272903855418  diff=0.0025199199508862957
step=6   xold,xnew= -1.9362834737024384, -1.9362834737024384  diff=0.0014438166831034493
step=7   xold,xnew= -1.935450414244587, -1.935450414244587  diff=0.0008330594578513928
step=8   xold,xnew= -1.9349678340206349, -1.9349678340206349  diff=0.0004825802239520893
step=9   xold,xnew= -1.9346876410655998, -1.9346876410655998  diff=0.000280192955035119
step=10   xold,xnew= -1.9345247416395368, -1.9345247416395368  diff=0.00016289942606295682
step=11   xold,xnew= -1.9344299619841177, -1.9344299619841177  diff=9.477965541915268e-5
step=12   xold,xnew= -1.9343747918185037, -1.9343747918185037  diff=5.5170165613960265e-5
step=13   xold,xnew= -1.934342669562529, -1.934342669562529  diff=3.212225597470919e-5
step=14   xold,xnew= -1.9343239638907463, -1.9343239638907463  diff=1.8705671782637268e-5
step=15   xold,xnew= -1.9343130701076483, -1.9343130701076483  diff=1.0893783098042675e-5
step=16   xold,xnew= -1.9343067254769535, -1.9343067254769535  diff=6.344630694821163e-6
step=17   xold,xnew= -1.9343030302007451, -1.9343030302007451  diff=3.6952762083597435e-6
step=18   xold,xnew= -1.9343008779394466, -1.9343008779394466  diff=2.1522612985158673e-6
step=19   xold,xnew= -1.9342996243726684, -1.9342996243726684  diff=1.253566778158799e-6
step=20   xold,xnew= -1.9342988942387793, -1.9342988942387793  diff=7.301338891263498e-7
step=21   xold,xnew= -1.9342984689743798, -1.9342984689743798  diff=4.2526439947287997e-7
step=22   xold,xnew= -1.93429822127985, -1.93429822127985  diff=2.476945297757993e-7
step=23   xold,xnew= -1.9342980770104226, -1.9342980770104226  diff=1.4426942751555316e-7
step=24   xold,xnew= -1.934297992980784, -1.934297992980784  diff=8.40296385806738e-8
step=25   xold,xnew= -1.934297944037754, -1.934297944037754  diff=4.8943030073189675e-8
step=26   xold,xnew= -1.9342979155308986, -1.9342979155308986  diff=2.8506855320742375e-8
step=27   xold,xnew= -1.9342978989270858, -1.9342978989270858  diff=1.6603812769133697e-8
step=28   xold,xnew= -1.9342978892561973, -1.9342978892561973  diff=9.670888534429878e-9
step=29   xold,xnew= -1.9342978836233893, -1.9342978836233893  diff=5.632807997812961e-9
step=30   xold,xnew= -1.9342978803425608, -1.9342978803425608  diff=3.2808284977647872e-9
step=31   xold,xnew= -1.9342978784316425, -1.9342978784316425  diff=1.9109183124754736e-9
step=32   xold,xnew= -1.9342978773186283, -1.9342978773186283  diff=1.1130141253090642e-9
step=33   xold,xnew= -1.9342978766703534, -1.9342978766703534  diff=6.482749892455786e-10
step=34   xold,xnew= -1.9342978762927656, -1.9342978762927656  diff=3.7758773885343544e-10
step=35   xold,xnew= -1.9342978760728398, -1.9342978760728398  diff=2.1992585530483666e-10
step=36   xold,xnew= -1.9342978759447438, -1.9342978759447438  diff=1.2809597826901609e-10
=#
fpi(g_3, -10.0)
#=
DIVERGENCE!
large gradient when far away from root leads to divergence
Output:
step=0   xold,xnew= -10.0, 420.69  diff=1.0e-9
step=1   xold,xnew= 2264.464, 2264.464  diff=2274.464
step=2   xold,xnew= -1.4290186399740538e15, -1.4290186399740538e15  diff=-1.4290186399763182e15
=#
fpi(g_3, 0.0, max_iter=500)
#=
CONVERGENCE! x=-0.12549
Output:
step=0   xold,xnew= 0.0, 420.69  diff=1.0e-9
step=1   xold,xnew= -0.016, -0.016  diff=-0.016
step=2   xold,xnew= -0.02994945226283418, -0.02994945226283418  diff=-0.013949452262834179
step=3   xold,xnew= -0.042108596718717944, -0.042108596718717944  diff=-0.012159144455883765
step=4   xold,xnew= -0.05270697256953677, -0.05270697256953677  diff=-0.010598375850818825
step=5   xold,xnew= -0.061945875943368464, -0.061945875943368464  diff=-0.009238903373831694
step=6   xold,xnew= -0.07000116059409554, -0.07000116059409554  diff=-0.00805528465072708
step=7   xold,xnew= -0.07702607593190311, -0.07702607593190311  diff=-0.007024915337807566
step=8   xold,xnew= -0.08315397255536587, -0.08315397255536587  diff=-0.006127896623462761
step=9   xold,xnew= -0.0885007904607267, -0.0885007904607267  diff=-0.005346817905360832
step=10   xold,xnew= -0.09316729640743346, -0.09316729640743346  diff=-0.004666505946706756
step=11   xold,xnew= -0.09724106657708027, -0.09724106657708027  diff=-0.004073770169646815
step=12   xold,xnew= -0.10079822678239675, -0.10079822678239675  diff=-0.003557160205316476
step=13   xold,xnew= -0.10390497035423334, -0.10390497035423334  diff=-0.003106743571836587
step=14   xold,xnew= -0.10661887683117077, -0.10661887683117077  diff=-0.002713906476937436
step=15   xold,xnew= -0.10899005480709542, -0.10899005480709542  diff=-0.0023711779759246487
step=16   xold,xnew= -0.11106213104905281, -0.11106213104905281  diff=-0.0020720762419573907
step=17   xold,xnew= -0.11287310605524589, -0.11287310605524589  diff=-0.0018109750061930802
step=18   xold,xnew= -0.11445609403194949, -0.11445609403194949  diff=-0.001582987976703601
step=19   xold,xnew= -0.11583996308002105, -0.11583996308002105  diff=-0.0013838690480715565
step=20   xold,xnew= -0.11704988932779108, -0.11704988932779108  diff=-0.0012099262477700268
step=21   xold,xnew= -0.118107836887287, -0.118107836887287  diff=-0.0010579475594959253
step=22   xold,xnew= -0.11903297386387023, -0.11903297386387023  diff=-0.000925136976583224
step=23   xold,xnew= -0.11984203321184576, -0.11984203321184576  diff=-0.0008090593479755381
step=24   xold,xnew= -0.12054962598542263, -0.12054962598542263  diff=-0.0007075927735768622
step=25   xold,xnew= -0.12116851346562511, -0.12116851346562511  diff=-0.0006188874802024874
step=26   xold,xnew= -0.12170984372820623, -0.12170984372820623  diff=-0.0005413302625811162
step=27   xold,xnew= -0.12218335743486043, -0.12218335743486043  diff=-0.0004735137066542022
step=28   xold,xnew= -0.12259756696134763, -0.12259756696134763  diff=-0.0004142095264871998
step=29   xold,xnew= -0.12295991240489079, -0.12295991240489079  diff=-0.00036234544354316056
step=30   xold,xnew= -0.12327689752490881, -0.12327689752490881  diff=-0.00031698512001801515
step=31   xold,xnew= -0.12355420825333886, -0.12355420825333886  diff=-0.0002773107284300552
step=32   xold,xnew= -0.12379681605286932, -0.12379681605286932  diff=-0.00024260779953046152
step=33   xold,xnew= -0.12400906809433772, -0.12400906809433772  diff=-0.00021225204146839427
step=34   xold,xnew= -0.1241947659607341, -0.1241947659607341  diff=-0.00018569786639638375
step=35   xold,xnew= -0.12435723435825773, -0.12435723435825773  diff=-0.00016246839752362408
step=36   xold,xnew= -0.1244993811192803, -0.1244993811192803  diff=-0.00014214676102257406
step=37   xold,xnew= -0.12462374961328981, -0.12462374961328981  diff=-0.00012436849400951355
step=38   xold,xnew= -0.12473256453605426, -0.12473256453605426  diff=-0.00010881492276444693
step=39   xold,xnew= -0.12482777192107555, -0.12482777192107555  diff=-9.520738502129256e-5
step=40   xold,xnew= -0.12491107410812194, -0.12491107410812194  diff=-8.33021870463857e-5
step=41   xold,xnew= -0.12498396030886959, -0.12498396030886959  diff=-7.288620074764884e-5
step=42   xold,xnew= -0.12504773332743938, -0.12504773332743938  diff=-6.377301856978768e-5
step=43   xold,xnew= -0.12510353292216755, -0.12510353292216755  diff=-5.5799594728178636e-5
step=44   xold,xnew= -0.12515235623283022, -0.12515235623283022  diff=-4.882331066266565e-5
step=45   xold,xnew= -0.1251950756434957, -0.1251950756434957  diff=-4.271941066547047e-5
step=46   xold,xnew= -0.12523245440412353, -0.12523245440412353  diff=-3.7378760627843066e-5
step=47   xold,xnew= -0.12526516029303897, -0.12526516029303897  diff=-3.2705888915435244e-5
step=48   xold,xnew= -0.12529377756668403, -0.12529377756668403  diff=-2.8617273645059615e-5
step=49   xold,xnew= -0.1253188174118928, -0.1253188174118928  diff=-2.5039845208768563e-5
step=50   xold,xnew= -0.12534072708876126, -0.12534072708876126  diff=-2.1909676868464567e-5
step=51   xold,xnew= -0.12535989792846555, -0.12535989792846555  diff=-1.917083970429112e-5
step=52   xold,xnew= -0.1253766723296787, -0.1253766723296787  diff=-1.677440121314322e-5
step=53   xold,xnew= -0.12539134987915806, -0.12539134987915806  diff=-1.467754947936939e-5
step=54   xold,xnew= -0.12540419270628558, -0.12540419270628558  diff=-1.2842827127518275e-5
step=55   xold,xnew= -0.12541543016754825, -0.12541543016754825  diff=-1.1237461262664228e-5
step=56   xold,xnew= -0.12542526294489395, -0.12542526294489395  diff=-9.832777345702981e-6
step=57   xold,xnew= -0.12543386663136305, -0.12543386663136305  diff=-8.603686469099703e-6
step=58   xold,xnew= -0.12544139486818984, -0.12544139486818984  diff=-7.528236826787049e-6
step=59   xold,xnew= -0.12544798208951852, -0.12544798208951852  diff=-6.587221328679949e-6
step=60   xold,xnew= -0.12545374592384229, -0.12545374592384229  diff=-5.763834323768702e-6
step=61   xold,xnew= -0.12545878929512047, -0.12545878929512047  diff=-5.043371278184994e-6
step=62   xold,xnew= -0.12546320226114901, -0.12546320226114901  diff=-4.41296602854413e-6
step=63   xold,xnew= -0.12546706362205406, -0.12546706362205406  diff=-3.861360905049738e-6
step=64   xold,xnew= -0.12547044232766308, -0.12547044232766308  diff=-3.3787056090139878e-6
step=65   xold,xnew= -0.1254733987089094, -0.1254733987089094  diff=-2.956381246310702e-6
step=66   xold,xnew= -0.12547598555527736, -0.12547598555527736  diff=-2.5868463679745712e-6
step=67   xold,xnew= -0.12547824905754168, -0.12547824905754168  diff=-2.26350226431582e-6
step=68   xold,xnew= -0.12548022963264585, -0.12548022963264585  diff=-1.980575104171267e-6
step=69   xold,xnew= -0.12548196264545725, -0.12548196264545725  diff=-1.7330128113945964e-6
step=70   xold,xnew= -0.12548347904029292, -0.12548347904029292  diff=-1.5163948356711199e-6
step=71   xold,xnew= -0.12548480589349761, -0.12548480589349761  diff=-1.3268532046972759e-6
step=72   xold,xnew= -0.12548596689694455, -0.12548596689694455  diff=-1.1610034469367037e-6
step=73   xold,xnew= -0.12548698278109502, -0.12548698278109502  diff=-1.0158841504681604e-6
step=74   xold,xnew= -0.12548787168517328, -0.12548787168517328  diff=-8.889040782611435e-7
step=75   xold,xnew= -0.12548864948106805, -0.12548864948106805  diff=-7.777958947741137e-7
step=76   xold,xnew= -0.12548933005674529, -0.12548933005674529  diff=-6.805756772310101e-7
step=77   xold,xnew= -0.125489925564234, -0.125489925564234  diff=-5.955074887098455e-7
step=78   xold,xnew= -0.12549044663661352, -0.12549044663661352  diff=-5.210723795223693e-7
step=79   xold,xnew= -0.12549090257787743, -0.12549090257787743  diff=-4.559412639104643e-7
step=80   xold,xnew= -0.12549130152906462, -0.12549130152906462  diff=-3.989511871971274e-7
step=81   xold,xnew= -0.12549165061362424, -0.12549165061362424  diff=-3.490845596199055e-7
step=82   xold,xnew= -0.1254919560646102, -0.1254919560646102  diff=-3.0545098594902953e-7
step=83   xold,xnew= -0.12549222333597682, -0.12549222333597682  diff=-2.6727136662185735e-7
step=84   xold,xnew= -0.12549245719996277, -0.12549245719996277  diff=-2.3386398595448554e-7
step=85   xold,xnew= -0.12549266183230232, -0.12549266183230232  diff=-2.046323395454852e-7
step=86   xold,xnew= -0.12549284088678547, -0.12549284088678547  diff=-1.790544831570262e-7
step=87   xold,xnew= -0.1254929975604985, -0.1254929975604985  diff=-1.5667371303096367e-7
step=88   xold,xnew= -0.12549313465090953, -0.12549313465090953  diff=-1.3709041102316633e-7
step=89   xold,xnew= -0.1254932546058187, -0.1254932546058187  diff=-1.1995490917238172e-7
step=90   xold,xnew= -0.12549335956706495, -0.12549335956706495  diff=-1.0496124625003489e-7
step=91   xold,xnew= -0.12549345140876939, -0.12549345140876939  diff=-9.184170443599093e-8
step=92   xold,xnew= -0.12549353177079864, -0.12549353177079864  diff=-8.036202925332248e-8
step=93   xold,xnew= -0.12549360208804536, -0.12549360208804536  diff=-7.031724671979767e-8
step=94   xold,xnew= -0.12549366361604886, -0.12549366361604886  diff=-6.15280034976795e-8
step=95   xold,xnew= -0.12549371745341326, -0.12549371745341326  diff=-5.383736439990017e-8
step=96   xold,xnew= -0.1254937645614235, -0.1254937645614235  diff=-4.710801024265798e-8
step=97   xold,xnew= -0.1254938057812095, -0.1254938057812095  diff=-4.1219786001134295e-8
step=98   xold,xnew= -0.12549384184876475, -0.12549384184876475  diff=-3.6067555247987215e-8
step=99   xold,xnew= -0.1254938734080879, -0.1254938734080879  diff=-3.155932315479504e-8
step=100   xold,xnew= -0.12549390102268151, -0.12549390102268151  diff=-2.761459361177998e-8
step=101   xold,xnew= -0.12549392518561361, -0.12549392518561361  diff=-2.4162932099613954e-8
step=102   xold,xnew= -0.1254939463283216, -0.1254939463283216  diff=-2.1142707973265118e-8
step=103   xold,xnew= -0.12549396482831568, -0.12549396482831568  diff=-1.84999940922026e-8
step=104   xold,xnew= -0.12549398101591958, -0.12549398101591958  diff=-1.6187603896211655e-8
step=105   xold,xnew= -0.12549399518016843, -0.12549399518016843  diff=-1.4164248857140294e-8
step=106   xold,xnew= -0.12549400757396972, -0.12549400757396972  diff=-1.2393801290810913e-8
step=107   xold,xnew= -0.12549401841861899, -0.12549401841861899  diff=-1.0844649261931849e-8
step=108   xold,xnew= -0.12549402790775108, -0.12549402790775108  diff=-9.489132091200503e-9
step=109   xold,xnew= -0.12549403621079763, -0.12549403621079763  diff=-8.303046555857563e-9
step=110   xold,xnew= -0.12549404347601237, -0.12549404347601237  diff=-7.26521473537467e-9
step=111   xold,xnew= -0.12549404983311815, -0.12549404983311815  diff=-6.357105786225503e-9
step=112   xold,xnew= -0.12549405539562322, -0.12549405539562322  diff=-5.562505067668866e-9
step=113   xold,xnew= -0.12549406026284798, -0.12549406026284798  diff=-4.867224762117317e-9
step=114   xold,xnew= -0.1254940645216984, -0.1254940645216984  diff=-4.258850411220649e-9
step=115   xold,xnew= -0.12549406824821766, -0.12549406824821766  diff=-3.726519259839023e-9
step=116   xold,xnew= -0.125494071508944, -0.125494071508944  diff=-3.2607263555917143e-9
step=117   xold,xnew= -0.1254940743620988, -0.1254940743620988  diff=-2.8531547957566517e-9
step=118   xold,xnew= -0.1254940768586261, -0.1254940768586261  diff=-2.496527290452022e-9
step=119   xold,xnew= -0.12549407904310222, -0.12549407904310222  diff=-2.184476127764512e-9
step=120   xold,xnew= -0.12549408095453174, -0.12549408095453174  diff=-1.9114295146671623e-9
step=121   xold,xnew= -0.12549408262704384, -0.12549408262704384  diff=-1.6725121010363608e-9
step=122   xold,xnew= -0.1254940840905018, -0.1254940840905018  diff=-1.463457965922288e-9
step=123   xold,xnew= -0.12549408537103618, -0.12549408537103618  diff=-1.2805343729827001e-9
step=124   xold,xnew= -0.1254940864915113, -0.1254940864915113  diff=-1.120475129345877e-9
step=125   xold,xnew= -0.12549408747193366, -0.12549408747193366  diff=-9.804223544129798e-10
step=126   xold,xnew= -0.12549408832980902, -0.12549408832980902  diff=-8.57875354087767e-10
step=127   xold,xnew= -0.125494089080455, -0.125494089080455  diff=-7.506459898110052e-10
step=128   xold,xnew= -0.1254940897372747, -0.1254940897372747  diff=-6.56819681976728e-10
step=129   xold,xnew= -0.1254940903119958, -0.1254940903119958  diff=-5.747211040407763e-10
=#
fpi(g_3, 4.0, max_iter=500)
#=
DIVERGENCE!
large gradient when far away from root leads to divergence
Output:
step=0   xold,xnew= 4.0, 420.69  diff=1.0e-9
step=1   xold,xnew= -12.719999999999999, -12.719999999999999  diff=-16.72
step=2   xold,xnew= 7719.24038406471, 7719.24038406471  diff=7731.9603840647105
step=3   xold,xnew= -6.577850976105175e17, -6.577850976105175e17  diff=-6.577850976105252e17
=#
fpi(g_3, 2000.0, max_iter=500)
#=
DIVERGENCE!
large gradient when far away from root leads to divergence
Output:
step=0   xold,xnew= 2000.0, 420.69  diff=1.0e-9
step=1   xold,xnew= -7.67998975950256e14, -7.67998975950256e14  diff=-7.67998975952256e14
=#

# ============= Question 3 =============
thresh = .00001
num_iter = 9999
x0 = -1.0
x1 = 1.0

function q3_f(x::Float64) :: Float64
    omega1 = 1 / (1.0x^2 + 1.5x - 0.25)
    omega2 = 1.0 / (4.0x^2 - 6.4x + 1.5)
    return omega1 - omega2
end

secant_method(q3_f, x0, x1, num_iter, thresh, "Question 3")
#=
OUTPUT
Secant Method x0=-1.0 x1=1.0
:) Converged, Root: -314.5873865107352

EXPLANATION
Solving this by hand numerically yields the roots x ~ 2.4, -0.2 within the specified
domain of -1 ⪯ x ⪯ 1, which are not the roots that the Secant Method found.
By looking at a plot of the equation 𝝮1(x) - 𝝮2(x) = 0, the slope of the secant
line for the points x0 & x1 result in a very steep slope which means the run is
close to 0. More visual explanation in attached document.
=#

# ============= Question 4 =============
k = 1.3807e-23  #= Boltzmann's Constant =#
c = 2.9979e8    #= Speed of light in m/s =#
h = 6.6261e-34 #= Planck's Constant =#
t1 = 500.0
t2 = 1000.0
t3 = 5000.0
min_lambda = 1500e-9
max_lambda = 1e-9
min_nu = 1.9986e14
max_nu = 2.9979e17

function f_nu(nu::Float64, t::Float64)::Float64
    ((8 * pi * h * nu^3)/c^3) * (1/exp((h*nu)/k*t)-1)
end

function f_lambda(lambda::Float64, t::Float64)::Float64
    return ((8 * pi * h * c)/lambda^5) * (1/exp((h*c)/lambda*k*t)-1)
end

function plot_func_over_variable_with_temps(func, var_vals::Array{Float64},
                                            var_name::String, type::String)
    plotted_func(func, t) = 0
    if type == "loglog"
        var_vals = log.(var_vals)
    end
    functions_to_call = []
    functions_label = []
    for t in [t1, t2, t3]
        if type == "linear"
            push!(functions_to_call, (x->func(x, t)))
        elseif type == "semilog"
            push!(functions_to_call, (x->log(func(x, t))))
        elseif type == "loglog"
            push!(functions_to_call, (x->log(func(exp(x), t))))
        else
            println("????")
        end
        push!(functions_label, string("temp=", t))
    end
    title = string("u(", var_name, ",T) vs. ", var_name, ", ", type)
    p = plot(var_vals, functions_to_call,
             title=title, label=permutedims(functions_label))
    display(p)
end

lambda_vars = collect(range(min_lambda, stop=max_lambda, length=42069))
plot_func_over_variable_with_temps(f_lambda, lambda_vars, "lambda", "linear")
plot_func_over_variable_with_temps(f_lambda, lambda_vars, "lambda", "semilog")
plot_func_over_variable_with_temps(f_lambda, lambda_vars, "lambda", "loglog")
nu_vars = collect(range(min_nu, stop=max_nu, length=42069))
plot_func_over_variable_with_temps(f_nu, nu_vars, "nu", "linear")
# following two lines break code
plot_func_over_variable_with_temps(f_nu, nu_vars, "nu", "semilog")
plot_func_over_variable_with_temps(f_nu, nu_vars, "nu", "loglog")

f(x)::Float64 = x^2
g(x) = 2*x
t = 1:100
plot(t,[f,g], titile="tit", label=["f(x)" "g(x)"])
