#====
Numerical Analysis
Workshop #4

Contributors:
Perry Deng
Clarissa Xue
Owen Miller
====#


# Question 1
FACC = 10.0^-10
MAXACC = 10.0^20
max_iter = 50

function fpi(func, xold)
    diff = 10*FACC
    step = 0

    while step < max_iter && abs(diff) > FACC && abs(diff) < MAXACC
        xnew = func(xold)
        diff = xnew-xold
        xold = xnew

        println("step=", step, "   xold,new= ", xold, ", ", xnew, "  diff=", diff)
        step += 1
    end
end

using Plots
plotly()
function fpi_plot(func, xold)
    diff = 10*FACC
    step = 0
    x = []
    y = []
    while step < max_iter && abs(diff) > FACC && abs(diff) < MAXACC
        xnew = func(xold)
        diff = xnew-xold
        xold = xnew
        println("step=", step, "   xold,new= ", xold, ", ", xnew, "  diff=", diff)
        step += 1

        push!(x, step)
        push!(y, abs(diff))
    end
    println(x)
    println(y)
    plot(x, y, yaxis=:log)
    xaxis!("steps")
    yaxis!("log| x_n - x_n+1 |")
end


# =========== Question 2 ============
function g(x)
    return (x^2 - x - 1) + x
end

fpi(g, 1)
#= - Output -
step=0   xold,new= 0, 0  diff=-1
step=1   xold,new= -1, -1  diff=-1
step=2   xold,new= 0, 0  diff=1
step=3   xold,new= -1, -1  diff=-1
step=4   xold,new= 0, 0  diff=1
step=5   xold,new= -1, -1  diff=-1
step=6   xold,new= 0, 0  diff=1
step=7   xold,new= -1, -1  diff=-1
step=8   xold,new= 0, 0  diff=1

Description:
We get a two cycle between -1 and 0.
=#

fpi(g, 1.5)
#= - Output -
step=0   xold,new= 1.25, 1.25  diff=-0.25
step=1   xold,new= 0.5625, 0.5625  diff=-0.6875
step=2   xold,new= -0.68359375, -0.68359375  diff=-1.24609375
step=3   xold,new= -0.5326995849609375, -0.5326995849609375  diff=0.1508941650390625
step=4   xold,new= -0.7162311521824449, -0.7162311521824449  diff=-0.18353156722150743
step=5   xold,new= -0.4870129366434073, -0.4870129366434073  diff=0.22921821553903765
step=6   xold,new= -0.7628183995419646, -0.7628183995419646  diff=-0.27580546289855734
step=7   xold,new= -0.4181080893202356, -0.4181080893202356  diff=0.34471031022172904
step=8   xold,new= -0.8251856256449819, -0.8251856256449819  diff=-0.4070775363247463
step=9   xold,new= -0.3190686832288998, -0.3190686832288998  diff=0.5061169424160821
step=10   xold,new= -0.898195175382576, -0.898195175382576  diff=-0.5791264921536762
step=11   xold,new= -0.19324542691946356, -0.19324542691946356  diff=0.7049497484631124
step=12   xold,new= -0.9626562049747143, -0.9626562049747143  diff=-0.7694107780552507
step=13   xold,new= -0.07329303102368079, -0.07329303102368079  diff=0.8893631739510335
step=14   xold,new= -0.9946281316033617, -0.9946281316033617  diff=-0.9213351005796809
step=15   xold,new= -0.010714879823205625, -0.010714879823205625  diff=0.9839132517801561
step=16   xold,new= -0.9998851913503742, -0.9998851913503742  diff=-0.9891703115271686
step=17   xold,new= -0.0002296041182253683, -0.0002296041182253683  diff=0.9996555872321489
step=18   xold,new= -0.9999999472819489, -0.9999999472819489  diff=-0.9997703431637235
step=19   xold,new= -1.0543609951785271e-7, -1.0543609951785271e-7  diff=0.9999998418458493
step=20   xold,new= -0.9999999999999889, -0.9999999999999889  diff=-0.9999998945638894
step=21   xold,new= -2.220446049250313e-14, -2.220446049250313e-14  diff=0.9999999999999667
step=22   xold,new= -1.0, -1.0  diff=-0.9999999999999778
step=23   xold,new= 0.0, 0.0  diff=1.0
step=24   xold,new= -1.0, -1.0  diff=-1.0
step=25   xold,new= 0.0, 0.0  diff=1.0

Description: Our values appear to randomly bounce around until it also loops
in another 2-cycle between -1 and 0.
=#

fpi(g, 2.0)
#= - Output -
step=0   xold,new= 3.0, 3.0  diff=1.0
step=1   xold,new= 8.0, 8.0  diff=5.0
step=2   xold,new= 63.0, 63.0  diff=55.0
step=3   xold,new= 3968.0, 3968.0  diff=3905.0
step=4   xold,new= 1.5745023e7, 1.5745023e7  diff=1.5741055e7
step=5   xold,new= 2.47905749270528e14, 2.47905749270528e14  diff=2.47905733525505e14
step=6   xold,new= 6.14572605213819e28, 6.14572605213819e28  diff=6.145726052138165e28

Description:
Our initial value of 2.0 cause the values to explode and go beyond our limits.
=#


# ========= Question 3 ===========
function g_alpha(alpha)
    return function g(x)
        return alpha*(x^2 - x - 1) + x
    end
end

fpi(g_alpha(.5), 1)
#= - Output -
step=0   xold,new= 0.5, 0.5  diff=-0.5
step=1   xold,new= -0.125, -0.125  diff=-0.625
step=2   xold,new= -0.5546875, -0.5546875  diff=-0.4296875
step=3   xold,new= -0.623504638671875, -0.623504638671875  diff=-0.068817138671875
step=4   xold,new= -0.6173733021132648, -0.6173733021132648  diff=0.006131336558610201
step=5   xold,new= -0.6181117539755141, -0.6181117539755141  diff=-0.0007384518622493008
step=6   xold,new= -0.6180248067864138, -0.6180248067864138  diff=8.694718910029486e-5
step=7   xold,new= -0.6180350724915149, -0.6180350724915149  diff=-1.0265705101097922e-5
step=8   xold,new= -0.6180338608309615, -0.6180338608309615  diff=1.2116605534462366e-6
step=9   xold,new= -0.6180340038486686, -0.6180340038486686  diff=-1.4301770712155104e-7
step=10   xold,new= -0.6180339869677263, -0.6180339869677263  diff=1.6880942310670832e-8
step=11   xold,new= -0.6180339889602513, -0.6180339889602513  diff=-1.992525033855941e-9
step=12   xold,new= -0.6180339887250657, -0.6180339887250657  diff=2.351856487337045e-10
step=13   xold,new= -0.6180339887528256, -0.6180339887528256  diff=-2.7759905485424952e-11

Description:
It looks like we found the other root!
=#

fpi(g_alpha(.5), 2)
fpi(g_alpha(.5), -1)
#= - Output -
2 :
step=0   xold,new= 2.5, 2.5  diff=0.5
step=1   xold,new= 3.875, 3.875  diff=1.375
step=2   xold,new= 8.9453125, 8.9453125  diff=5.0703125
step=3   xold,new= 43.981964111328125, 43.981964111328125  diff=35.036651611328125
step=4   xold,new= 988.6975656007417, 988.6975656007417  diff=944.7156014894135
step=5   xold,new= 489255.2868952168, 489255.2868952168  diff=488266.58932961605
step=6   xold,new= 1.196856125046039e11, 1.196856125046039e11  diff=1.19685123249317e11
step=7   xold,new= 7.16232292036094e21, 7.16232292036094e21  diff=7.162322920241255e21

-1 :
step=0   xold,new= -0.5, -0.5  diff=0.5
step=1   xold,new= -0.625, -0.625  diff=-0.125
step=2   xold,new= -0.6171875, -0.6171875  diff=0.0078125
step=3   xold,new= -0.618133544921875, -0.618133544921875  diff=-0.000946044921875
step=4   xold,new= -0.6180222327820957, -0.6180222327820957  diff=0.0001113121397793293
step=5   xold,new= -0.6180353762845644, -0.6180353762845644  diff=-1.3143502468726531e-5
step=6   xold,new= -0.6180338249726807, -0.6180338249726807  diff=1.5513118837295892e-6
step=7   xold,new= -0.6180340080811593, -0.6180340080811593  diff=-1.8310847860192325e-7
step=8   xold,new= -0.6180339864681484, -0.6180339864681484  diff=2.161301082548306e-8
step=9   xold,new= -0.6180339890192185, -0.6180339890192185  diff=-2.551070021894475e-9
step=10   xold,new= -0.6180339887181056, -0.6180339887181056  diff=3.0111291238199556e-10
step=11   xold,new= -0.618033988753647, -0.618033988753647  diff=-3.5541458665022674e-11

Description:
It looks as though if our initial guess is close enough to the negative root
we will get convergence to it.
=#

fpi_plot(g_alpha(-0.5), 1.0)
title!("FPI: α=-.05 & x_0 = 1")
#= - Output -
step=0   xold,new= 1.5, 1.5  diff=0.5
step=1   xold,new= 1.625, 1.625  diff=0.125
step=2   xold,new= 1.6171875, 1.6171875  diff=-0.0078125
step=3   xold,new= 1.618133544921875, 1.618133544921875  diff=0.000946044921875
step=4   xold,new= 1.6180222327820957, 1.6180222327820957  diff=-0.0001113121397793293
step=5   xold,new= 1.6180353762845645, 1.6180353762845645  diff=1.3143502468837553e-5
step=6   xold,new= 1.6180338249726804, 1.6180338249726804  diff=-1.551311884062656e-6
step=7   xold,new= 1.6180340080811593, 1.6180340080811593  diff=1.8310847882396786e-7
step=8   xold,new= 1.6180339864681486, 1.6180339864681486  diff=-2.1613010714460756e-8
step=9   xold,new= 1.6180339890192184, 1.6180339890192184  diff=2.55106979984987e-9
step=10   xold,new= 1.6180339887181057, 1.6180339887181057  diff=-3.0111269033739063e-10
step=11   xold,new= 1.618033988753647, 1.618033988753647  diff=3.554134764272021e-11

Description:
Houston... we have convergence!
The convergence is linear, as can be seen in the plot fpi_1.png
=#

#= Part 3
Using f'(r) = 2ϕ - 1 ≈ 2.24 we shoudl get the fastest convergence at
α = -1/f'(r) ≈ -.447
=#



# ========= Question 4 ==========
function g_beta(β)
    return function g(x)
        return (β * x + (x^2 - 1)) / (β + 1)
    end
end

fpi(g_beta(1), 1)
fpi(g_beta(1), 2)
fpi(g_beta(1), -1)
#= - Output -
1 -
step=0   xold,new= 0.5, 0.5  diff=-0.5
step=1   xold,new= -0.125, -0.125  diff=-0.625
step=2   xold,new= -0.5546875, -0.5546875  diff=-0.4296875
step=3   xold,new= -0.623504638671875, -0.623504638671875  diff=-0.068817138671875
step=4   xold,new= -0.6173733021132648, -0.6173733021132648  diff=0.006131336558610201
step=5   xold,new= -0.6181117539755141, -0.6181117539755141  diff=-0.0007384518622493008
step=6   xold,new= -0.6180248067864138, -0.6180248067864138  diff=8.694718910029486e-5
step=7   xold,new= -0.6180350724915149, -0.6180350724915149  diff=-1.0265705101097922e-5
step=8   xold,new= -0.6180338608309615, -0.6180338608309615  diff=1.2116605534462366e-6
step=9   xold,new= -0.6180340038486686, -0.6180340038486686  diff=-1.4301770712155104e-7
step=10   xold,new= -0.6180339869677263, -0.6180339869677263  diff=1.6880942310670832e-8
2 -
step=0   xold,new= 2.5, 2.5  diff=0.5
step=1   xold,new= 3.875, 3.875  diff=1.375
step=2   xold,new= 8.9453125, 8.9453125  diff=5.0703125
step=3   xold,new= 43.981964111328125, 43.981964111328125  diff=35.036651611328125
step=4   xold,new= 988.6975656007417, 988.6975656007417  diff=944.7156014894135
step=5   xold,new= 489255.2868952168, 489255.2868952168  diff=488266.58932961605
step=6   xold,new= 1.196856125046039e11, 1.196856125046039e11  diff=1.19685123249317e11
step=7   xold,new= 7.162322920360941e21, 7.162322920360941e21  diff=7.162322920241256e21
-1 -
step=0   xold,new= -0.5, -0.5  diff=0.5
step=1   xold,new= -0.625, -0.625  diff=-0.125
step=2   xold,new= -0.6171875, -0.6171875  diff=0.0078125
step=3   xold,new= -0.618133544921875, -0.618133544921875  diff=-0.000946044921875
step=4   xold,new= -0.6180222327820957, -0.6180222327820957  diff=0.0001113121397793293
step=5   xold,new= -0.6180353762845644, -0.6180353762845644  diff=-1.3143502468726531e-5
step=6   xold,new= -0.6180338249726807, -0.6180338249726807  diff=1.5513118837295892e-6
step=7   xold,new= -0.6180340080811593, -0.6180340080811593  diff=-1.8310847860192325e-7
step=8   xold,new= -0.6180339864681484, -0.6180339864681484  diff=2.161301082548306e-8
step=9   xold,new= -0.6180339890192185, -0.6180339890192185  diff=-2.551070021894475e-9
step=10   xold,new= -0.6180339887181054, -0.6180339887181054  diff=3.01113023404298e-10
step=11   xold,new= -0.618033988753647, -0.618033988753647  diff=-3.5541569687325136e-11

Description:
It look like we found the other root to our function with x = 1. For other values
it looks like if we are close enough to the root it will converge.
=#

fpi_plot(g_beta(-2.5), 1)
title!("FPI: β=-2.5 & x=1")
#= - Output -
step=0   xold,new= 1.6666666666666667, 1.6666666666666667  diff=0.6666666666666667
step=1   xold,new= 1.5925925925925926, 1.5925925925925926  diff=-0.07407407407407418
step=2   xold,new= 1.6300868770004573, 1.6300868770004573  diff=0.03749428440786473
step=3   xold,new= 1.612022643954693, 1.612022643954693  diff=-0.018064233045764322
step=4   xold,new= 1.6209597368427024, 1.6209597368427024  diff=0.008937092888009479
step=5   xold,new= 1.616592582427729, 1.616592582427729  diff=-0.004367154414973484
step=6   xold,new= 1.6187399190059788, 1.6187399190059788  diff=0.0021473365782498366
step=7   xold,new= 1.617687248087643, 1.617687248087643  diff=-0.0010526709183358296
step=8   xold,new= 1.618204058395824, 1.618204058395824  diff=0.0005168103081809594
step=9   xold,new= 1.6179505142538968, 1.6179505142538968  diff=-0.00025354414192713115
step=10   xold,new= 1.6180749460401953, 1.6180749460401953  diff=0.0001244317862985067
step=11   xold,new= 1.6180138893983378, 1.6180138893983378  diff=-6.105664185751536e-5
step=12   xold,new= 1.6180438514732722, 1.6180438514732722  diff=2.996207493444203e-5
step=13   xold,new= 1.61802914892848, 1.61802914892848  diff=-1.4702544792299932e-5
step=14   xold,new= 1.6180363636926522, 1.6180363636926522  diff=7.2147641723141476e-6
step=15   xold,new= 1.6180328233332597, 1.6180328233332597  diff=-3.540359392539827e-6
step=16   xold,new= 1.6180345606328996, 1.6180345606328996  diff=1.737299639925638e-6
step=17   xold,new= 1.618033708119832, 1.618033708119832  diff=-8.52513067695071e-7

Description:
It looks like we have linear convergence based on fpi_2.png.
=#

#= Part 3
We shoudlget β = -2ϕ ≈ -3.236
=#
