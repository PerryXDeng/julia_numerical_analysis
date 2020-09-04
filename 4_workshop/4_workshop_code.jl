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


# Question 2
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
