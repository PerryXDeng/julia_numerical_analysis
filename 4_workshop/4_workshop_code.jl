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
MAXACC = 10.0^50
max_iter = 50

function fpi(f, xold)
    diff = 10*FACC
    step = 0

    while step < max_iter && abs(diff) > FACC && abs(diff) < MAXACC
        xnew = f(xold)
        diff = xnew-xold
        xold = xnew

        println("step=", step, "   xold,new= ", xold, ", ", xnew, "  diff=", diff)
        step += 1
    end
end

function g(x)
    return (x^2 - x - 1) + x
end

fpi(g, 1.5)
