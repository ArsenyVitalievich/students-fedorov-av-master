using ModelingToolkit
using OrdinaryDiffEq
using Plots
using DifferentialEquations

function solve_sirsde()
    @parameters t β γ
    @variables s(t) i(t) r(t)
    D = Differential(t)
    N=s+i+r

    eqs = [D(s) ~ -(β*s*i)/N,
        D(i) ~ (β*s*i)/N - γ*i ,
        D(r) ~ γ*i]

    noiseeqs = [sqrt(abs( ((β*s*i)/N*γ*i) )) ,
                sqrt(abs( ((β*s*i)/N)^2 + ((β*s*i)/N + γ*i)^2 + (γ*i)^2 )),
                sqrt(abs( ((β*s*i)/N)*γ*i ))
                ]

    sde = SDESystem(eqs, noiseeqs, t, [s, i, r], [β, γ])

    u0 = [s => 2000.0, #количество восприимчивых
        i => 10.0, #количество инфицированных
        r => 3.0] #количество переболевших

    p  = [β => 0.1, #коэффициент заболевания
        γ => 0.01, #коэффициент выздоровления
        ]

    tspan = (0.0, 1000.0)

    prob = SDEProblem(sde, u0, tspan, p)
    sol = solve(prob, EM(), dt = 1)

    plotly()

    p = plot(sol,vars=(s,i,r)) #SIR
    savefig(p, "C:\\Users\\arsen\\Desktop\\p.html") #mt_sde_phase.png
    println("Phase graph was builded and is saved to C:\\Users\\arsen\\Desktop\\p.html")

    s = plot(sol)
    savefig(s,"C:\\Users\\arsen\\Desktop\\s.html")
    println("Solution graph was builded and is saved to C:\\Users\\arsen\\Desktop\\s.html")

    return sol
end

sol = solve_sirsde();
println("Solution is $(sol)")