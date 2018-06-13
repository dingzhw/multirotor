# this file is used for GA module applied into ffw method in trimming

include(pwd()*"//gamodule//gamd.jl")

@everywhere using Plots; gr()
using GAMD

# 确定操纵量的取值范围；用于遗传算法寻优
@everywhere varmin = [0., -20., -20., 0., -20., -20.]/(180/π)
@everywhere varmax = [42., 20., 20., 42., 20., 20.]/(180/π)

# functions
@everywhere function fitness(ct::Creature)
    # 适应度函数
    vars = ct.vars
    soltmp = sol_tdm(vars)

    if soltmp[3][1] > 0
        fits = ξ/2*abs(T-soltmp[3][1])/10 + ξ/2*abs(0.-soltmp[3][2])/2
            (1-ξ)/2*(abs(0.-soltmp[1][3])+abs(0.-soltmp[1][4]))/(π/180)+
            (1-ξ)/2*(abs(0.-soltmp[2][3])+abs(0.-soltmp[2][4]))/(π/180)
    else
        fits = lossξ/2*abs(T-soltmp[3][1])/10 + lossξ/2*abs(0.-soltmp[3][2])/2
            (1-lossξ)/2*(abs(0.-soltmp[1][3])+abs(0.-soltmp[1][4]))/(π/180)+
            (1-lossξ)/2*(abs(0.-soltmp[2][3])+abs(0.-soltmp[2][4]))/(π/180)
    end

    return fits
end

@everywhere function constraint(vars::Array{Float64})
    return true
end

# run ffw_ga
pop = Population(ncre, nmax, nmin, npar, mr, rpr)
gatmp = garun(pop, varmin, varmax, cfit)
thetmp = gatmp[3]
the = thetmp.vars
trresult = trfm(the, sol_tdm, ntiter)
