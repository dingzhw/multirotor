# this file is used for GA module applied into ffw method in trimming

include(pwd()*"//gamodule//gamd.jl")

@everywhere using Plots; gr()
using GAMD

# functions
@everywhere function fitness(ct::Creature)
    vars = ct.vars
    fittmp = sol_ffw(vars[1], vars[2], vars[3])
    fits = ξ*abs(T-fittmp[1])/10 +(1-ξ)/2*abs(0.-fittmp[3])+(1-ξ)/2*abs(0.-fittmp[4])
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
trresult = trfm(the[1], the[2], the[3], sol_ffw, ntiter)
