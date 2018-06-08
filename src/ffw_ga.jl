# this file is used for GA module applied into ffw method in trimming

include(pwd()*"//gamodule//gamd.jl")
include(pwd()*"//src//ffw_tr.jl")

using Plots; gr()
using GAMD

# set paremeters for ga module
const ncre = 20
const nmax = 20*20
const nmin = 20
const npar = 6
const mr = 0.5
const rpr = 0.8

const varmin = [-22.,-22.,-22.]*π/180
const varmax = [22.,22.,22.]*π/180

const ξ = 0.6 # triming weighting coefficient
const cfit = 1. # if minfits is less than cfit, then convergence is reached too

# functions
@everywhere function fitness(ct::Creature)
    vars = ct.vars
    fittmp = sol_ffw(vars[1], vars[2], vars[3])
    fits = ξ*abs(T-fittmp[1])/10 +(1-ξ)/2*abs(0.-fittmp[3]/(π/180))+(1-ξ)/2*abs(0.-fittmp[4])/(π/180)
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
