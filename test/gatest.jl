# file for ga test

include(pwd()*"//gamodule//gamd.jl")
include(pwd()*"//gamodule//gafunc.jl")

using Plots;gr()
using GAMD

const ncre = 200
const nmax = 200*200
const nmin = 20
const npar = 6
const mr = 0.5
const rpr = 0.8

const varmin = [-42.,-42.,-42.,-42.,-42.]
const varmax = [42.,42.,42.,42.,42.]

pop = Population(ncre, nmax, nmin, npar, mr, rpr)
garun(pop, varmin, varmax)
