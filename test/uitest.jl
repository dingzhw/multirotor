# test file for uniform-inflow

using Plots;gr()
# test start here
tic();
include(pwd()*"//src//00_const.jl")
include(pwd()*"//src//00_mathfunctions.jl")
include(pwd()*"//src//01_uniforminflow.jl")
include(pwd()*"//src//02_bladeflap.jl")
include(pwd()*"//src//clcd.jl")
include(pwd()*"//src//rotorforce.jl")
include(pwd()*"//src//solfunction.jl")
include(pwd()*"//src//trim.jl")

ttainp = [18./57.3, 0., 0.]
trresult = trfm(ttainp[1], ttainp[2], ttainp[3], sol_ui, ntiter)

toc()
