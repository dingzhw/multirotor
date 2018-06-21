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

ttainp = [0.3068026697055875, 0.010809323652109164, -0.03433742775135728]
trresult = trfm(ttainp[1], ttainp[2], ttainp[3], sol_ui, ntiter)
show(trresult)

toc()

# cost extra 1139.1348 seconds 
# ---> just make [0.3068026697055875, 0.010809323652109164, -0.03433742775135728]
# ---> become [0.3071093557471181, -5.913775258831376e-5, -0.0359532799927037]
# ---> it is not necessary and effecitve for fast analysis
