# main test

include(pwd()*"//src//00_const.jl");
include(pwd()*"//src//00_mathfunctions.jl")
include(pwd()*"//src//ffw.jl")
include(pwd()*"//src//ffw_clcd.jl")
include(pwd()*"//src//ffw_bf.jl")
include(pwd()*"//src//ffw_rf.jl")
include(pwd()*"//src//ffw_sl.jl")
include(pwd()*"//src//ffw_tr.jl")

tic(); # at this point the project runs
trresult = trfm(15.0/57.3, 0., 0., sol_ffw, ntiter)

toc(); # at this points the project runs out
