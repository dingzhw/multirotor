# main test

include(pwd()*"//src//00_const.jl")
include(pwd()*"//src//00_mathfunctions.jl")
include(pwd()*"//src//ffw.jl")
include(pwd()*"//src//ffw_clcd.jl")
include(pwd()*"//src//ffw_bf.jl")
include(pwd()*"//src//ffw_rf.jl")
include(pwd()*"//src//tdm_sol.jl")
include(pwd()*"//src//tdm_tr.jl")

tic(); # at this point the project runs
θtdm1 = [11/57.3, -0.0209714, 0.0154191]
θtdm2 = [9/57.3, -0.0268583, 0.0326203]
trresult = trfm(θtdm1, θtdm2, sol_tdm, ntiter)

toc(); # at this points the project runs out
