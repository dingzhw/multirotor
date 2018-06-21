# main test

# tic(); # at this point the project runs
#
# # parallel computating settings
# np = nprocs()
# print("*** Input your cores number and press 'Enter': \n")
# np_ = readline()
# np_ = parse(Int64, np_)
# print("*** Your cores number is $(np_) .\n\n")
# addprocs(np_)

@everywhere include(pwd()*"//src//00_const.jl");
include(pwd()*"//src//00_mathfunctions.jl")
include(pwd()*"//src//01_uniforminflow.jl")
include(pwd()*"//src//ffw.jl")
include(pwd()*"//src//ffw_clcd.jl")
include(pwd()*"//src//ffw_bf.jl")
include(pwd()*"//src//ffw_rf.jl")
include(pwd()*"//src//ffw_sl.jl")
include(pwd()*"//src//ffw_tr.jl")
# include(pwd()*"//src//ffw_ga.jl")

using Plots;gr()
tic();
trresult = trfm( 0.35371398448316366, 0.010398093373873027, -0.04081833065142815, sol_ffw, ntiter)

toc(); # at this points the project runs out
