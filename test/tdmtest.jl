# main test

tic(); # at this point the project runs

# parallel computating settings
np = nprocs()
print("*** Input your cores number and press 'Enter': \n")
np_ = "4" # readline()
np_ = parse(Int64, np_)
print("*** Your cores number is $(np_) .\n\n")
addprocs(np_)

@everywhere include(pwd()*"//src//00_const.jl")
include(pwd()*"//src//00_mathfunctions.jl")
include(pwd()*"//src//01_uniforminflow.jl")
include(pwd()*"//src//ffw.jl")
include(pwd()*"//src//ffw_clcd.jl")
include(pwd()*"//src//ffw_bf.jl")
include(pwd()*"//src//ffw_rf.jl")
include(pwd()*"//src//tdm_sol.jl")
include(pwd()*"//src//tdm_tr.jl")
include(pwd()*"//src//tdm_ga.jl")

toc(); # at this points the project runs out
