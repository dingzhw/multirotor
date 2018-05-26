# this file is written for test several methods

################################################################################

# # --- the block below is the test for uniforminflow ---
# include(pwd()*"//src//00_const.jl");
# include(pwd()*"//src//00_mathfunctions.jl")
# include(pwd()*"//src//01_uniforminflow.jl")
# include(pwd()*"//src//02_bladeflap.jl")
# include(pwd()*"//src//clcd.jl")
# include(pwd()*"//src//rotorforce.jl")
# include(pwd()*"//src//solfunction.jl")
# # --- the test for uniforminflow is end ---


################################################################################

# --- the block below is the test for FFW method ---
# Γ = 10
# A = [1.,0,0]
# B = [0.,1,0]
# C = [-1.,0,0]
# D = [0.,-1,0]
# P = [[1.,1.,0],[0.99,0.99,0]]
#
# vr1 = vortexring1(Γ,A,B,C,D,croc,vrtosys,systovr,vringvind)
# vr2 = vortexring1(-20,A,B,C,D,croc,vrtosys,systovr,vringvind)
# --- the test of FFW method is end ---


################################################################################

# --- the block below is the test for fast free wake
include(pwd()*"//src//00_const.jl");
include(pwd()*"//src//00_mathfunctions.jl")
include(pwd()*"//src//ffw.jl")
include(pwd()*"//src//ffw_clcd.jl")
include(pwd()*"//src//ffw_bf.jl")
include(pwd()*"//src//ffw_rf.jl")
include(pwd()*"//src//ffw_sl.jl")
# --- the test for ffw is end ---
