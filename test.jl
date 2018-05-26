# include(pwd()*"\\src\\00_const.jl");
# include(pwd()*"\\src\\00_mathfunctions.jl")
# include(pwd()*"\\src\\01_uniforminflow.jl")
# include(pwd()*"\\src\\02_bladeflap.jl")
# include(pwd()*"\\src\\clcd.jl")
# include(pwd()*"\\src\\rotorforce.jl")
# include(pwd()*"\\src\\solfunction.jl")
#
# # testtmp = solfunc_ui()
# # print(show(testtmp))

Γ = 10
A = [1.,0,0]
B = [0.,1,0]
C = [-1.,0,0]
D = [0.,-1,0]
P = [[1.,1,1]]

vr1 = vortexring1(Γ,A,B,C,D,P,croc,vrtosys,systovr,vringvind)
vr2 = vortexring1(-20,A,B,C,D,P,croc,vrtosys,systovr,vringvind)
