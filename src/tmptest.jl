# 变量初始化
θcp = 12.0*π/180
θ_lat = 0.0
θ_lon = 0.0
twistr = 0.5*R
twist1 = 0.0
twist2 = 0.0
chord = zeros(Nbe)
rb = zeros(Nbe) # position of blade control points
dr = zeros(Nbe) # distance of each two blade control points
cbe = 0.8   # Range[0,1] for adjust the denisty of rd[] in tip, lager for more dense
βlon = 0.0
βlat = 0.0

for i in 1:Nbe
    chord[i] = chroot*taper*(i-1)/Nbe
    dr[i] = (R*(1-ecut)*(sin(i/Nbe*cbe*π/2)-
                sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2))
    rb[i] = (ecut*R+R*(1-ecut)*(sin(i/Nbe*cbe*π/2)+
                sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2)/2)
end

θ0 = the0(θcp, twistr, twist1, twist2, rb) # calculate the install angle

β = zeros(npsi+1)
dβ = zeros(npsi+1)
ddβ = zeros(npsi+1)

vr = vortexring1[]
vdiskind = zeros(npsi,Nbe)

pdisk = Array{Vector}(npsi,Nbe) # disk control points initilization
for i in 1:npsi
    ψ = (i-1)*dψ
    for j in 1:Nbe
        pdisk[i,j] = [rb[j]*cos(ψ),rb[j]*sin(ψ),(rb[j]-eflap)*sin(β[i])]
    end
end

vindj = Array{Any}(2)
iternum = 1
t = 0.
# 变量初始化完成
