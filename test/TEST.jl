# while true
ggtmp = gamget(T0)
Γ = ggtmp[1]
dτ = ggtmp[2]
# new vortex ring generate
push!(vr,vortexring1(Γ,A0,B0,C0,D0,croc,vrtosys,systovr,vringvind))
# calculate the induced velocity
vdiskind = zeros(npsi,Nbe)
for i in 1:length(vr)
    # calculate the disk induced velocity distribution
    p = reshape(pdisk,1,npsi*Nbe)
    vindtmp = vr[i].vrtosys(vr[i], p)
    vind = reshape(vindtmp,npsi,Nbe)
    # vind = tstrans(vind)
    vdiskind += vind
end

# record the No. i-1 and i induced velocity for compare
if iternum == 1
    vindj[1] = mean(vdiskind)[3]/(Ω*R)
    Tj[1]     = T0
elseif iternum == 2
    vindj[2] = mean(vdiskind)[3]/(Ω*R)
    Tj[2]     = T0
else
    vindj[1] = vindj[2]
    Tj[1]     = Tj[2]
    vindj[2] = mean(vdiskind)[3]/(Ω*R)
    Tj[2]     = T0
end

# judge if the indeuce velocity is trimmed
# rmsind = trvind(vindj[1],vindj[2])
if iternum>10
    print("=== Induced Velocity can not be converaged ===\n")
    print("=== Thrust is $(T0) ===\n")
    # # print("=== The Total Power need is $(power)===\n")
    # return fz_s, β0, βlon, βlat, fy_s, power, vind_
    # break
end

if iternum>2 && abs(vindj[2]-vindj[1])<=1e-3# && abs(Tj[2]-Tj[1])<= 5
    print("=== Induced Velocity is converaged ===\n")
    print("=== Thrust is $(T0) ===\n")
    # # print("=== The Total Power need is $(power)===\n")
    # return fz_s, β0, βlon, βlat, fy_s, power, vind_
    # break
end

# calculate the total velocity distribution in disk
vbertmp = vbe(vdiskind, β, dβ, rb)
vber    = vbertmp[1]

# # calculate the blade flap using 精细数值方法
# bftmp = bladeflap(β, dβ, ddβ, vber, chord, θ0, θ_lat, θ_lon, dr, rb)
# if bftmp[1]   # judge if the flap iteration converaged
#     β    = bftmp[2]
#     dβ   = bftmp[3]
#     ddβ  = bftmp[4]
#     β0   = bftmp[5]
#     βlon = bftmp[6]
#     βlat = bftmp[7]
# end

# 使用经验公式求解挥舞
vind_ = mean(vdiskind)[3]/(Ω*R)
βtmp = staticbf(θ75, (twist1+twist2), θ_lon, θ_lat, μ_air, abs.(λ_air+vind_))
β = βtmp[1]
dβ = βtmp[2]
ddβ = βtmp[3]
β0 = βtmp[4]
βlon = βtmp[5]
βlat = βtmp[6]

# calculate force, moment and power
rftmp = rotoraero(vber, chord, β, dβ, ddβ, θ0, θ_lat, θ_lon, dr, rb)
fx_s = rftmp[1]
fy_s = rftmp[2]
fz_s = rftmp[3]
MQ   = rftmp[4]
power= rftmp[5]

# vortex rings move
vr = vrdel(vr, cut)
vr = vrmove(vr, dτ)
vr = vrdel(vr, cut)

# T0 change
T0 = fz_s
# print("=== Thrust is $(fz_s) ===\n")
# print("=== Remain is $(Tj[2]-Tj[1]) ===\n")
# print("=== No is $(iternum) ===\n\n\n")

iternum += 1
t += dτ
# if iternum>100
#     return fz_s, β0/π*180, βlon/π*180, βlat/π*180, fy_s, power
#     break
# end
# end
