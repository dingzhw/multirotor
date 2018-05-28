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
    vindj[1] = vdiskind
elseif iternum == 2
    vindj[2] = vdiskind
else
    vindj[1] = vindj[2]
    vindj[2] = vdiskind
end

# # judge if the indeuce velocity is trimmed
# if iternum>1 && trvind(vindj[1],vindj[2])
#     print("=== Induced Velocity is converaged ===\n")
#     return fz_s, β0/π*180, βlon/π*180, βlat/π*180, fy_s, power
#     break
# end

# calculate the total velocity distribution in disk
vbertmp = vbe(vdiskind, β, dβ, rb)
vber    = vbertmp[1]

# calculate the blade flap
bftmp = bladeflap(β, dβ, ddβ, vber, chord, θ0, θ_lat, θ_lon, dr, rb)
if bftmp[1]   # judge if the flap iteration converaged
    β    = bftmp[2]
    dβ   = bftmp[3]
    ddβ  = bftmp[4]
    β0   = bftmp[5]
    βlon = bftmp[6]
    βlat = bftmp[7]
end

# calculate force, moment and power
rftmp = rotoraero(vber, chord, β, dβ, ddβ, θ0, θ_lat, θ_lon, dr, rb)
fx_s = rftmp[1]
fy_s = rftmp[2]
fz_s = rftmp[3]
MQ   = rftmp[4]
power= rftmp[5]

# vortex rings move
vr = vrmove(vr, dτ)

# T change
T0 = fz_s
print("=== T is $(fz_s) ===\n")
print("=== No is $(iternum) ===\n")
print("=== vind is $(vdiskind[1]) ===\n\n\n")
iternum += 1
t += dτ
