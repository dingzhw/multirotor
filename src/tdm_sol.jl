# this solution function is used for tandem rotors calculation

# 针对单组操纵输入的完整求解函数
# --- that means in the function we will keep the control input
# --- and calculate the force, moment, flap angle and so on

@everywhere function sol_tdm(θtdm1::Array{Float64}, θtdm2::Array{Float64})
    # fast free wake solution function

    # 变量初始化
    θcp1 = θtdm1[1]
    θ_lat1 = θtdm1[2]
    θ_lon1 = θtdm1[3]
    θcp2 = θtdm2[1]
    θ_lat2 = θtdm2[2]
    θ_lon2 = θtdm2[3]
    T0 = T
    T01 = T0/2
    T02 = T0/2
    twistr = 0.8*R
    twist1 = -0.0*π/180
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

    θ0_r1 = the0(θcp1, twistr, twist1, twist2, rb) # calculate the install angle
    θ0_r2 = the0(θcp2, twistr, twist1, twist2, rb)

    β_r1 = zeros(npsi+1)
    dβ_r1 = zeros(npsi+1)
    ddβ_r1 = zeros(npsi+1)

    β_r2 = zeros(npsi+1)
    dβ_r2 = zeros(npsi+1)
    ddβ_r2 = zeros(npsi+1)

    vr = vortexring1[]
    vdiskind_r1 = zeros(npsi,Nbe)
    vdiskind_r2 = zeros(npsi,Nbe)

    pdisk_r1 = Array{Vector}(npsi,Nbe) # disk control points initilization
    pdisk_r2 = Array{Vector}(npsi,Nbe)
    for i in 1:npsi
        ψ1 = (i-1)*dψ
        ψ2 = -(i-1)*dψ
        for j in 1:Nbe
            pdisk_r1[i,j] = [rb[j]*cos(ψ1),rb[j]*sin(ψ1),(rb[j]-eflap)*sin(β_r1[i])]
            pdisk_r2[i,j] = [rb[j]*cos(ψ2),rb[j]*sin(ψ2),(rb[j]-eflap)*sin(β_r2[i])]+[0., disr, hr]
        end
    end
    p1 = reshape(pdisk_r1,1,npsi*Nbe)
    p2 = reshape(pdisk_r2,1,npsi*Nbe)

    vindj = Array{Any}(2)
    Tj    = zeros(2)
    iternum = 1
    t1 = 0.
    t2 = 0.
    # 变量初始化完成

    while true
        ggtmp1 = gamget(T01)
        Γ1 = ggtmp1[1]
        dτ = ggtmp1[2]
        ggtmp2 = gamget(T02)
        Γ2 = ggtmp2[1]
        # dτ2 = ggtmp2[2]
        # new vortex ring generate
        push!(vr,vortexring1(Γ1,A0,B0,C0,D0,croc,vrtosys,systovr,vringvind)) # rotor 1
        push!(vr,vortexring1(Γ2,A01,B01,C01,D01,croc,vrtosys,systovr,vringvind)) # rotor 2
        # calculate the induced velocity
        vdiskind_r1 = zeros(npsi,Nbe)
        vdiskind_r2 = zeros(npsi,Nbe)
        for i in 1:length(vr)
            # calculate the disk induced velocity distribution
            vindtmp1 = vr[i].vrtosys(vr[i], p1)
            vindtmp2 = vr[i].vrtosys(vr[i], p2)
            vind1 = reshape(vindtmp1,npsi,Nbe)
            vind2 = reshape(vindtmp2,npsi,Nbe)
            vdiskind_r1 += vind1
            vdiskind_r2 += vind2
        end

        # record the No. i-1 and i induced velocity for compare
        if iternum == 1
            vindj[1] = vdiskind_r1
            Tj[1]    = T0
        elseif iternum == 2
            vindj[2] = vdiskind_r1
            Tj[2]    = T0
        else
            vindj[1] = vindj[2]
            Tj[1]    = Tj[2]
            vindj[2] = vdiskind_r1
            Tj[2]    = T0
        end

        # judge if the indeuce velocity is trimmed
        # rmsind = trvind(vindj[1],vindj[2])
        if iternum>2 && abs(Tj[2]-Tj[1])<=10 #trvind(vindj[1],vindj[2])[1]
            print("=== Induced Velocity is converaged ===\n")
            print("=== Rotor1 Thrust is $(T01) ===\n")
            print("=== Rotor2 Thrust is $(T02) ===\n")
            print("=== Delta of MQ is $(abs(abs(MQ1)-abs(MQ2))) ===\n")
            # print("=== The Total Power need is $(power)===\n")
            return [[fz_s1, β0_r1/π*180, βlon_r1/π*180, βlat_r1/π*180, fy_s1, power1],
                    [fz_s2, β0_r2/π*180, βlon_r2/π*180, βlat_r2/π*180, fy_s2, power2],
                    [T0, abs(abs(MQ1)-abs(MQ2))]]
            break
        end

        # calculate the total velocity distribution in disk
        vbertmp1 = vbe(vdiskind_r1, β_r1, dβ_r2, rb)
        vbertmp2 = vbe(vdiskind_r2, β_r2, dβ_r2, rb, 2)
        vber1    = vbertmp1[1]
        vber2    = vbertmp2[1]

        # calculate the blade flap
        bftmp1 = bladeflap(β_r1, dβ_r1, ddβ_r1, vber1, chord, θ0_r1, θ_lat1, θ_lon1, dr, rb)
        bftmp2 = bladeflap(β_r2, dβ_r2, ddβ_r2, vber2, chord, θ0_r2, θ_lat2, θ_lon2, dr, rb, 2)

        if bftmp1[1] && bftmp2[1]   # judge if the flap iteration converaged
            β_r1    = bftmp1[2]
            dβ_r1   = bftmp1[3]
            ddβ_r1  = bftmp1[4]
            β0_r1   = bftmp1[5]
            βlon_r1 = bftmp1[6]
            βlat_r1 = bftmp1[7]
            β_r2    = bftmp2[2]
            dβ_r2   = bftmp2[3]
            ddβ_r2  = bftmp2[4]
            β0_r2   = bftmp2[5]
            βlon_r2 = bftmp2[6]
            βlat_r2 = bftmp2[7]
        end

        # calculate force, moment and power
        rftmp1 = rotoraero(vber1, chord, β_r1, dβ_r1, ddβ_r1,
                            θ0_r1, θ_lat1, θ_lon1, dr, rb)
        rftmp2 = rotoraero(vber2, chord, β_r2, dβ_r2, ddβ_r2,
                            θ0_r2, θ_lat2, θ_lon2, dr, rb, 2)
        fx_s1 = rftmp1[1]
        fy_s1 = rftmp1[2]
        fz_s1 = rftmp1[3]
        MQ1   = rftmp1[4]
        power1= rftmp1[5]

        fx_s2 = rftmp2[1]
        fy_s2 = rftmp2[2]
        fz_s2 = rftmp2[3]
        MQ2   = rftmp2[4]
        power2= rftmp2[5]


        # vortex rings move
        vr = vrdel(vr, cut)
        vr = vrmove(vr, dτ)
        vr = vrdel(vr, cut)

        # T0 change
        T01 = fz_s1
        T02 = fz_s2
        T0 = T01+T02
        # print("=== Thrust is $(T0) ===\n")
        # print("=== Remain is $(Tj[2]-Tj[1]) ===\n")
        # print("=== No is $(iternum) ===\n\n\n")

        iternum += 1
        t1 += dτ
        t2 += dτ
    end


end
