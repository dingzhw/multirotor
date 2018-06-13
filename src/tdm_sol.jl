# this solution function is used for tandem rotors calculation

# 针对单组操纵输入的完整求解函数
# --- that means in the function we will keep the control input
# --- and calculate the force, moment, flap angle and so on

@everywhere function sol_tdm(θtdm::Array{Float64})
    # fast free wake solution function

    # 变量初始化
    θtdm1 = θtdm[1:3]
    θtdm2 = θtdm[4:6]

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
    twist1 = -6.5*π/180
    twist2 = -6.5*π/180
    chord = zeros(Nbe)
    rb = zeros(Nbe) # position of blade control points
    dr = zeros(Nbe) # distance of each two blade control points
    cbe = 0.8   # Range[0,1] for adjust the denisty of rd[] in tip, lager for more dense

    for i in 1:Nbe
        chord[i] = chroot*taper*(i-1)/Nbe
        dr[i] = (R*(1-ecut)*(sin(i/Nbe*cbe*π/2)-
                    sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2))
        rb[i] = (ecut*R+R*(1-ecut)*(sin(i/Nbe*cbe*π/2)+
                    sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2)/2)
    end

    θ0tmp1 = the0(θcp1, twistr, twist1, twist2, rb) # calculate the install angle
    θ0tmp2 = the0(θcp2, twistr, twist1, twist2, rb)
    θ0_r1 = θ0tmp1[1]
    θ0_r2 = θ0tmp2[1]
    θ75_r1 = θ0tmp1[2]
    θ75_r2 = θ0tmp2[2]

    # 旋翼 1 的挥舞初始值；由挥舞经验公式解得
    βtmp1 = staticbf(θ75_r1, (twist1+twist2), θ_lon1, θ_lat1, μ_air, abs(λ_air+uniforminflow(T01)[1]))
    β1 = βtmp1[1] # zeros(npsi+1)
    dβ1 = βtmp1[2] # zeros(npsi+1)
    ddβ1 = βtmp1[3] # zeros(npsi+1)
    β01 = βtmp1[4]
    βlon1 = βtmp1[5]
    βlat1 = βtmp1[6]

    # 旋翼 2 的挥舞初始值；由挥舞经验公式解得
    βtmp2 = staticbf(θ75_r2, (twist1+twist2), θ_lon2, θ_lat2, μ_air, abs(λ_air+uniforminflow(T02)[1]), 2)
    β2 = βtmp2[1] # zeros(npsi+1)
    dβ2 = βtmp2[2] # zeros(npsi+1)
    ddβ2 = βtmp2[3] # zeros(npsi+1)
    β02 = βtmp2[4]
    βlon2 = βtmp2[5]
    βlat2 = βtmp2[6]

    vr = vortexring1[]
    vdiskind_r1 = zeros(npsi,Nbe)
    vdiskind_r2 = zeros(npsi,Nbe)

    pdisk_r1 = Array{Vector}(npsi,Nbe) # disk control points initilization
    pdisk_r2 = Array{Vector}(npsi,Nbe)
    for i in 1:npsi
        ψ1 = (i-1)*dψ
        ψ2 = -(i-1)*dψ
        for j in 1:Nbe
            pdisk_r1[i,j] = [rb[j]*cos(ψ1),rb[j]*sin(ψ1),(rb[j]-eflap)*sin(β1[i])]
            pdisk_r2[i,j] = [rb[j]*cos(ψ2),rb[j]*sin(ψ2),(rb[j]-eflap)*sin(β2[i])]+[0., disr, hr]
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
            vindj[1] = mean(vdiskind_r1+vdiskind_r2)[3]/(Ω*R)
            Tj[1]    = T0
        elseif iternum == 2
            vindj[2] = mean(vdiskind_r1+vdiskind_r2)[3]/(Ω*R)
            Tj[2]    = T0
        else
            vindj[1] = vindj[2]
            Tj[1]    = Tj[2]
            vindj[2] = mean(vdiskind_r1+vdiskind_r2)[3]/(Ω*R)
            Tj[2]    = T0
        end

        # judge if the indeuce velocity is trimmed
        # rmsind = trvind(vindj[1],vindj[2])
        if iternum>2 && abs(vindj[2]-vindj[1]) <= 1e-3# abs(Tj[2]-Tj[1])<=10
            print("++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
            print("=== Induced Velocity is converaged ===\n")
            print("=== Rotor1 Thrust is $(T01) ===\n")
            print("=== Rotor2 Thrust is $(T02) ===\n")
            print("=== Delta of MQ is $(abs(abs(MQ1)-abs(MQ2))) ===\n")
            print("=== Rotor1 Power is $(power1) ===\n")
            print("=== Rotor2 Power is $(power2) ===\n")
            print("++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
            # print("=== The Total Power need is $(power)===\n")
            return [[fz_s1, β01, βlon1, βlat1, fy_s1, power1, vind1_],
                    [fz_s2, β02, βlon2, βlat2, fy_s2, power2, vind2_],
                    [T0, abs(abs(MQ1)-abs(MQ2))]]
            break
        end

        # calculate the total velocity distribution in disk
        vbertmp1 = vbe(vdiskind_r1, β1, dβ2, rb)
        vbertmp2 = vbe(vdiskind_r2, β2, dβ2, rb, 2)
        vber1    = vbertmp1[1]
        vber2    = vbertmp2[1]

        # <--- BEGIN OF BLADE FLAP CALCULATION --->
        # 经验公式求解旋翼挥舞；速度快，效率高，对于稳定飞行，置信度足矣
        # 求解各旋翼桨盘平均诱导速度
        vind1_ = mean(vdiskind_r1)[3]/(Ω*R)
        vind2_ = mean(vdiskind_r2)[3]/(Ω*R)

        # 旋翼 1 的挥舞初始值；由挥舞经验公式解得
        βtmp1 = staticbf(θ75_r1, (twist1+twist2), θ_lon1, θ_lat1, μ_air, abs(λ_air+vind1_))
        β1 = βtmp1[1] # zeros(npsi+1)
        dβ1 = βtmp1[2] # zeros(npsi+1)
        ddβ1 = βtmp1[3] # zeros(npsi+1)
        β01 = βtmp1[4]
        βlon1 = βtmp1[5]
        βlat1 = βtmp1[6]

        # 旋翼 2 的挥舞初始值；由挥舞经验公式解得
        βtmp2 = staticbf(θ75_r2, (twist1+twist2), θ_lon2, θ_lat2, μ_air, abs(λ_air+vind2_), 2)
        β2 = βtmp2[1] # zeros(npsi+1)
        dβ2 = βtmp2[2] # zeros(npsi+1)
        ddβ2 = βtmp2[3] # zeros(npsi+1)
        β02 = βtmp2[4]
        βlon2 = βtmp2[5]
        βlat2 = βtmp2[6]

        # 经验公式求解挥舞完成
        # <--- END OF BLADE FLAP CALCULATION --->

        # # calculate the blade flap; 高阶数值方法，求解效率低，置信度更高，可用于机动飞行
        # bftmp1 = bladeflap(β1, dβ1, ddβ1, vber1, chord, θ0_r1, θ_lat1, θ_lon1, dr, rb)
        # bftmp2 = bladeflap(β2, dβ2, ddβ2, vber2, chord, θ0_r2, θ_lat2, θ_lon2, dr, rb, 2)
        #
        # if bftmp1[1] && bftmp2[1]   # judge if the flap iteration converaged
        #     β1    = bftmp1[2]
        #     dβ1   = bftmp1[3]
        #     ddβ1  = bftmp1[4]
        #     β01   = bftmp1[5]
        #     βlon1 = bftmp1[6]
        #     βlat1 = bftmp1[7]
        #     β2    = bftmp2[2]
        #     dβ2   = bftmp2[3]
        #     ddβ2  = bftmp2[4]
        #     β02   = bftmp2[5]
        #     βlon2 = bftmp2[6]
        #     βlat2 = bftmp2[7]
        # end

        # calculate force, moment and power
        rftmp1 = rotoraero(vber1, chord, β1, dβ1, ddβ1,
                            θ0_r1, θ_lat1, θ_lon1, dr, rb)
        rftmp2 = rotoraero(vber2, chord, β2, dβ2, ddβ2,
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
