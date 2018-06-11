# 针对单组操纵输入的完整求解函数
# --- that means in the function we will keep the control input
# --- and calculate the force, moment, flap angle and so on

@everywhere function sol_ffw(θcp, θ_lat, θ_lon)
    # fast free wake solution function

    # 变量初始化
    # θcp = 2.0*π/180
    # θ_lat = 0.0
    # θ_lon = 0.0
    T0 = T
    twistr = 0.8*R
    twist1 = -6.50*π/180
    twist2 = -6.50*π/180
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

    θ0tmp = the0(θcp, twistr, twist1, twist2, rb) # calculate the install angle
    θ0 = θ0tmp[1]
    θ75 = θ0tmp[2]

    βtmp = staticbf(θ75, (twist1+twist2), θ_lon, θ_lat, μ_air, abs(λ_air+uniforminflow(T)[1]))
    β = βtmp[1] # zeros(npsi+1)
    dβ = βtmp[2] # zeros(npsi+1)
    ddβ = βtmp[3] # zeros(npsi+1)
    β0 = βtmp[4]
    βlon = βtmp[5]
    βlat = βtmp[6]
    # print("β is $(β)++++++\n
    # dβ is $(dβ)++++++\n
    # ddβ is $(ddβ)+++++\n")

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
    Tj    = zeros(2)
    iternum = 1
    t = 0.
    # 变量初始化完成

    while true
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
            Tj[1]     = T0
        elseif iternum == 2
            vindj[2] = vdiskind
            Tj[2]     = T0
        else
            vindj[1] = vindj[2]
            Tj[1]     = Tj[2]
            vindj[2] = vdiskind
            Tj[2]     = T0
        end

        # judge if the indeuce velocity is trimmed
        # rmsind = trvind(vindj[1],vindj[2])
        if iternum>100
            print("=== Induced Velocity can not be converaged ===\n")
            print("=== Thrust is $(T0) ===\n")
            # print("=== The Total Power need is $(power)===\n")
            return fz_s, β0/π*180, βlon/π*180, βlat/π*180, fy_s, power
            break
        end

        if iternum>2 && abs(Tj[2]-Tj[1])<= 5 #trvind(vindj[1],vindj[2])[1]
            print("=== Induced Velocity is converaged ===\n")
            print("=== Thrust is $(T0) ===\n")
            # print("=== The Total Power need is $(power)===\n")
            return fz_s, β0/π*180, βlon/π*180, βlat/π*180, fy_s, power
            break
        end

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
    end
end
