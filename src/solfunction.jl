# 针对单组操纵输入的完整求解函数

@everywhere function solfunc(wakemethod="uniforminflow")
    if wakemethod == "Uniform Inflow"
        return solfunc_ui()
    elseif wakemethod == "Dynamic Inflow"
        return solfunc_di()
    elseif wakemethod == "Fast Free Wake"
        return solfunc_ffw()
    elseif wakemethod == "Typical Free Wake"
        return solfunc_tfw()
    elseif wakemethod == "Vortex Particle Method"
        return solfunc_vpm()
    else
        print("The method is not contained of the scripts, please
                contact the author.")
        return false
    end
end

# --- that means in the function we will keep the control input
# --- and calculate the force, moment, flap angle and so on

@everywhere function sol_ui(θcp, θ_lat, θ_lon)
    # fast free wake solution function

    # 变量初始化
    # θcp = 2.0*π/180
    # θ_lat = 0.0
    # θ_lon = 0.0
    T0 = T
    twistr = 0.5*R
    twist1 = -6.50*π/180
    twist2 = -6.50*π/180
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

    vindj = Array{Any}(2)
    Tj    = zeros(2)
    iternum = 1
    # 变量初始化完成

    while true
        uitmp = uniforminflow(T0) # 均匀入流求解诱导速度
        lmdui = uitmp[1]
        vall_s = uitmp[2]
        vdiskind = uitmp[3]

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
        if iternum>100
            print("=== Induced Velocity can not be converaged ===\n")
            print("=== Thrust is $(T0) ===\n")
            # print("=== The Total Power need is $(power)===\n")
            return fz_s, β0, βlon, βlat, fy_s, power, lmdui
            break
        end

        if iternum>2 && abs(vindj[2]-vindj[1])<=1e-3# && abs(Tj[2]-Tj[1])<= 5
            print("=== Induced Velocity is converaged ===\n")
            print("=== Thrust is $(T0) ===\n")
            # print("=== The Total Power need is $(power)===\n")
            return fz_s, β0, βlon, βlat, fy_s, power, lmdui
            break
        end

        # calculate the total velocity distribution in disk
        vbertmp = vbeui(vall_s, β, dβ, rb)
        vber    = vbertmp

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
        vind_ = lmdui # mean(vdiskind)[3]/(Ω*R)
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

        # T0 change
        T0 = fz_s
        # print("=== Thrust is $(fz_s) ===\n")
        # print("=== Remain is $(Tj[2]-Tj[1]) ===\n")
        # print("=== No is $(iternum) ===\n\n\n")

        iternum += 1
    end
end
