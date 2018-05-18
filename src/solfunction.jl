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

@everywhere function solfunc_ui(θcp, θ_lat, θ_lon, T)
    # 均匀入流求解函数

    # 变量初始化
    # θcp = 2.0*π/180
    # θ_lat = 0.0
    # θ_lon = 0.0
    twsitr = 0.5*R
    twist1 = 0.0
    twist2 = 0.0
    chord = zeros(Nb,Nbe)
    rb = zeros(Nb,Nbe)
    dr = zeros(Nb,Nbe)
    cbe = 0.8
    β1c = 0.0
    β1s = 0.0

    for k in 1:Nb
        for i in 1:Nbe
            chord[k,i] = chroot*taper*(i-1)/Nbe
            dr[k,i] = (R*(1-ecut)*(sin(i/Nbe*cbe*π/2)-
                        sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2))
            rb[k,i] = (ecut*R+R*(1-ecut)*(sin(i/Nbe*cbe*π/2)+
                        sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2)/2)
        end
    end
    β = zeros(Nb,npsi+1)
    dβ = zeros(Nb,npsi+1)
    ddβ = zeros(Nb,npsi+1)
    # 变量初始化完成

    uitmp = uniforminflow(T)
    λind = uitmp[1]
    vall_s = uitmp[2]

    θ0 = the0(θcp, twsitr, rb)

    betatmp = bladeflap(β, dβ, ddβ, vall_s, chord, θ0, dr)
    if betatmp[1]
        β = betatmp[2]
        dβ = betatmp[3]
        ddβ = betatmp[4]
        β1c = betatmp[5]
        β1s = betatmp[6]
    end

    rftmp = rotoraero(vall_s, chord, β, ddβ, θ0, dr, rb)
    fx_s = rftmp[1]
    fy_s = rftmp[2]
    fz_s = rftmp[3]
    MQ   = rftmp[4]
    power= rftmp[5]

    return fz_s, β1c/π*180, β1s/π*180, fy_s, power
end
