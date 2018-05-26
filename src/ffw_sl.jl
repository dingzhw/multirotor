# 针对单组操纵输入的完整求解函数
# --- that means in the function we will keep the control input
# --- and calculate the force, moment, flap angle and so on

@everywhere function sol_ffw(θcp, θ_lat, θ_lon, T)
    # fast free wake solution function

    # 变量初始化
    # θcp = 2.0*π/180
    # θ_lat = 0.0
    # θ_lon = 0.0
    twsitr = 0.5*R
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
    # 变量初始化完成

    t = 0
    ggtmp = gamget(T)
    Γ = ggtmp[1]
    dτ = ggtmp[2]
    push!(vr,vortexring1(Γ,A0,B0,C0,D0,croc,vrtosys,systovr,vringvind))
    for i in 1:length(vr)
        # calculate the disk induced velocity distribution
        p = reshape(pdisk,1,npsi*Nbe)
        vindtmp = vr[i].vrtosys(vr[i], p)
        vind = reshape(vindtmp,Nbe,npsi)
        transpose!(vind)
        vdiskind += vind
    end

    # calculate the total velocity distribution in disk
    vber = vbe(vdiskind, β, dβ)

    # calculate the blade flap
    bftmp = bladeflap(β, dβ, ddβ, vber, chord, θ0, θ_lat, θ_lon, dr, rb)
    if betatmp[1]   # judge if the flap iteration converaged
        β    = betatmp[2]
        dβ   = betatmp[3]
        ddβ  = betatmp[4]
        β0   = betatmp[5]
        βlon = betatmp[6]
        βlat = betatmp[7]
    end

    # calculate force, moment and power
    rftmp = rotoraero(vber, chord, β, dβ, ddβ, θ0, θ_lat, θ_lon, dr, rb)
    fx_s = rftmp[1]
    fy_s = rftmp[2]
    fz_s = rftmp[3]
    MQ   = rftmp[4]
    power= rftmp[5]

    return fz_s, β0/π*180, βlon/π*180, βlat/π*180, fy_s, power, λind
end
