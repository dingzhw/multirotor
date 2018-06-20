# The aerodynamics force of rotor

@everywhere function beaero(vbe_r, chord, α, β, θ, dr, rotor=1)
    # 叶素气动力
    # α的范围为0~360°

    # vbe = rotobe(vbe_r, β, θ) # 将速度从旋转坐标系转换到叶素当地坐标系
    # blade elements local lift unit vector
    if fcl(α)>=0.0
        lbe = [0, cos(α+π/2), sin(α+π/2)]
    else
        lbe = [0, cos(α-π/2), sin(α-π/2)]
    end

    dbe = [0, cos(α), sin(α)] # blade elements local drag unit vector

    # !!! how to get the VECTOR of airfoil lift correcttly is important
    # $$$ 目前该问题通过上述判断语句分开解决，后需考虑是否有更优方法
    # --- for now, the better method to judge Cl is use "Total Vector"
    # --- as shown above
    lift = lbe*abs(1/2*ρ*norm(vbe_r)^2*chord*fcl(α))
    drag = dbe*abs(1/2*ρ*norm(vbe_r)^2*chord*fcd(α))
    lift_beta = betobeta(lift, β)
    drag_beta = betobeta(drag, β)
    fz_β      = (lift_beta[3]+drag_beta[3])*dr # the blade element z-force in flap coordination
    lift_ro = betoro(lift, β, θ, rotor)
    drag_ro = betoro(drag, β, θ, rotor)
    fz_be = (lift_ro[3]+drag_ro[3])*dr # the blade element aerodynamic force in rotation coordination
    fy_be = (lift_ro[2]+drag_ro[2])*dr
    return fy_be, fz_be, fz_β
end

@everywhere function bladeaero(vall_r, chord, α, β, ddβ, θ, dr, rb, rotor=1)
    # vall_r here is the column value of the total vber ```

    fy_blade = 0.0
    fz_blade = 0.0
    # Mβ       = 0.0
    MQ       = 0.0
    Maero    = 0.0

    for j in 1:Nbe
        fbe = beaero(vall_r[j], chord[j], α[j], β, θ[j], dr[j], rotor)
        fy_blade += fbe[1]
        fz_blade += fbe[2]
        MQ       += fbe[1]*rb[j]
        Maero    += fbe[3]*(rb[j]-eflap)

        # Mβ       += (m_*g*dr[j]*(rb[j]-eflap)*cos(β)+
        #              m_*dr[j]*ddβ*(rb[j]-eflap)^2+
        #              m_*dr[j]*Ω^2*(eflap+(rb[j]-eflap)*cos(β))*(rb[j]-eflap)*sin(β)-
        #              fbe[3]*(rb[j]-eflap))
    end
    return fy_blade, fz_blade, MQ, Maero
end

@everywhere function rotoraero(vber, chord, β, dβ, ddβ, θ0, θ_lat, θ_lon, dr, rb, rotor=1)
    # summary all the force from all blades
    fy_r = 0.0
    fz_r = 0.0
    fx_s = 0.0
    fy_s = 0.0
    fz_s = 0.0
    MQ   = 0.0
    frec = Float64[]

    # vber = vbe(vind, β, dβ)
    θ  = theget(θ0, θ_lat, θ_lon, rotor)
    α = aoaget(vber, β, θ, rotor)
    for i in 1:npsi
        # print("============\n")
        # print("PSI IS: $(ψ) ;\n Vro is : $(vall_r) ;\n\n Alpha is $(α) ;\n\n
        #         BETA is $(β[i]) ;")
        # print("\n============\n")
        ψ = (-1)^(rotor-1)*(i-1)*dψ
        fblade = bladeaero(vber[i,:], chord, α[i,:], β[i], ddβ[i], θ[i,:], dr, rb, rotor)
        push!(frec, fblade[2])
        fy_r += fblade[1]
        fz_r += fblade[2]
        fx_s += -fy_r*sin(ψ)
        fy_s += fy_r*cos(ψ)
        MQ   += fblade[3]
    end

    x = 1:npsi
    fplot = plot(x, frec, title="force change with azimuth", xlabel="Azimuth Angle", ylabel="Lift")
    display(fplot)

    fy_r = fy_r/npsi*Nb
    fz_r = fz_r/npsi*Nb
    MQ   = MQ/npsi*Nb
    power = (-1)^rotor*MQ*Ω/1e3

    fx_s = fx_s/npsi
    fy_s = fy_s/npsi
    fz_s = fz_r
    return fx_s, fy_s, fz_s, MQ, power
end
