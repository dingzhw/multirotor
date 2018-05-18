# The aerodynamics force of rotor

@everywhere function beaero(vbe_r, chord, α, β, θ, dr)
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
    lift = lbe*abs(1/2*ρ*norm(vbe_r)^2*chord*fcl(α))
    drag = dbe*abs(1/2*ρ*norm(vbe_r)^2*chord*fcd(α))
    lift_ro = betoro(lift, β, θ)
    drag_ro = betoro(drag, β, θ)
    fz_be = (lift_ro[3]+drag_ro[3])*dr
    fy_be = (lift_ro[2]+drag_ro[2])*dr
    return fy_be, fz_be
end

@everywhere function bladeaero(vall_r, chord, α, β, ddβ, θ, dr, rb, k=1)
    fy_blade = 0.0
    fz_blade = 0.0
    Mβ       = 0.0
    MQ       = 0.0

    for j in 1:Nbe
        fbe = beaero(vall_r[k,j], chord[k,j], α[k,j], β, θ[k,j], dr[k,j])
        fy_blade += fbe[1]
        fz_blade += fbe[2]
        MQ       += fbe[1]*rb[k,j]
        Mβ       += (m_*g*dr[k,j]*(rb[k,j]-eflap)*cos(β)+m_*dr[k,j]*ddβ*(rb[k,j]-eflap)^2+
                     m_*dr[k,j]*Ω^2*(rb[k,j]-eflap)*cos(β)*(rb[k,j]-eflap)*sin(β)-
                     fbe[2]*dr[k,j]*(rb[k,j]-eflap))
    end
    return fy_blade, fz_blade, MQ, Mβ
end

@everywhere function rotoraero(vall_s, chord, β, ddβ, θ0, dr, rb)
    # summary all the force from all blades
    fy_r = 0.0
    fz_r = 0.0
    fx_s = 0.0
    fy_s = 0.0
    fz_s = 0.0
    MQ   = 0.0

    for i in 1:npsi
        for k in 1:Nb
            ψ = dψ*i+2*π/Nb*(k-1)
            vall_r = Array{Vector}(Nb,Nbe)
            vall_r = vallr(vall_s, ψ, β[k,i], dβ[k,i])
            θ  = theget(ψ, θ0, θ_lat, θ_lon)
            α = aoaget(vall_r, β[k,i], θ)
            fblade = bladeaero(vall_r, chord, α, β[k,i], ddβ[k,i], θ, dr, rb, k)
            fy_r += fblade[1]
            fz_r += fblade[2]
            fx_s += -fy_r*sin(ψ)
            fy_s += fy_r*cos(ψ)
            MQ   += fblade[3]
        end
    end

    fy_r = fy_r/npsi
    fz_r = fz_r/npsi
    MQ   = MQ/npsi
    power = -MQ*Ω/1e3

    fx_s = fx_s/npsi
    fy_s = fy_s/npsi
    fz_s = fz_r
    return fx_s, fy_s, fz_s, MQ, power
end
