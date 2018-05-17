# The aerodynamics force of rotor

@everywhere function beaero(vall_r, chord, α, β, θ, dr)
    # 叶素气动力
    vall_be = rotobe(vall_r, β, θ) # 将速度从旋转坐标系转换到叶素当地坐标系
    if vall_be[2]>=0 # blade elements local lift unit vector
        if fcl(α)>=0.0
            lbe = [0, cos(α+π/2), sin(α+π/2)]
        else
            lbe = [0, cos(α-π/2), sin(α-π/2)]
        end
    else
        if fcl(α)>=0.0
            lbe = [0, cos(α-π/2), sin(α-π/2)]
        else
            lbe = [0, cos(α+π/2), sin(α+π/2)]
        end
    end
    dbe = [0, cos(α), sin(α)] # blade elements local drag unit vector

    # !!! how to get the VECTOR of airfoil lift correcttly is important
    # $$$ 目前该问题通过上述判断语句分开解决，后需考虑是否有更优方法
    lift = lbe*abs(1/2*ρ*norm(vall_r)^2*chord*fcl(α))
    drag = dbe*abs(1/2*ρ*norm(vall_r)^2*chord*fcd(α))
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
        fbe = beaero(vall_r[k,j], chord, α, β, θ, dr[j])
        fy_blade += fbe[1]
        fz_blade += fbe[2]
        MQ       += fbe[1]*rb[j]
        Mβ       += (m_*g*dr[j]*(rb[j]-eflap)*cos(β)+m_*dr[j]*ddβ*(rb[j]-eflap)^2+
                     m_*dr[j]*Ω^2*(rb[j]-eflap)*cos(β)*(rb[j]-eflap)*sin(β)-
                     fbe[2]*dr[j]*(rb[j]-eflap))
    end
    return fy_blade, fz_blade, MQ, Mβ
end

@everywhere function rotoraero(vall_s, chord, α, β, ddβ, θ, dr, rb)
    # summary all the force from all blades
    fy_r = 0.0
    fz_r = 0.0
    MQ   = 0.0

    for i in 1:npsi
        ψ = dψ*i+2*π/Nb*(k-1)
        vall_r = vallr(vall_s, ψ)
        for k in 1:Nb
            fblade = bladeaero(vall_r, chord, α, β, ddβ, θ, dr, k)
            fy_r += fblade[1]
            fz_r += fblade[2]
            MQ   += fblade[3]
        end
    end

    fy_r = fy_r/npsi
    fz_r = fz_r/npsi
    MQ   = MQ/npsi
    power = MQ*Ω

    fx_s = -fy_r*sin(ψ)
    fy_s = fy_r*cos(ψ)
    fz_s = fz_r
    return fx_s, fy_s, fz_s, mx, my, MQ, power
end
