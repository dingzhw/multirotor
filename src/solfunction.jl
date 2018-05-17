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

@everywhere function solfunc_ui()
    # 均匀入流求解函数

    # 变量初始化

    # 变量初始化完成

    uitmp = uniforminflow(T)
    λind = uitmp[1]
    vall_s = uitmp[2]

    θ0 = the0(θcp, twsitr)
    θ  = theget(ψ, θ_lat, θ_lon)
    α  = aoaget(vall_r, β, θ)

    clcdtmp = clcdget(vall_r, α)
    cl = clcdtmp[1]
    cd = clcdtmp[2]

    rftmp = rotoraero(vall_s, chord, α, β, ddβ, θ, dr, rb)
    betatmp = bladeflap(β, dβ, ddβ, vall_s, chord, α, θ, dr)
