# 求解aoa以及对应的cl,cd
@everywhere function the0get(θcp, twistr)
    θ0 = Array{Float64}(Nb,Nbe)
    for k in Nb
        for i in Nbe
            if rb[i]<=twistr
                θ[k,i] = θcp[k]+twist1*(rb[k,i]/twistr-1)
            else
                θ[k,i] = θcp[k]+twist2*((rb[k,i]-twistr)/(R-twistr))
            end
        end
    end
    return θ0
end

@everywhere function theget(ψ, θ_lat,θ_lon)
    θ = Array{Float64}(Nb,Nbe)
    for k in 1:Nb
        ψk = ψ+(k-1)*2*π/Nb
        for i in 1:Nbe
            θ[k,i] = θ0[k,i]+θ_lat*cos(ψk)+θ_lon*sin(ψk)
        end
    end
    return θ
end

@everywhere function aoaget(vall_r, β, θ)
    α = Array{Float64}(Nb, Nbe)
    for k in 1:Nb
        for i in 1:Nbe
            windbe = rotobe(vall_r, β, θ[k,i])
            windyz = [windbe[2],windbe[3]]
            α      = aoaang(windyz)
        end
    end
    return α
end

@everywhere function fcl(α, ma, Re=1e6)
    cl = 0.65
    return cl
end

@everywhere function fcd(α, ma, Re=1e6)
    cd = 0.5
    return cd
end
