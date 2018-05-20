# 求解aoa以及对应的cl,cd
@everywhere function the0(θcp, twistr, twist1, twist2, rb)
    θ0 = Array{Float64}(Nb,Nbe)
    for k in 1:Nb
        for i in 1:Nbe
            if rb[i]<=twistr
                θ0[k,i] = θcp+twist1*rb[k,i]/twistr
            else
                θ0[k,i] = θcp+twist1+twist2*((rb[k,i]-twistr)/(R-twistr))
            end
        end
    end
    return θ0
end

@everywhere function theget(ψ, θ0, θ_lat, θ_lon)
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
    # 直接在叶素坐标系下通过矢量关系求解气动迎角
    # 并且通过aoaang函数来确保迎角范围为0~360°（见"\\src\\00_mathfunctions.jl"）
    α = Array{Float64}(Nb, Nbe)
    for k in 1:Nb
        for i in 1:Nbe
            windbe = rotobe(vall_r[k,i], β, θ[k,i]) # 叶素坐标系下来流的速度矢量
            windyz = [windbe[2],windbe[3]]
            α[k,i] = aoaang(windyz)
        end
    end
    return α
end

@everywhere function fcl(α, ma=0.0, Re=1e6, itp=clitp_na12)
    deg = α/π*180
    index_ma = ma*10+1
    index_deg = deg/5+1
    cl = itp[index_deg, index_ma]
    return cl
end

@everywhere function fcd(α, ma=0.0, Re=1e6, itp=cditp_na12)
    deg = α/π*180
    index_ma = ma*10+1
    index_deg = deg/5+1
    cd = itp[index_deg, index_ma]
    return cd
end

# @everywhere function clcdget(vall_r, α, Re=1e6)
#     # 本函数意图将cl和cd以数组的形式列举出来，以此来减少函数运算，从而提高性能
#     # 从目前Interpolations插值包的计算速度来看，本函数的意义不大
#     # 因此暂时弃用本函数
#     ma = Array{Float64}(Nb,Nbe)
#     for k in 1:Nb # 计算叶素微段马赫数
#         for i in 1:Nbe
#             ma[k,i] = norm(vall_r[k,i])/v_sound
#         end
#     end
#     clift = Array{Float64}(Nb,Nbe)
#     cdrag = Array{Float64}(Nb,Nbe)
#     for k in 1:Nb
#         for i in 1:Nbe
#             clift[k,i] = fcl(α[k,i], ma[k,i], Re)
#             cdrag[k,i] = fcd(α[k,i], ma[k,i], Re)
#         end
#     end
#     return clfit, cdrag, ma
# end
