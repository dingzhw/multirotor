# NACA 0012 airfoil data import and interpolate
@everywhere using Interpolations;
# @everywhere function caeroimport(filename=pwd()*"\\input\\naca0012_cl")
@everywhere function caeroimport(filename=pwd()*"//input//naca0012_cl")
    # 本函数的作用是读取相应的数据文件，并由此生成插值序列
    cafile = open(filename,"r")
    calines = readlines(cafile)
    rownum = length(calines)
    catmp = Float64[]
    for i in 1:rownum
        if calines[i] == ""
            break
        end
        cama = split(calines[i],"\t")
        for j in 2:length(cama)-1
            camaf = parse(Float64,cama[j])
            append!(catmp, camaf)
        end
    end
    colnum = Int64(length(catmp)/rownum)
    cacof = reshape(catmp, colnum, rownum)
    caitp = interpolate(cacof, BSpline(Linear()), OnGrid())
    return caitp
end

@everywhere clitp_na12 = caeroimport() # 生成升力系数插值序列
# cditp_na12 = caeroimport(pwd()*"\\input\\naca0012_cd") # 生成阻力系数插值序列
@everywhere cditp_na12 = caeroimport(pwd()*"//input//naca0012_cd") # 生成阻力系数插值序列

# 求解aoa以及对应的cl,cd
@everywhere function the0(θcp, twistr, twist1, twist2, rb)
    # get the collective angle in the blade root
    # --- θcp = θ_install + θ_collectivepitch

    θ0 = Array{Float64}(Nbe)

    for i in 1:Nbe
        if rb[i]<=twistr
            θ0[i] = θcp+twist1*rb[i]/twistr
        else
            θ0[i] = θcp+twist1+twist2*((rb[i]-twistr)/(R-twistr))
        end
    end

    θ75 = θ0[trunc(Int,Nbe*0.75)]

    return θ0, θ75
end

@everywhere function theget(θ0, θ_lat, θ_lon, rotor=1)
    θ = Array{Float64}(npsi,Nbe)
    for i in 1:npsi
        ψ = (-1)^(rotor-1)*(i-1)*dψ
        for j in 1:Nbe
            θ[i,j] = θ0[j]+θ_lat*cos(ψ)+θ_lon*sin(ψ)
        end
    end
    return θ
end

@everywhere function aoaget(vber, β, θ, rotor=1)
    # 直接在叶素坐标系下通过矢量关系求解气动迎角
    # 并且通过aoaang函数来确保迎角范围为0~360°（见"\\src\\00_mathfunctions.jl"）
    α = Array{Float64}(npsi,Nbe)

    for i in 1:npsi
        for j in 1:Nbe
            windbe = rotobe(vber[i,j], β[i], θ[i,j], rotor) # 叶素坐标系下来流的速度矢量
            windyz = [windbe[2],windbe[3]]
            α[i,j] = aoaang(windyz)
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
#             ma[i] = norm(vall_r[i])/v_sound
#         end
#     end
#     clift = Array{Float64}(Nb,Nbe)
#     cdrag = Array{Float64}(Nb,Nbe)
#     for k in 1:Nb
#         for i in 1:Nbe
#             clift[i] = fcl(α[i], ma[i], Re)
#             cdrag[i] = fcd(α[i], ma[i], Re)
#         end
#     end
#     return clfit, cdrag, ma
# end
