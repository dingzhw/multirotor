# This script is for blade flap calculation
# 给出旋翼挥舞模型

@everywhere function staticbf(θ75, θtw, θlon, θlat, μ, λtpp, rotor=1)
    # static solution for blade flap
    # ---> use as the initilization values for Runge-Kutta Method

    β = zeros(npsi+1)
    dβ = zeros(npsi+1)
    ddβ = zeros(npsi+1)

    # parameters
    a = 2*π
    c = chroot
    γ = ρ*a*c*R^4/Iβ

    # simplified solution
    βlon = -θlon-((8/3)*μ*(θ75-(3/4)*λtpp))/(1+(3/2)*μ^2)
    β0 = γ*(θ75/8*(1+μ^2)-μ^2/60*θtw-λtpp/6+μ/6*(βlon+θlon))
    βlat = θlat-((4/3)*μ*β0)/(1+(1/2)*μ^2)

    for i in 1:npsi+1
        ψ = (-1)^(rotor-1)*(i-1)*dψ
        β[i] = β0+βlon*cos(ψ)+βlat*sin(ψ)
        dβ[i] = (-βlon*sin(ψ)+βlat*cos(ψ))*Ω
        ddβ[i] = (-βlon*cos(ψ)+βlat*sin(ψ))*Ω^2
    end

    β1   = β[Int64(1)]
    β2   = β[Int64(npsi/4)]
    β3   = β[Int64(npsi/2)]
    β4   = β[Int64(npsi/4*3)]
    β0   = ((β1+β3)/2+(β2+β4)/2)/2
    βlon = -(β1-β3)/2
    βlat = -(β2-β4)/2

    return β, dβ, ddβ, β0, βlon, βlat
end

@everywhere function bladeflap(β, dβ, ddβ, vber, chord, θ0, θ_lat, θ_lon, dr, rb, rotor=1)
    # 当k只取1的时候，认为所有的桨叶挥舞特性都是一致的
    # Initilize the value of β

    # rms recorder
    inum = 1
    rmsbeta = Float64[]
    θ  = theget(θ0, θ_lat, θ_lon, rotor)

    while true
        if inum>1000
            print("挥舞迭代次数过多！失败！\n")
            return false, rmsbeta
        end

        # vber = vbe(vind, β, dβ)
        α  = aoaget(vber, β, θ, rotor)

        # Runge Kutta 4th order method
        for i in 1:npsi
            t = (i-1)*dt
            # ddβ[i] = ddbeta(β[i], t, dt)
            f1 = ddbeta(β[i], t, dt, vber[i,:], chord, α[i,:], ddβ[i], θ[i,:], dr, rb)
            f2 = ddbeta(β[i]+1/2*dt*dβ[i]+1/8*dt^2*f1, t+1/2*dt, dt, vber[i,:], chord, α[i,:], ddβ[i], θ[i,:], dr, rb)
            f3 = ddbeta(β[i]+1/2*dt*dβ[i]+1/8*dt^2*f2, t+1/2*dt, dt, vber[i,:], chord, α[i,:], ddβ[i], θ[i,:], dr, rb)
            f4 = ddbeta(β[i]+dt*dβ[i]+1/2*dt^2*f3, dt, t+dt, vber[i,:], chord, α[i,:], ddβ[i], θ[i,:], dr, rb)
            β[i+1] = β[i]+dt*(dβ[i]+dt/6*(f1+f2+f3))
            dβ[i+1] = dβ[i]+dt/6*(f1+2*f2+2*f3+f4)
            ddβ[i+1] = f1
        end
        # print("=== β is $(β) ===\n
        # === dβ is $(dβ) ===\n
        # === ddβ is $(ddβ) ===\n\n")
        #
        β[1] 	= (β[1]+β[npsi+1])/2
        dβ[1]	= (dβ[1]+dβ[npsi+1])/2
        ddβ[1]	= (ddβ[1]+ddβ[npsi+1])/2

        # judge if flap converaged
        rmsβ = abs(β[1]-β[npsi+1]) +abs(dβ[1]-dβ[npsi+1])+
                    abs(ddβ[1]-ddβ[npsi+1])
        # print("+++++++++++++++++++++++++++++++++++\n")
        # print("+++++++++ rmsβ is $(rmsβ) +++++++++\n")
        # print("+++++++++++++++++++++++++++++++++++\n")
        push!(rmsbeta, rmsβ)
        if rmsβ<1e-1
            # print("%%%%%%%%%%%%%%%%%CONVERAGED%%%%%%%%%%%%%%%%\n")
            break
        end

        inum += 1
    end

    # 求出等效的纵横向挥舞角（用于配平）
    β1   = β[Int64(1)]
    β2   = β[Int64(npsi/4)]
    β3   = β[Int64(npsi/2)]
    β4   = β[Int64(npsi/4*3)]
    β0   = ((β1+β3)/2+(β2+β4)/2)/2
    βlon = -(β1-β3)/2
    βlat = -(β2-β4)/2

    return true, β, dβ, ddβ, β0, βlon, βlat, rmsbeta
end

@everywhere function ddbeta(β, t, dt, vber, α, θ, ddβ, chord, dr, rb)
    Maero = bladeaero(vber, chord, α, β, ddβ, θ, dr, rb)[4]
    ddbeta = -1/Iβ*(Sm*g*cos(β)+Sm*Ω^2*eflap*sin(β)+Iβ*Ω^2*cos(β)*sin(β)-Maero)

    return ddbeta
end
