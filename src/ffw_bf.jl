# This script is for blade flap calculation
# 给出旋翼挥舞模型

@everywhere function bladeflap(β, dβ, ddβ, vber, chord, θ0, θ_lat, θ_lon, dr, rb)
    # 当k只取1的时候，认为所有的桨叶挥舞特性都是一致的
    # Initilize the value of β

    # 挥舞迭代计数器
    inum = Int64(0)

    while true
        inum += 1

        if inum>50
            print("挥舞迭代次数过多！失败！")
            return false
        end

        # vber = vbe(vind, β, dβ)
        θ  = theget(θ0, θ_lat, θ_lon)
        α  = aoaget(vber, β, θ)
        for i in 2:(npsi+1) # 一周挥舞变化
            ψ = (i-1)*dψ
            ddβ[i] = ddβ[i-1]
            dβ[i] = dβ[i-1]+ddβ[i-1]*dt
            β[i] = β[i-1]+dβ[i-1]*dt+ddβ[i-1]*dt^2
            while true
                Mβ = bladeaero(vber[i-1,:], chord, α[i-1,:], β[i], ddβ[i], θ[i-1,:], dr, rb)[4]
                if abs(Mβ)<1 # 判断挥舞铰处总力矩是否为零
                    break
                end
                ddβ[i] += -Mβ/Iβ
            end
        end

        rmsβ = abs(β[1]-β[npsi+1])# +abs(dβ[1]-dβ[npsi+1])+
                #abs(ddβ[1]-ddβ[npsi+1])
        if rmsβ<1e-2
            break
        end
		β[1] 	= (β[1]+β[npsi+1])/2
		dβ[1]	= (dβ[1]+dβ[npsi+1])/2
		ddβ[1]	= (ddβ[1]+ddβ[npsi+1])/2
    end

    # 求出等效的纵横向挥舞角（用于配平）
    β1   = β[Int64(1)]
    β2   = β[Int64(npsi/4)]
    β3   = β[Int64(npsi/2)]
    β4   = β[Int64(npsi/4*3)]
    β0   = ((β1+β3)/2+(β2+β4)/2)/2
    βlon = -(β1-β3)/2
    βlat = -(β2-β4)/2

    return true, β, dβ, ddβ, β0, βlon, βlat
end
