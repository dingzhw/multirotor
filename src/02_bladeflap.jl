# This script is for blade flap calculation
# 给出旋翼挥舞模型

@everywhere function bladeflap(β, dβ, ddβ)
    # Initilize the value of β

    # 挥舞迭代计数器
    inum = Int8(0)

    while true
        inum += 1

        if inum>50
            print("挥舞迭代次数过多！中断！")
            return false
        end

        for i in 2:(npsi+1) # 一周挥舞变化
            ψ = i*dψ
            dβ[i] = dβ[i-1]+ddβ[i-1]*dt
            β[i] = β[i-1]+dβ[i-1]*dt+ddβ[i-1]*dt^2
            ddβ[i] = ddβ[i-1]
            while true
                Mβ = rotorforce()[5]
                if abs(Mβ)<1e-3 # 判断挥舞铰处总力矩是否为零
                    break
                end
                ddβ += Mβ/Iβ
            end
        end

        rmsβ = abs(β[1]-β[npsi+1])+abs(dβ[1]-dβ[npsi+1])+abs(
                    ddβ[1]-ddβ[npsi+1])
        if rmsβ<1e-3
            break
        end
    end
    return true, β, dβ, ddβ
end
