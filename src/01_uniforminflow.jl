# The file calculate the uniform induced velocity
# 给出均匀入流模型

@everywhere function uniforminflow(T)
    vind_s = [0.0,0.0,-sqrt(abs(T/(2*ρ*A)))] # 固定坐标系诱导速度初值
    vall_s = v_air + vind_s # 固定坐标系总的滑流速度初值
    # lmdaui = -T/(2*ρ*A*norm(vall_s))/(Ω*R) # 前飞均匀诱导速度系数
    vindui_s = [0.0,0.0,-T/(2*ρ*A*norm(vall_s))]  # 第一步迭代的固定坐标系均匀入流值

    while true  # 迭代求解均匀入流
    if norm(vindui_s-vind_s)<=1e-3
      break
    else
      vind_s = vindui_s
      vall_s = v_air + vind_s
      vindui_s = [0.0,0.0,-T/(2*ρ*A*norm(vall_s))]
    end
    end
    lmdaui = vindui_s[3]/(Ω*R) # 收敛前飞均匀诱导速度系数

    vdiskind = Array{Vector}(npsi, Nbe)
    for i in 1:npsi
    	for j in 1:Nbe
    		vdiskind[i,j] = vall_s
    	end
    end

	return lmdaui, vall_s, vdiskind
end

@everywhere function vbeui(vall_s, β, dβ, rb, rotor=1)
    # 叶素当地位于旋转坐标系中的速度

    vall_r = Array{Vector}(npsi,Nbe)
    for i in 1:npsi
        ψ = (-1)^(rotor-1)*(i-1)*dψ
        for j in 1:Nbe
            vall_r[i,j] = (systoro(vall_s, ψ)+[0.0,(-1)^rotor*Ω*rb[j],0.0]+
                            betatoro([0,0,dβ[i]*rb[j]],β[i]))
        end
    end

    return vall_r
end
