# The file calculate the uniform induced velocity
# 给出均匀入流模型

@everywhere function uniforminflow(T)
    vind_s = [0.0,0.0,-T/(2*ρ*A*norm(v_air))] # 固定坐标系诱导速度初值
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

	return lmdaui, vall_s
end

@everywhere function vallr(vall_s, ψ, β, dβ)
    vall_r = Array{Vector}(Nb,Nbe)
    for k in 1:Nb # 叶素当地来流速度（包含诱导速度、前方来流、挥舞流动以及旋转来流）
        ψ = ψ+(k-1)*2*π/Nb
        for i in 1:Nbe
            vall_r[k,i] = (systoro(vall_s, ψ)+[0.0,-Ω*rb[k,i],0.0]+
                            betatoro([0,0,dβ*rb[k,i]],β))
        end
    end
    return vall_r
end
