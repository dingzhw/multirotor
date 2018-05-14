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
