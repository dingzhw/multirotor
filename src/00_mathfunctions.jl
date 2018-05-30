# 计算中需要用到的数学函数

# 坐标系转换矩阵函数
@everywhere function Mstor(ψ) # system to rotation
    Msr =  [cos(ψ) sin(ψ) 0;
            -sin(ψ) cos(ψ) 0;
            0   0   1]
    return Msr
end

@everywhere function Mrtobeta(β) # rotation to beta
    Mrbt = [cos(β) 0 sin(β);
            0   1   0;
            -sin(β) 0 cos(β)]
    return Mrbt
end

@everywhere function Mbetatothe(θ, rotor=1) # beta to theta
    if rotor == 1
        Mbthe = [-1 0 0;
                  0 -cos(θ) -sin(θ);
                  0 -sin(θ) cos(θ)]
    else
        Mbthe = [1 0 0;
                 0 cos(θ) -sin(θ);
                 0 sin(θ) cos(θ)]
    end
    return Mbthe
end

# 坐标系转换的具体实现
@everywhere function systoro(vec::Array,ψ,M=Mstor) # sys to rotation
    mro = zeros(Float64,3)
    mro = M(ψ)*vec
    return mro
end

@everywhere function rotosys(vec::Array,ψ,M=Mstor) # ratation to sys
    msys = zeros(Float64,3)
    msys = inv(M(ψ))*vec
    return msys
end

@everywhere function rotobeta(vec::Array,β,M=Mrtobeta) # rotation to beta
    mbeta = zeros(Float64,3)
    mbeta = M(β)*vec
    return mbeta
end

@everywhere function betatoro(vec::Array,β,M=Mrtobeta) # beta to rotation
    mro = zeros(Float64,3)
    mro = inv(M(β))*vec
    return mro
end

@everywhere function betatobe(vec::Array, θ, rotor=1, M=Mbetatothe) # beta to theta
    mthe = zeros(Float64,3)
    mthe = M(θ, rotor)*vec
    return mthe
end

@everywhere function betobeta(vec::Array, θ, rotor=1, M=Mbetatothe) # theta to beta
    mbeta = zeros(Float64,3)
    mbeta = inv(M(θ, rotor))*vec
    return mbeta
end

@everywhere function rotobe(vec::Array, β, θ, rotor=1, M1=Mrtobeta, M2=Mbetatothe)
    # rotation to theta
    mthe = zeros(Float64,3)
    mthe = M2(θ, rotor)*(M1(β)*vec)
    return mthe
end

@everywhere function betoro(vec::Array, β, θ, rotor=1, M1=Mrtobeta, M2=Mbetatothe)
    # theta to rotation
    mro = zeros(Float64,3)
    mro = inv(M1(β))*(inv(M2(θ, rotor))*vec)
    return mro
end

# 矢量夹角函数
@everywhere function aoaang(vec_::Array{Float64,1},x_=[1.0,0.0]) # 二维坐标系下，任意矢量与轴夹角(-180deg~180deg)
  if vec_[2]>=x_[2] # 判断矢量与x轴夹角方向
    return acos(dot(vec_,x_)/(norm(vec_)*norm(x_)))
  else
    return 2*π-acos(dot(vec_,x_)/(norm(vec_)*norm(x_)))
  end
end

@everywhere function cosvecang(a,b) #计算夹角余弦值
  if norm(a)<=1e-3||norm(b)<=1e-3
    return 0.0
  elseif norm(a+b)<1e-3||norm(a-b)<1e-3
    return cos(π)
  elseif norm(a-b)<1e-3
    return cos(0.0)
  else
    return dot(a,b)/(norm(a)*norm(b))
  end
end

# Tensor operation functions

@everywhere function tstrans(ts)
    row = length(ts[:,1])
    col = length(ts[1,:])

    tst = Array{Any}(col, row)

    for i in 1:row
        for j in 1:col
            tst[j,i] = ts[i,j]
        end
    end

    return tst
end
