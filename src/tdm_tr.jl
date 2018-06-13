# the script is for trim of the whole aircraft

@everywhere function trfm(θtdm::Array{Float64}, solf::Function, trfunc::Function)
    # make the force and moment of the rotor have converged

	θtdm1 = θtdm[1:3]
	θtdm2 = θtdm[4:6]

    itrim = 1
    while true
        soltmp = solf(θtdm)
        res    = [T-soltmp[3][1];
                  0.-soltmp[3][2];
                  0.-soltmp[1][3];
                  0.-soltmp[1][4];
                  0.-soltmp[2][3];
                  0.-soltmp[2][4]] # trim values difference
        if itrim%100 == 0 # in order to get rid of too many iteration steps
            print("Hard to Trim!\n")
            print("If the calculation shoulbe continue? (Type 'y' or 'n')\n")
            resp = readline()
            if resp =="n"
                print("Trim Failed!\n")
                return soltmp, θtdm
                break
            end
        end

        if abs(res[1])<=10 && abs(res[2])<=1 && mean(abs.(res[3:6]))<=1e-2
            print("Trim Succeeded!\n")
            return soltmp, θtdm
            break
        else
            print("\n*********第$(itrim)次尝试配平，很遗憾，暂时还没成功*********\n\n")
            mat = trfunc(θtdm, solf)
            # print("|--------------------------------|\n")
            # print("====== The Jocobi Matrix ======\n")
            # for i in 1:length(mat[:,1])
            #     print("$(mat[i,:])\n")
            # end
            # print("|--------------------------------|\n")
            optmp = inv(mat)*res

            # ensure all the control varibles are less than 90°
            for i in 1:length(θtdm1)
                θtdm1[i] += optmp[i]
                if abs(θtdm1[i])>=π/2
                    θtdm1[i] = -π/2+π*rand()
                end
                θtdm2[i] += optmp[i+length(θtdm1)]
                if abs(θtdm2[i])>=π/2
                    θtdm2[i] = -π/2+π*rand()
                end
            end
            itrim += 1
        end
    end
end

@everywhere function ntiter(θtdm, solf::Function, ϵ=0.5*π/180)
	# Newton iteration method

	θtdm1 = θtdm[1:3] # 旋翼 1 操纵量
	θtdm2 = θtdm[4:6] # 旋翼 2 操纵量

	# 参数初始化
    mat  = zeros(6,6)
    dT   = zeros(1,6)
    dQ   = zeros(1,6)
    dβlon1 = zeros(1,6)
    dβlat1 = zeros(1,6)
    dβlon2 = zeros(1,6)
    dβlat2 = zeros(1,6)
    cba    = Array{Any}(3)
    cfo    = Array{Any}(3)
    cba2   = Array{Any}(3)
    cfo2    = Array{Any}(3)

    for i in 1:3
		# 为防止操纵量数组被指针修改，用copy函数复制一个位于不同地址的相同数组
        θ1  = copy(θtdm1)
        θ1_ = copy(θtdm1)
        θ2  = copy(θtdm2)
        θ2_ = copy(θtdm2)
        # print("|--------------------------------|\n")
        # print("====== The Control Values ======\n")
        # print("Rotor1 : $(θ1)\n")
        # print("Rotor2 : $(θ2)\n")
        # print("|--------------------------------|\n")
        θ1[i] += ϵ
        θ1_[i] -= ϵ
        θ2[i] += ϵ
        θ2_[i] -= ϵ
        cba[i] = solf([θ1;θtdm2])
        cfo[i] = solf([θ1_;θtdm2])
        cba2[i] = solf([θtdm1;θ2])
        cfo2[i] = solf([θtdm1;θ2_])
        dT[i]  = cba[i][3][1] - cfo[i][3][1]
        dQ[i]  = cba[i][3][2] - cfo[i][3][2]
        dβlon1[i] = cba[i][1][3] - cfo[i][1][3]
        dβlat1[i] = cba[i][1][4] - cfo[i][1][4]
        dβlon2[i] = cba[i][2][3] - cfo[i][2][3]
        dβlat2[i] = cba[i][2][4] - cfo[i][2][4]
        dT[i+length(θ1)]  = cba2[i][3][1] - cfo2[i][3][1]
        dQ[i+length(θ1)]  = cba2[i][3][2] - cfo2[i][3][2]
        dβlon1[i+length(θ1)] = cba2[i][1][3] - cfo2[i][1][3]
        dβlat1[i+length(θ1)] = cba2[i][1][4] - cfo2[i][1][4]
        dβlon2[i+length(θ1)] = cba2[i][2][3] - cfo2[i][2][3]
        dβlat2[i+length(θ1)] = cba2[i][2][4] - cfo2[i][2][4]
    end

    mat = [dT;
           dQ;
           dβlon1;
           dβlat1;
           dβlon2;
           dβlat2]/(2*ϵ)

    return mat
end
