# the script is for trim of the whole aircraft

@everywhere function trvind(vind1, vind2)
    # check if the calculation of vdiskind has converged

    remain = 0.0
    for i in 1:npsi
        for j in 1:Nbe
            remain += norm(vind2[i,j]-vind1[i,j])
        end
    end
    re = remain/npsi
    if re<=1e-1
        return true, re
    else
        return false, re
    end
end

@everywhere function trfm(θtdm1, θtdm2, solf::Function, trfunc::Function)
    # make the force and moment of the rotor have converged

    itrim = 1
    while true
        soltmp = solf(θtdm1, θtdm2)
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
                return soltmp, θtdm1, θtdm2
                break
            end
        end

        if abs(res[1])<=20 && abs(res[2])<=5 # || mean(abs.(res[3:6]))<=1e-1
            print("Trim Succeeded!\n")
            return soltmp, θtdm1, θtdm2
            break
        else
            mat = trfunc(θtdm1, θtdm2, solf)
            print("|--------------------------------|\n")
            print("====== The Jocobi Matrix ======\n")
            for i in 1:length(mat[:,1])
                print("$(mat[i,:])\n")
            end
            print("|--------------------------------|\n")
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

@everywhere function ntiter(θtdm1, θtdm2, solf::Function, ϵ=0.1*π/180)
    # Newton iteration method

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
        θ1  = copy(θtdm1)
        θ1_ = copy(θtdm1)
        θ2  = copy(θtdm2)
        θ2_ = copy(θtdm2)
        print("|--------------------------------|\n")
        print("====== The Control Values ======\n")
        print("Rotor1 : $(θ1)\n")
        print("Rotor2 : $(θ2)\n")
        print("|--------------------------------|\n")
        θ1[i] += ϵ
        θ1_[i] -= ϵ
        θ2[i] += ϵ
        θ2_[i] -= ϵ
        cba[i] = solf(θ1, θtdm2)
        cfo[i] = solf(θ1_, θtdm2)
        cba2[i] = solf(θtdm1, θ2)
        cfo2[i] = solf(θtdm1, θ2_)
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
