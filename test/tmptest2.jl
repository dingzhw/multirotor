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
        # return soltmp, θtdm1, θtdm2
        # # break
    end
end

if abs(res[1])<=20 && abs(res[2])<=10 # || mean(abs.(res[3:6]))<=1e-1
    print("Trim Succeeded!\n")
    # return soltmp, θtdm1, θtdm2
    # break
else
    mat = trfunc(θtdm1, θtdm2, solf)
    print("|--------------------------------|\n")
    print("====== The Jocobi Matrix ======\n")
    for i in 1:length(mat[:,1])
        print("$(mat[1,:])\n")
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
