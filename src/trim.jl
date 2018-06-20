# the script is for trim of the whole rotor using "wind trim rules"

@everywhere function trfm(θcp, θ_lat, θ_lon, solf::Function, trfunc::Function)
    # make the force and moment of the rotor have converged

    itrim = 1
    while true
        soltmp = solf(θcp, θ_lat, θ_lon)
        res    = [T-soltmp[1];
                  0.-soltmp[3];
                  0.-soltmp[4]] # trim values difference
        if itrim%100 == 0 # in order to get rid of too many iteration steps
            print("Hard to Trim!\n")
            print("If the calculation shoulbe continue? (Type 'y' or 'n')\n")
            resp = readline()
            if resp =="n"
                print("Trim Failed!\n")
                return soltmp, θcp, θ_lat, θ_lon, optmp
                break
            end
        end

        if abs(res[1])<=1 && abs(res[2])<=2e-2 && abs(res[3])<=2e-2
            print("Trim Succeeded!\n")
            return soltmp, θcp, θ_lat, θ_lon, optmp
            break
        else
            # print("+++++++++++++++++++++++++= test line = ++++++++++++++\n\n\n")
            mat = trfunc(θcp, θ_lat, θ_lon, solf)
            dmat = det(mat)
            # print("奇异矩阵？？？ $(dmat) ------------\n\n\n")
            if abs(dmat) == 0.0
                θcp += 0.05*rand()
                θ_lat = 0.05*rand()
                θ_lon = 0.05*rand()
            else
                optmp = inv(mat)*res
                # print("++++++ 优化步进值为：$(optmp) ++++++\n")

                # ensure all the control varibles are less than 90°
                θcp += optmp[1]
                if abs(θcp) >= π/2
                    θcp = -π/2+π*rand()
                end
                θ_lat += optmp[2]
                if abs(θ_lat) >= π/2
                    θ_lat = -π/2+π*rand()
                end
                θ_lon += optmp[3]
                if abs(θ_lon) >= π/2
                    θ_lon = -π/2+π*rand()
                end

                # # 配平稳定收敛方法
                # itrim = 1
                # while itrim<=10
                #     # trim acceleration
                #     resba = res[1]^2+res[2]^2+res[3]^2
                #     soltmp = solf(θcp, θ_lat, θ_lon)
                #     res    = [T-soltmp[1];
                #                 0.-soltmp[3];
                #                 0.-soltmp[4]]
                #     resfo = res[1]^2+res[2]^2+res[3]^2
                #
                #     if resfo<resba
                #         break
                #     else
                #         optmp = 0.5*optmp
                #         θcp += optmp[1]
                #         if abs(θcp) >= π/2
                #             θcp = -π/2+π*rand()
                #         end
                #         θ_lat += optmp[2]
                #         if abs(θ_lat) >= π/2
                #             θ_lat = -π/2+π*rand()
                #         end
                #         θ_lon += optmp[3]
                #         if abs(θ_lon) >= π/2
                #             θ_lon = -π/2+π*rand()
                #         end
                #     end
                #     itrim += 1
                # end

            end
        end

        print("=== θcp is $(θcp) ====\n
        === θlat is $(θ_lat) ===\n
        === θlon is $(θ_lon) ===\n\n")
        itrim += 1
    end
end

@everywhere function ntiter(θcp, θ_lat, θ_lon, solf::Function, ϵ=0.1*π/180)
    # Newton iteration method

    mat  = zeros(3,3)

    # change the collective pitch
    cpback = solf(θcp+ϵ, θ_lat, θ_lon)
    cpforw = solf(θcp-ϵ, θ_lat, θ_lon)
    dT     = cpback[1] - cpforw[1]
    dβlon  = cpback[3] - cpforw[3]
    dβlat  = cpback[4] - cpforw[4]
    mat[:,1] = [dT dβlon dβlat]/(2*ϵ)

    # change the lateral cyclic pitch
    clatba = solf(θcp, θ_lat+ϵ, θ_lon)
    clatfo = solf(θcp, θ_lat-ϵ, θ_lon)
    dT     = clatba[1] - clatfo[1]
    dβlon  = clatba[3] - clatfo[3]
    dβlat  = clatba[4] - clatfo[4]
    mat[:,2] = [dT dβlon dβlat]/(2*ϵ)

    # change the longitudinal cyclic pitch
    clonba = solf(θcp, θ_lat, θ_lon+ϵ)
    clonfo = solf(θcp, θ_lat, θ_lon-ϵ)
    dT     = clonba[1] - clonfo[1]
    dβlon  = clonba[3] - clonfo[3]
    dβlat  = clonba[4] - clonfo[4]
    mat[:,3] = [dT dβlon dβlat]/(2*ϵ)

    return mat
end
