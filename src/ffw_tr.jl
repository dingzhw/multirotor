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
                return soltmp, θcp, θ_lat, θ_lon
                break
            end
        end

        if abs(res[1])<=5 && abs(res[2])<=1e-2 && abs(res[3])<=1e-2
            print("Trim Succeeded!\n")
            return soltmp, θcp, θ_lat, θ_lon
            break
        else
            mat = trfunc(θcp, θ_lat, θ_lon, solf)
            optmp = inv(mat)*res

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

            itrim += 1
        end
    end
end

@everywhere function ntiter(θcp, θ_lat, θ_lon, solf::Function, ϵ=0.1*π/180)
    # Newton iteration method

    res  = zeros(3)
    mat  = zeros(3,3)

    # change the collective pitch
    cpback = solf(θcp+ϵ, θ_lat, θ_lon)
    cpforw = solf(θcp-ϵ, θ_lat, θ_lon)
    dT     = cpback[1] - cpforw[1]
    dβlon  = cpback[3] - cpforw[3]
    dβlat  = cpback[4] - cpforw[4]
    mat[:, 1] = [dT dβlon dβlat]/(2*ϵ)

    # change the lateral cyclic pitch
    clatba = solf(θcp, θ_lat+ϵ, θ_lon)
    clatfo = solf(θcp, θ_lat-ϵ, θ_lon)
    dT     = clatba[1] - clatfo[1]
    dβlon  = clatba[3] - clatfo[3]
    dβlat  = clatba[4] - clatfo[4]
    mat[:, 2] = [dT dβlon dβlat]/(2*ϵ)

    # change the longitudinal cyclic pitch
    clonba = solf(θcp, θ_lat, θ_lon+ϵ)
    clonfo = solf(θcp, θ_lat, θ_lon-ϵ)
    dT     = clonba[1] - clonfo[1]
    dβlon  = clonba[3] - clonfo[3]
    dβlat  = clonba[4] - clonfo[4]
    mat[:, 3] = [dT dβlon dβlat]/(2*ϵ)

    return mat
end
