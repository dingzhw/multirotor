# the script is for trim of the whole aircraft

@everywhere function trvind(vind1, vind2)
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

@everywhere function trfm()
    return true
end
