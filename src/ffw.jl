# this script is for Fast Free Wake Calculation
# the result of this script should be the induced
# --- velocity of the whole rotor pane

# using packages
using Elliptic

# define calculation vars
const ϵr = 0.7
const A0 = [1.,0.,0]*0.7*R
const B0 = [0.,1.,0]*0.7*R
const C0 = [-1,0.,0]*0.7*R
const D0 = [0.,-1,0]*0.7*R

################################################################################

@everywhere function gamget(T, kΓ=1.2, vrnum=4)
    # Calculation of Γ of vortex ring
    # Γ of vortex ring is calculate at the instant of its release
    vad = sqrt(T/(2*ρ*A))
    kp  = R/(vad+abs(vair))
    dτ  = kp/vrnum
    return -4*kΓ*kp*T/(ρ*Vtip*A*σ), dτ
end

@everywhere type vortexring1
    # define a TYPE of vortex ring
    # what is confuesd here is that the control points may
    # --- transfer with the local velocity and become sth
    # --- but not a square, should it be handeled like this?

    # vortex strength
    Γ::Float64

    # four control points
    A::Vector
    B::Vector
    C::Vector
    D::Vector

    # # effects
    # P::Array{Vector}
    # # vind::Vector

    # functions
    croc::Function
    vrtosys::Function
    systovr::Function
    vringvind::Function
end

@everywhere function croc(vr::vortexring1)
    # calculate the center and radius of the circle



    A  = vr.A
    B  = vr.B
    C  = vr.C
    D  = vr.D
    O  = 1/4*(A+B+C+D)
    OA = A-O
    OB = B-O
    OC = C-O
    OD = D-O
    R  = 1/(4*ϵr)*(norm(OA)+norm(OB)+norm(OC)+norm(OD))
    return O, R
end

@everywhere function vringvind(vr::vortexring1,P)
    # calculation of the induced velocity of vortex ring
    # --- in its local coordination and in system coordination

    vind_array = Vector[]
    pvr = vr.systovr(vr,P)
    for i in 1:length(P)
        R   = vr.croc(vr)[2]
        p   = pvr[i]
        pxy = [p[1],p[2]]
        z   = p[3]
        r   = norm(pxy)
        δ   = 0.1*R # vortex core radius
        Γ   = vr.Γ

        # r and R is different
        # should fix r
        A = (r-R)^2+z^2+δ^2
        a = sqrt((r+R)^2+z^2+δ^2)
        m = 4*r*R/a^2

        ur = Γ/(2*π*a)*z/r*((r^2+R^2+z^2+δ^2)/A*Elliptic.E(m)-Elliptic.K(m))
        uz = Γ/(2*π*a)*(-(r^2-R^2+z^2+δ^2)/A*Elliptic.E(m)+Elliptic.K(m))

        p_ = pxy/r
        vind = [ur*p_[1],ur*p_[2],uz]
        push!(vind_array,vind)
    end

    return vind_array
end

@everywhere function systovr(vr::vortexring1,P)
    # transfer vector in system coordination to vortex ring local coordination

    a = vr.A
    b = vr.B
    c = vr.C
    d = vr.D

    ca = a-c
    db = b-d

    xvr = ca/norm(ca)
    yvr = db/norm(db)
    zvr = cross(xvr, yvr)

    m = [   xvr[1] xvr[2] xvr[3];
            yvr[1] yvr[2] yvr[3];
            zvr[1] zvr[2] zvr[3]
            ]

    pvr_array = Vector[]
    for i in 1:length(P)
        pvr = m*P[i]
        push!(pvr_array,pvr)
    end

    return pvr_array
end

@everywhere function vrtosys(vr::vortexring1,P)
    # transfer vector in vortex ring coordination to system coordination

    a = vr.A
    b = vr.B
    c = vr.C
    d = vr.D

    ca = a-c
    db = b-d

    xvr = ca/norm(ca)
    yvr = db/norm(db)
    zvr = cross(xvr, yvr)

    m = [   xvr[1] xvr[2] xvr[3];
            yvr[1] yvr[2] yvr[3];
            zvr[1] zvr[2] zvr[3]
            ]

    mt = inv(m)
    vind = vr.vringvind(vr,P)
    vindsys_array = Vector[]
    for i in 1:length(vind)
        vindsys = mt*vind[i]
        push!(vindsys_array, vindsys)
    end

    return vindsys_array
end

@everywhere function vbe(vind, β, dβ, rb)
    # calculate the blade elements velocity

    # vinds = reshape(vind, Nbe, npsi)
    # transpose!(vinds)

    vall_s = Array{Vector}(npsi,Nbe)
    for i in 1:npsi
        for j in 1:Nbe
            vall_s[i,j] = v_air+vind[i,j]
        end
    end

    vbe = Array{Vector}(npsi,Nbe)
    for i in 1:npsi
        ψ = (i-1)*dψ
        for j in 1:Nbe
            vbe[i,j] = (systoro(vall_s[i,j], ψ)+[0.0,-Ω*rb[j],0.0]+
                            betatoro([0,0,dβ[i]*rb[j]],β[i]))
        end
    end

    return vbe, vall_s
end
