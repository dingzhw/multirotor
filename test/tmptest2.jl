θcp = 2.0*π/180
θ_lat = 0.0
θ_lon = 0.0
# T0 = T
twistr = 0.8*R
twist1 = -12.0*π/180
twist2 = 0.0
chord = zeros(Nbe)
rb = zeros(Nbe) # position of blade control points
dr = zeros(Nbe) # distance of each two blade control points
cbe = 0.8   # Range[0,1] for adjust the denisty of rd[] in tip, lager for more dense
βlon = 0.0
βlat = 0.0

for i in 1:Nbe
    chord[i] = chroot*taper*(i-1)/Nbe
    dr[i] = (R*(1-ecut)*(sin(i/Nbe*cbe*π/2)-
                sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2))
    rb[i] = (ecut*R+R*(1-ecut)*(sin(i/Nbe*cbe*π/2)+
                sin((i-1)/Nbe*cbe*π/2))/sin(cbe*π/2)/2)
end

θ0 = the0(θcp, twistr, twist1, twist2, rb)
