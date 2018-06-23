# This file defines the const varibles

const ρ = 1.225 #空气密度 量纲kg*m^-3
const ν = 1.46e-4 #运动黏性系数 无量纲
const g = 9.8 #重力加速度  量纲kg*m*s^-2
const v_sound = 343.2 #音速 量纲m/s
const Nro = Int64(2) #Number of rotors
const Nb = Int64(4)  #Number of blades
const R = 0.8255 #旋翼半径 量纲m
const chroot = 0.095*R # 桨根弦长
const airfoil = "NACA0012" # Airfoil
const ecut = 0.24 #桨叶根切比例 无量纲
const eflap = 0.05 #桨叶挥舞铰偏置量，量纲m
const m_ = 0.222 #桨叶质量密度 量纲kg/m
const Sm = 1/2*(R-eflap)^2*m_ # mass inertia
const Ω = 230  #旋翼转速 量纲rad/s
const Vtip = Ω*R # blade tip velocity
const αs = -0.0*π/180.0  # 旋翼轴倾角  量纲rad
const Kβ = 0.0  # 桨叶根部挥舞弹簧刚度 量纲？？？
const vair = 0.0 # 来流速度 量纲 m/s
const T = 1000/cos(αs) #飞行器重量 (量纲为kg*m*s^-2)
const dpsideg = 10.0  # 方位角步进长度（量纲为deg）
const betap = 0.0/180*π # Precone
# const βang0 = 0.0/180*π # 挥舞角初值
# const dβ0 = 0.0/180*π # 挥舞角速度初值
# const θ7 = 2.0/180*π # 旋翼总距
# const θtw = 0.0/180*π # 桨叶扭转角
# const thelat = -0.0/180*π # 横向周期变距
# const thelon = 0.0/180*π # 纵向周期变距
const Nbe = Int64(10) #旋翼分段数量
const taper = 1.0   # 桨叶尖削比
const Iβ = m_/3*R^3(1-eflap/R)^3  #桨叶挥舞惯量 量纲kg*m^2
const v_air = [vair*cos(αs),0.0,vair*sin(αs)] #forward wind speed  量纲m/s
const μ_air = v_air[1]/(Ω*R)  #来流在桨盘edgewise方向分量  无量纲
const λ_air = v_air[3]/(Ω*R)  #来流在桨盘轴向分量  无量纲
const A = π*R^2 #参考面积 量纲为m^2
const σ = Nb*R*chroot/A
const fnonc = ρ*A*Ω*R^2*Ω #力的无量纲化参数 量纲kg*m*s^-2
const mnonc = ρ*A*Ω^2*R^3 #力矩的无量纲化参数 量纲kg*m^2/s^2
# const CT = T/fnonc  #无量纲重量系数
const dψ = dpsideg*π/180 #方位角步进步长 (量纲为rad)
const dt = dψ/Ω # 方位角步进时间 （量纲为s）
const npsi = Int64(360/dpsideg) # 周向分割步数

# 双旋翼参数设置
const cut = 10*R # wake cutoff distance
const disr = 1.5*R # the distance between two rotors in Y coordination
const hr = 0.2*R # the distance between two rotors in Z coordination

# define calculation vars for ffw
const ϵr = 0.7 # the distance from rotor center to control points cycle
const A0 = [1.,0.,0]*ϵr*R # new released vortex ring control points
const B0 = [0.,1.,0]*ϵr*R
const C0 = [-1,0.,0]*ϵr*R
const D0 = [0.,-1,0]*ϵr*R
const A01 = A0+[0., disr, hr] # control points of the other tandem rotor
const B01 = B0+[0., disr, hr]
const C01 = C0+[0., disr, hr]
const D01 = D0+[0., disr, hr]

# set paremeters for ga module
const ncre = 20
const nmax = 20*20
const nmin = 20
const npar = 6
const mr = 0.5
const rpr = 0.8

const ξ = 0.6 # triming weighting coefficient
const lossξ = 6 # 针对劣势个体的惩罚因子
const cfit = 0.01 # if minfits is less than cfit, then convergence is reached too
