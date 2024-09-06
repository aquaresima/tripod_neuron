
# Copyright (c) 2022 Alessio Quaresima
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT


import Unitful: μm, cm, m, Ω, GΩ, F, μF, pF
module MyUnits
using Unitful
@unit Sim "Sim" Siemens 1u"1/Ω" true
end

import Unitful: @u_str
using Unitful
Unitful.register(MyUnits)

struct Physiology
    Ri::typeof(1.0 * Ω * cm)
    Rd::typeof(1.0 * Ω * cm^2)
    Cd::typeof(1.0 * μF / cm^2)
end

function G_axial(; Ri = Ri, d = d, l = l)
    l_ = uconvert(cm, d)
    d_ = uconvert(cm, l)
    R_ = Ri * l / (π * d * d / 4)
    return uconvert(u"nSim", 1 / R_)
end

function G_mem(; Rd = Rd, d = d, l = l)
    l_ = uconvert(cm, d)
    d_ = uconvert(cm, l)
    R_ = Rd / l_ / d_ / π
    return uconvert(u"nSim", 1 / R_)
end

function C_mem(; Cd = Cd, d = d, l = l)
    l_ = uconvert(cm, d)
    d_ = uconvert(cm, l)
    # return C = uconvert(pF*cm*cm,Cd)
    C_ = Cd * π * d * l
    return uconvert(pF, C_)
end
