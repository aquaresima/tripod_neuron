
include("../../base.jl")


t = TN.Tripod()

t.d[1].pm.g_ax
t.d[2].pm.g_ax
AdEx.gl

##
AMPA = TN.Esyn_dend.AMPA
NMDA = TN.Esyn_dend.NMDA
GABAa = TN.Esyn_dend.GABAa
GABAb = TN.Esyn_dend.GABAb

g(rec::TN.AbstractReceptor) = rec.gsyn * (1 / rec.τr⁻ + 1 / rec.τd⁻)
g(GABAa)
g(GABAb)
GABAa.E_rev
##

ui(u) = -(GABAa.E_rev * g(GABAa) + GABAb.E_rev * g(GABAb) - u * (g(GABAa) + g(GABAb)))

ue(u) = (-u) * (g(AMPA) + g(NMDA) * TN.NMDA_nonlinear(TN.Esyn_dend.NMDA, u))

ue(-44.0)
ui(-44.0)


pm = t.d[1].pm
gL(pm) =
    ((TN.AdEx.gl + pm.g_ax) + 1 / pm.Rm * (TN.AdEx.gl + 2 * pm.g_ax)) /
    (TN.AdEx.gl + 2 * pm.g_ax)
gL(t.d[1].pm)
gL(t.d[2].pm)

λi(λe, u, pm) = (λe * ue(u) + (AdEx.Er - u) * gL(pm)) / ui(u)

λi(1, -55.0, t.d[1].pm)
