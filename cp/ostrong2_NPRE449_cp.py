# Input (optional)
import sys
args = sys.argv

# Solving
import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import newton

# Post-processing/plotting
from matplotlib import pyplot as plt

## Fundamental Constants
PI = 3.14159265358979323
WATER = 'IF97::Water'
P_CRIT_WATER = 22.064 * 1e6 # MPa -> Pa
M_WATER = 18.05 # g / mol
K = lambda deg_C: deg_C + 273.15
g = 9.81 # m / s²

#### ---- User Configuration ---- ####
PWR = 1
BWR = 2
option = PWR # 1 for PWR, 2 for BWR
if (len(args) > 1):
    if args[1] in ['1', 'PWR']:
        option = PWR
    elif args[2] in ['2', 'BWR']:
        option = BWR
#### -END User Configuration ---- ####


## Shared setup
q_lin = lambda q0_lin, H, z: q0_lin*np.sin(PI*(z/H))

def calc_Re(D, G, mu):
    Re = G*D / mu
    return Re

def calc_f_1ph(Re):
    return 0.316*(Re**-0.25)

def calc_dP(G, drρdz, f, A, Pwet, ρ, g):
    return -((G**2)*drρdz + 0.5*f*(G**2)*Pwet/(ρ*A) + ρ*g)

def calc_dχedz(hfg, qlin, Af, G, dhfdP, χe, dhfgdP, dPdz):
    return (qlin/(Af*G) - ((dhfdP + χe*dhfgdP)*dPdz)) / hfg 

def propPT(key, P, T):
    return PropsSI(key, 'P', P, 'T', T, WATER)

def fgSatPropP(key, P):
    fProp = PropsSI(key, 'P', P, 'Q', 0, WATER)
    gProp = PropsSI(key, 'P', P, 'Q', 1, WATER)
    return fProp, gProp

def satPropP(key, P, χ):
    fP, gP = fgSatPropP(key, P)
    return χ*gP + (1-χ)*fP

def calc_Tf(χe, P, cp):
    hf = PropsSI('H', 'P', P, 'Q', 0, WATER)
    hg = PropsSI('H', 'P', P, 'Q', 1, WATER)
    Tsat = PropsSI('T', 'P', P, 'Q', 0, WATER)
    hfg = hg - hf
    Tf = χe*hfg/cp + Tsat
    return Tf

'''
Calcuates Nu_D using the Dittus-Bolter correlation.
This is designed for flow within a smooth circular tube, but in this case we assume it may apply to the flow around repeating tubes.

ASSUMES HEATING CONDITIONS, SUBCOOLED FLUID.
'''
def calc_Nu_DB(Re, Pr):
    n = 0.4 # HEATING
    return 0.023 * (Re**0.8) * (Pr**n)

'''
Assumes subcooled fluid
'''
def calc_h_DB(Re, L, P, T):
    Pr = propPT('PRANDTL', P, T)
    Nu = calc_Nu_DB(Re, Pr)
    k = propPT('CONDUCTIVITY', P, T)
    # Nu = hL/k
    # -> h = Nu k / L
    return Nu*k / L



def calc_Tw_LW(Re, L, P, Tf, χ, qpp, Tw_guess):
    Prl = satPropP('PRANDTL', P, 0)
    ρl, ρv = fgSatPropP('D', P)

    F = (1 + χ*Prl*(ρl/ρv - 1))**0.35
    S = 1/(1 + 0.055*(F**0.1)*(Re**0.16))
    # print(f"S: {S}")

    Tsat = PropsSI('T', 'P', P, 'Q', 0, WATER)
    # NuL = calc_Nu_DB(Re, Prl)
    # k = PropsSI('CONDUCTIVITY', 'P', P, 'Q', 0, WATER)
    # hL = NuL*k / L
    hL = calc_h_DB(Re, L, P, Tf)
    pr = P / P_CRIT_WATER

    hPool = 55 * (pr**0.12) * (qpp**(2/3)) * ((-np.log10(pr))**(-0.55)) * (M_WATER**(-0.5))
    # S = 0
    zero_func = lambda Tw: np.sqrt((F*hL*(Tw-Tf))**2 + (S*hPool*(Tw-Tsat))**2) - qpp
    Tw = newton(zero_func, Tw_guess)
    return Tw






## Variation

if option == 1:
    ## Givens for PWR
    H = 4.1 #m
    Drod = 0.94 / 100 # m, originally cm
    Pitch = 1.25 / 100 # m, originally cm
    Dfuel = 0.81 / 100 # ''
    tgap = 0.005 / 100 # ''
    Dgap = Dfuel + 2*tgap # '', diameter to where gap meets cladding
    kgap = 0.25 # W/mK
    kfuel = 2.6 # W/mK
    kclad = 21.5 # W/mK

    G = 3800 # kg / m²s
    qp0 = 380 * 100 # W / m, originally W/cm
    P0 = 15 * 1e6 # Pa, originally MPa
    Tf0 = K(277) # K, orignally ⁰C

    Af = Pitch**2 - (0.25 * PI * Drod**2)
    Pwet = PI*Drod
    Pheat = Pwet
    Lchar = Drod
    Lchar = 4 * Af / Pwet
    print(f"Fluid area Af: {Af} m²")
    print(f"Characteristic length: {Lchar} m")

    ## Properties constant 
    prop00 = lambda name: propPT(name, P0, Tf0)
    μ = prop00('viscosity')
    cp = prop00('CPMASS')
    ρ = prop00('D')
    print(f"Dynamic viscosity: {μ} Pa s")


    Re = calc_Re(Lchar, G, μ)
    print(f"Reynolds number: {Re} >> 10,000 -> Turbulent flow")
    print(f"Which leads to friction factor f = {calc_f_1ph(Re)}")

    n = 50000
    dz = H / n

    z = np.linspace(0, H, n)
    hfhg = lambda P: (
                PropsSI('H', 'P', P, 'Q', 0, WATER),
                PropsSI('H', 'P', P, 'Q', 1, WATER))

    ## P
    P = np.zeros(n)
    P[0] = P0
    f = calc_f_1ph(Re)
    Tf = np.zeros(n)
    Tf[0] = Tf0
    χe = np.zeros(n)
    hf0, hg0 = hfhg(P0)
    hfg0 = hg0-hf0
    Tsat0 = PropsSI('T', 'P', P0, 'Q', 0.5, WATER)
    χe[0] = cp*(Tf0-Tsat0)/hfg0
    print(f"T at inlet: {Tf0} K")
    print(f"cp at inlet: {cp} J/kgK")
    print(f'h_fg at inlet: {hfg0} J/kg')
    print(f"χe0: {χe[0]}")

    for i in range(n-1):
        dPdz = calc_dP(G, 0, f, G, Pwet, ρ, g)
        dP = dPdz * dz
        P[i+1] = P[i] + dP

        # Properties
        
        
        hfi, hgi = hfhg(P[i])
        hfgi = hgi - hfi

        hfn, hgn = hfhg(P[i])
        hfgn = hgn - hfn

        dhfdP = (hfn-hfi) / dP
        dhfgdP = (hfgn-hfgi) / dP

        qp = q_lin(qp0, H, z[i])
        dχedz = calc_dχedz(hfgi, qp, Af, G, dhfdP, χe[i], dhfgdP, dPdz)
        χe[i+1] = χe[i] + dχedz*dz
        Tf[i+1] = calc_Tf(χe[i+1], P[i+1], cp)

    ## Secondary temperature evaluations

    Tw = np.zeros(n)
    # Tw[0]
    Tgap = np.zeros(n) # At edge of gap and inner of clad
    Tfuel = np.zeros(n) # At outer of fuel and inner of gap
    Tcl = np.zeros(n) # At center of fuel
    Tsat = np.zeros(n)

    # hf = np.zeros(n) # heat transfer 

    ONE_PHASE = 1
    TWO_PHASE = 2
    ht_status = ONE_PHASE    
    for i in range(n):
        qp = q_lin(qp0, H, z[i])
        Tsat[i] = PropsSI('T', 'P', P[i], 'Q', 0, WATER)

        qpp = qp / Pheat
        ## Heat transfer, then 
        if ht_status == ONE_PHASE:
            h = calc_h_DB(Re, Lchar, P[i], Tf[i])
            # hf[i] = h
            Tw[i] = Tf[i] + qpp / h
        elif ht_status == TWO_PHASE:
            Tw[i] = calc_Tw_LW(Re, Lchar, P[i], Tf[i], 0, qpp, Tw[i-1]) # Calls previous Tw as a starting guess, but this is never called on z=0
        if Tw[i] >= Tsat[i]:
            ht_status = TWO_PHASE
        elif Tw[i] < Tsat[i]:
            ht_status = ONE_PHASE
        ## within clad
        # Tin = Tout + q'/2kπ * ln(Rout/Rin) for annulars
        Tgap[i] = Tw[i] + (qp / (2*kclad*PI)) * np.log(Drod / Dgap)
        ## Edge of fuel
        Tfuel[i] = Tgap[i] + (qp / (2*kgap*PI)) * np.log(Dgap / Dfuel)

        #Within fuel
        qpp = qp / (0.25 * PI * (Dfuel**2))
        Tcl[i] = Tfuel[i] + (qpp/(4*kfuel))*(0.25*(Dfuel**2))




    Tfd = np.zeros(n)
    Tfd[0] = Tf0
    ## T pass
    for i in range(n-1):
        qp = q_lin(qp0, H, z[i])
        dT = qp*dz / (G*Af*cp)
        Tfd[i+1] = Tfd[i] + dT

    ## Sanity check
    qtot_an = qp0*2*H/PI
    # qtot_nu = sum(q_lin(qp0, H, z))*dz
    # print(f"Qtot analytic: {qtot_an} W")
    # print(f"Qtot numeric:  {qtot_nu} W")
    m = G*Af
    hin = propPT('H', P0, Tf0)  
    hout = hin + qtot_an/m
    print(f"h goes from\n{hin} J/kg to\n{hout} J/kg")
    # Tfout = PropsSI('T', 'P', P[-1], 'H', hout)
    Tfout = Tf0 + qtot_an/(m*cp)
    print(f"Which translates to final Tf of {Tfout-K(0)}⁰C")






    plt.plot(Tf-K(0), z*100, label="Fluid T from χe")
    plt.plot(Tfd-K(0), z*100, linestyle='--', label="Fuid T (alternate calculation)")
    plt.plot(Tw-K(0), z*100, label="Outer wall T")
    plt.plot(Tgap-K(0), z*100, label="Inside clad T")
    plt.plot(Tfuel-K(0), z*100, label="Outside fuel T")
    plt.plot(Tcl-K(0), z*100, label="Fuel Centerline T")

    plt.plot(Tsat-K(0), z*100, linestyle='--', label="T$_{sat}$")
    plt.xlabel("T (⁰C)")
    plt.ylabel("z (cm)")
    plt.legend()
    plt.show()

    # print(f"χe: {χe}")
    # plt.plot(χe, z*100, label="Equilibrium quality χe")
    # plt.xlabel("χe")
    # plt.ylabel("z (cm)")
    # plt.legend()
    # plt.show()


    # plt.plot(P, z*100, label="Pressure")
    # plt.xlabel("P (Pa)")
    # plt.ylabel("z (cm)")
    # plt.legend()
    # plt.show()

    # plt.plot(hf, z*100, label="h (W/m²K)")
    # plt.xlabel("h (W/m²K)")
    # plt.ylabel("z (cm)")
    # plt.legend()
    # plt.show()
    
    pass
elif option == 2:


    pass
else:
    raise ValueError(f"'{option}' is not a valid case option (1 for PWR, 2 for BWR)")