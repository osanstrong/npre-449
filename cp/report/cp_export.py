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
K = lambda deg_C: deg_C + 237.15
g = 9.81 # m / s2

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

def calc_dP(G, drrhodz, f, A, Pwet, rho, g):
    return -((G**2)*drrhodz + 0.5*f*(G**2)*Pwet/(rho*A) + rho*g)

def calc_dchiedz(hfg, qlin, Af, G, dhfdP, chie, dhfgdP, dPdz):
    return (qlin/(Af*G) - ((dhfdP + chie*dhfgdP)*dPdz)) / hfg 

def propPT(key, P, T):
    return PropsSI(key, 'P', P, 'T', T, WATER)

def calc_Tf(chie, P, cp):
    hf = PropsSI('H', 'P', P, 'Q', 0, WATER)
    hg = PropsSI('H', 'P', P, 'Q', 1, WATER)
    Tsat = PropsSI('T', 'P', P, 'Q', 0, WATER)
    hfg = hg - hf
    Tf = chie*hfg/cp + Tsat
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

    G = 3800 # kg / m2s
    qp0 = 380 * 100 # W / m, originally W/cm
    P0 = 15 * 1e6 # Pa, originally MPa
    Tf0 = K(277) # K, orignally degC

    Af = Pitch**2 - (0.25 * PI * Drod**2)
    Pwet = PI*Drod
    Pheat = Pwet
    Lchar = Drod
    Lchar = 4 * Af / Pwet
    print(f"Characteristic length: {Lchar} m")

    ## Properties constant 
    prop00 = lambda name: propPT(name, P0, Tf0)
    mu = prop00('viscosity')
    cp = prop00('CPMASS')
    rho = prop00('D')
    print(f"Dynamic viscosity: {mu} Pa s")


    Re = calc_Re(Lchar, G, mu)
    print(f"Reynolds number: {Re} >> 10,000 -> Turbulent flow")
    print(f"Which leads to friction factor f = {calc_f_1ph(Re)}")

    n = 500
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
    chie = np.zeros(n)
    hf0, hg0 = hfhg(P0)
    hfg0 = hg0-hf0
    Tsat0 = PropsSI('T', 'P', P0, 'Q', 0.5, WATER)
    chie[0] = cp*(Tf0-Tsat0)/hfg0

    for i in range(n-1):
        dPdz = calc_dP(G, 0, f, G, Pwet, rho, g)
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
        dchiedz = calc_dchiedz(hfgi, qp, Af, G, dhfdP, chie[i], dhfgdP, dPdz)
        chie[i+1] = chie[i] + dchiedz*dz
        Tf[i+1] = calc_Tf(chie[i+1], P[i+1], cp)

    ## Secondary temperature evaluations

    Tw = np.zeros(n)
    Tgap = np.zeros(n) # At edge of gap and inner of clad
    Tfuel = np.zeros(n) # At outer of fuel and inner of gap
    Tcl = np.zeros(n) # At center of fuel
    Tsat = np.zeros(n)

    ONE_PHASE = 1
    TWO_PHASE = 2
    ht_status = ONE_PHASE    
    for i in range(n):
        qp = q_lin(qp0, H, z[i])
        Tsat[i] = PropsSI('T', 'P', P[i], 'Q', 0, WATER)


        ## Heat transfer, then 
        if ht_status == ONE_PHASE:
            h = calc_h_DB(Re, Lchar, P[i], Tf[i])
            # q'' = h (Tw - Tf)
            # -> Tw = Tf + q''/h

            # q'' = q' / Pheat
            qpp = qp / Pheat
            Tw[i] = Tf[i] + qpp / h

            if Tw[i] >= Tsat[i]:
                ht_status = TWO_PHASE
        elif ht_status == TWO_PHASE:

            pass

        ## within clad
        # Tin = Tout + q'/2kpi * ln(Rout/Rin) for annulars
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
    print(f"Which translates to final Tf of {Tfout-K(0)}degC")






    plt.plot(Tf-K(0), z*100, label="Fluid T from chie")
    plt.plot(Tfd-K(0), z*100, linestyle='--', label="Fuid T from direct calculation")
    plt.plot(Tw-K(0), z*100, label="Outer wall T")
    plt.plot(Tgap-K(0), z*100, label="Inside clad T")
    plt.plot(Tfuel-K(0), z*100, label="Outside fuel T")
    plt.plot(Tcl-K(0), z*100, label="Fuel Centerline T")

    plt.plot(Tsat-K(0), z*100, label="T$_{sat}$")
    plt.xlabel("T (degC)")
    plt.ylabel("z (cm)")
    plt.legend()
    plt.show()
    
    plt.plot(chie, z*100, label="Equilibrium quality chie")
    plt.xlabel("chie")
    plt.ylabel("z (cm)")
    plt.legend()
    plt.show()


    plt.plot(P, z*100, label="Pressure")
    plt.xlabel("P (Pa)")
    plt.ylabel("z (cm)")
    plt.legend()
    plt.show()
    
    pass
elif option == 2:


    pass
else:
    raise ValueError(f"'{option}' is not a valid case option (1 for PWR, 2 for BWR)")