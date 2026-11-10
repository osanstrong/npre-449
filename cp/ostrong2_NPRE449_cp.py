# Input (optional)
import argparse

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
g = 9.81 # m / s²

#### ---- User Configuration ---- ####
option = 1 # 1 for PWR, 2 for BWR
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

def calc_Tf(χe, P, cp):
    hf = PropsSI('H', 'P', P, 'Q', 0, WATER)
    hg = PropsSI('H', 'P', P, 'Q', 1, WATER)
    Tsat = PropsSI('T', 'P', P, 'Q', 0, WATER)
    hfg = hg - hf
    Tf = χe*hfg/cp + Tsat
    return Tf


## Variation

if option == 1:
    ## Givens for PWR
    H = 4.1 #m
    Drod = 0.94 / 100 # m, originally cm
    Pitch = 1.25 / 100 # m, originally cm
    Dfuel = 0.81 / 100 # ''
    tgap = 0.005 / 100 # ''
    kgap = 0.25 # W/mK
    kfuel = 2.6 # W/mK
    kclad = 21.5 # W/mK

    G = 3800 # kg / m²s
    qp0 = 380 * 100 # W / m, originally W/cm
    P0 = 15 * 1e6 # Pa, originally MPa
    Tf0 = K(277) # K, orignally ⁰C

    Af = Pitch**2 - (0.25 * PI * Drod**2)
    Pwet = PI*Drod

    ## Properties constant 
    prop00 = lambda name: propPT(name, P0, Tf0)
    μ = prop00('viscosity')
    cp = prop00('CPMASS')
    ρ = prop00('D')


    Re = calc_Re(Drod, G, μ)
    print(f"Reynolds number: {Re} >> 10,000 -> Turbulent flow")

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
    χe = np.zeros(n)
    hf0, hg0 = hfhg(P0)
    hfg0 = hg0-hf0
    Tsat0 = PropsSI('T', 'P', P0, 'Q', 0.5, WATER)
    χe[0] = cp*(Tf0-Tsat0)/hfg0


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




    Tfd = np.zeros(n)
    Tfd[0] = Tf0
    ## T pass
    for i in range(n-1):
        qp = q_lin(qp0, H, z[i])
        dT = qp*dz / (G*Af*cp)
        Tfd[i+1] = Tfd[i] + dT


    plt.plot(Tf-K(0), z*100, label="Fuel centerline T")
    plt.plot(Tfd-K(0), z*100, linestyle='--', label="T from χe")
    plt.xlabel("T (⁰C)")
    plt.ylabel("z (cm)")
    plt.legend()
    plt.show()
    
    plt.plot(χe, z*100, label="Equilibrium quality χe")
    plt.xlabel("χe")
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