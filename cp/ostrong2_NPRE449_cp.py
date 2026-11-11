# Input (optional)
import sys
args = sys.argv

# Solving
import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import newton

# Post-processing/plotting
from matplotlib import pyplot as plt


## Debugging
PRINT_DEBUG = True
def dbprint(str):
    if PRINT_DEBUG:
        print(str)

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
    # return -(0 + ρ*g)

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
    hL = calc_h_DB(Re, L, P, Tf)
    pr = P / P_CRIT_WATER

    hPool = 55 * (pr**0.12) * (qpp**(2/3)) * ((-np.log10(pr))**(-0.55)) * (M_WATER**(-0.5))

    zero_func = lambda Tw: np.sqrt((F*hL*(Tw-Tf))**2 + (S*hPool*(Tw-Tsat))**2) - qpp
    Tw = newton(zero_func, Tw_guess, maxiter=150)
    dbprint(f"pr: {pr}, ")
    dbprint(f"Tf: {Tf} P: {P} q'': {qpp} W/m²")
    dbprint(f"hL: {hL} hPool: {hPool}")
    dbprint(f"F: {F} S: {S}")
    dbprint(f" -> Tw: {Tw} K")
    return Tw



class Config:
    def prop00(self, name):
        return propPT(name, self.P0, self.Tf0)
    
    def __init__(self, rtype):
        self.rtype = rtype
        if rtype == PWR:
            # PWR params
            self.P0 = 15 * 1e6
            ## Givens for PWR
            self.H = 4.1 #m
            self.Drod = 0.94 / 100 # m, originally cm
            self.Pitch = 1.25 / 100 # m, originally cm
            self.Dfuel = 0.81 / 100 # ''
            self.tgap = 0.005 / 100 # ''
            self.kgap = 0.25 # W/mK
            self.kfuel = 3.6 # W/mK
            self.kclad = 21.5 # W/mK

            self.G = 3800 # kg / m²s
            self.qp0 = 380 * 100 # W / m, originally W/cm
            self.P0 = 15 * 1e6 # Pa, originally MPa
            self.Tf0 = K(277) # K, orignally ⁰C


            


            pass
        elif rtype == BWR:
            # BWR params
            pass
        else:
            raise ValueError(f"Type {rtype} does not correspond to either PWR ({PWR}) or BWR ({BWR})")
        self.recalc_secondaries()

    def recalc_secondaries(self):
        self.Dgap = self.Dfuel + 2*self.tgap # '', diameter to where gap meets cladding

        self.Af = self.Pitch**2 - (0.25 * PI * self.Drod**2)
        self.Pwet = PI*self.Drod
        self.Pheat = self.Pwet
        self.Lchar = self.Drod
        self.Lchar = 4 * self.Af / self.Pwet
        dbprint(f"Pwet and Pheat: {self.Pwet} m")
        dbprint(f"Fluid area Af: {self.Af} m²")
        dbprint(f"Characteristic length: {self.Lchar} m")

        ## Properties constant 
        self.μ = self.prop00('viscosity')
        self.cp = self.prop00('CPMASS')
        self.ρ = self.prop00('D')
        dbprint(f"Dynamic viscosity: {self.μ} Pa s")


        self.Re = calc_Re(self.Lchar, self.G, self.μ)
        dbprint(f"Reynolds number: {self.Re} >> 10,000 -> Turbulent flow")
        dbprint(f"Which leads to friction factor f = {calc_f_1ph(self.Re)}")
        pass
        

def solve_distributions(c, n):



    dz = c.H / n

    z = np.linspace(0, c.H, n)
    hfhg = lambda P: (
                PropsSI('H', 'P', P, 'Q', 0, WATER),
                PropsSI('H', 'P', P, 'Q', 1, WATER))

    ## P
    P = np.zeros(n)
    P[0] = c.P0
    f = calc_f_1ph(c.Re)    
    Tf = np.zeros(n)
    Tf[0] = c.Tf0
    χe = np.zeros(n)
    hf0, hg0 = hfhg(c.P0)
    hfg0 = hg0-hf0
    Tsat0 = PropsSI('T', 'P', c.P0, 'Q', 0.5, WATER)
    χe[0] = c.cp*(c.Tf0-Tsat0)/hfg0

    dbprint(f"Tf at inlet: {c.Tf0} K")
    dbprint(f"Tsat at inlet: {Tsat0} K")
    dbprint(f"cp at inlet: {c.cp} J/kgK")
    dbprint(f'h_fg at inlet: {hfg0} J/kg')
    dbprint(f"χe0: {χe[0]}")

    for i in range(n-1):
        dPdz = calc_dP(c.G, 0, f, c.Af, c.Pwet, c.ρ, g)
        dP = dPdz * dz
        P[i+1] = P[i] + dP

        # Properties
        
        
        hfi, hgi = hfhg(P[i])
        hfgi = hgi - hfi

        hfn, hgn = hfhg(P[i])
        hfgn = hgn - hfn

        dhfdP = (hfn-hfi) / dP
        dhfgdP = (hfgn-hfgi) / dP

        qp = q_lin(c.qp0, c.H, z[i])
        qpp = qp / c.Pheat
        dbprint(f"q': {qp} W/m -> q'': {qpp}")
        dχedz = calc_dχedz(hfgi, qp, c.Af, c.G, dhfdP, χe[i], dhfgdP, dPdz)
        χe[i+1] = χe[i] + dχedz*dz
        Tf[i+1] = calc_Tf(χe[i+1], P[i+1], c.cp)

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
        qp = q_lin(c.qp0, c.H, z[i])
        Tsat[i] = PropsSI('T', 'P', P[i], 'Q', 0, WATER)

        qpp = qp / c.Pheat
        ## Heat transfer, then 
        if ht_status == ONE_PHASE:
            h = calc_h_DB(c.Re, c.Lchar, P[i], Tf[i])
            # hf[i] = h
            Tw[i] = Tf[i] + qpp / h
            if (Tw[i] >= Tsat[i]):
                dbprint(f"h right before LW: {h} W/m²K")
        elif ht_status == TWO_PHASE:
            Tw[i] = calc_Tw_LW(c.Re, c.Lchar, P[i], Tf[i], 0, qpp, Tw[i-1]) # Calls previous Tw as a starting guess, but this is never called on z=0
        if Tw[i] >= Tsat[i]:
            ht_status = TWO_PHASE
        elif Tw[i] < Tsat[i]:
            ht_status = ONE_PHASE
        ## within clad
        # Tin = Tout + q'/2kπ * ln(Rout/Rin) for annulars
        Rc = c.Drod / 2
        Rf = c.Dfuel / 2
        Rg = Rf + c.tgap

        Tgap[i] = Tw[i] + (qp / (2*c.kclad*PI)) * np.log(Rc / Rg)
        ## Edge of fuel
        Tfuel[i] = Tgap[i] + (qp / (2*c.kgap*PI)) * np.log(Rg / Rf)



        #Within fuel
        qppp = qp / (PI * (Rf**2))
        Tcl[i] = Tfuel[i] + ((Rf**2)*qppp/(4*c.kfuel))


    ## Sanity check
    qtot_an = c.qp0*2*c.H/PI
    # qtot_nu = sum(q_lin(qp0, H, z))*dz
    # print(f"Qtot analytic: {qtot_an} W")
    # print(f"Qtot numeric:  {qtot_nu} W")
    m = c.G*c.Af
    hin = propPT('H', c.P0, c.Tf0)  
    hout = hin + qtot_an
    dbprint(f"h goes from\n{hin} J/kg to\n{hout} J/kg")
    # Tfout = PropsSI('T', 'P', P[-1], 'H', hout)
    Tfout = c.Tf0 + qtot_an/(m*c.cp)
    dbprint(f"Which translates to (analytic) final Tf of {Tfout-K(0)}⁰C")


    return z, Tf, Tw, Tgap, Tfuel, Tcl, Tsat, χe, P
    pass


## Variation

if option == 1:
    base_config = Config(PWR)
    n = 500
    dz = base_config.H / n

    ## Main solve
    #0, 1,  2,    3,     4,   5,    6,  7, 8
    z, Tf, Tw, Tgap, Tfuel, Tcl, Tsat, χe, P = solve_distributions(base_config, n)
    IDX_TCL = 5
    IDX_χE = 7

    ## Parameter sweeps
    P0_sweep = [P_MPa*1e6 for P_MPa in [14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18]]
    dT_sweep = [-50, -35, -25, -15, -5, 0, 5, 15, 25, 35, 50]
    dT_sweep = [-50, -35, -25, -15, -5, 0, 5, 7, 8, 8.7]
    Tf0_sweep = [base_config.Tf0 + dT for dT in dT_sweep]
    clad_sweep = [base_config.Drod * mod for mod in [0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05]]

    def plot_χe_Tcl(params, param_label, results):
        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel(param_label)
        ax1.set_ylabel('χe at outlet (-)', color=color)

        
        ax1.scatter(params, [res[IDX_χE][-1] for res in results], color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

        color = 'tab:blue'
        ax2.set_ylabel('Peak $T_{cl}$ (⁰C)', color=color)  # we already handled the x-label with ax1
        ax2.scatter(params, [max(res[IDX_TCL]) - K(0) for res in results], color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()
        return fig, ax1, ax2

    def conf_custom_P0(P0):
        conf = Config(PWR)
        conf.P0 = P0
        conf.recalc_secondaries()
        return conf

    # P0_results = [solve_distributions(conf_custom_P0(P0), n) for P0 in P0_sweep]
    # plot_χe_Tcl([P0/1e6 for P0 in P0_sweep], 'Inlet P (MPa)', P0_results)


    def conf_custom_Tf0(Tf0):
        conf = Config(PWR)
        conf.Tf0 = Tf0
        conf.recalc_secondaries()
        return conf

    Tf0_results = [solve_distributions(conf_custom_Tf0(Tf0), n) for Tf0 in Tf0_sweep]
    plot_χe_Tcl([Tf0-K(0) for Tf0 in Tf0_sweep], 'Inlet Tf (⁰C)', Tf0_results)

    def conf_custom_Drod(Drod):
        conf = Config(PWR)
        conf.Drod = Drod
        conf.recalc_secondaries()
        return conf
    # Drod_results = [solve_distributions(conf_custom_Drod(Drod), n) for Drod in clad_sweep]
    # plot_χe_Tcl([(Drod - base_config.Dgap) * 100 / 2 for Drod in clad_sweep], 'Clad thickness (cm)', Drod_results)



    # Tfd = np.zeros(n)
    # Tfd[0] = base_config.Tf0
    # ## T pass
    # for i in range(n-1):
    #     qp = q_lin(base_config.qp0, base_config.H, z[i])
    #     dT = qp*dz / (base_config.G*base_config.Af*base_config.cp)
    #     Tfd[i+1] = Tfd[i] + dT


    # plt.plot(Tf-K(0), z*100, label="Fluid T from χe")
    # plt.plot(Tfd-K(0), z*100, linestyle='--', label="Fuid T (alternate calculation)")
    # plt.plot(Tw-K(0), z*100, label="Outer wall T")
    # plt.plot(Tgap-K(0), z*100, label="Inside clad T")
    # plt.plot(Tfuel-K(0), z*100, label="Outside fuel T")
    # plt.plot(Tcl-K(0), z*100, label="Fuel Centerline T")

    # plt.plot(Tsat-K(0), z*100, linestyle='--', label="T$_{sat}$")
    # plt.xlabel("T (⁰C)")
    # plt.ylabel("z (cm)")
    # plt.legend()
    # plt.show()

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