
from CoolProp.CoolProp import PropsSI
import numpy as np
from matplotlib import pyplot as plt


WATER = 'IF97::Water'

PI = 3.141592653589
K0C = 273.15
Pa_1ATM = 101_325

def K(C):
    return C + K0C
Tin = K(25)
Tout = K(75)
Ts = K(100)

m_base = 2 #kg / s
mu = PropsSI('VISCOSITY', 'T', K(50), 'P', Pa_1ATM, WATER) # Pa s -> kg/m/s
# print(f"Viscosity: {mu} kg/mg/s")
rho = PropsSI('D', 'T', K(50), 'P', Pa_1ATM, WATER) # kg / m³
# print(f"Density: {rho} kg / m³")
Pr = PropsSI('PRANDTL', 'T', K(50), 'P', Pa_1ATM, WATER)
cp = PropsSI('CPMASS', 'T', K(50), 'P', Pa_1ATM, WATER) # J / kg K

def calc_Re(D, m):
    # D = 0.04 #m
    A = PI*(0.5*D)**2

    Re = m*D / (A*mu)
    return Re
    # print(f'Re: {Re}')


def calc_NuD(D, m):
    # Eq. 8.60, Dittus-Boelter
    NuD = 0.023*(calc_Re(D, m)**0.8)*(Pr**0.4)
    return NuD

print(f"NuD for original D: {calc_NuD(0.04, m_base)}")

m_range = np.array([1, 1.5, 2, 2.5, 3]) # kg / s
D_range = np.linspace(0.03, 0.05) # m

k = PropsSI('CONDUCTIVITY', 'T', K(50), 'P', Pa_1ATM, WATER)

def calc_h(D, m):
    NuD = calc_NuD(D, m)
    return NuD * k / D

print(f'default h: {calc_h(0.04, m_base)} W/m²K')


def get_L(D, m):
    P = D*PI
    L = (m*cp/(P*calc_h(D, m))) * np.log((Ts-Tin)/(Ts-Tout))
    return L

def calc_f(D, m):
    Re = calc_Re(D, m)
    f = (0.79*np.log(Re) - 1.64)**-2
    return f

def calc_dp(D, m):
    f = calc_f(D, m)
    A = PI*(D**2)/4
    u = m / (rho*A)
    L = get_L(D, m)
    dp = (f*u*u / (2*D)) * L
    return dp

print(f"Pipe length base params: {get_L(0.04, m_base)} m")
print(f"Friction factor base params: {calc_f(0.04, m_base)}")
print(f"Pressure drop base params: {calc_dp(0.04, m_base)}")



for m in m_range:
    Re_range = calc_Re(D_range, m)
    idx = Re_range.argsort()
    f_range = calc_f(D_range, m)
    print(f"m: {m} kg/s -> Re = ({Re_range[idx[0]]} - {Re_range[idx[-1]]}) > 10,000 : turbulent!")
    L_range = get_L(D_range, m)
    plt.plot(D_range, L_range, label='$\cdot{m}$'+f' = {m} kg/s')
    # plt.plot(D_range, Re_range, label='$\cdot{m}$'+f' = {m} kg/s')

plt.legend()
plt.xlabel("Diameter (m)")
plt.ylabel("Required length (m)")
plt.show()



for m in m_range:
    Re_range = calc_Re(D_range, m)
    idx = Re_range.argsort()
    f_range = calc_f(D_range, m)
    print(f"m: {m} kg/s -> Re = ({Re_range[idx[0]]}, {Re_range[idx[-1]]}) within (3000, 5e+6) : friction dp valid!")
    # L_range = get_L(D_range, m)
    dp_range = calc_dp(D_range, m)
    plt.plot(D_range, dp_range, label='$\cdot{m}$'+f' = {m} kg/s')
    # plt.plot(D_range, Re_range, label='$\cdot{m}$'+f' = {m} kg/s')

plt.legend()
plt.xlabel("Diameter (m)")
plt.ylabel("Pressure drop (Pa)")
plt.show()