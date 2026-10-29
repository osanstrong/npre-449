from CoolProp.CoolProp import PropsSI
# import numpy as np
# from matplotlib import pyplot as plt


# WATER = 'IF97::Water'
AIR = 'AIR.PPF'
Pa_1ATM = 101_325

PI = 3.141592653589

def K(C):
    return C+273.15
m = 0.05 #kg/s
cp = 1011 # J / kg K
Tout = K(77) 
Tin = K(103)
T0 = K(0)
h0 = 6 # W / m²K
D = 0.15
L = 5
Pr = 0.7 # @ 350k
k = PropsSI('conductivity','P', Pa_1ATM, 'T', Tout, AIR)

print(f"Q into pipe: {m*cp*(Tout-Tin)}")

# mu = PropsSI('V', 'P', Pa_1ATM, 'T', Tout, AIR) # Realized for this one getting from the table was easier
# print(f"Viscosity: {mu} Pa s")
mu = 2.08e-5 # N s / m², at 350K
print(f"Viscosity: {mu} Pa s")
# print(f"Viscosity: {208.2e-7} Pa s")

Re = 4*m / (PI*D*mu)
print(f"Re: {Re}")

Nu = 0.023*(Re**0.8)*(Pr**0.3)
print(f"Nu: {Nu}")

hL = Nu*k/D
print(f"hL: {hL} W/m²K")

Ts = (hL*Tout + h0*T0) / (hL+h0)
print(f"Ts: {Ts-K(0)}⁰C")

qii = h0*(Ts-T0)
qii2 = hL*(Tout-Ts)
print(f"q''(L): {qii} or {qii2} W/m²")
