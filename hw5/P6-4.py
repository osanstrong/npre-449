from CoolProp.CoolProp import PropsSI
WATER = 'IF97::Water'

p1 = 10e6 # 10MPa
p2 = 0.5e6
h1 = 3000e3 # 3000 kJ/kg
h2 = 2600e3
v1 = 150 # m/s
v2 = 0
rho1 = PropsSI('D', 'P', p1, 'H', h1, WATER)
rho2 = PropsSI('D', 'P', p2, 'H', h2, WATER)
print(f"ρ1: {rho1} kg/m³")
print(f"ρ2: {rho2} kg/m³")

print(f'e from dh: {(h1-h2)/1000} kJ/kg')
e_dv = 0.5 * (v1**2 - v2**2)
print(f"e from dv: {e_dv/1000} kJ/kg")

w12 = 384.45e3

q12 = w12 + h2 - h1 - 0.5*v1**2
print(f"q12: {q12/1000} kJ/kg")

print("______")
h2l = PropsSI("H", 'P', p2, 'Q', 0, WATER)
h2v = PropsSI('H', 'P', p2, 'Q', 1, WATER)
print(f"h2l: {h2l/1000} kJ/kg")
print(f"h2v: {h2v/1000} kJ/kg")
x2 = (h2-h2l) / (h2v-h2l)
print(f"x2: {x2}")

