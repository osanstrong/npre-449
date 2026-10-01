from CoolProp.CoolProp import PropsSI
WATER = 'IF97::Water'

p1 = 10e6 # 10MPa
p2 = 0.5e6
h1 = 3000e3 # 3000 kJ/kg
h2 = 2600e3
v1 = 150 # m/s
v2 = 0
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

'''
e from dh: 400.0 kJ/kg
e from dv: 11.25 kJ/kg
q12: -26.8 kJ/kg
______
h2l: 640.1853353633849 kJ/kg
h2v: 2748.1076146579694 kJ/kg
x2: 0.9297376302187321
'''

