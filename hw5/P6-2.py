from CoolProp.CoolProp import PropsSI
WATER = 'IF97::Water'



# Return quality corresponding to the given property 'key' with value 'val', assuming it lies between that of saturated liquid and vapor
# Defined using temperature, as this problem gives isotherms
def get_qual(temperature, key, val) -> float:
    liq_val = PropsSI(key, 'T', temperature, 'Q', 0, WATER)
    vap_val = PropsSI(key, 'T', temperature, 'Q', 1, WATER)

    return (val - liq_val) / (vap_val - liq_val)

# Converted to K
T_hi = 293 + 273.15
T_lo =  33 + 273.15

### Cycle A
print("_______\nCycle A")

h1 = PropsSI('H', "T", T_hi, "Q", 1, WATER)
h3 = PropsSI('H', "T", T_lo, "Q", 0, WATER)
h4 = PropsSI('H', "T", T_hi, "Q", 0, WATER)
print(f"h1: {h1/1000} kJ/kg")
print(f"h3: {h3/1000} kJ/kg")
print(f"h4: {h4/1000} kJ/kg")

# Interpolate for x2 such that s2=s1
s1 = PropsSI('S', "T", T_hi, 'Q', 1, WATER)
s33v = PropsSI('S', 'T', T_lo, 'Q', 1, WATER)
s33l = PropsSI('S', 'T', T_lo, 'Q', 0, WATER)
x2 = (s1 - s33l) / (s33v - s33l)
h2 = PropsSI('H', 'T', T_lo, 'Q', x2, WATER) # i.e. x2 * h33v + (1-x2) * h33l
print(f"x2: {x2} -> h2: {h2/1000} kJ/kg")

# Interpolate for h3p
s3 = s33l
p4 = PropsSI('P', 'T', T_hi, 'Q', 0, WATER)
h3p = PropsSI("H", 'P', p4, 'S', s3, WATER)
print(f"h3p: {h3p/1000} kJ/kg")

# Thermal efficiency
eff = ((h1-h2) - (h3p-h3)) / (h1-h3p)
print(f"Efficiency: {eff}")

# Steam rate (kg / We-hr):
sr = 3600 / (h1-h2 - (h3p-h3))
print(f"Steam rate: {sr*1000} kg/kWe-hr")
print(f"Added 3'->4: {(h4-h3p)/1000} kJ/kg")
print(f"Added 4 ->1: {(h1-h4) /1000} kJ/kg")

'''
Cycle A
h1: 2761.949736266131 kJ/kg
h3: 138.28554911597698 kJ/kg
h4: 1306.077914527653 kJ/kg
x2: 0.667517832005123 -> h2: 1755.483410547582 kJ/kg
h3p: 146.069372134698 kJ/kg
Efficiency: 0.38177682603287816
Steam rate: 3.5768707884288555 kg/kWe-hr
'''

### Cycle B
print("_______\nCycle B")

s4 = PropsSI('S', 'T', T_hi, 'Q', 0, WATER)
x3 = (s4 - s33l) / (s33v - s33l)
h3 = PropsSI("H", "T", T_lo, 'Q', x3, WATER)
print(f"x3 = {x3} -> h3: {h3 / 1000} kJ/kg")

eff = (h1-h2 - (h4-h3)) / (h1-h4)
print(f"Efficiency: {eff}")

# Steam rate (kg / We-hr):
sr = 3600 / (h1-h2 - (h4-h3))
print(f"Steam rate: {sr*1000} kg/kWe-hr")
print(f"Added 4 ->1: {(h1-h4) /1000} kJ/kg")

'''
Cycle B
x3 = 0.34256175421728796 -> h3: 968.2111477540751 kJ/kg
Efficiency: 0.4592434230553254
Steam rate: 5.384388834597327 kg/kWe-hr
'''

### Cycle C
p = 5e6 # 5MPa
print("_______\nCycle C")

h1 = PropsSI("H", "T", T_hi, 'P', p, WATER)
h3 = PropsSI('H', "T", T_lo, "Q", 0, WATER)
h4 = PropsSI('H', "T", T_hi, "Q", 0, WATER)
h5 = PropsSI('H', "T", T_hi, "Q", 1, WATER)
print(f"h1: {h1/1000} kJ/kg")
print(f"h3: {h3/1000} kJ/kg")
print(f"h4: {h4/1000} kJ/kg")
print(f"h5: {h5/1000} kJ/kg")

s1 = PropsSI('S', 'T', T_hi, 'P', p, WATER)
x2 = (s1 - s33l) / (s33v - s33l)
h2 = PropsSI('H', 'T', T_lo, 'Q', x2, WATER) # i.e. x2 * h33v + (1-x2) * h33l
print(f"x2: {x2} -> h2: {h2/1000} kJ/kg")

s3 = PropsSI('S', 'T', T_lo, 'Q', 0, WATER)
h3p = PropsSI("H", "P", p, 'S', s3, WATER)
print(f"h3p: {h3p/1000} kJ/kg")

# Thermal efficiency
eff = ((h1-h2) - (h3p-h3)) / (h1-h3p)
print(f"Efficiency: {eff}")

# Steam rate (kg / We-hr):
sr = 3600 / (h1-h2 - (h3p-h3))
print(f"Steam rate: {sr*1000} kg/kWe-hr")
print(f"Added 3'->4: {(h4-h3p)/1000} kJ/kg")
print(f"Added 4 ->1: {(h1-h4) /1000} kJ/kg")

'''
Cycle C
h1: 2903.0028084402525 kJ/kg
h3: 138.28554911597698 kJ/kg
h4: 1306.077914527653 kJ/kg
h5: 2761.949736266131 kJ/kg
x2: 0.7194361740011986 -> h2: 1881.266165133955 kJ/kg
h3p: 143.2899698793983 kJ/kg
Efficiency: 0.3684195718975912
Steam rate: 3.5407552944435037 kg/kWe-hr
Added 3'->4: 1162.7879446482546 kJ/kg
Added 4 ->1: 1596.9248939125996 kJ/kg
'''

### Misc
print("________\nMisc")
eff_max = 1 - (T_lo / T_hi)
print(f"Carnot efficiency: {eff_max}")

'''
Misc
Carnot efficiency: 0.4592422502870265
'''


