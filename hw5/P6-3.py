from CoolProp.CoolProp import PropsSI
WATER = 'IF97::Water'


# Return quality corresponding to the given property 'key' with value 'val', assuming it lies between that of saturated liquid and vapor
# Defined using pressure, since specifying isobars seems more common
def get_qual(pressure, key, val) -> float:
    liq_val = PropsSI(key, 'P', pressure, 'Q', 0, WATER)
    vap_val = PropsSI(key, 'P', pressure, 'Q', 1, WATER)

    return (val - liq_val) / (vap_val - liq_val)

# Sanity Test quality function
P_hi = 6.89e6 # Pa
h89 = PropsSI('H', 'P', P_hi, 'Q', 0.89, WATER)
print(f"Sanity test: quality for 6.89MPa, h={h89} J/kg: {get_qual(P_hi, 'H', h89)} (should be 0.89)")

#### Identify states

# Find the h resulting from moving from P1, h1 to P2 without changing entropy
def get_h_isentropic(P1, h1, P2) -> float:
    s12 = PropsSI('S', 'P', P1, 'H', h1, WATER)
    h2s = PropsSI('H', 'P', P2, 'S', s12, WATER)
    return h2s

# Find h out of a turbine with the given isentropic efficiency
def get_h_out_turbine(P1, h1, P2, efficiency) -> float:
    h2s = get_h_isentropic(P1, h1, P2)
    h2 = h1 - efficiency*(h1 - h2s)
    return h2

# Find h out of a pump with given isentropic efficiency
def get_h_out_pump(P1, h1, P2, efficiency) -> float:
    h2s = get_h_isentropic(P1, h1, P2)
    h2 = h1 + (h2s - h1) / efficiency
    return h2


eff_pumps = 0.85
eff_turbs = 0.9

P_hi = 6.89e6
P_int = 1.38e6
P_low = 6.89e3

h1 = PropsSI('H', 'P', P_hi, 'Q', 1 , WATER)
h3 = PropsSI("H", 'P', P_int, 'Q', 1, WATER)
h4 = PropsSI("H", 'P', P_int, 'Q', 0, WATER)
h6 = PropsSI("H", 'P', P_low, 'Q', 0, WATER)
print(f"Units of h: kJ / kg")
print(f"h1: {h1/1000}")
print(f"h3: {h3/1000}")
print(f"h4: {h4/1000}")
print(f"h6: {h6/1000}")
# This is where I learned that apparently the corresponding enthalpy can decrease w/ pressure for sat. steam at too high pressures
# But I guess actually that makes sense, because eventually the latent energy between liquid and steam decreases (faster than pressure adds work) until you get to the triple point
# The more you know!

## Solve for the first degree unknowns

h2s = get_h_isentropic(P_hi, h1, P_int)
h2r = get_h_out_turbine(P_hi, h1, P_int, eff_turbs)
print(f"h2s: {h2s/1000}\nh2r: {h2r/1000}")
h5s = get_h_isentropic(P_int, h3, P_low)
h5r = get_h_out_turbine(P_int, h3, P_low, eff_turbs)
print(f"h5s: {h5s/1000}\nh5r: {h5r/1000}")

h7s = get_h_isentropic(P_low, h6, P_int)
h7r = get_h_out_pump(P_low, h6, P_int, eff_pumps)
print(f"h7s: {h7s/1000}\nh7r: {h7r/1000}")

## Second/third/fourth degree unknowns: x2, h8, h9 (x5 would also be one but we don't actually really need that since we're not separating/searching for other properties)
x2s = get_qual(P_int, 'H', h2s)
x2r = get_qual(P_int, 'H', h2r)
print(f"x2s: {x2s}\nx2r: {x2r}")

h8s = x2s*h7s + (1-x2s)*h4
h8r = x2r*h7r + (1-x2r)*h4
print(f"h8s: {h8s/1000}\nh8r: {h8r/1000}")

h9s = get_h_isentropic(P_int, h8s, P_hi)
h9r = get_h_out_pump(P_int, h8r, P_hi, eff_pumps)
print(f"h9s: {h9s/1000}\nh9r: {h9r/1000}")

## ALL WORK/HEAT RATES NORMALIZED OVER TOTAL MASS FLOW RATE
WT1s = h1-h2s
WT1r = h1-h2r

WT2s = x2s * (h3-h5s)
WT2r = x2r * (h3-h5r)

WP1s = x2s * (h7s - h6)
WP1r = x2s * (h7r - h6)

WP2s = h9s - h8s
WP2r = h9r - h8r

Qins = h1 - h9s
Qinr = h1 - h9r

eff_s = (sum([WT1s, WT2s, -WP1s, -WP2s]) / Qins)
eff_r = (sum([WT1r, WT2r, -WP1r, -WP2r]) / Qinr)

print(f"Efficiency (s): ({WT1s} + {WT2s} - {WP1s} - {WP1s}) / {Qins}")
print(f"Efficiency (r): ({WT1r} + {WT2r} - {WP1r} - {WP1r}) / {Qinr}")
print("Across the board r gets less turbine work out, takes more work to pump, and more Q in")

print('_________')
print(f"eff_s: {eff_s}\neff_r: {eff_r}")



    