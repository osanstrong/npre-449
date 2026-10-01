from CoolProp.CoolProp import PropsSI

#Specific heat capacity of Water at 500 K and 1 atm
print(PropsSI('C','T',500,'P',101325,'IF97::Water'))

#Density of Water at 500 K and 1 atm.
print(PropsSI('D','T',500,'P',101325,'IF97::Water'))