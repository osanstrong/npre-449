from CoolProp.CoolProp import PropsSI

WATER = 'IF97::Water'

print(PropsSI('T', 'Q', 1, 'P', 17_200_000, WATER))
print(PropsSI('T', 'Q', 1, 'P', 101_325, WATER))
