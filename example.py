import OpenFoamFields as OF


elbowFields=OF.OpenFoamFields()
# Introduce the path to your OpenFoam case and the number of dimensions
elbowFields.initialize('elbow-test',2)
# Now, here are some examples of what you can do. However, you may read through the comments
# of OpenFoamFields.py for more information.

# 1. Find the value of the velocity and pressure fields at a given point (x,y,z)
#    and time t.
# If t is beyond the timespan of the simulation steady state is assumed and the
# last saved field is extrapolated.
t=3.26 
# If the point is outside the domain nan is returned
xyz=[3,4,0]
txyz2uvwp=elbowFields.interpolate(t,xyz[0],xyz[1],xyz[2])
print(("At time %s in point (%s,%s,%s):") %(t,xyz[0],xyz[1],xyz[2]))
print(("Ux=%s") %(txyz2uvwp[0]))
print(("Uy=%s") %(txyz2uvwp[1]))
print(("Uz=%s") %(txyz2uvwp[2]))
print(("P=%s") %(txyz2uvwp[3]))
print('\n')

# 2. Find the maximum and minimum value of a field at any given time t:
t=2.12345
print(("Minimum and maximum values at time: %s") %t)
# x-Velocity
Ux_MinAndMax=elbowFields.minmax('Ux',t)
print(("min(Ux)=%s   max(Ux)=%s") %(Ux_MinAndMax[0],Ux_MinAndMax[1]))
# y-Velocity
Uy_MinAndMax=elbowFields.minmax('Uy',t)
print(("min(Uy)=%s   max(Uy)=%s") %(Uy_MinAndMax[0],Uy_MinAndMax[1]))
# z-Velocity
Uz_MinAndMax=elbowFields.minmax('Uz',t)
print(("min(Uz)=%s   max(Uz)=%s") %(Uz_MinAndMax[0],Uz_MinAndMax[1]))
# Velocity magnitude
Umag_MinAndMax=elbowFields.minmax('Umag',t)
print(("min(Umag)=%s   max(Umag)=%s") %(Umag_MinAndMax[0],Umag_MinAndMax[1]))
# Pressure
P_MinAndMax=elbowFields.minmax('P',t)
print(("min(P)=%s   max(P)=%s") %(P_MinAndMax[0],P_MinAndMax[1]))
print('\n')

# 3. Find the minimum and maximum values of the pressure and velocity fields 
#    throughout the whole timespan of the simulation.
print("All-time Minimum and maximum values:")
# x-Velocity
Ux_MinAndMax=elbowFields.allTimeMinmax('Ux')
print(("min(Ux)=%s   max(Ux)=%s") %(Ux_MinAndMax[0],Ux_MinAndMax[1]))
# y-Velocity
Uy_MinAndMax=elbowFields.allTimeMinmax('Uy')
print(("min(Uy)=%s   max(Uy)=%s") %(Uy_MinAndMax[0],Uy_MinAndMax[1]))
# z-Velocity
Uz_MinAndMax=elbowFields.allTimeMinmax('Uz')
print(("min(Uz)=%s   max(Uz)=%s") %(Uz_MinAndMax[0],Uz_MinAndMax[1]))
# Velocity magnitude
Umag_MinAndMax=elbowFields.allTimeMinmax('Umag')
print(("min(Umag)=%s   max(Umag)=%s") %(Umag_MinAndMax[0],Umag_MinAndMax[1]))
# Pressure
P_MinAndMax=elbowFields.allTimeMinmax('P')
print(("min(P)=%s   max(P)=%s") %(P_MinAndMax[0],P_MinAndMax[1]))

