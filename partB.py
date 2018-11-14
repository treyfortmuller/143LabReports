# Plots accompanying the theoretical caluclations in Lab Report 1

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import special # for the complimentary error function
# import models # custom script with theoretical models regarding processing physics



# # PHOSPHORUS AFTER PRE-DIFFUSION
# plt.figure(2) 
# x = np.linspace(0, 5e-7, 100)

# No = 1.25e21 # solid solubility limited phosphorous concentration
# D = 9.207e-14 # diffusion coefficient
# t = 0.17 # hours

# C = 1.0/(2*np.sqrt(D*t))

# plt.plot(x, No*special.erfc(C*x))

# plt.title("After Pre-diffusion Phosphorous Distribution")
# plt.xlabel("Wafer Depth")
# plt.ylabel("Dopant Concentration")



# # PHOSPHORUS AFTER DRIVE-IN
# plt.figure(3)
# x = np.linspace(0, 5e-7, 100)

# Q = 1.765e14 # total dose from the pre-diffusion step above
# D = 9.207e-14 # diffusion coefficient
# t = 0.62 # hours, both the wet oxidation and anneal time

# plt.title("After Drive-in Phosphorous Distribution")
# plt.xlabel("Wafer Depth")
# plt.ylabel("Dopant Concentration")
# plt.plot(x, (Q / np.sqrt(np.pi*D*t))*np.exp(-(x / (2*np.sqrt(D*t)))**2))

###########################################################

# AFTER PRE-DIFFUSION PLOT
plt.figure(4)
plt.title("Dopant Distributions After Pre-diffusion")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
plt.ylim(1e-40, 10e25)
plt.xlim(2e-5, 3.5e-5)
x = np.linspace(1e-5, 5e-5, 10000) # centimeters

# the phosphorus distribution
No = 1.25e21 # solid solubility limited phosphorous concentration
D = 9.207e-14 # diffusion coefficient
t = 0.17 # hours
C = 1.0/(2*np.sqrt(D*t))
phosphorous = No*special.erfc(C*(x-2.638e-5))
plt.loglog(x, phosphorous, label='phosphorous distribution')

# the gate oxide depth
oxide_depth = 2.638e-5
plt.axvline(x=oxide_depth, color='g', label='oxide interface')

# the boron distribution after pre-diffusion
Q = 6e12 # total dose from the pre-diffusion step above
Dt = 2.379e-13
Rp = 2e-5
boron = (Q / np.sqrt(np.pi*Dt))*np.exp(-((x - Rp) / (2*np.sqrt(Dt)))**2)
plt.loglog(x, boron, label='boron distribution')

# the intersection points of the dopants
dopant_intersect_x = np.argwhere(np.diff(np.sign(phosphorous - boron))).flatten()
plt.plot(x[dopant_intersect_x], phosphorous[dopant_intersect_x], 'ro')

# # print the junction depth
# print("intersection depth after pre-diffusion: " + str(x[dopant_intersect_x]))

# # surface concentration of phosphorus
# surface_conc = No*special.erfc(C*(oxide_depth-2.638e-5))
# print("surface concentration is: "+ str(surface_conc))

plt.legend()
plt.grid(True,which="major",ls="-")



# AFTER DRIVE-IN PLOT
plt.figure(5)
plt.title("Dopant Distributions After Drive-in")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
x = np.linspace(1e-5, 10e-5, 10000)

# the phosphorus distribution
Q = 1.765e14 # total dose from the pre-diffusion step above
D = 9.207e-14 # diffusion coefficient
t = 0.79 # hours, both the wet oxidation and anneal time
phosphorous = (Q / np.sqrt(np.pi*D*t))*np.exp(-((x-2.638e-5) / (2*np.sqrt(D*t)))**2)
plt.plot(x, phosphorous, label='phosphorous distribution')

# the gate oxide depth
oxide_depth = 3.471e-5
plt.axvline(x=oxide_depth, color='g', label='oxide interface')
plt.axvline(x=2.638e-5, color='k', linestyle='--', label='previous interface') # old oxide depth


# the boron distribution after pre-diffusion
Q = 6e12 # total dose from the pre-diffusion step above
Dt = 2.949e-13
Rp = 2e-5
boron = (Q / np.sqrt(np.pi*Dt))*np.exp(-((x - Rp) / (2*np.sqrt(Dt)))**2)
plt.loglog(x, boron, label='boron distribution')

# the intersection points of the dopants
dopant_intersect_x = np.argwhere(np.diff(np.sign(phosphorous - boron)))
plt.plot(x[dopant_intersect_x], phosphorous[dopant_intersect_x], 'ro')

# # print the junction depth
# print("intersection depth after drive-in: " + str(x[dopant_intersect_x]))

# # surface concentration of phosphorus
# surface_conc = (Q / np.sqrt(np.pi*D*t))*np.exp(-((3.471e-5-3.471e-5) / (2*np.sqrt(D*t)))**2)
# print("surface concentration is: "+ str(surface_conc))

plt.legend()
plt.grid(True,which="major",ls="-")


plt.show()


# add a vertical line somewhere
# plt.axvline(x=0.22058956)


# plot the intersection of two curves f and g.
# x = np.arange(0, 1000)
# f = np.arange(0, 1000)
# g = np.sin(np.arange(0, 10, 0.01) * 2) * 1000

# plt.plot(x, f, '-')
# plt.plot(x, g, '-')

# idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
# plt.plot(x[idx], f[idx], 'ro')
# plt.show()


# def diffusivity(temp):
# 	# for calculating the diffusivity of Boron as a function of temperature in Kelvin
# 	Do = 10.5 # diffusion coefficient [cm^2 / s]
# 	Ea = 3.69 # activation energy [eV]
# 	k = 8.617e-5 # boltzmann constant [eV / K]

# 	D = Do*np.exp(-Ea/(k*temp))
# 	return D

# # get the diffusivity*time for each of the five thermal step in the process
# steps = ["FOX", "GOX", "Poly", "Depo", "Drive"]
# temps = [1273, 1373, 923, 1323, 1323] # in Kelvin
# times = [1.42, 0.62, 4.0, 0.17, 0.62] # in hours

# Dt_products = []
# print("D*t products for each thermal step:")
# for i in range(len(temps)):
# 	D = diffusivity(temps[i])
# 	t = times[i]

# 	print(steps[i] + ": " + str(D*t))
# 	Dt_products.append(D*t)

# print("sum of D*t after pre-diffusion:")
# print(sum(Dt_products[:4])) # sum of the first four terms

# print("sum of D*t after drive-in")
# print(sum(Dt_products)) # sum of all of the terms


