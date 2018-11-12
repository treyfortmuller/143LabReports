# Plots accompanying the theoretical caluclations in Lab Report 1

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import special # for the complimentary error function
import models # custom script with theoretical models regarding processing physics


# initial ion implantation concentration (step 1)
x = np.linspace(0, 4e-5, 100)

# add units here
Np = 4.35e17
delta_Rp = 5.5e-6
Rp = 2e-5

# Np * np.exp(-(x-delta_Rp)**2/(2*delta_Rp**2))

plt.figure(1)
plt.plot(x, Np * np.exp((-(x-delta_Rp)**2)/(2*delta_Rp**2)))

plt.title("Initial Ion Implantation Distribution")
plt.xlabel("Wafer Depth")
plt.ylabel("Ion Implantation Concentration")

# constant source, solid solubility limited phosphorus diffusion
plt.figure(2) 
x = np.linspace(0, 5e-7, 100)

No = 1.25e21 # solid solubility limited phosphorous concentration
D = 9.207e-14 # diffusion coefficient
t = 0.17 # hours

C = 1.0/(2*np.sqrt(D*t))

plt.plot(x, No*special.erfc(C*x))

plt.title("Constant Source Pre-diffusion Distribution")
plt.xlabel("Wafer Depth")
plt.ylabel("Diffused Dopant Concentration")

# limited source drive-in distribution
plt.figure(3)
x = np.linspace(0, 5e-7, 100)

Q = 1.765e14 # total dose from the pre-diffusion step above
D = 9.207e-14 # diffusion coefficient
t = 0.62 # hours, both the wet oxidation and anneal time

plt.title("Limited Source Drive-in Distribution")
plt.xlabel("Wafer Depth")
plt.ylabel("Diffused Dopant Concentration")

plt.plot(x, (Q / np.sqrt(np.pi*D*t))*np.exp(-(x / (2*np.sqrt(D*t)))**2))

plt.show()









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


