# Plots accompanying the theoretical caluclations in Lab Report 1

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import special # for the complimentary error function

def diffusivity(temp):
	# for calculating the diffusivity of Boron as a function of temperature in Kelvin
	Do = 10.5 # diffusion coefficient [cm^2 / s]
	Ea = 3.69 # activation energy [eV]
	k = 8.617e-5 # boltzmann constant [eV / K]

	D = Do*np.exp(-Ea/(k*temp))
	return D

def limited_source_diffusion(x_peak, Q, Dt, x):
	# the gaussian limited source diffusion model
	# takes the location of the peak, the dose, diffusivity times time, and x linspace
	# returns the dopant distrubtion
	return (Q/np.sqrt(np.pi*Dt))*np.exp(-((x - x_peak)/(2*np.sqrt(Dt)))**2)

def oxidation(A, B, tau, time):
	# the deal-grove oxidation model for a <100> Si wafer
	# takes the linear and parabolic diffusivity constants, the inital oxidation time, and time in hours
	# returns oxide thickness in microns
	return (A/2.0)*(np.sqrt(1+((4*B)/(A**2)*(time+tau))) - 1)

print("field oxide thickness: " + str(oxidation(0.4245, 0.3151, 0, 1.42)))
print("gate oxide thickness: " + str(oxidation(0.1396, 0.0236, 0.1744, 0.62)))
print("intermediate oxide thickness: " + str(oxidation(0.274, 0.412, 0, 0.2)))

# get the diffusivity*time for each of the five thermal step in the process
steps = ["FOX", "GOX", "Poly", "Depo", "Drive", "Sinter"]
temps = [1273, 1373, 923, 1323, 1323, 673] # in Kelvin
times = [1.42, 0.62, 4.0, 0.17, 0.62, 0.33] # in hours

Dt_products = []
print("D*t products for each thermal step:")
for i in range(len(temps)):
	D = diffusivity(temps[i])
	t = times[i]

	print(steps[i] + ": " + str(D*t))
	Dt_products.append(D*t)

print("sum of D*t after drive-in")
print(sum(Dt_products)) # sum of all of the terms

