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

def ion_implantation(Dt_term, x):
	Q = 6e12 # total dose from the pre-diffusion step above
	Dt = Dt_term
	Rp = 2e-5
	return (Q / np.sqrt(np.pi*Dt))*np.exp(-((x - Rp) / (2*np.sqrt(Dt)))**2)

# print("field oxide thickness: " + str(oxidation(0.4245, 0.3151, 0, 1.42)))
# print("gate oxide thickness: " + str(oxidation(0.1396, 0.0236, 0.1744, 0.62)))
# print("intermediate oxide thickness: " + str(oxidation(0.274, 0.412, 0, 0.2)))

# get the diffusivity*time for each of the five thermal step in the process
steps = ["FOX", "GOX", "Poly", "Depo", "Drive", "Sinter"]
temps = [1273, 1373, 923, 1323, 1323, 673] # in Kelvin
times = [1.42, 0.62, 4.0, 0.17, 0.62, 0.33] # in hours

oxides = [4.895e-5, 8.39e-6, 1.811e-5] # FOX, GOX, IOX growth in cm
oxides_sum = []
for i in range(len(oxides)):
	oxides_sum.append(sum(oxides[:i+1]))

Dt_products = []
# print("D*t products for each thermal step:")
for i in range(len(temps)):
	D = diffusivity(temps[i])
	t = times[i]

	# print(steps[i] + ": " + str(D*t))
	Dt_products.append(D*t)

# print(Dt_products) # sum of all of the terms

# get the successive sums of the DT terms for each thermal step
Dt_subsequent_sums = []
for i in range(len(Dt_products)):
	Dt_subsequent_sums.append(sum(Dt_products[:i+1]))

####################################

# PLOT EVERYTHING SUPERIMPOSED
plt.figure(1)

plt.title("All Steps Dopant Distributions")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")

plt.xlim(1e-6, 3.5e-4)
# plt.ylim(1e-10, 1e20)

colors = ['r', 'g', 'b']

x = np.linspace(0, 10e-5, 10000)

for i in range(len(oxides_sum)):
	plt.axvline(x=0.46*oxides_sum[i], color=colors[i], linestyle='--', label='oxide interface'+str(i))

for i in range(len(Dt_subsequent_sums)):
	plt.loglog(x, ion_implantation(Dt_subsequent_sums[i], x), label='boron distribution'+str(i))

# the phosphorus distribution after pre-diffusion
No = 1.25e21 # solid solubility limited phosphorous concentration
D = 9.207e-14 # diffusion coefficient
t = 0.17 # hours
C = 1.0/(2*np.sqrt(D*t))
phosphorous = No*special.erfc(C*(x-2.638e-5))
plt.loglog(x, phosphorous, label='phosphorous distribution 1')

# the phosphorus distribution after drive-in
Q = 1.765e14 # total dose from the pre-diffusion step above
D = 9.207e-14 # diffusion coefficient
t = 0.79 # hours, both the wet oxidation and anneal time
phosphorous = (Q / np.sqrt(np.pi*D*t))*np.exp(-((x-3.471e-5) / (2*np.sqrt(D*t)))**2)
plt.plot(x, phosphorous, label='phosphorous distribution 2')

# the phosphorus distribution after sintering
Q = 1.765e14 # total dose from the pre-diffusion step above
D = 9.207e-14 # diffusion coefficient
t = 1.12 # hours, both the wet oxidation and anneal time and the sintering
phosphorous = (Q / np.sqrt(np.pi*D*t))*np.exp(-((x-3.471e-5) / (2*np.sqrt(D*t)))**2)
plt.plot(x, phosphorous, label='phosphorous distribution 3')

####################################

# PLOT INDIVIDUAL STEPS
plt.figure(2)
plt.grid(True,which="both",ls="-")

plt.title("Dopant Distribution After Field Oxidation")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
plt.xlim(1e-6, 3.5e-4)
plt.axvline(x=0.46*oxides_sum[0], color='r', linestyle='--', label='FOX interface')
plt.loglog(x, ion_implantation(Dt_subsequent_sums[0], x), label='boron distribution')
plt.legend()


plt.figure(3)
plt.grid(True,which="both",ls="-")

plt.title("Dopant Distribution After Gate Oxidation")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
plt.xlim(1e-6, 3.5e-4)
plt.axvline(x=0.46*oxides_sum[0], color='k', linestyle='--', label='FOX interface')
plt.axvline(x=0.46*oxides_sum[1], color='r', linestyle='--', label='GOX interface')
plt.loglog(x, ion_implantation(Dt_subsequent_sums[1], x), label='boron distribution')
plt.legend()


plt.figure(4)
plt.grid(True,which="both",ls="-")

plt.title("Dopant Distribution After Poly-Deposition")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
plt.xlim(1e-6, 3.5e-4)
plt.axvline(x=0.46*oxides_sum[0], color='k', linestyle='--', label='FOX interface')
plt.axvline(x=0.46*oxides_sum[1], color='r', linestyle='--', label='GOX interface')
plt.loglog(x, ion_implantation(Dt_subsequent_sums[2], x), label='boron distribution')
plt.legend()


plt.figure(5)
plt.grid(True,which="both",ls="-")

plt.title("Dopant Distribution After Pre-Diffusion")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
plt.xlim(1e-6, 3.5e-4)
plt.axvline(x=0.46*oxides_sum[0], color='k', linestyle='--', label='FOX interface')
plt.axvline(x=0.46*oxides_sum[1], color='r', linestyle='--', label='GOX interface')
plt.loglog(x, ion_implantation(Dt_subsequent_sums[3], x), label='boron distribution')
# the phosphorus distribution after pre-diffusion
No = 1.25e21 # solid solubility limited phosphorous concentration
D = 9.207e-14 # diffusion coefficient
t = 0.17 # hours
C = 1.0/(2*np.sqrt(D*t))
phosphorous = No*special.erfc(C*(x-2.638e-5))
plt.loglog(x, phosphorous, label='phosphorous distribution 1')
plt.legend()


plt.figure(6)
plt.grid(True,which="both",ls="-")

plt.title("Dopant Distribution After Drive-in")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
plt.xlim(1e-6, 3.5e-4)
plt.axvline(x=0.46*oxides_sum[0], color='0.85', linestyle='--', label='FOX interface')
plt.axvline(x=0.46*oxides_sum[1], color='k', linestyle='--', label='GOX interface')
plt.axvline(x=0.46*oxides_sum[2], color='r', linestyle='--', label='IOX interface')
plt.loglog(x, ion_implantation(Dt_subsequent_sums[4], x), label='boron distribution')
# the phosphorus distribution after drive-in
Q = 1.765e14 # total dose from the pre-diffusion step above
D = 9.207e-14 # diffusion coefficient
t = 0.79 # hours, both the wet oxidation and anneal time
phosphorous = (Q / np.sqrt(np.pi*D*t))*np.exp(-((x-3.471e-5) / (2*np.sqrt(D*t)))**2)
plt.plot(x, phosphorous, label='phosphorous distribution 2')
plt.legend()


plt.figure(7)
plt.grid(True,which="both",ls="-")

plt.title("Dopant Distribution After Sintering")
plt.xlabel("Wafer Depth ($cm$)")
plt.ylabel("Dopant Concentration (atoms/$cm^2$)")
plt.xlim(1e-6, 3.5e-4)
plt.axvline(x=0.46*oxides_sum[0], color='0.85', linestyle='--', label='FOX interface')
plt.axvline(x=0.46*oxides_sum[1], color='k', linestyle='--', label='GOX interface')
plt.axvline(x=0.46*oxides_sum[2], color='r', linestyle='--', label='IOX interface')
plt.loglog(x, ion_implantation(Dt_subsequent_sums[5], x), label='boron distribution')
# the phosphorus distribution after sintering
Q = 1.765e14 # total dose from the pre-diffusion step above
D = 9.207e-14 # diffusion coefficient
t = 1.12 # hours, both the wet oxidation and anneal time and the sintering
phosphorous = (Q / np.sqrt(np.pi*D*t))*np.exp(-((x-3.471e-5) / (2*np.sqrt(D*t)))**2)
plt.plot(x, phosphorous, label='phosphorous distribution 3')
plt.legend()

plt.show()

