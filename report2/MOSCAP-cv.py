# Section IV Part C of EE143 Lab Report 2
# plotting the CV curve for a Gate Oxide MOS cap, Device 4

# acc. | V_fb | depl. | V_t | inv.

import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np

def cv_curve_fn_hf(gate_volts):
	# accumulation region
	if(gate_volts <= V_fb):
		return C_acc

	# depletion region
	if(V_fb < gate_volts <= V_t):
		return depletion_capacitance(gate_volts)

	# inversion region
	else:
		return C_acc

def cv_curve_fn_lf(gate_volts):
	# accumulation region
	if(gate_volts <= V_fb):
		return C_acc

	# depletion region
	if(V_fb < gate_volts <= V_t):
		return depletion_capacitance(gate_volts)

	# inversion region
	else:
		return C_lf

def depletion_capacitance(gate_volts):
	return 1 / (np.sqrt((1/(C_ox**2)) + (2*(gate_volts - V_fb))/(4.41e-16)))

# Farads / cm^2
C_ox = 3.354e-8 
C_acc = C_ox
# C_inv = 2.645e-8 # calculated from first principals
C_inv = C_acc # assuming HF conditions

# Volts
V_fb = -0.923
V_t = 1.41

# low frequency capacitance
C_lf = depletion_capacitance(V_t)

cv_curve_hf = np.vectorize(cv_curve_fn_hf)
cv_curve_lf = np.vectorize(cv_curve_fn_lf)

style.use("bmh")

v = np.arange(-10, 10, 0.02)
plt.plot(v, cv_curve_hf(v), label='High Frequency')
plt.plot(v, cv_curve_lf(v), label='Low Frequency', linestyle='--')
plt.axvline(x=V_fb, label='Flatband Voltage', linestyle=':', color='green')
plt.axvline(x=V_t, label='Threshold Voltage', linestyle=':', color='orange')

plt.title("Gate Oxide MOS Capacitor CV Curve")
plt.xlabel("Gate Voltage [V]")
plt.ylabel("Capacitance Per Area [F/cm^2]")
plt.ylim(0, 5e-8)

plt.grid(True)
plt.legend()

plt.show()
