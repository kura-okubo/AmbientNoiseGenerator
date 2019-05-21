import numpy as np

#To plot dispersion curve, we used python wrapper for surf96.
#For instration, see https://github.com/miili/pysurf96

from pysurf96 import surf96
import matplotlib.pyplot as plt

# Define the velocity model in km and km/s
#pick up numbers from Haskell(1953) Case II

thickness = np.array([28.38, 0])
vs = np.array([3.39, 4.65])
vp = np.array([6.14, 8.26])
rho = np.array([2.7, 3.0])

# Periods we are interested in
periods = np.linspace(.1, 100., 50)

phase_velocities = surf96(thickness, vp, vs, rho, periods,
                    wave='rayleigh', mode=1, velocity='phase', flat_earth=True)

group_velocities = surf96(thickness, vp, vs, rho, periods,
                    wave='rayleigh', mode=1, velocity='group', flat_earth=True)

#nondimensionalize
n_periods = [x*vs[0]/thickness[0] for x in periods]
n_phase_velocities = [x/vs[0] for x in phase_velocities]
n_group_velocities = [x/vs[0] for x in group_velocities]

plt.plot(n_periods, n_phase_velocities, label='Phase velocity')
plt.plot(n_periods, n_group_velocities, label='Group velocity')
#plt.show()
plt.legend()
plt.xlabel('Nondimensionalized Period T cs/H')
plt.ylabel('Nondimensionalized phase velocity cr/cs')
plt.savefig('dispersion_curve.png')

#save profiles
with open('./dispersion_profile.in','w') as fo:
	
	fo.write("#nondimensinoalized period, nondimensinoalized velocity\n")

	for i in range(len(n_periods)):
		fo.write("%12.8f, %12.8f\n"%(n_periods[i], n_phase_velocities[i]))

with open('./dispersion_group_velocity_profile.in','w') as fo:
	
	fo.write("#nondimensinoalized period, nondimensinoalized group velocity\n")

	for i in range(len(n_periods)):
		fo.write("%12.8f, %12.8f\n"%(n_periods[i], n_group_velocities[i]))