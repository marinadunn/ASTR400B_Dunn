# Created by Marina Dunn on 3/29/18.
# ASTR400B Spring 2018, Prof. Gurtina Besla


# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

M31Orbit = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True)
M33Orbit = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True)
M33Analytic = np.genfromtxt('M33Analytical_orbit.txt',dtype=None,names=True)

M33Simulation_separation = np.sqrt((M33Orbit['x']-M31Orbit['x'])**2+(M33Orbit['y']-M31Orbit['y'])**2+(M33Orbit['z']-M31Orbit['z'])**2)
M33Analytic_separation = np.sqrt((M33Analytic['x'])**2+(M33Analytic['y'])**2+M33Analytic['z']**2)

M33Simulation_velocity = np.sqrt((M33Orbit['vx']-M31Orbit['vx'])**2+(M33Orbit['vy']-M31Orbit['vy'])**2+(M33Orbit['vz']-M31Orbit['vz'])**2)
M33Analytic_velocity = np.sqrt((M33Analytic['vx'])**2+(M33Analytic['vy'])**2+M33Analytic['vz']**2)

#First plot M33-M31 position comparison

#plot figure
figure1 = plt.figure(1,figsize=(10,10))
ax = plt.subplot(111)
#plot the simulation vs. analytic predictions
plt.plot(M31Orbit['t']/1e6, M33Simulation_separation, c='b', label='Besla Simulation')
plt.plot(M33Analytic['t'], M33Analytic_separation, c='y', label='Analytic Solution')

#add plot title
plt.title('Separation of M33 and M31')

#add plot axis labels
plt.xlabel('Time (Gyr)', fontsize=18)
plt.ylabel('Separation (kpc)', fontsize=18)

#Add plot legend
plt.legend(loc='best',fontsize=9)

#save plot
plt.savefig('M33-M31_Separation.jpg')
plt.close()


#Next plot M33-M31 velocity comparison

#plot figure
figure2 = plt.figure(1,figsize=(10,10))
ax = plt.subplot(111)
#plot the simulation vs. analytic predictions
plt.plot(M31Orbit['t']/1e6, M33Simulation_velocity, c='r', label='Besla Simulation')
plt.plot(M33Analytic['t'], M33Analytic_velocity, c='g', label='Analytic Solution')

#add plot title
plt.title('Velocity Comparison of M33 and M31')

#add plot axis labels
plt.xlabel('Time (Gyr)', fontsize=18)
plt.ylabel('Velocity (km/s)', fontsize=18)

#Add plot legend
plt.legend(loc='best',fontsize=9)

#save plot
plt.savefig('M33-M31_Velocity.jpg')
plt.close()
