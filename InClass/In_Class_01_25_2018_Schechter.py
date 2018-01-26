#Created by Marina Dunn on January 25, 2018

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt


#Parameters used from Smith+2009 for K band
h = 70.4/100
psi_star = 1.66*10**-2
M_star = (-23.19 - 5*np.log(h))
alpha = -0.81

#First step of this program is to plot the Schechter Function
def Schechter(alpha, M_star, psi_star, M):
    
    psi = (0.4*np.log(10))*psi_star*(10**(0.4*(M_star - M)*(alpha+1)))*np.exp(-10**(0.4*(M_star - M)))

    return psi

M = np.arange(-26,-17,0.1) #M is an array with magnitude range of -17 to -26
results = Schechter(alpha, M_star, psi_star, M)

#Set plot parameters
plt.semilogy(M,results)
plt.xlim(-17,-26) #set x limits
plt.show()
