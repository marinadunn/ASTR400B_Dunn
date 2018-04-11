#Created by Marina Dunn on January 25, 2018
#InClassLab1
#Last edited: 4/10/18

#import modules
import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt

#Creating a Schechter luminosity function (In terms of magnitude) that takes in an array of absolute magnitudes
#Parameters used from Smith+2009 for K band
#h = 70.4/100
psi_star = 0.0166
M_star = -23.19 #- 5*np.log(h)
#alpha = -0.81

#First step of this program is to plot the Schechter Function
def Schechter(alpha, M_star, psi_star, M):
    
    return 0.4*np.log(10)*psi_star*10**(0.4*(M_star - M)*(alpha+1.0))*np.exp(-10**(0.4*(M_star - M)))

#Creating an aray that stores K band magnitudes
MK = np.arange(-26,-17,0.1)

#Plot Schechter Function
figure = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

plt.semilogy(MK,Schechter(-0.81, M_star, psi_star, MK), color='blue', linewidth=5, label='Smith+2009')
plt.semilogy(MK,Schechter(-0.6, M_star, psi_star, MK), color='pink', linewidth=3, label=r'low $\alpha$')
plt.semilogy(MK,Schechter(-1.35, M_star, psi_star, MK), color='yellow', linewidth=3, label=r'high $\alpha$')

#Add axis labels
plt.xlabel(r'M$_k$)', fontsize=22)
plt.ylabel(r'$\Phi$ (Mpc$^{-3}h^3$/mag)', fontsize=22)
           
#set axis limits
#plt.xlim(-17,-26)
           
#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
           
# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')
           
# Save to a file
ax.set_rasterized(True)
plt.savefig('Schechter.png', rasterized=True, dpi=350)


###IMF
from scipy.integrate import quad
#creating Salpeter function that defines the Salpeter IMF
#It takes in an array of stellar masses, M, and returns Z(M), the fractional number of stars
#Must determine the normaliztion, Zo, by integrating over mass from 0.1 to 120 Msun and setting equal to 1
def Salpeter(M,Mmin,Mmax):
    normaliztion = quad(lambda M: M**(-2.35), Mmin, Mmax) #finding magnitude of integral
    A = 1./normaliztion[0]
    return A*M**(-2.35) #normalized function

Z = quad(lambda M: Salpeter(M,0.1,120),1.0,120) #fraction of stars more massive than Sun
print("Fraction of stars more massive than Sun:", np.around(Z[0],3))
