#Marina Dunn, ASTR400B Spring 2018, Dr. Gurtina Besla
#Last edited: 10 April, 2018
# In Class Lab 4 
# Surface Brightness Profile

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
#get_ipython().magic('matplotlib inline')

# Read in the provided file
# mass profile of the bulge  M(r) 
# Headers:   'Bulge' for mass profiles
#            'R' for radius 
# usage :   Mass = BulgeMass['Bulge']   Radius = BulgeMass['R']
BulgeMass = np.genfromtxt('BulgeMass.dat',dtype=None,names=True) 

# storing data from file 
Mass = BulgeMass['Bulge']
R = BulgeMass['R']

# Total Bulge mass 
Btot = 0.96e10 #Msun

# set a mass to light ratio, say 1.5
ML = 1.5

# Determine Re: The Half Mass Radius
def HalfMassRadius(Bmass,R,tot):
    # input, Bulge mass profile, Radius, Total Mass of Bulge
    # returns: Radius where mass is half the total mass 
  
    # half the total mass in units of 1e10
    HalfMass= tot/2.0/1e10
    
    # find where mass profile yields half the total mass 
    # note the use of "np.logical_and"   
    index = np.where( np.logical_and(Bmass/1e10 < (HalfMass+0.1), Bmass/1e10 > (HalfMass-0.1)))
 
    return R[index]

#Determine Half Mass Radius for Bulge 
Re = HalfMassRadius(BulgeMass['Bulge'],BulgeMass['R'],Btot)
print(Re)

#####Start InClass edit here######
#defining a function that takes in Re, radius, sersic index, M/L ratio, and total mass of bulge and returns the function
n=4
def Sersic(Re,R,n,ML,Btot):
    L = Btot/ML
    Ie = L/7.2/np.pi/(Re)**2.
    Sersic = Ie*np.exp(-7.67*((R/Re)**(1.0/n)-1.0))
    return Sersic

Sersic = Sersic(Re,R,n,ML,Btot)

# Plot the Bulge density profile vs
# the Sersic profile
####################################


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)


# bulge luminositiy density :   bulge mass profile/ Volume
BulgeI = BulgeMass['Bulge']/4.0*3.0/R**3/np.pi/ML

# plot the bulge luminosity density as a proxy for surface brighntess
plt.semilogy(R,BulgeI, color='green', linestyle="--",linewidth=3, label='Bulge Density')

#Plotting Sersic fit to the surface brightness
# plot the sersic bulge density profile
plt.semilogy(R,Sersic, color='pink', linestyle="--",linewidth=3, label='Sersic n=4')


# Add axis labels
plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Log(I)', fontsize=22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


# Save to a file
ax.set_rasterized(True)
plt.savefig('BulgeSersicProfile.eps', rasterized=True, dpi=350)
