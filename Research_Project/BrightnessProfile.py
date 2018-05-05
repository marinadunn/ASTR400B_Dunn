#Created by Marina Dunn on 4/4/18
#Part of Research Project for ASTR400B, Spring 2018, Dr. Gurtina Besla
#Last edited on 5/5/2018

#This code will build upon code used throughout past Homeworks and InClassLabs. Please refer to these
#directories in my ASTR400B repo for further details.

#Goal of this code: Find the resultant brightness profile post-merger of the remnant, and compare profile to expected
#predictions (whether it obeys elliptical Sersic index)

#First, I will import relevant modules
import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterofMass import CenterOfMass
from MassProfile import *
#import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
#get_ipython().magic('matplotlib inline')

#Parameters for M31 and MW
ML = 1.5 #setting a mass-to-light ratio of about 1.5
def MassEnclosed(ptype, radii):
    #Finding positions of particles in COM frame
    COM = CenterOfMass(filename,2)
    
    #Store COM position of galaxy
    COMPgal = COM.COM_P(1.0,2.0)
    
    #For all particles of a specific ptype, creating an array that stores an index of particles
    index = np.where(self.data['type'] == ptype)
    
    #Now we will store the positions and mass of the particles for certain ptype
    mG = self.m[index]
    xG = self.x[index] - COMPgal[0]
    yG = self.y[index] - COMPgal[1]
    zG = self.z[index] - COMPgal[2]
    #Computing magnitude of 3D radius
    radius = np.sqrt(xG**2. + yG**2. + zG**2. )
    
    #initialize the array using np.zeros, storing the mass
    mass = np.zeros(np.size(radii))
    
    #We want to loop over the radius array in order to define particles within a given radius for every
    #array element
    for i in range(np.size(radii)):
        #Find the mass within the radius
        indexR = np.where(radius <= radii[i]*u.kpc)
        mass[i] = np.sum(mG[indexR])*1e10

        #Want to return an array containing masses in units of solar masses
        return mass*u.Msun

#radius array out to 30kpc
R = np.arange(0.5,30,0.5)
DMass = MassEnclosed(2, R) #ptype = 2

#use mass info in order to find the luminosity info.
#I will do this out until approximately 20 kpc. I expect the plot of this to have a negative slope because the disk particles
#are denser towards the center of the galaxies

#Below values taken from Homework 3
MWhalo = 1.975e12 #Msun
MWdisk = 0.075e12 #Msun
MWbulge = 0.01e12 #Msun
MWtot = 2.06e12 #Msun
M31halo = 1.921e12 #Msun
M31disk = 0.12e12 #Msun
M31bulge = 0.019e12 #Msun
M31tot = 2.06e12 #Msun

#defining a function to calculate the Half Mass radius, which will be used to find the Half Light radius of
#each galaxy at snapshot 0
#Half mass radius is the radius at where half the mass of the galaxy is contained
#inputs are disk mass profile, radius, total mass of disk
def HalfMassRadius(self,Dmass,R,Mtot):
    
  
    #Half the total mass in units of 1e10
    HalfMass = Mtot/2.0/1e10
                        
    #finding where mass profile yields half of total mass, using "np.logical_and"
    index = np.where(np.logical_and(Dmass/1e10 < (HalfMass+0.1), Dmass/1e10 > (HalfMass-0.1)))
    return R[index]

###Finding the initial Sersic indexes for each galaxy
def MWSersic(MW_HMR,r,n,ML,MWtot):
##For MW
                        
    #Luminosity
    L = MWtot/ML
    Ie = L/7.2/np.pi/MW_HMR**2
    MWsersic = Ie*np.exp(-7.67*((r/self.MW_HMR)**(1.0/n)-1.0))
    print MWsersic
    return MWsersic
                        
def M31Sersic(M31_HMR,r,n,ML,M31tot):
##For M31
                        
    #Luminosity
    L = M31tot/ML
    Ie = L/7.2/np.pi/M31_HMR**2
    M31sersic = Ie*np.exp(-7.67*((r/self.M31_HMR)**(1.0/n)-1.0))
    print M31sersic
    return M31sersic
                        
                        
 #----------------------------------------------------------------------------------------------------------------

#Determine Half Mass Radius for MW Disk
MW_HMR = HalfMassRadius(Dmass,R,MWtot)
print MW_HMR
                        
#MW Disk Luminosity density: disk mass profile/volume
MWDiskI = Dmass/4.0*3.0/R**3/np.pi/ML
                        
#Determine Half Mass Radius for M31 Disk
M31_HMR = HalfMassRadius(Dmass,R,M31tot)
print M31_HMR
                        
#M31 Disk Luminosity density: disk mass profile/volume
M31DiskI = Dmass/4.0*3.0/R**3/np.pi/ML
                        
"""
# I want to find the initial surface brightness profiles for M31 and Milky Way at snapshot 0 and compare to Sersic
#indices.
#----------------------------------------------------------------------------------------------------------------
#Creating an array to store radii, from 0 to 30 kpc, in intervals of 0.5 kpc
Radii = np.arange(0.5,30,0.5)
                        
###Plot Disk MW density profile vs sersic profile at snapshot 0
###Also trying to overplot the initial brightness profiles of M31 and Milky Way and compare to Sersic indexes
#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
                        
#plot disk luminosity density as a proxy for surface brightness
plt.semilogy(Radii,MWDiskI, color='green',linewidth=2,label='MW Disk Density')
#add sersic fit to surface brightess sersic fit
plt.semilogy(Radii,Sersic(MW_HMR,Radii,4,ML,MWdisk), color='red',linewidth=2,label='MW Initial Sersic')
                        
#plot disk luminosity density as a proxy for surface brightness
plt.semilogy(Radii,M31DiskI, color='pink',linewidth=2,label='M31 Disk Density')
#Could also use scipy.optimize.curve_fit instead of astropy.modeling for this step if desired
#add sersic fit to surface brightess sersic fit
plt.semilogy(Radii,Sersic(M31_HMR,Radii,4,ML,M31disk), color='blue',linewidth=2,label='M31 Initial Sersic')
                        
s1 = Sersic1D(amp=1, r_eff=5)
for n in range (1,10):
    s1.n = n
plt.plot = (radii, s1(radii))
plt.axis([1e-1, 30, 1e-2, 1e3])
plt.text(.25, 1.5, 'n=1')
plt.text(.25, 300, 'n=10')
                        
# Adding the axis labels and title
plt.xlabel('Log Radius (kpc)', fontsize=16)
plt.ylabel('Log Surface Brightness', fontsize=16)
                        
#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
                        
# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')
                        
plt.show()
                        
# Save to a file
ax.set_rasterized(True)
plt.savefig('MW/M31_Initial_Disk_Sersic_Profile.eps',rasterized=True, dpi=350)
plt.close()
                        
#----------------------------------------------------------------------------------------------------------------

"""
                        
