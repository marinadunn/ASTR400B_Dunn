#Created by Marina Dunn on 4/3/18
#Part of Research Project for ASTR400B, Spring 2018, Dr. Gurtina Besla
#Last edited on 5/4/2018
####Goal of this code:
#Part 1: find the initial brightness profiles for M31 and Milky Way at Snapshot 0 and compare results to
#current literature

#This code will build upon code used throughout past Homeworks and InClassLabs. Please refer to these
#directories in my ASTR400B repo for further details.


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

def InitialMassProfile(galaxy, snap, radii):
    #Next, for a given filename, we only want the first characters that specify which galaxy we are talking about
    #We will add a string of the filenumber to the value "000"
    ilbl = '000' + str(snap)
    #Then, we want the part of the filename that specifies the snapnumber
    #remove everything except the last 3 digits
    ilbl = ilbl[-3:]
    filename = "%s_"%(galaxy) + ilbl + '.txt'
    #i.e. if I wanted "MW_010.txt", I would input "MW" and "10"
    #optional print statement to see if filename is printed correctly
    #print filename
    
    #This calls the function to read in these 3 data types
    time, total_particles, data = Read(filename)
    
    #This stores the mass, positions, velocities of only the particles of the given type
    m = data['m']*1e10
    x = data['x']
    y = data['y']
    z = data['z']
    #vx = data['vx'][index]
    #vy = data['vy'][index]
    #vz = data['vz'][index]

    #initializing mass components
    #Below values taken from Homework 3
    MWhalo = 1.975e12 #Msun
    MWdisk = 0.075e12 #Msun
    MWbulge = 0.01e12 #Msun
    MWtot = 2.06e12 #Msun
    M31halo = 1.921e12 #Msun
    M31disk = 0.12e12 #Msun
    M31bulge = 0.019e12 #Msun
    M31tot = 2.06e12 #Msun

    #Finding positions of particles in COM frame
    COM = CenterOfMass(filename,2)

    #Store COM position of galaxy
    COMPgal = COM.COM_P(1.0,2.0) #Choosing a tolerance and VolDec for COM
    
    #mass profile of the objects
    mass_profile = MassProfile(galaxy,snap)
    DiskMass = mass_profile.MassEnclosed(2,radii)
    
    #For all particles of a specific ptype, creating an array that stores an index of particles of desired Ptype
    #In this case, we want ptype=2 for disk particles
    index = np.where(data['type'] == 2)
    
    #Now we will store the positions and mass of the particles for certain ptype
    mG = m[index]
    xG = x[index] - COMPgal[0]
    yG = y[index] - COMPgal[1]
    zG = z[index] - COMPgal[2]
    #Computing magnitude of 3D radius
    r = np.sqrt(xG**2. + yG**2. + zG**2.)
    print len(r)
    print r

    #initialize the array using np.zeros, storing the mass
    mass = np.zeros(np.size(radii))
    
    #We want to loop over the radius array in order to define particles within a given radius for every
    #array element
    for i in range(np.size(radii)):
        #Find the mass within a thin shell with volume 4/3*pi*(r1^3-r2^3), where r1 and r2 define the width of the shell
        V = (4./3.)*(np.pi)*(radii[i]**3. - radii[i - 1]**3.)
        index = np.where((r<radii[i]) & (r>radii[i - 1.]))
        #Find the mass within the radius
        print radii[i]
        print np.shape(index)
        mass[i - 1] = np.sum(mG[index])
            
    #Want to return an array containing masses in units of solar masses
    return mass#*u.Msun

#radius array out to 30kpc
R = np.arange(0.01,30,0.5)

MW = InitialMassProfile("MW",0,R)
M31 = InitialMassProfile("M31",0,R)


#Mass profiles for both galaxies at different key snapshots
MW_Mass1 = InitialMassProfile('MW', 0, R)
MW_Mass2 =  InitialMassProfile('MW', 70, R)
MW_Mass3 =  InitialMassProfile('MW', 480, R)
MW_Mass4 =  InitialMassProfile('MW', 710, R)
MW_Mass5 =  InitialMassProfile('MW', 800, R)
        
M31_Mass1 =  InitialMassProfile('M31', 0, R)
M31_Mass2 =  InitialMassProfile('M31', 70, R)
M31_Mass3 =  InitialMassProfile('M31', 480, R)
M31_Mass4 =  InitialMassProfile('M31', 710, R)
M31_Mass5 = InitialMassProfile('M31', 800, R)
        
#Plot mass profiles
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, xscale='log')
        
plt.semilogy(R,MW_Mass1,color='red',linewidth=2,label='MW 0 Gyr')
plt.semilogy(R,MW_Mass2,color='orange',linewidth=2,label='MW 1 Gyr')
plt.semilogy(R,MW_Mass3,color='yellow',linewidth=2,label='MW 6.86 Gyr')
plt.semilogy(R,MW_Mass4,color='green',linewidth=2,label='MW 10.14 Gyr')
plt.semilogy(R,MW_Mass5,color='blue',linewidth=2,label='MW 12 Gyr')
        
plt.semilogy(R,M31_Mass1,color='purple',linewidth=2,label='M31 0 Gyr')
plt.semilogy(R,M31_Mass2,color='pink',linewidth=2,label='M31 1 Gyr')
plt.semilogy(R,M31_Mass3,color='brown',linewidth=2,label='M31 6.86 Gyr')
plt.semilogy(R,M31_Mass4,color='firebrick',linewidth=2,label='M31 10.14 Gyr')
plt.semilogy(R,M31_Mass5,color='aqua',linewidth=2,label='M31 12 Gyr')
        
# Adding the axis labels and title
plt.xlabel('Radius (kpc)', fontsize=16)
plt.ylabel('Mass(Msun)', fontsize=16)
        
#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
        
# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')
        
# Save to a file
ax.set_rasterized(True)
plt.savefig('MW/M31_Mass_Profiles.eps',rasterized=True, dpi=350)
plt.close()


