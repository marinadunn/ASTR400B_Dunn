#Created by Marina Dunn on 4/3/18
#Part of Research Project for ASTR400B, Spring 2018, Dr. Gurtina Besla
#Last edited on 4/10/2018
####Goal of this code:
#Part 1: find the initial brightness profiles for M31 and Milky Way at Snapshot 0 and compare results to
#current literature
#Part 2: Find the resultant brightness profile post-merger of the remnant, and compare profile to expected
#predictions (whether it obeys elliptical Sersic index)

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
get_ipython().magic('matplotlib inline')


#Parameters for M31 and MW
ptype = 2 #p_type refers to which particle we want, which in this case is 2, since we only want disk particles
delta = 0.3 #Choosing a tolerance for COM


#Creating a class that will have a set of functions designed to
class MergerRemnantSersic:

    def __init__(self, galaxy, snap):
        #Next, for a given filename, we only want the first characters that specify which galaxy we are talking about
        #We will add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        #Then, we want the part of the filename that specifies the snapnumber
        #remove everything except the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt' #stores filename as global property
        #i.e. if I wanted "MW_010.txt", I would input "MW" and "10"
        #optional print statement to see if filename is printed correctly
        #print filename
        
        #This calls the function to read in these 3 data types
        self.time, self.total_particles, self.data = Read(filename)
        
        #We want to store the galaxy name as a 'global property' self.gname
        self.gname = galaxy
        
        #Insert gravitational constant, need to adjust units; using the G below, in units of kpc^3/Gyr^2/Msun
        #Stores G as global property
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #Using 'self' refers to the quantities common to an object; the values are stored so that data does not
        #need to be read in each time
        
        #This creates an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
        
        #This stores the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        self.filename = filename
    
        #initializing mass components
        #Below values taken from Homework 3
        self.MWhalo = 1.975e12 #Msun
        self.MWdisk = 0.075e12 #Msun
        self.MWbulge = 0.01e12 #Msun
        self.MWtot = 2.06e12 #Msun
        self.M31halo = 1.921e12 #Msun
        self.M31disk = 0.12e12 #Msun
        self.M31bulge = 0.019e12 #Msun
        self.M31tot = 2.06e12 #Msun
    
    
    def MassEnclosed(self, ptype, radii):
        #Finding positions of particles in COM frame
        COM = CenterOfMass(self.filename,2)

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
            
    def DiskMass(self, radii):
        DiskMass = self.MassEnclosed(2, radii) #ptype = 2 for disk particles
        return DiskMass


###Part 1### I want to find the initial surface brightness profiles for M31 and Milky Way at snapshot 0 and compare to Sersic
#indices.
#----------------------------------------------------------------------------------------------------------------

###Plot Disk MW density profile vs sersic profile at snapshot 0
###Also trying to overplot the initial brightness profiles of M31 and Milky Way and compare to Sersic indexes
#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#MW Disk Luminosity density: disk mass profile/volume
MWDiskI = MWDiskMass['Disk']/4.0*3.0/MWR**3/np.pi/ML

#plot disk luminosity density as a proxy for surface brightness
plt.semilogy(R,MWDiskI, color='pink',linewidth=2,label='MW Disk Density')
#add sersic fit to surface brightess sersic fit
plt.semilogy(R,Sersic(Re,MWR,4,ML,MWDiskMass), color='red',linewidth=2,label='Initial Sersic')

#MW Disk Luminosity density: disk mass profile/volume
M31DiskI = M31DiskMass['Disk']/4.0*3.0/M31R**3/np.pi/ML

#plot disk luminosity density as a proxy for surface brightness
plt.semilogy(R,M31DiskI, color='pink',linewidth=2,label='M31 Disk Density')
#I decided to use astropy.modeling instead of scipy.optimize.curve_fit
#add sersic fit to surface brightess sersic fit
plt.semilogy(R,Sersic(Re,M31R,4,ML,MWDiskMass), color='red',linewidth=2,label='Initial Sersic')

s1 = Sersic1D(amp=1, r_eff=5)
for n in range (1,10):
    s1.n = n
    plt.plot = (R, s1(R), color = str(float(n)/15))
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

# Save to a file
ax.set_rasterized(True)
plt.savefig('MW/M31_Initial_Disk_Sersic_Profile.eps',rasterized=True, dpi=350)
plt.close()

#----------------------------------------------------------------------------------------------------------------

###Part2###
#I want to define very close, concentric shells and calculate the mass enclosed in each in order to find the flux.
#I will do this out until approximately 20 kpc. I expect the plot of this to have a negative slope because the disk particles
#are denser towards the center of the galaxies

while (r < 20 ): #set up while loop
    delta_r = 0.2
    np.where((r + delta_r) & (r - delta_r)): #setting an upper and lower limits for the shell width
    r = r + 2.*delta_r #Want delta_r to be small but don't want to overlap with other shells so that particles are counted twice
    index = np.where(self.data['type'] == ptype)
        
        mass = np.zeros(
#set a mass to light ratio, say 1.5
ML = 1.5
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

Re = HalfMassRadius(DiskMass['Disk'],DiskMass['R'],Dtot)
print(Re)
    
    ###Finding the initial Sersic indexes for each galaxy
    def MWSersic(Re,r,n,ML,MWtot)
        ##For MW
                         
        #Luminosity
        L = self.MWtot/ML
        Ie = L/7.2/np.pi/Re**2
        
        return Ie*np.exp(-7.67*((r/Re)**(1.0/n)-1.0))
       
    def M31Sersic(Re,r,n,ML,M31tot)
        ##For M31
                         
        #Luminosity
        L = self.MWtot/ML
        Ie = L/7.2/np.pi/Re**2
                         
        return Ie*np.exp(-7.67*((r/Re)**(1.0/n)-1.0))


#----------------------------------------------------------------------------------------------------------------







