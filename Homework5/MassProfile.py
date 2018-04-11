#Created by Marina Dunn on February 7, 2018
#ASTR400B, Spring 2018, Dr. Gurtina Besla, Homework 5

#The goal of this program is to calculate the mass distribution of a galaxy at Snapnumber 0, and subsequently
#use this to find each
#galaxy's rotation curve

#Import relevant modules first
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterofMass import CenterOfMass
#import plotting modules, used in InClassLab solution 1
import matplotlib
import matplotlib.pyplot as plt

#inputs for this class will be the galaxy name (string) and snapshot number
class MassProfile:
    def __init__(self, galaxy, snap):
        #Next, for a given filename, we only want the first characters that specify which galaxy we are talking
        #about
        #We will add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        #Then, we want the part of the filename that specifies the snapnumber
        #remove everything except the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt' #stores filename as global property
        #i.e. if I wanted "MW_010.txt", I would input "MW" and "10"
        #optional print statement to see if filename is printed correctly
        #print filename
        
        #Insert gravitational constant, need to adjust units; using the G below, in units of kpc^3/Gyr^2/Msun
        #Stores G as global property
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
    
        #Read in the filename, and data for particles of a certain ptype
        self.time, self.total_particles, self.data = Read(self.filename)
        
        #We want to store the galaxy name as a 'global property' self.gname
        self.gname = galaxy
    
        #Now we will store the positions and mass of the particles
        self.m = self.data['m']#*u.Msun
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
    
    #Next we will define a function that, with given a radius (kpc) fora galaxy's COM position and a vector
    #component, will calculate its mass
    #Inputs will be particle type, and an array containing radii for COM
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
    
    #Now we will create a function that will take in the radii array and return a mass array (in solar masses)
    #containing the total mass of all the
    #galactic components(halo+disk+bulge(if applicable))
    def MassEnclosedTotal(self, radii):
        #recall ptype = 1 is for halo, 2 is for disk, 3 is for bulge
        #M33 does not have a bulge component, which will have to be accounted for
        halomass = self.MassEnclosed(1, radii) #ptype = 1
        diskmass = self.MassEnclosed(2, radii) #ptype = 2
        bulgemass = self.MassEnclosed(3, radii) #ptype = 3
        Mtot = halomass + diskmass + bulgemass
        
        if (self.gname == 'M33'): #M33 doesn't have bulge component
            Mtot = halomass + diskmass
    
        return Mtot

    #Next, we will create a function that can compute the enclosed mass within a given radius (kpc), using the theoretical Hernquist profile
    #Recall that the Hernquist sphere model has a finite-density core and falls off as r^-4 at large radii, and has a cusp, where the profile steepens
    #in the core (original paper can be found on ApJ: DOI 10.1086/168845)
    #Consult InClassLab2 for similar problem for finding halo mass using Hernquist profile
    #input: radius array (kpc), scale factor 'a', and halo mass Mhalo
    def HernquistMass(self, radii, a, Mhalo):
        HernquistMass = Mhalo*radii**2./(radii+a)**2.*u.Msun
        #return halo mass in solar masses
        return HernquistMass
                      
    #Creating a function that takes in the radii array and ptype, and returns an array containing circular speeeds in [km/s]
    def CircularVelocity(self, ptype, radii):
        #Find enclosed mass of certain ptype
        mass = self.MassEnclosed(ptype, radii)
        
        #Calculate circular velocity
        Vcirc = np.around(np.sqrt(self.G*mass/radii/u.kpc),2)
        return Vcirc
    
    #Defining a function that takes in the radius array, and gives back an array of circular velocities in [km/s], using total mass enclosed
    def CircularVelocityTotal(self, radii):
         #Find enclosed mass of given radii
        mass = self.MassEnclosedTotal(radii)
        
        #Calculate circular velocity
        Vcirc_tot = np.around(np.sqrt((self.G*mass)/(radii/u.kpc)),2)
                      
        return Vcirc_tot
                      
    #Defining a function that takes in radius array, scale factor 'a', and halo mass Mhalo, and computes circular speed based on the Hernquist mass
    #profile, and returns an array of circular velocities in [km/s]
    def HernquistVCirc(self, radii, Mhalo, a):
        #Find enclosed mass of given radii
        mass = self.HernquistMass(radii, Mhalo, a)*u.Msun
                      
        #Calculate circular velocity
        HernquistVCirc = np.around(np.sqrt((self.G*mass)/(radii/u.kpc)),2)
                      
        return HernquistVCirc
#test
RR = 30
testR = np.arange(0.1,RR+1,1.0)

###Calculations & Plotting###

##Consulting InClassLab4
                      
#Creating arrays of galaxy names, radii, and scale factor 'a' based on Mass Profile plots
                      
MW = MassProfile("MW",0)
M31 = MassProfile("M31",0)
M33 = MassProfile("M33",0)

#Creating an array to store radii, from 0 to 30 kpc, in intervals of 0.5 kpc
Radii = np.arange(0.5,30,0.5)
a = [62,62,25]

##Milky Way Mass Profile
MW_Mtot = 1.97e12 #taken from Homework 3
MW_scale = 61.0 #MW Hernquist scale length
                      
#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
    
#Plot the Masses vs. Radius
plt.semilogy(Radii, MW.MassEnclosed(1,Radii), color='pink', label='Halo Mass')
plt.semilogy(Radii, MW.MassEnclosed(2,Radii), color='purple', label='Disk Mass')
plt.semilogy(Radii, MW.MassEnclosed(3,Radii), color='red', label='Bulge Mass')
plt.semilogy(Radii, MW.MassEnclosedTotal(Radii), color='yellow', label='Total Mass')
plt.semilogy(Radii, MW.HernquistMass(Radii,MW_scale,MW_Mtot), color='blue', label='Hernquist Mass, a =61 kpc')
                      
# Adding the axis labels and title
plt.xlabel('Radius (kpc)', fontsize=16)
plt.ylabel(r'Log(Mass Enclosed(M$_\odot$))', fontsize=16)
plt.title('MW Mass Profile')

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
ax.set_rasterized(True)
plt.savefig('MW_Mass_Profile.eps',rasterized=True, dpi=350)
plt.close()

#------------------------------------------------------------
##Plot M31 Mass Profile

M31_Mtot = 1.921e12 #taken from Homework 3
M31_scale = 62.0 #M31 Hernquist scale length

#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#Plot the Masses vs. Radius
plt.semilogy(Radii, M31.MassEnclosed(1,Radii), color='pink', label='Halo Mass')
plt.semilogy(Radii, M31.MassEnclosed(2,Radii), color='purple', label='Disk Mass')
plt.semilogy(Radii, M31.MassEnclosed(3,Radii), color='red', label='Bulge Mass')
plt.semilogy(Radii, M31.MassEnclosedTotal(Radii), color='yellow', label='Total Mass')
plt.semilogy(Radii, M31.HernquistMass(Radii,M31_scale,M31_Mtot), color='blue', label='Hernquist Mass, a =62 kpc')
 
# Adding the axis labels and title
plt.xlabel('Radius (kpc)', fontsize=16)
plt.ylabel(r'Log(Mass Enclosed(M$_\odot$))', fontsize=16)
plt.title('M31 Mass Profile')

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
ax.set_rasterized(True)
plt.savefig('M31_Mass_Profile.eps',rasterized=True, dpi=350)
plt.close()
    
#------------------------------------------------------------
##Plot M33 Mass Profile
M33_Mtot = 0.187e12 #taken from Homework 3
M33_scale = 25.0 #M33 Hernquist scale length

#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#Plot the Masses vs. Radius
plt.semilogy(Radii, M33.MassEnclosed(1,Radii), color='pink', label='Halo Mass')
plt.semilogy(Radii, M33.MassEnclosed(2,Radii), color='purple', label='Disk Mass')
plt.semilogy(Radii, M33.MassEnclosedTotal(Radii), color='yellow', label='Total Mass')
plt.semilogy(Radii, M33.HernquistMass(Radii,M33_scale,M33_Mtot), color='blue', label='Hernquist Mass, a =25 kpc')

# Adding the axis labels and title
plt.xlabel('Radius (kpc)', fontsize=16)
plt.ylabel(r'Log(Mass Enclosed(M$_\odot$))', fontsize=16)
plt.title('M33 Mass Profile')

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
ax.set_rasterized(True)
plt.savefig('M33_Mass_Profile.eps',rasterized=True, dpi=350)
plt.close()
    
#------------------------------------------------------------
##Plot MW Circular Velocity

#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#Plot the Masses vs. Radius
plt.semilogy(Radii, MW.CircularVelocity(1,Radii), color='pink', label='Halo Mass')
plt.semilogy(Radii, MW.CircularVelocity(2,Radii), color='purple', label='Disk Mass')
plt.semilogy(Radii, MW.CircularVelocity(3,Radii), color='purple', label='Bulge Mass')
plt.semilogy(Radii, MW.CircularVelocityTotal(Radii), color='yellow', label='Total Mass')
plt.semilogy(Radii, MW.HernquistVCirc(Radii,MW_scale,MW_Mtot), color='blue', label='Hernquist Mass, a =61 kpc')

# Adding the axis labels and title
plt.xlabel('Radius (kpc)', fontsize=16)
plt.ylabel('Circular Velocity (km/s)', fontsize=16)
plt.title('MW Circular Velocity')

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
ax.set_rasterized(True)
plt.savefig('MW_Circular_Velocity.eps',rasterized=True, dpi=350)
plt.close()
    
#------------------------------------------------------------
##Plot M31 Circular Velocity

#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#Plot the Masses vs. Radius
plt.semilogy(Radii, M31.CircularVelocity(1,Radii), color='pink', label='Halo Mass')
plt.semilogy(Radii, M31.CircularVelocity(2,Radii), color='purple', label='Disk Mass')
plt.semilogy(Radii, M31.CircularVelocity(3,Radii), color='purple', label='Bulge Mass')
plt.semilogy(Radii, M31.CircularVelocityTotal(Radii), color='yellow', label='Total Mass')
plt.semilogy(Radii, M31.HernquistVCirc(Radii,MW_scale,MW_Mtot), color='blue', label='Hernquist Mass, a =62 kpc')

# Adding the axis labels and title
plt.xlabel('Radius (kpc)', fontsize=16)
plt.ylabel('Circular Velocity (km/s)', fontsize=16)
plt.title('M31 Circular Velocity')

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
ax.set_rasterized(True)
plt.savefig('M31_Circular_Velocity.eps',rasterized=True, dpi=350)
plt.close()
    
#------------------------------------------------------------
##Plot M33 Circular Velocity

#Plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#Plot the Masses vs. Radius
plt.semilogy(Radii, M33.CircularVelocity(1,Radii), color='pink', label='Halo Mass')
plt.semilogy(Radii, M33.CircularVelocity(2,Radii), color='purple', label='Disk Mass')
plt.semilogy(Radii, M33.CircularVelocityTotal(Radii), color='yellow', label='Total Mass')
plt.semilogy(Radii, M33.HernquistVCirc(Radii,MW_scale,MW_Mtot), color='blue', label='Hernquist Mass, a =25 kpc')

# Adding the axis labels and title
plt.xlabel('Radius (kpc)', fontsize=16)
plt.ylabel('Circular Velocity (km/s)', fontsize=16)
plt.title('M33 Circular Velocity')

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
ax.set_rasterized(True)
plt.savefig('M33_Circular_Velocity.eps',rasterized=True, dpi=350)
plt.close()

