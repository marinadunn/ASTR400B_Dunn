#Created by Marina Dunn on February 7, 2018
#ASTR400B, Spring 2018, Dr. Gurtina Besla, Homework 5

#The goal of this program is to calculate the mass distribution of a galaxy at Snapnumber 0, and subsequently use this to find each
#galaxy's rotation curve

#Import relevant modules first
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterofMass import CenterofMass
#import plotting modules, used in InClassLab solution 1
from scipy.integrate import quad
import matplotlib
import matplotlib.pyplot as plt

#inputs for this class will be the galaxy name (string) and snapshot number
class MassProfile:
    def _init_(self, galaxy, snap):
        #Next, for a given filename, we only want the first characters that specify which galaxy we are talking about
        #We will add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        #Then, we want the part of the filename that specifies the snapnumber
        #remove everything except the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = '/home/marinadunn'+ "%s_"%(galaxy) + ilbl + '.txt' #stores filename as global property
        #i.e. if I wanted "MW_010.txt", I would input "MW" and "10"
        #optional print statement to see if filename is printed correctly
        #print filename
    
        self.delta = 0.3 #set a value for the tolerance
    
        #Insert gravitational constant, need to adjust units; using the G below, in units of kpc^3/Gyr^2/Msun
        #Stores G as global property
        self.G = G.to(u.kpc*u.km**2/u.s**2/u/Msun)
    
        #We want to store the galaxy name as a 'global property' self.gname
        galaxy = self.gname
    
        #Now we will store the positions and mass of the particles
        self.m = self.data['m']
        self.x = self.data['x']
        self.y = self.data['y']
        self.z = self.data['z']
    
        #Read in the filename, and data for particles of a certain ptype
        self.time, self.total_particles, self.data = Read(self.filename)

    #Next we will define a function that, with given a radius (kpc) fora galaxy's COM position and a vector component, will calculate its mass
    #Inputs will be particle type, and an array containing radii for COM
    def MassEnclosed(self, ptype, radii):
        #For all particles of a specific ptype
        index = np.where(self.data['type'] == ptype)
        
        
        #Now we will store the positions and mass of the particles for certain ptype
        m = self.m[index]
        x = self.x[index]
        y = self.y[index]
        z = self.z[index]
        
        #Finding positions of particles in COM frame
        COMgal = CenterofMass(self.filename, ptype)
        COMgal_Xcom, COMgal_Ycom, COMgal_Zcom = COMgal.COM_P(self.delta)
        radius = np.sqrt((self.x[index]-COMgal_Xcom)**2 + (self.y[index]-COMgal_Ycom)**2 + (self.z[index]-COMgal_Zcom)**2)

        #initialize the array using np.zeros, storing the mass
        mass = np.zeros(len(radii))*1e10*u.Msun
        
        pmass = self.m[index]
      
        #We want to loop over the radius array in order to define particles within a given radius for every array element
        for i in range(len(radii)):
            #Find the mass within the radius
            j = np.where(radius<=radii[i])[0]
            mass[i] = np.sum(m[j])*1e10
        
        #Want to return an array containing masses in units of solar masses
        return Mass
    
    #Now we will create a function that will take in the radii array and return a mass array (in solar masses) containing the total mass of all the
    #galactic components(halo+disk+bulge(if applicable))
    def MassEnclosedTotal(self, radii):
        #recall ptype = 1 is for halo, 2 is for disk, 3 is for bulge
        #M33 does not have a bulge component, which will have to be accounted for
        halomass = self.MassEnclosed(1, radii) #ptype = 1
        diskmass = self.MassEnclosed(2, radii) #ptype = 2
                      
        if self.gname != 'M33':
            Mtot = halomass + diskmass
       
        else:
            bulgemass = self.MassEnclosed(3, radii) #ptype = 3
            Mtot = halomass + diskmass + bulgemass
        
        return Mtot

    #Next, we will create a function that can compute the enclosed mass within a given radius (kpc), using the theoretical Hernquist profile
    #Recall that the Hernquist sphere model has a finite-density core and falls off as r^-4 at large radii, and has a cusp, where the profile steepens
    #in the core (original paper can be found on ApJ: DOI 10.1086/168845)
    #Consult InClassLab2 for similar problem for finding halo mass using Hernquist profile
    #input: radius array (kpc), scale factor 'a', and halo mass Mhalo
    def HernquistMass(self, radii, a, Mhalo):
        HernquistMass = (Mhalo*(radii**2))/((radii+a)**2)/1e12*u.Msun
        #return halo mass in solar masses
        return HernquistMass
                      
    #Creating a function that takes in the radii array and ptype, and returns an array containing circular speeeds in [km/s]
    def CircularVelocity(self, ptype, radii):
        #Find enclosed mass of certain ptype
        Mass = self.MassEnclosed(ptype, radii)*u.Msun
        
        #Calculate circular velocity
        circular_velocity = np.around(np.sqrt((self.G*Mass)/(radii*u.kpc)),2)
        return circular_velocity
    
    #Defining a function that takes in the radius array, and gives back an array of circular velocities in [km/s], using total mass enclosed
    def CircularVelocityTotal(self, radii):
         #Find enclosed mass of given radii
        Mass = MassEnclosedTotal(radii)*u.Msun
        
        #Calculate circular velocity
        circular_velocity_total = np.around(np.sqrt((self.G*Mass)/(radii*u.kpc)),2)
                      
        return circular_velocity_total
                      
    #Defining a function that takes in radius array, scale factor 'a', and halo mass Mhalo, and computes circular speed based on the Hernquist mass
    #profile, and returns an array of circular velocities in [km/s]
    def HernquistVCirc(self, radii, Mhalo, a):
        #Find enclosed mass of given radii
        Mass = self. HernquistMass(radii, Mhalo, a)*u.Msun
                      
        #Calculate circular velocity
        HernquistVCirc = np.around(np.sqrt((self.G*Mass)/(radii*u.kpc)),2)
                      
        return HernquistVCirc


###Calculations & Plotting###

##Consulting InClassLab4
                      
#Creating arrays of galaxy names, radii, and scale factor 'a' based on Mass Profile plots
                      
galaxies = ['MW', 'M31', 'M33']
Radii = np.arange(0.1, 30, 0.5)
a = [62, 62, 25]
import matplotlib
import matplotlib.pyplot as plt

for ii in range(len(galaxies)):
    #Initialize galaxy
    galaxy = MassProfile(galaxies[ii], 00)
                      
    #Find halo and diak masses within radii
    Mhalo = galaxy.MassEnclosed(1, Radii)
    Mdisk = galaxy.MassEnclosed(2, Radii)
              
    #Find bulge masses within radii, except for M33
    if galaxy.gname != 'M33':
        Mbulge = galaxy.MassEnclosed(3, Radii)
                      
    #Find total masses within radii
    Mtot = galaxy.MassEnclosedTotal(Radii)
        
    #Find total halo mass
    index = np.where(galaxy.data['type'] == 1)
    m = galaxy.data['m'[index]]
    MhaloTot = np.sum(m)*1e10
                      
    #Find Hernquist Mass within radii
    Hernquistmass = galaxy.HernquistMass(Radii, a[ii], MhaloTot)
                      
    #Plot
    fig = plt.figure()
                      
    #Plot the Masses vs. Radius
    plt.semilogy(Radii, Mhalo, color='pink', label='Halo Mass')
    plt.semilogy(Radii, Mdisk, color='purple', label='Disk Mass')
    if galaxy.gname != 'M33':
        plt.semilogy(Radii, Mbulge, color='red', label='Bulge Mass')
    plt.semilogy(Radii, Mtot, color='yellow', label='Total Mass')
    plt.semilogy(Radii, Hernquistmass, color='blue', label='Hernquist Mass, a =' +str(a[ii])
                      
    # Adding the axis labels and title
    plt.xlabel('Radius (kpc)', fontsize=22)
    plt.ylabel('Mass (Msun)', fontsize=22)
    plt.title(galaxies[ii] + ' Mass Profile')

    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size
    matplotlib.rcParams['ytick.labelsize'] = label_size

    # add a legend with some customizations.
    plt.legend(loc='best',fontsize='x-large')

    # Save to a file
    plt.savefig(galaxies[ii] + '_Mass_Profile.eps')

for jj in range(len(galaxies)):
    #Initialize galaxy
    galaxy = MassProfile(galaxies[jj], 00)
                     
    #Find halo and diak circular velocities within radii
    VCircHalo = galaxy.CircularVelocity(1, Radii)
    VCircDisk = galaxy.CircularVelocity(2, Radii)
                     
    #Find bulge circular velocity within radii, except for M33
    if galaxy.gname != 'M33':
        VCircBulge = galaxy.CircularVelocity(3, Radii)
                     
    #Find total circular velocities within radii
    VCircTot = galaxy.CircularVelocityTotal(Radii)
                     
    #Find total halo mass
    index = np.where(galaxy.data['type'] == 1)
    m = galaxy.data['m'[index]]
    MhaloTot = np.sum(m)*1e10
                     
    #Find Hernquist Mass within radii
    VCircHernquist = galaxy.HernquistVCirc(Radii, a[ii], MhaloTot)
                     
    #Plot
    fig = plt.figure()
                     
    #Plot the Masses vs. Radius
    plt.semilogy(Radii, VCircHalo, color='pink', label='Halo Circular Velocity')
    plt.semilogy(Radii, VCircDisk, color='purple', label='Disk Circular Velocity')
    if galaxy.gname = 'M33':
        plt.semilogy(Radii, VCircBulge, color='red', label='Bulge Circular Velocity')
    plt.semilogy(Radii, VCircTot, color='yellow', label='Total Circular Velocity')
    plt.semilogy(Radii, VCircHernquist, color='blue', label='Hernquist Circular Velocity, a =' +str(a[ii])
                                  
    # Adding the axis labels and title
    plt.xlabel('Radius (kpc)', fontsize=22)
    plt.ylabel('Circular Velocity (km/s)', fontsize=22)
    plt.title(galaxies[jj] + ' Velocity Profile')
                                  
                                  
    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size
    matplotlib.rcParams['ytick.labelsize'] = label_size
                                  
    # add a legend with some customizations.
    plt.legend(loc='best',fontsize='x-large')
                                  
    # Save to a file
    plt.savefig(galaxies[jj] + '_Velocity_Profile.eps')
