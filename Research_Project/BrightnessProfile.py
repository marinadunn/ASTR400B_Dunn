#Created by Marina Dunn on 4/4/18
#Part of Research Project for ASTR400B, Spring 2018, Dr. Gurtina Besla
#Last edited on 5/4/2018

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
get_ipython().magic('matplotlib inline')


###Part2###
#I want to define very close, concentric shells and calculate the mass enclosed in each in order to find the flux.
#I will do this out until approximately 20 kpc. I expect the plot of this to have a negative slope because the disk particles
#are denser towards the center of the galaxies

while (r < 20 ): #set up while loop
    delta_r = 0.2
    #Find the mass within a thin shell with volume 4/3*pi*(r1^3-r2^3), where r1 and r2 define the width of the shell
    V = (4./3.)*(np.pi)*(radii[radius]**3. - radii[radius - 1.]**3.)
    
    index = np.where((radius[R]) & (radii - delta_r)): #setting an upper and lower limits for the shell width
    radii = radii + 2.*delta_r #Want delta_r to be small but don't want to overlap with other shells so that particles are
    #counted twice
    index = np.where(self.data['type'] == ptype)
        
        mass = np.zeros(
                    
                      
                        
                        """
                            #defining a function to calculate the Half Mass radius, which will be used to find the Half Light radius of
                            #each galaxy at snapshot 0
                            #Half mass radius is the radius at where half the mass of the galaxy is contained
                            #inputs are disk mass profile, radius, total mass of disk
                            def HalfMassRadius(self,Dmass,radii,Mtot):
                            
                            #Half the total mass in units of 1e10
                            HalfMass = Mtot/2.0/1e10
                            
                            #finding where mass profile yields half of total mass, using "np.logical_and"
                            index = np.where(np.logical_and(Dmass/1e10 < (HalfMass+0.1), Dmass/1e10 > (HalfMass-0.1)))
                            return radii[index]
                            
                            ###Finding the initial Sersic indexes for each galaxy
                            def MWSersic(self,MW_HMR,r,n,ML,MWtot):
                            ##For MW
                            
                            #Luminosity
                            L = self.MWtot/ML
                            Ie = L/7.2/np.pi/MW_HMR**2
                            
                            return Ie*np.exp(-7.67*((r/self.MW_HMR)**(1.0/n)-1.0))
                            
                            def M31Sersic(self,M31_HMR,r,n,ML,M31tot):
                            ##For M31
                            
                            #Luminosity
                            L = self.M31tot/ML
                            Ie = L/7.2/np.pi/M31_HMR**2
                            
                            return Ie*np.exp(-7.67*((r/self.M31_HMR)**(1.0/n)-1.0))
                            """
                        
                        #----------------------------------------------------------------------------------------------------------------
                        
                        
                        
                        
                        
                        MWI = InitialMass("MW",0)
                        M31I = InitialMass("M31",0)
                        
                        #Creating an array to store radii, from 0 to 30 kpc, in intervals of 0.5 kpc
                        Radii = np.arange(0.5,30,0.5)
                        a = [62,62,25]
                        
                        #Determine Half Mass Radius for MW Disk
                        MW_HMR = MWI.HalfMassRadius(self.MWdisk,radii,self.MWtot)
                        print MW_HMR
                        
                        #MW Disk Luminosity density: disk mass profile/volume
                        MWDiskI = self.MWdisk/4.0*3.0/radii**3/np.pi/ML
                        
                        #Determine Half Mass Radius for M31 Disk
                        M31_HMR = M31I.HalfMassRadius(self.M31disk,radii,self.M31tot)
                        print M31_HMR
                        
                        #M31 Disk Luminosity density: disk mass profile/volume
                        M31DiskI = self.M31disk/4.0*3.0/radii**3/np.pi/ML
                        
                        
                        # I want to find the initial surface brightness profiles for M31 and Milky Way at snapshot 0 and compare to Sersic
                        #indices.
                        #----------------------------------------------------------------------------------------------------------------
                        
                        
                        ###Plot Disk MW density profile vs sersic profile at snapshot 0
                        ###Also trying to overplot the initial brightness profiles of M31 and Milky Way and compare to Sersic indexes
                        #Plot
                        fig = plt.figure(figsize=(10,10))
                        ax = plt.subplot(111)
                        
                        #plot disk luminosity density as a proxy for surface brightness
                        plt.semilogy(radii,MWDiskI, color='green',linewidth=2,label='MW Disk Density')
                        #add sersic fit to surface brightess sersic fit
                        plt.semilogy(radii,Sersic(MW_HMR,radii,4,ML,self.MWdisk), color='red',linewidth=2,label='MW Initial Sersic')
                        
                        #plot disk luminosity density as a proxy for surface brightness
                        plt.semilogy(radii,M31DiskI, color='pink',linewidth=2,label='M31 Disk Density')
                        #I decided to use astropy.modeling instead of scipy.optimize.curve_fit
                        #add sersic fit to surface brightess sersic fit
                        plt.semilogy(radii,Sersic(M31_HMR,radii,4,ML,self.M31disk), color='blue',linewidth=2,label='M31 Initial Sersic')
                        
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

                        
                        
