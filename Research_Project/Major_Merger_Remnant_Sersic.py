#Created by Marina Dunn on 4/3/18
#Part of Research Project for ASTR400B, Spring 2018, Dr. Gurtina Besla
#Last edited on 4/3/2018
####Goal of this code:
#Part 1: find the initial brightness profiles for M31 and Milky Way at Snapshot 0 and compare results to
#current literature
#Part 2: Find the resultant brightness profile post-merger of the remnant, and compare profile to expected
#predictions (whether it obeys elliptical Sersic index)

#This code will build upon code used throughout past Homeworks and InClassLabs. Please refer to these
#directories in my ASTR400B repo for further details.

###Part 1###

#First, I will import relevant modules
import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterofMass import CenterOfMass
from MassProfile import MassEnclosed
#import plotting modules
import matplotlib
import matplotlib.pyplot as plt

#Parameters for M31 and MW
ptype = 2 #p_type refers to which particle we want, which in this case is 2, since we only want disk particles
delta = 0.3 #Choosing a tolerance for COM

#---------------------------------------------------
#Working on reading in files for each galaxy to create new files with just a list of radii and masses for the disks
#and the use those for the HalfMassRadius and Sersic Profile
#---------------------------------------------------

#Creating a class that will have a set of functions designed to
class Sersic:

    def __init__(self, filename):
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
    
    ###Later Step###
    #defining a function to calculate the Half Mass radius, which will be used to find the Half Light radius of
    #each galaxy at snapshot 0
    #Half mass radius is the radius at where half the mass of the galaxy is contained
    def HalfMassRadius(self,Mtot,R,Mdisk):
  
        #set a mass to light ratio, say 1.5
        ML = 1.5

        #Half the total mass in units of 1e10
        HalfMass = Mtot/2.0/1e10
        
        #finding where mass profile yields half of total mass, using "np.logical_and"
        index = np.where(np.logical_and(Mdisk/1e10 < (HalfMass+0.1), Mdisk/1e10 > (HalfMass-0.1)))
        return R[index]
                         
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
                         
                         
