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
#import plotting modules
import matplotlib
import matplotlib.pyplot as plt

#Parameters for M31 and MW
ptype = 2 #p_type refers to which particle we want, which in this case is 2, since we only want disk particles
delta = 0.3 #Choosing a tolerance for COM

#Creating a class that will have a set of functions designed to
class

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
        self.MWdisk = 0.075e12 #Msun, taken from Homework 3
        self.M31disk = 0.12e12 #Msun, taken from Homework 3

     def MassEnclosed(self, ptype, radii):


