#Created by Marina M. Dunn on January 16, 2018
#This script reads in and returns the time elapsed, the total number of particles,
#the particles' type and masses,and their velocities in the
#x, y, and z directions

#First, import Numpy and Astropy modules
import numpy as np
import astropy.units as u

#We now define a function that will read in data from a file
def Read(filename):
    
    file = open(filename,'r') #opens the file
    
    #Reads in first line of file, which is the time elapsed
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Gyr
    #print time
    #Reads in the second line of the file, which is the total number of particles
    line2 = file.readline()
    label, value = line2.split()
    total_particles = float(value)
    #print total_particles
    file.close() #closes the file
    
    #We want to store the rest of the file in an array,
    #but skip the first 3 lines of the file, which are separated by white spaces
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    #print data
    #We want to get back the time, total number of particles, and rest of the data array
    return time, total_particles, data

#print Read('MW_000.txt')
#print Read('MW_000.txt')[2]
#The program 'ParticleProperties' will use the particle type and
#number from 'ReadFile' as inputs to calculate the 3D distance and
#velocity, as well as mass, for any given particle type (disk, halo, etc).
