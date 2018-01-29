#Created by Marina M. Dunn on January 20, 2018
##This program is designed to take in the particle type and
#total number of particles from 'ReadFile' as inputs to calculate
#the 3D distance and velocity, as well as mass, for any given
#particle type (disk, halo, etc).

#First, import Numpy and Astropy modules
import numpy as np
import astropy.units as u
from ReadFile import Read

#defining a function that will read in inputs and calculate the 3D distance, velocity, and mass of particle
def ParticlePropertyInfo(ptype, pnum, filename):
    #p_type refers to which particle we want, p_num is the partical number
    #filename is the file from which we are importing data, in this case, MW_000.txt

    #This calls the function to read in these 3 data types
    time, total_particles, data = Read(filename)

    #optional print statement to see if inputs were read in properly
    print "time in Gyr = %s" %time
    print "total number of particles = %s" %total_particles
    
    #We want to create an array that stores the index of particle types
    index = np.where(data['type']==ptype)
    
    #Here we are assigning mass, position, and velocity for the Nth particle
    mnew = data['m'][index]
    xnew = data['x'][index]
    ynew = data['y'][index]
    znew = data['z'][index]
    vxnew = data['vx'][index]
    vynew = data['vy'][index]
    vznew = data['vz'][index]
    
    #calculates the 3D distance and velocity
    #For particle p_num, the array index is p_num -1
    distance = np.sqrt(xnew[pnum-1]**2 + ynew[pnum-1]**2 + znew[pnum-1]**2)*u.kpc
    velocity = np.sqrt(vxnew[pnum-1]**2 + vynew[pnum-1]**2 + vznew[pnum-1]**2)*u.km/u.s
    distance_new = np.around(distance,3)
    velocity_new = np.around(velocity,3)
    m = mnew[pnum-1]*1e10*u.solMass

    return distance_new, velocity_new, m

#Print information for the 100th disk particle with SnapNumber 0
pnum = 100
ptype = 2
results = ParticlePropertyInfo(ptype,pnum,'MW_000.txt')

print "The 3D distance of the particle %s of type %s in kpc is: " %(pnum,ptype), results[0]
print "The 3D velocity of the particle %s of type %s in km/s is: " %(pnum,ptype), results[1]
print "The mass of the particle %s of type in Solar Masses is:", results[2]

distance_Ly = np.around(results[0].to(u.lyr),3)
velocity_Ly = np.around(results[1],3)

print "The 3D distance of the particle %s of type %s in ly is: " %(pnum,ptype), distance_Ly
print "The 3D velocity of the particle %s of type %s in km/s is: " %(pnum,ptype), velocity_Ly
