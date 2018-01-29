#Created by Marina Dunn on January 23, 2018
#Goal of this program is to compute the mass of 3 galaxies in the Local Group (SnapNumber 0),
#the Milky Way itself, Andromeda(M31), and the Triangulum Galaxy(M33), the most massive members
#We do this by first finding the mass of the disk, halo, and bulge (if present) components, then summing them


#First, import Numpy and Astropy modules
import numpy as np
import astropy.units as u

#Call ReadFile program to use Read function
from ReadFile import Read

#The purpose of this function is to calculate the total mass of a galaxy
def ComponentMass(filename, ptype):
    #ptype can either refer to a halo (ptype=1), a disk (ptype=2), or a bulge (ptype=3)
    #filename is the file from which we are importing data
    
    #This calls the function to read in these 3 data types
    time, total_particles, data = Read(filename)
    
    #We want to create an array that stores the index of particle types
    index = np.where(data['type']==ptype)
    
    #Assigning a mass for the Nth particle
    mnew = data['m'][index]*1e10*u.solMass
    mass = np.sum(mnew) #summing over all the entries
    total_mass = np.around(mass,3)/1e12 #We want to round to 3 decimal places
    return total_mass


#execute the function for filename = MW_000.txt for Milky Way, filename = M31_000.txt for Andromeda, and
#filename = M33_000.txt for Triangulum
MW_halo = np.around(ComponentMass('MW_000.txt',ptype=1),3)
MW_disk = np.around(ComponentMass('MW_000.txt',ptype=2),3)
MW_bulge = np.around(ComponentMass('MW_000.txt',ptype=3),3)

M31_halo = np.around(ComponentMass('M31_000.txt',ptype=1),3)
M31_disk = np.around(ComponentMass('M31_000.txt',ptype=2),3)
M31_bulge = np.around(ComponentMass('M31_000.txt',ptype=3),3)

M33_halo = np.around(ComponentMass('M33_000.txt',ptype=1),3)
M33_disk = np.around(ComponentMass('M33_000.txt',ptype=2),3)
#M33 does not have a bulge, so its total mass will consist of only the halo and disk components
#sum up the masses of all the galactic components

MW_total = np.around((MW_halo + MW_disk + MW_bulge),3)
M31_total = np.around((M31_halo + M31_disk + M31_bulge),3)
M33_total = np.around((M33_halo + M33_disk),3)

#Finally, we want to calculate fbar, or the baryon fraction, for each of the galaxies, and the entire Local
#Group. This can be done by taking the ratio of total stellar mass to the total overall mass, including both
#dark and stellar.
fbar_MW = np.around(((MW_disk + MW_bulge) /MW_total),3)
fbar_M31 = np.around(((M31_disk + M31_bulge) /M31_total),3)
fbar_M33 = np.around(((M33_disk)/M33_total),3)
fbar_Local_Group = np.around(((MW_disk + MW_bulge + M31_disk + M31_bulge + M33_disk + M33_halo)/(MW_total + M31_total + M33_total)),3)

#prints the mass and number of particles for each of the three galaxies
#masses have the units of 10^12 Solar Masses, so the outputs for the masses are actually multiplied by 10^12
#solar masses
print "The mass of the MW halo is:", MW_halo
print "The mass of the MW disk is:", MW_disk
print "The mass of the MW bulge is:", MW_bulge
print "The mass of the entire MW is:", MW_total
                        
print "The mass of the M31 halo is:", M31_halo
print "The mass of the M31 disk is:", M31_disk
print "The mass of the M31 bulge is:", M31_bulge
print "The total mass of M31 is:", M31_total
                        
print "The mass of the M33 halo is:", M33_halo
print "The mass of the M33 disk is:", M33_disk
print "The total mass of M33 is:", M33_total

print "The baryon fraction for the Milky Way is:", fbar_MW
print "The baryon fraction for the Andromeda Galaxy is:", fbar_M31
print "The baryon fraction for the Triangulum Galaxy is:", fbar_M33
print "The baryon fraction for the Local Group is:", fbar_Local_Group



