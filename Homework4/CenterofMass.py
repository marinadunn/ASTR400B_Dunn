# Edited original file on February 4, 2018
#Marina Dunn, Homework 4, ASTR400B
# Center of Mass Position and Velocity

# First, import relevant modules
import numpy as np
import astropy.units as u
from ReadFile import Read

#A class creates a new type of object, and has attributes that characterize it. CenterOfMass is the class, and the functions
#we will create will describe it
class CenterOfMass:
    #Here, we will initialize the class so that for each object, the data from the simulation file will be stored, depending on particle type
    def __init__(self, filename, ptype):
        #This reads in the filename and particle type
        self.time, self.total, self.data = Read(filename)
            
        #This creates an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        #This stores the mass, positions, velocities of only the particles of the given type
        #Using 'self' refers to the quantities common to an object; the values are stored so that data does not need to be read in each time
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
      
    def total_mass(m):
          #Note: you can add other keyword arguments into the function, but 'self' must be first
          return np.sum(m)
    
    #Defining a function that will calculate the 3D coordinates for the center of mass (based on position and velocity) for any galaxy
    #inputs will be random variables that will change later
    #Calculate the center of mass position
    def COMdefine(self,x,y,z,m):
        Xcom = np.sum(x*m)/self.total_mass(m)
        Ycom = np.sum(y*m)/self.total_mass(m)
        Zcom = np.sum(z*m)/self.total_mass(m)

        return Xcom, Ycom, Zcom
    
    #Defining a function that will call COMdefine and take in the position and velocity vectors for a
    #given particle type, and return the center of mass for a galaxy
    def COM_P(self, delta):
        #Call COMdefine to compute a first estimate for COM position, and use "self" again since it has been defined inside of the class
        Xcom, Ycom, Zcom = self.COMdefine(self.x,self.y,self.z,self.m)
        
        #We will compute a first estimate for the magnitude of the Center of Mass vector for a particular ptype
        RCOM = np.sqrt((Xcom)**2+(Ycom)**2+(Zcom)**2)
        
        #Want to transfer to the Center of Mass reference frame, so we must take a difference between the
        #first guess and the position vectors
        XNEW = self.x - Xcom
        YNEW = self.y - Ycom
        ZNEW = self.z - Zcom
        
        #We create an array that stores the magnitude of the new position vectors in the COM frame
        RNEW = np.sqrt((XNEW)**2+(YNEW)**2+(ZNEW)**2)
        
        #We want to find the max 3D separation between the COM coordinates and the reference frame, then divide
        #that by half. We will continue to do this in order to check if the position is converging.
        RMAX = np.amax(RNEW)/2.0
        
        #Setting an initial value for the difference
        diff = 100
        #Creating iteration to continue if the separation in COM positions is bigger than delta
        while(diff > delta):
                #Defining particles within a certain radius
                index = np.where(RNEW<RMAX)
                
                #Calculates COM position for particles inside RMAX
                XNEW, YNEW, ZNEW = self.COMdefine(self.x[index],self.y[index],self.z[index],self.m[index])
                
                #Difference in particle and new COM positions
                XNEW = self.x - Xcom
                YNEW = self.y - Ycom
                ZNEW = self.z - Zcom
                RNEW = np.sqrt(XNEW**2+YNEW**2+ZNEW**2)
                
                #Max 3D separation for new COM coordinates, so a smaller radius
                RMAX = np.amax(RNEW)/2.0
                #Magnitude of distance for new COM
                RCOM2 = np.sqrt(Xcom**2+Ycom**2+Zcom**2)
                
                #difference in old and new COM for each component in position vector
                diff = np.abs(RCOM-RCOM2)
                RCOM=RCOM2
              

        return XNEW, YNEW, ZNEW

    def COM_V(self, delta):
        #Defining particles within a certain radius
        i, j, k = self.COM_P(0.3)
        
        #Separation to COM frame
        XNEW = self.x - i
        YNEW = self.y - j
        ZNEW = self.z - k
        r = np.sqrt(XNEW**2+YNEW**2+ZNEW**2)
        
        #Want velocities for particles within a radius of 15 kpc from COM position
        index = np.where(np.abs(r)<15)
                
        #Calculates COM velocity for particles inside RMAX, and use "self" again since it has been defined inside of the class
        VXNEW, VYNEW, VZNEW = self.COMdefine(self.vx[index],self.vy[index],self.vz[index],self.m[index])

        return VXNEW, VYNEW, VZNEW
#####
###Test if the code works
#####

#To find the COM position and velocity for disk particles in MW, M31, and M33, use ptype = 2
#For MW, filename = MW_0000.txt, for M31, filename = M31_000.txt, for M33, filename = M33_000.txt

#Choosing a tolerance for COM
delta = 0.3

print("Answer to Question 1:")
MWCOM = CenterOfMass("MW_000.txt", 2)
# Calculate quantities for MW data
MW_Xcom, MW_Ycom, MW_Zcom = MWCOM.COM_P(delta)
#Print COM position vector components for Milky Way
print("MW COM Position Vector Components:"), MW_Xcom, MW_Ycom, MW_Zcom
MW_VXcom, MW_VYcom, MW_VZcom = MWCOM.COM_V(delta)
#Print COM velocity vector for Milky Way
print("MW COM Velocity Vector Components:"),MW_VXcom, MW_VYcom, MW_VZcom

M31COM = CenterOfMass("M31_000.txt", 2)
# Calculate quantities for M31 data
M31_Xcom, M31_Ycom, M31_Zcom = M31COM.COM_P(delta)
#Print COM position vector for M31
print("M31 COM Position Vector Components:"), M31_Xcom, M31_Ycom, M31_Zcom
M31_VXcom, M31_VYcom, M31_VZcom = M31COM.COM_V(delta)
#Print COM velocity vector for M33
print("M31 COM Velocity Vector Components:"),M31_VXcom, M31_VYcom, M31_VZcom

M33COM = CenterOfMass("M33_000.txt", 2)
# Calculate quantities for M33 data
M33_Xcom, M33_Ycom, M33_Zcom = M33COM.COM_P(delta)
#Print COM position vector for M33
print("M33 COM Position Vector Components:"), M33_Xcom, M33_Ycom, M33_Zcom
M33_VXcom, M33_VYcom, M33_VZcom = M33COM.COM_V(delta)
#Print COM velocity vector for M33
print("M33 COM Velocity Vector Components:"),M33_VXcom, M33_VYcom, M33_VZcom

print("The components above do not seem quite correct, but I cannot find a place in the code that changes the output correctly")

print ("Answer to Question 2:")
#Print the magnitude of current separation between MW and M31
print("Magnitude of current separation between MW and M31:"), np.sqrt((M31_Xcom - MW_Xcom)**2+(M31_Ycom - MW_Ycom)**2+(M31_Zcom - MW_Zcom)**2)
print("Magnitude of velocity between MW and M31:"), np.sqrt((M31_VXcom - MW_VXcom)**2+(M31_VYcom - MW_VYcom)**2+(M31_VZcom - MW_VZcom)**2)

print ("Answer to Question 3:")
#Print the magnitude of current separation between M31 and M33
print("Magnitude of current separation between M31 and M33:"), np.sqrt((M31_Xcom - M33_Xcom)**2+(M31_Ycom - M33_Ycom)**2+(M31_Zcom - M33_Zcom)**2)
print("Magnitude of velocity between M31 and M33:"), np.sqrt((M31_VXcom - M33_VXcom)**2+(M31_VYcom - M33_VYcom)**2+(M31_VZcom - M33_VZcom)**2)

print("Answer to Question 4:")
print("As the galaxies become closer, the total separation decreases and velocity increases, thus changing the COM, especially after they start interacting, so we must use an iterative process as time evolves to update the COM calculation.")
