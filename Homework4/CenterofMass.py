# Last Edited: May 28, 2020
#Marina Dunn, Homework 4, ASTR400B
# Center of Mass Position and Velocity

# First, import relevant modules
import numpy as np
import astropy.units as u
from ReadFile import Read

#A class creates a new type of object, and has attributes that characterize it. CenterOfMass is the class, and
#the functions we will create describe it
class CenterOfMass:
    #Here, we will initialize the class so that for each object, the data from the simulation file will be stored,
    #depending on particle type
    def __init__(self, filename, ptype):
        #This reads in the filename and particle type
        self.time, self.total, self.data = Read(filename)
            
        #This creates an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        #This stores the mass, positions, velocities of only the particles of the given type
        #Using 'self' refers to the quantities common to an object; the values are stored so that data does not
        #need to be read in each time
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
      
    def total_mass(self,m):
          #Note: you can add other keyword arguments into the function, but 'self' must be first
          return np.sum(self.m)
    
    #Defining a function that will calculate the 3D coordinates for the center of mass (based on position and
    #velocity) for any galaxy
    #inputs will be random variables that will change later
    #Calculate the center of mass position
    def COMdefine(self,x,y,z,m):
        Xcom = np.sum(x*m)/np.sum(m)
        Ycom = np.sum(y*m)/np.sum(m)
        Zcom = np.sum(z*m)/np.sum(m)

        return Xcom, Ycom, Zcom
    
    #Defining a function that will call COMdefine and take in the position and velocity vectors for a
    #given particle type, and return the center of mass for a galaxy
    def COM_P(self, delta):
        #Call COMdefine to compute a first estimate for COM position, and use "self" again since it has been defined inside of the class
        Xcom, Ycom, Zcom = self.COMdefine(self.x,self.y,self.z,self.m)
        
        #We will compute a first estimate for the magnitude of the Center of Mass vector for a particular ptype
        RCOM = np.sqrt((Xcom)**2.+(Ycom)**2.+(Zcom)**2.)
        
        #Want to transfer to the Center of Mass reference frame, so we must take a difference between the
        #first guess and the position vectors
        XNEW = self.x - Xcom
        YNEW = self.y - Ycom
        ZNEW = self.z - Zcom
        
        #We create an array that stores the magnitude of the new position vectors in the COM frame
        RNEW = np.sqrt((XNEW)**2.+(YNEW)**2.+(ZNEW)**2.)
        
        #We want to find the max 3D separation between the COM coordinates and the reference frame, then divide
        #that by half. We will continue to do this in order to check if the position is converging.
        #print RNEW
        RMAX = np.max(RNEW)/2.
        
        #Setting an initial value for the difference
        diff = 100.0
        #Creating iteration to continue if the separation in COM positions is bigger than delta
        while(diff > delta):
                #Defining particles within a certain radius
                index2 = np.where(RNEW<RMAX)

                #Difference in particle and new COM positions
                X2 = self.x[index2]
                Y2 = self.y[index2]
                Z2 = self.z[index2]
                m2 = self.m[index2]
                
                #Calculates COM position for particles inside smaller radius
                Xcom2, Ycom2, Zcom2 = self.COMdefine(X2,Y2,Z2,m2)
                
                #Magnitude of distance for new COM
                RCOM2 = np.sqrt(Xcom2**2 + Ycom2**2 + Zcom2**2)
                
                #difference in old and new COM for each component in position vector
                diff = np.abs(RCOM - RCOM2)
              
                #Max 3D separation for new COM coordinates, so a smaller radius
                RMAX = RMAX/2.0

                #Changes reference frame to new COM
                XNEW = self.x - Xcom2
                YNEW = self.y - Ycom2
                ZNEW = self.z - Zcom2
                RNEW = np.sqrt(XNEW**2 + YNEW**2 + ZNEW**2)
    
                #Setting COM positions to refined values
                Xcom = Xcom2
                Ycom = Ycom2
                Zcom = Zcom2
                RCOM = RCOM2
        
                #Creating vector that stores rounded COM positions
                COMP = [Xcom, Ycom, Zcom]
                    
        return np.round(COMP,2)*u.kpc #returns COM position vector

    def COM_V(self, Xcom, Ycom, Zcom):
        #Defining particles within a certain radius
        RVMAX = 15.0*u.kpc
        
        #Separation to COM frame
        VX = self.x[:]*u.kpc - Xcom
        VY = self.y[:]*u.kpc - Ycom
        VZ = self.z[:]*u.kpc - Zcom
        RV = np.sqrt(VX**2 + VY**2 +VZ**2)
        
        #Want velocities for particles within a radius of 15 kpc from COM position
        index = np.where(RV < RVMAX)
                
        #Calculates COM velocity for particles inside RVMAX, and use "self" again since it has been defined inside
        #of the class
        VXcomNEW = self.vx[index]
        VYcomNEW = self.vy[index]
        VZcomNEW = self.vz[index]
        mnew = self.m[index]

        #Computing the COM velocity
        VXcom, VYcom, VZcom = self.COMdefine(VXcomNEW,VYcomNEW,VZcomNEW,mnew)

        #Creating vector that stores rounded COM velocities
        COMV = [VXcom, VYcom, VZcom]
        return np.round(COMV,2)*u.km/u.s
#####
###Test if the code works
#####

#To find the COM position and velocity for disk particles in MW, M31, and M33, use ptype = 2
#For MW, filename = MW_0000.txt, for M31, filename = M31_000.txt, for M33, filename = M33_000.txt

#Choosing a tolerance for COM
delta = 0.3

print("Answer to Question 1:")
MWCOM = CenterOfMass("MW_000.txt", 2.)
# Calculate quantities for MW data
MW_COMP = MWCOM.COM_P(0.3)
#Print COM position vector for Milky Way
print("MW COM Position Vector:"), MW_COMP
MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])
#Print COM velocity vector for Milky Way
print("MW COM Velocity Vector:"), MW_COMV

M31COM = CenterOfMass("M31_000.txt", 2.)
# Calculate quantities for M31 data
M31_COMP = M31COM.COM_P(0.3)
#Print COM position vector for M31
print("M31 COM Position Vector:"), M31_COMP
M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])
#Print COM velocity vector for M33
print("M31 COM Velocity Vector:"), M31_COMV

M33COM = CenterOfMass("M33_000.txt", 2.)
# Calculate quantities for M33 data
M33_COMP = M33COM.COM_P(0.3)
#Print COM position vector for M33
print("M33 COM Position Vector:"), M33_COMP
M33_COMV = M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2])
#Print COM velocity vector for M33
print("M33 COM Velocity Vector:"), M33_COMV


print ("Answer to Question 2:")
#Print the magnitude of current separation between MW and M31
MW_M31_sep = np.sqrt((M31_COMP[0] - MW_COMP[0])**2 + (M31_COMP[1] - MW_COMP[1])**2 + (M31_COMP[2] - MW_COMP[2])**2)
print("Magnitude of current separation between MW and M31:"), np.round(MW_M31_sep)
MW_M31_vel = np.sqrt((M31_COMV[0] - MW_COMV[0])**2 + (M31_COMV[1] - MW_COMV[1])**2 + (M31_COMV[2] - MW_COMV[2])**2)
print("Magnitude of velocity between MW and M31:"), np.round(MW_M31_vel)

print ("Answer to Question 3:")
#Print the magnitude of current separation between M31 and M33
M33_M31_sep = np.sqrt((M33_COMP[0] - M31_COMP[0])**2 + (M33_COMP[1] - M31_COMP[1])**2 + (M33_COMP[2] - M31_COMP[2])**2)
print("Magnitude of current separation between M31 and M33:"), np.round(M33_M31_sep)
M33_M31_vel = np.sqrt((M33_COMV[0] - M31_COMV[0])**2 + (M33_COMV[1] - M31_COMV[1])**2 + (M33_COMV[2] - M31_COMV[2])**2)
print("Magnitude of velocity between M31 and M33:"), np.round(M33_M31_vel)

print("Answer to Question 4:")
print("As the galaxies become closer, the total separation decreases and velocity increases, thus changing the COM, especially after they start interacting and the stars all change radii, so we must use an iterative process as time evolves to update the COM calculation.")
