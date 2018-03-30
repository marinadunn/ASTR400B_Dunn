# Created by Marina Dunn on 3/22/18.
# ASTR400B Spring 2018, Prof. Gurtina Besla
# Last edited on 3/29/18
# First, import relevant modules
import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterofMass import CenterOfMass

#Creating a class that will have a set of functions designed to calculate M33's acceleration felt by M31,
#and integrate over its current position and velocity as time moves forward.

###Step 1
class M33AnalyticOrbit:
    #Here, we will initialize the class, where input is the name of the file where the integrated orbit is stored,
    def __init__(self, filename):
        #Define G, acceleration, as a global variable, but not the units
        self.G = 4.498768e-6 #in kpc^3/Msun/Gyr
            
        #This stores the COM positions and velocities of M33 relative to M31
        #Using 'self' refers to the quantities common to an object; the values are stored so that data does not need to
        #be read in each time
        self.filename = filename
        #Using final correct (rounded) values from HW4 at Snapshot 0
        self.x = -99
        self.y = -120
        self.z = -128
        self.vx = -28
        self.vy = 174
        self.vz = 93
        #Defining the mass of M31's Bulge, Disk, and Halo with data from HW3
        #disk scale length
        self.rd = 5. # u.kpc
        self.Mdisk = 0.12e12 # u.Msun #taken from HW3
        #bulge scale length
        self.rbulge = 1. #kpc
        self.Mbulge = 0.019e12 # u.Msun #taken from HW3
        #halo scale length, use Hernquist scale length 'a' from HW5
        self.rhalo = 62
        self.Mhalo = 1.921e12 # u.Msun #taken from HW3
        
        return None
        
###Step 2
                 
#################################
#Defining the Acceleration Terms#
#################################
# This will be used for the Halo, Disk, and Bulge components for M31. We want to find the gravitational acceleration
#from each component in the galaxy

###Halo and Bulge Acceleration
#Defining a function that takes in the scale length (Ra), total halo or bulge mass (M), self, and the position
#coordinates (x,y,z) and a dummy variable used to define which accerelation component we want, and outputs the
#acceleration from the Hernquist potential in the direction of said variable
#Function only used for Halo and Bulge of M31
    def HernquistAccel(self, M, Ra, x, y, z, n):
        if n=='x':
            k = x
        if n=='y':
            k = y
        if n=='z':
            k = z

        r = np.sqrt(x**2+y**2+z**2) #position vector
        a_k1 = -((self.G*M)/(r*(Ra+r)**2))*k #Hernquist acceleration formula
        return a_k1

###Disk Acceleration
#We will now define a function for the approximation for exponential disk profile when very far away from the
#disk. It is formally known as the Miyamoto-Nagai 1975 profile. It has the following potential:

# phi(r) = (-self.G*self.Mdisk)/(np.sqrt(R**2+(Rd+np.sqrt(z**2+Zd**2)**2)))

#The function will take in self, M, the disk scale length (rd), a dummy variable, and the position coordinates (x,y,z)
#It returns the acceleration from the Miyamoto-Nagai 1975 potential in the direction of the dummy variable
    def MiyamotoNagaiAccel(self,M,rd,x,y,z,n):
    #Its acceleration for the x or y component is the following
        #Zd is the disk scale height
        Zd = self.rd/5.
        R = np.sqrt(x**2+y**2)
        B = self.rd+(np.sqrt(z**2+Zd**2)) #Zd is the disk scale height
        if n=='x':
            a_k2 = (-(self.G*M)/((R**2+B**2)**1.5))*x
        if n=='y':
            a_k2 = (-(self.G*M)/((R**2+B**2)**1.5))*y
        if n=='z':
            a_k2 = (-(self.G*M*B)/((R**2+B**2)**1.5)*np.sqrt(z**2+Zd**2))*z
        return a_k2
###M31Acceleration
#Defining a function that takes in self, the position components (x,y,z), and a dummy variable, and outputs the sum of all
#the acceleration terms from each galactic component in the direction of dummy variable
    def M31Accel(self,x,y,z,n):
        a_bulge = self.HernquistAccel(self.Mbulge, self.rbulge, x, y, z, n)
        a_halo = self.HernquistAccel(self.Mhalo, self.rhalo, x, y, z, n)
        a_disk = self.MiyamotoNagaiAccel(self.Mdisk, self.rd, x, y, z, n)
        accel_sum = a_bulge + a_halo + a_disk
        return accel_sum

###Step 3

#####################
#Build An Integrator#
#####################

#Want to solve M33's orbit by integrating equation of motion as time progresses. We will do this by defining a function
#called 'LeapFrog' that will take in the time interval for integration (delta_t), the M33 COM starting position vector
#(x,y,z), and the starting velocity vector (vx,vy,vz). These were initialized in part 1. We will be treating M33 as a
#point mass in order to do this integration.

#We will update the positions and velocities using standard kinematic equations
    def LeapFrog(self, delta_t, x, y, z, vx, vy, vz):
        #First, predict M33's COM 3D position at the middle of the timestep (delta_t) using the current COM velocity
        #and position for each component of the position (x,y,z) using the corresponding velocity component.

        x_half = x+(vx*(delta_t/2.))
        y_half = y+(vy*(delta_t/2.))
        z_half = z+(vz*(delta_t/2.))

        #Second, the COM position and velocity both go through another time step using the acceleration at another
        #timestep

        M31Accel_X = self.M31Accel(x_half, y_half, z_half, 'x')
        M31Accel_Y = self.M31Accel(x_half, y_half, z_half, 'y')
        M31Accel_Z = self.M31Accel(x_half, y_half, z_half, 'z')

        #delta_t can technically be positive or negative so the simulation can run both forward or backward in time.
        #This is because LeapFrog integrators are 'symplectic.' For our purposes, we want it to be positive because
        #we are calculating future orbits.
        vx_full = vx+M31Accel_X*delta_t
        vy_full = vy+M31Accel_Y*delta_t
        vz_full = vz*M31Accel_Z*delta_t

        x_full = x+(vx+vx_full)/2.*delta_t
        y_full = y+(vy+vy_full)/2.*delta_t
        z_full = z+(vz+vz_full)/2.*delta_t

        return x_full, y_full, z_full, vx_full, vy_full, vz_full

###Step 4

#####################
#Integrate the Orbit#
#####################

#We will now define a function that takes in slef, the start time (t0), the time interval (delta_t), and the final
#time (t_max), and will loop over the LeapFrog integrator in order to solve the equations of motion and calculate
#the orbit of M33 10 Gyr in the future.
    def OrbitIntegrator(self, t0, delta_t, t_max):
        #Defining the initial COM position and velocity components for M33 relative to M31 again. Variables will update
        #as we integrate.
        x = self.x
        y = self.y
        z = self.z
        vx = self.vx
        vy = self.vy
        vz = self.vz

        #Now define variable t and integration start at t0, then continue looping with a while loop over LeapFrog until
        #reaching t_max = 10 Gyr; update time, positions, and velocities while doing so. Return results in an array,
        #initialized outside while loop.
        Orbit = np.zeros(shape=(int((t_max-t0)/delta_t)+2,7))
        #store initial values in array
        t0 = Orbit[0][0]
        x = Orbit[0][1]
        y = Orbit[0][2]
        z = Orbit[0][3]
        vx = Orbit[0][4]
        vy = Orbit[0][5]
        vz = Orbit[0][6]
    
        #loop from t0 to t_max
        t = t0
        i = 1 #counter
        while t <= t_max:
                integrator = self.LeapFrog(delta_t, x, y, z, vx, vy, vz) #integrates over 1 step

                #Optional print statement
                print i
                #store position and velocity vectors in array
                Orbit[i,0] = t + delta_t
                Orbit[i,1] = integrator[0]
                Orbit[i,2] = integrator[1]
                Orbit[i,3] = integrator[2]
                Orbit[i,4] = integrator[3]
                Orbit[i,5] = integrator[4]
                Orbit[i,6] = integrator[5]

                #reinitialize variable for new orbit
                i = i + 1
                t = t + delta_t
                x = integrator[0]
                y = integrator[1]
                z = integrator[2]
                vx = integrator[3]
                vy = integrator[4]
                vz = integrator[5]

        #save output in text file
        np.savetxt('M33Analytical_orbit.txt', Orbit, header='t x y z vx vy vz', comments='#', fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])
        return i

AnalyticalObject = M33AnalyticOrbit('M33Analytical_orbit.txt')
M33Analytic_file = AnalyticalObject.OrbitIntegrator(0,0.5,10)


print "Question 1: I cannot seem to get the analyrical prediction to plot in addition to the simulation. I would expect the separation to sharply increase around 6 Gyr and continue for the analytical solution. I would also expect M33's velocity to sharply increase around 2 Gttr and 6 Gyr as the galaxy make the passbys."
print "Question 2: Analytic prediction assumes M33 is a point mass, which is not true; it is a large collection of particles that are also pulling gravitationally on each other, so taking this into account would be extremely difficult, but more accurate."
print "Question 3: Milky Way's gravitational pull is relatively unsignificant until the galaxies start to become very near each other and make that first passby. In that case, it's not just M31's mass pulling on M33, but the combined mass of M31 and Milky Way, creating an even stronger attraction."
