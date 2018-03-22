#  Created by Marina Dunn on 2/17/18.
# ASTR400B Homework 6, spring 2018
#Note on creating symbolic link: I ssh'ed into Nimoy and used ln -s ../astr400b/VLowRes ./VLowRes to copy the files into my home directory.
#I then logged out and used  scp -r mdunn@nimoy.as.arizona.edu:VLowRes/ . to copy the files to my personal laptop

# First, import relevant modules
import numpy as np
import astropy.units as u
from astropy.constants import G

import matplotlib
import matplotlib.pyplot as plt
#Read in CenterOfMass, which has been updated to
from CenterofMass import *
#Read in ReadFile, which has been updated not include the original factor of 10 Myr for time
from ReadFile import Read


#We will create a function that takes in the galaxy name, snap number,
#For this code, we want n=5, and orbits up to Snapshot 800, which corresponds to ~12 Gyr
def OrbitCOM(galaxy, start, end, n):

    fileout = 'Orbit_%s.txt'
    
    #We will create an array that will store the positions nd velocities of the COM of the galaxy at each snapshot, so there are 7 columns
    Orbit = np.zeros((int(end/n)+1,7))

    #Set a tolerance and VolDec for each CenterOfMass object (MW, M33, M31)
    #For MW and M31:
    delta = 0.3
    VolDec = 2.0

    #For M33:
    if (galaxy == "M33"):
    delta = 0.5
    VolDec = 4.

    #Creating a for loop that will go from 'start' to 'end + 1' in intervals of n
    for i in np.arange(start, end+n, n):
        #Next, for a given filename, we only want the first characters that specify which galaxy we are talking about
        #We will add a string of the filenumber to the value "000"
        ilbl = '000' + str(i)
        #Then, we want the part of the filename that specifies the snapnumber
        #remove everything except the last 3 digits
        ilbl = ilbl[-3:]
        filename = '%s_'%(galaxy) + ilbl + '.txt' #stores filename as global property
        #i.e. if I wanted "MW_010.txt", I would input "MW" and "10"
        #optional print statement to see if filename is printed correctly
        #print filename

        #Creating a COM object for disk particles
        COM = CenterOfMass(filename, 2) #use ptype = 2 for disk component

        #Need to turn quantities back into floats by dividing by units in order to use them in an array
        Orbit[int(i/n), 0] = float(COM.time/u.Myr/1000.) #stores time in Gyr

        #Store the COM position in the Orbit array and divide out the units
        XcomNEW, YcomNEW, ZcomNEW = COM.COM_P(delta, VolDec)
        Orbit[int(i/n), 1] = XcomNEW
        Orbit[int(i/n), 2] = YcomNEW
        Orbit[int(i/n), 3] = ZcomNEW

        #Store the COM velocity in the Orbit array and divide out the units
        VXcomNEW, VYcomNEW, VZcomNEW = COM.COM_V(delta, VolDec)
        Orbit[int(i/n), 4] = VXcomNEW
        Orbit[int(i/n), 5] = VYcomNEW
        Orbit[int(i/n), 6] = VZcomNEW
            
            
        #Prints counter to know where code is
        print "Counter:", i


        #Save Orbit array to a file
        fileout = 'M31_orbit.txt' #Change to alternate name when running on another galaxy with different tolerance and VolDec
        np.savetxt(fileout, Orbit, header='t x y z vx vy vz', comments='#', fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])

    #return Orbit
#Returns the Orbit array and the text file specified above, which uses the array



####test code with small range of snapshots###




#First off, in order to compute 3D COM position and velocity vectors, we will need to go from SnapNumber 0 to 800 in n=5 intervals. Each time, we will
#make 3 files that store these COM properies for each galaxy
print "MW Orbit:", OrbitCOM("MW",0,800,5)
print "MW Orbit:", OrbitCOM("M31",0,800,5)
print "MW Orbit:", OrbitCOM("M33",0,800,5)


#Next, we will read in these newly-created COM files, and we will need to skip the lines beginning with '#'
MW_data = np.genfromtxt('MW_Orbit.txt',dtype=None,names=True)
M31_data = np.genfromtxt('M31_Orbit.txt',dtype=None,names=True)
M33_data = np.genfromtxt('M33_Orbit.txt',dtype=None,names=True)

#Now subtracting the 3D coordinates of the Milky Way from M31, and also M33 from M31, we will be able to calculate the magnitude of the separation
#between the galaxies, then plot that over time

#Separation of MW and M31
MW_M31_separation = np.sqrt((MW_data['x']-M31_data['x'])**2+(MW_data['y']-M31_data['y'])**2+(MW_data['z']-M31_data['z'])**2)
#Velocity of MW and M31
MW_M31_velocity = np.sqrt((MW_data['vx']-M31_data['vx'])**2+(MW_data['vy']-M31_data['vy'])**2+(MW_data['vz']-M31_data['vz'])**2)

#Separation of M33 and M31
M33_M31_separation = np.sqrt((M33_data['x']-M31_data['x'])**2+(M33_data['y']-M31_data['y'])**2+(M33_data['z']-M31_data['z'])**2)
#Velocity of M33 and M31
M33_M31_velocity = np.sqrt((M33_data['vx']-M31_data['vx'])**2+(M33_data['vy']-M31_data['vy'])**2+(M33_data['vz']-M31_data['vz'])**2)


##########
#Plotting#
##########



#Plot a zoom of MW and M31 separation
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(MW_data['t'], MW_M31_separation, color='red', linewidth=5, label='MW-M31')
#Plot a zoom of MW and M31 separation also
fig = plt.figure(figsize=(10,10))
plt.plot(M33_data['t'], M33_M31_separation, color='purple', linewidth=5, label='M33-M31')


 # Adding the axis labels and title
plt.xlabel('Time (Gyr)', fontsize=16)
plt.ylabel('Separation (kpc)', fontsize=16)
plt.title("Center of Mass Separation", fontsize=20)

#set axis limits
plt.ylim(0,50)
plt.xlim(6,6.5)

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
ax.set_rasterized(True)
plt.savefig('MW-M31_Separation.eps',rasterized=True, dpi=350)
plt.close()

#Finally, instead of distance magnitude plotted against time, we will find the magnitude of the relative velocity and plot that against time, once again
#for MW and 31, and M33 and M31. This is done by subtracting the velocity components, VX, VY, and VZ rather than X, Y, and Z.



#Separation of MW and M31
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(MW_data['t'], MW_M31_velocity, color='blue', linewidth=5, label='MW-M31')

#Separation of M33 and M31
fig = plt.figure(figsize=(10,10))
plt.plot(M33_data['t'], M33_M31_velocity, color='pink', linewidth=5, label='M33-M31')

# Adding the axis labels and title
plt.xlabel('Time (Gyr)', fontsize=16)
plt.ylabel('Separation (kpc)', fontsize=16)
plt.title("Center of Mass Velocity", fontsize=20)

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
plt.savefig('MW_Separation.eps')
plt.close()
# Adding the axis labels and title
plt.xlabel('Time (Gyr)', fontsize=16)
plt.ylabel('Relative Velocity (km/s)', fontsize=16)
plt.title("Center of Mass Relative Velocity Between MW and M31", fontsize=20)

#adjust tick label font size
label_size = 16
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
plt.legend(loc='best',fontsize='medium')

# Save to a file
plt.savefig('MW_Relative_Velocity.eps')
plt.close()




#Question 1: How many close encounters will the MW and M31 experience in the future?
print "Answer to Question 1:"
print "3. They will get close about 2 times before merging completely into one. \n"

#Question 2: How is the time evolution of the separation and relative velocity related?
print "Answer to Question 2:"
print "Relative Velocity maxima are when the galaxies are closest, so as separation decreases, relative velocity increases.\n"

#Question 3: When do M31 and MW merge (zoom in in plot)? What happens to M33's orbit when they do?
print "Answer to Question 3:"
print "They merge at around 6.3 Gyr, and M33 becomes the closest it ever will to M31, then swings out and back to orbit the M31-MW merger remnant, so it stops around 7 Gyr.\n"

#Question 4: What is the approximate decay rate of M33's orbit after 6 Gyr? If this rate is constant, how long will it take M33 to merge with the combined
#MW+M31 remnant if it is at a distance of 75 kpc?
print "Answer to Question 4 (Bonus):"
print "I would say maybe about 4.5 kpc/Gyr. \n"

