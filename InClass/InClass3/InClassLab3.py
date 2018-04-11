#Marina Dunn, ASTR400B, Spring 2018, Dr. Gurtina Besla
#In Class Lab 3
# Created Feb 6 2018, Last updated 4/10/18

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
#get_ipython().magic('matplotlib inline')


# LOAD DATA
#****************
 
# Comes from   http://stellar.dartmouth.edu/models/isolf_new.html

# files have been modified from download.  ( M/Mo --> M;   Log L/Lo --> L)
# removed #'s from all lines except column heading

# NOTE SETTINGS USED:  Y = 0.245 default   [Fe/H] = -2.0  alpha/Fe = -0.2
# These could all be changed and it would generate a different isochrone


# Filename for data with Isochrone fit for 1-13 Gyr
filename1="Isochrone1.txt"
filename2="Isochrone2.txt"
filename3="Isochrone3.txt"
filename4="Isochrone4.txt"
filename5="Isochrone5.txt"
filename6="Isochrone6.txt"
filename7="Isochrone7.txt"
filename8="Isochrone8.txt"
filename9="Isochrone9.txt"
filename10="Isochrone10.txt"
filename11="Isochrone11.txt"
filename12="Isochrone12.txt"
filename13="Isochrone13.txt"


# READ IN DATA
# "dtype=None" means line is split using white spaces
# "skip_header=8"  skipping the first 8 lines 
# the flag "names=True" creates arrays to store the date
#       with the column headers given in line 8 

# Read in data for 1-13 Gyr
data1 = np.genfromtxt(filename1,dtype=None,names=True,skip_header=8)
data2 = np.genfromtxt(filename2,dtype=None,names=True,skip_header=8)
data3 = np.genfromtxt(filename3,dtype=None,names=True,skip_header=8)
data4 = np.genfromtxt(filename4,dtype=None,names=True,skip_header=8)
data5 = np.genfromtxt(filename5,dtype=None,names=True,skip_header=8)
data6 = np.genfromtxt(filename6,dtype=None,names=True,skip_header=8)
data7 = np.genfromtxt(filename7,dtype=None,names=True,skip_header=8)
data8 = np.genfromtxt(filename8,dtype=None,names=True,skip_header=8)
data9 = np.genfromtxt(filename9,dtype=None,names=True,skip_header=8)
data10 = np.genfromtxt(filename10,dtype=None,names=True,skip_header=8)
data11 = np.genfromtxt(filename11,dtype=None,names=True,skip_header=8)
data12 = np.genfromtxt(filename12,dtype=None,names=True,skip_header=8)
data13 = np.genfromtxt(filename13,dtype=None,names=True,skip_header=8)

# Plot Isochrones 
# For Carina
################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# PLOT ISOCHRONES

### EDIT THIS ! 
plt.plot(data1['B']-data1['R'], data1['R'], color='red', linewidth=5, label='1 Gyr')
plt.plot(data2['B']-data2['R'], data2['R'], color='orange', linewidth=5, label='2 Gyr')
plt.plot(data3['B']-data3['R'], data3['R'], color='yellow', linewidth=5, label='3 Gyr')
plt.plot(data4['B']-data4['R'], data4['R'], color='green', linewidth=5, label='4 Gyr')
plt.plot(data5['B']-data5['R'], data5['R'], color='blue', linewidth=5, label='5 Gyr')
plt.plot(data6['B']-data6['R'], data6['R'], color='purple', linewidth=5, label='6 Gyr')
plt.plot(data7['B']-data7['R'], data7['R'], color='black', linewidth=5, label='7 Gyr')
plt.plot(data8['B']-data8['R'], data8['R'], color='brown', linewidth=5, label='8 Gyr')
plt.plot(data9['B']-data9['R'], data9['R'], color='pink', linewidth=5, label='9 Gyr')
plt.plot(data10['B']-data10['R'], data10['R'], color='grey', linewidth=5, label='10 Gyr')
plt.plot(data11['B']-data11['R'], data11['R'], color='teal', linewidth=5, label='11 Gyr')
plt.plot(data12['B']-data12['R'], data12['R'], color='gold', linewidth=5, label='12 Gyr')
plt.plot(data13['B']-data13['R'], data13['R'], color='silver', linewidth=5, label='13 Gyr')

# Add axis labels
plt.xlabel('B-R', fontsize=22)
plt.ylabel(r'R', fontsize=22)

#set axis limits
plt.xlim(-0.5,2)
plt.ylim(5,-2.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.6, 0.15, 'Isochrone Carina', fontsize=22)


# Save to a file
ax.set_rasterized(True)
plt.savefig('IsochroneLabCarina.eps', rasterized=True, dpi=350)





