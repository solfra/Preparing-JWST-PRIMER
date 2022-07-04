############################################################
# Author : Frank Soldano                                   #
#                                                          #
# This code convert image create by Skymaker               #
# + egg-postskymaker in ADU/s to an image in ADU.          #
#                                                          #
# For this, you need exposure time                         #
#                                                          #
############################################################

# ----- necessary import -----
from astropy.io import fits
import numpy as np

# ----- ask for original files to read -----
file = input("file name : ")
a = fits.open(file)
d=a[0].data
h=a[0].header

# ----- ask for exposure time -----
expos=float(input("exposure time : ")) #need to be float for SExtractor

# ----- modify image -----
h['EXPTIME']=expos #change the exposure time
d2=np.round(d*expos) #multiply data by exposure time and round the value
h['BUNIT']='DN' #change DN/s to DN

# ----- write image -----
fits.writeto(file[:-5]+"_adu.fits",d2,h,overwrite=True)