############################################################
# Author : Frank Soldano                                   #
#                                                          #
# This code modify RE of the catalogue create by           #           
# egg-2skyamaker for create galaxies 2 time lower          #
# The code can also create image                           #
# associated to this catalog in ADU/s                      #
# (with egg-postskymake).                                  # 
# You can create only cat or cat and image.                #
#                                                          # 
# Before using this code, create your universe whith       #
# egg (egg-gencat) and compute cataloge whith              #
# egg-2skymaker. For create image, you need to             #
# have installed EGG & Skymaker software.                  #
#                                                          #
############################################################

# ----- necessary import -----
from astropy.io import ascii
import numpy as np
import os

# ----- ask for original files to read -----
files = input("cat file : ")
conf = input("conf file : ")

# ----- ask for img creation or only cat -----
todo = int(input("create only cat (1) or create cat + image (2) ? "))

# ----- read files -----
data = ascii.read(files)
data_conf = ascii.read(conf)

# ----- modification radius -----
data["bulge_radius"]=data["bulge_radius"]/2
data["disk_radius"]=data["disk_radius"]/2

# ----- create RE/2 sky -----
ascii.write(data, files[:-4]+"_reL2.cat", overwrite=True)

data_conf["SKY"][-3] = data_conf["SKY"][-3][:-5]+"_reL2.fits"
ascii.write(data_conf, conf[:-5]+"_reL2.conf", overwrite=True)

if todo == 2 :
    os.system("sky {} -c {}".format(files[:-4]+"_reL2.cat", conf[:-5]+"_reL2.conf"))
    os.system("egg-postskymaker conf={}".format(conf[:-5]+"_reL2.conf"))

print("Sky whith galaxies smaller by a factor 2 created")