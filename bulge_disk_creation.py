############################################################
# Author : Frank Soldano                                   #
#                                                          #
# This code modify bulge to total ratio of the catalog     #
# create by egg-2skyamaker. The code can also create       #
# image associated to this catalog in ADU/s                #
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

# ----- modification bulge to total ratio -----
bt_ratio = data["bt"]
n = len(bt_ratio)

bt_bulge = np.full(n,1)
bt_disk = np.full(n,0)

# ----- create bulge sky -----
data["bt"] = bt_bulge
ascii.write(data, files[:-4]+"_bulge.cat", overwrite=True)

data_conf["SKY"][-3] = data_conf["SKY"][-3][:-5]+"_bulge.fits"
ascii.write(data_conf, conf[:-5]+"_bulge.conf", overwrite=True)

if todo == 2 :
    os.system("sky {} -c {}".format(files[:-4]+"_bulge.cat", conf[:-5]+"_bulge.conf"))
    os.system("egg-postskymaker conf={}".format(conf[:-5]+"_bulge.conf"))
    print("bulge sky created")

# ----- create disk sky -----
data["bt"] = bt_disk
ascii.write(data, files[:-4]+"_disk.cat", overwrite=True)

data_conf["SKY"][-3] = data_conf["SKY"][-3][:-10]+"disk.fits" #-10 for rmv bulge added previouscly
ascii.write(data_conf, conf[:-5]+"_disk.conf", overwrite=True)

if todo == 2 :
    os.system("sky {} -c {}".format(files[:-4]+"_disk.cat", conf[:-5]+"_disk.conf"))
    os.system("egg-postskymaker conf={}".format(conf[:-5]+"_disk.conf"))
    print("disk sky created")