import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from tqdm import trange

def flag_results(data):
    """
    Flag Galfit results
    Flag from van der Wel + 2012

    input : catalog : data from the fits catalogue from Galapagos-c
    output : good / suspicious / bad / nan
    """
    flag_nan, = np.where((np.isnan(data['SERSIC_GALFIT'])) | (np.isnan(data['ELLIPTICITYERR_GALFIT'])) | (np.isnan(data['ARERR_GALFIT'])) 
                         |(np.isnan(data['SERSICERR_GALFIT'])) | (np.isnan(data['REERR_GALFIT'])) | (np.isnan(data['MAGERR_GALFIT'])))
    mask_nan=np.ones(len(data['FIT_DONE']),bool)
    mask_nan[flag_nan]=0
    flag_bad,=np.where((data['SERSIC_GALFIT']==0.2) | (data['SERSIC_GALFIT']==20.) | (data['REERR_GALFIT']==99999) | (data['RE_GALFIT']==400) ) 
    mask_bad=np.ones(len(data['FIT_DONE']),bool)
    mask_bad[flag_bad]=0
    flag_good,=np.where(data['FLAGS_GALFIT']=='GOOD')
    mask_good=np.ones(len(data['FIT_DONE']),bool)
    mask_good[flag_good]=0
    
    good,=np.where((mask_bad&mask_nan&np.logical_not(mask_good))==True)
    suspisious,=np.where((mask_bad&mask_nan&mask_good)==True)
    bad,=np.where((np.logical_not(mask_bad)&mask_nan)==True)
    non_exist,=np.where(np.logical_not(mask_nan))
    return good, suspisious, bad, non_exist

def flag_results_nF(data):
    """
    Flag Galfit results when sersic index is fixed in galfit
    Flag from van der Wel + 2012

    input : catalog : data from the fits catalogue from Galapagos-c
    output : good / suspicious / bad / nan
    """
    flag_nan, = np.where( (np.isnan(data['ELLIPTICITYERR_GALFIT'])) | (np.isnan(data['ARERR_GALFIT'])) 
                         | (np.isnan(data['REERR_GALFIT'])) | (np.isnan(data['MAGERR_GALFIT'])))
    mask_nan=np.ones(len(data['FIT_DONE']),bool)
    mask_nan[flag_nan]=0
    flag_bad,=np.where((data['SERSIC_GALFIT']==0.2) | (data['SERSIC_GALFIT']==20.) | (data['REERR_GALFIT']==99999) | (data['RE_GALFIT']==400) ) 
    mask_bad=np.ones(len(data['FIT_DONE']),bool)
    mask_bad[flag_bad]=0
    flag_good,=np.where(data['FLAGS_GALFIT']=='GOOD')
    mask_good=np.ones(len(data['FIT_DONE']),bool)
    mask_good[flag_good]=0
    
    good,=np.where((mask_bad&mask_nan&np.logical_not(mask_good))==True)
    suspisious,=np.where((mask_bad&mask_nan&mask_good)==True)
    bad,=np.where((np.logical_not(mask_bad)&mask_nan)==True)
    non_exist,=np.where(np.logical_not(mask_nan))
    return good, suspisious, bad, non_exist

def flag_results_mLim(data,magL):
    """
    Flag Galfit results
    Flag from van der Wel + 2012 and adding maglim 50 / 80

    input : catalog : data from the fits catalogue from Galapagos-c / maglim
    output : good / suspicious / bad / nan
    """
    flag_nan, = np.where((np.isnan(data['SERSIC_GALFIT'])) | (np.isnan(data['ELLIPTICITYERR_GALFIT'])) | (np.isnan(data['ARERR_GALFIT'])) 
                        |(np.isnan(data['SERSICERR_GALFIT'])) | (np.isnan(data['REERR_GALFIT'])) | (np.isnan(data['MAGERR_GALFIT'])))
    mask_nan=np.ones(len(data['FIT_DONE']),bool)
    mask_nan[flag_nan]=0
    flag_bad,=np.where((data['SERSIC_GALFIT']==0.2) | (data['SERSIC_GALFIT']==20.) | (data['REERR_GALFIT']==99999) | (data['RE_GALFIT']==400) | (data['MAG_BEST']>magL) )
    mask_bad=np.ones(len(data['FIT_DONE']),bool)
    mask_bad[flag_bad]=0
    flag_good,=np.where(data['FLAGS_GALFIT']=='GOOD')
    mask_good=np.ones(len(data['FIT_DONE']),bool)
    mask_good[flag_good]=0
    
    good,=np.where((mask_bad&mask_nan&np.logical_not(mask_good))==True)
    suspisious,=np.where((mask_bad&mask_nan&mask_good)==True)
    bad,=np.where((np.logical_not(mask_bad)&mask_nan)==True)
    non_exist,=np.where(np.logical_not(mask_nan))
    return good, suspisious, bad, non_exist

def combine_EGG_Galfit(data,ra_pix,dec_pix,mag_detect):
    """
    Create corespondence between EGG and Galfit catalogue

    input : 
    * data : data from the fits catalogue from Galapagos-c 
    * ra_pix / dec_pix : position (in pixel) of galaxies in egg catalogue 
    * mag_detect : np array - position in the EGG catalogue of galaxies whith a mag_lim under the detection limite

    output :
    * cores_out_in : array of conrespondance betwenn egg and galfit catalogue
    * ok : position in galfit gatalogue of galaxies withe corespondence
    * err_cores : position in galfit catalogue of galaxies with no corespondence
    """
    cores_out_in = []
    err_cores = []
    ok=[]
    ce = 0
    error=0
    for gal in trange(len(data['X_IMAGE'])):
        pos, = np.where( ( (data['X_IMAGE'][gal]-data['RE_GALFIT'][gal]) < ra_pix[mag_detect] ) & ((data['X_IMAGE'][gal]+data['RE_GALFIT'][gal]) > ra_pix[mag_detect]) & 
                  ((data['Y_IMAGE'][gal]-data['RE_GALFIT'][gal]) < dec_pix[mag_detect]) & ((data['Y_IMAGE'][gal]+data['RE_GALFIT'][gal]) > dec_pix[mag_detect]))
    
        if pos.shape[0]==1 :
            if int(pos[0]) not in cores_out_in : #check if corespondence fond is already assigned
                pos=np.array([pos[0]])
                ok.append(gal)
            elif int(pos[0]) in cores_out_in : 
                err_cores.append(gal)
                error+=1
                pos=np.array([0])
    
        elif pos.shape[0]>1 : #more than one galaxies on catalogue corespond
            ce+=1
            i=0
            while pos.shape[0]!=1 :
                pos, = np.where( ( (data['X_IMAGE'][gal]-i) < ra_pix[mag_detect] ) & ((data['X_IMAGE'][gal]+i) > ra_pix[mag_detect]) & 
                            ((data['Y_IMAGE'][gal]-i) < dec_pix[mag_detect]) & ((data['Y_IMAGE'][gal]+i) > dec_pix[mag_detect]))
                if pos.shape[0]==1 :
                    if int(pos[0]) not in cores_out_in : #check if corespondence fond is already assigned
                        pos=np.array([pos[0]])
                        ok.append(gal)
                    elif int(pos[0]) in cores_out_in : 
                        err_cores.append(gal)
                        error+=1
                        pos=np.array([0])
            
                if i>1*data['RE_GALFIT'][gal] : #security for infinite loop
                    error+=1
                    err_cores.append(gal)
                    pos=np.array([0])
                i+=0.1   
                
        elif pos.shape[0]==0 : #no corespondence
            ce+=1
            i=0
            while pos.shape[0]!=1 :
                pos, = np.where( ( (data['X_IMAGE'][gal]-i) < ra_pix[mag_detect] ) & ((data['X_IMAGE'][gal]+i) > ra_pix[mag_detect]) & 
                            ((data['Y_IMAGE'][gal]-i) < dec_pix[mag_detect]) & ((data['Y_IMAGE'][gal]+i) > dec_pix[mag_detect]))
                if pos.shape[0]==1 :
                    if int(pos[0]) not in cores_out_in :
                        pos=np.array([pos[0]])
                        ok.append(gal)
                    elif int(pos[0]) in cores_out_in : 
                        err_cores.append(gal)
                        error+=1
                        pos=np.array([0])
            
                if i>1*data['RE_GALFIT'][gal]: #security for infinite loop
                    error+=1
                    err_cores.append(gal)
                    pos=np.array([0])
                i+=0.1

                if (i>400) & (data['FIT_DONE'][gal]==0) : #security for infinite loop
                    error+=1
                    err_cores.append(gal)
                    pos=np.array([0])
    
        cores_out_in.append(int(pos[0]))
    
    cores_out_in=np.array(cores_out_in)
    ok=np.array(ok)
    print("number of error : ",error)
    return cores_out_in, ok, err_cores
