# -*- coding: utf-8 -*-

"""
This module implements a heavily cleaned up version of Jane Lin's GUESS code, as used by GALAH.

Original code taken from <https://github.com/jlin0504/GUESS>.
"""

import scipy.signal as sig
import numpy as np
from scipy.optimize import leastsq
import glob
import logging
import csv
from scipy.interpolate import UnivariateSpline
from numpy.linalg import lstsq
import pickle
import os
from os import path as os_path
from scipy.stats import sigmaclip

import fourgp_speclib
import fourgp_degrade

logger = logging.getLogger(__name__)

class RvInstanceGuess(object):
    """
    A class which is adapted from Jane Lin's GUESS code, as used by GALAH.
    """

    def __init__(self, spectrum_library):
        """
        Instantiate the RV code, and read from disk the library of template spectra.

        :param spectrum_library:
            A SpectrumLibrary containing the template spectra we use for modelling.

        :type spectrum_library:
            SpectrumLibrary
        """

        assert isinstance(spectrum_library, fourgp_speclib.SpectrumLibrary), \
            "Argument to RvInstance should be a SpectrumLibrary."
        self._spectrum_library = spectrum_library

        # Load template spectra
        spectrum_list = self._spectrum_library.search()

        # Load library of template spectra
        self._template_spectra = self._spectrum_library.open(ids=[i["specId"] for i in spectrum_list],
                                                             shared_memory=True)

        # Identifying 4MOST wavelength arms
        spectrograph_arms = ('blu', 'grn', 'red')
        all_wavelength_arms = {}

        wavelength_arms = fourgp_degrade.SpectrumProperties(
            wavelength_raster=self._template_spectra.wavelengths
        ).wavelength_arms()

        for arm_name, arm in zip(spectrograph_arms, wavelength_arms["wavelength_arms"]):
            arm_raster, mean_pixel_width = arm
            all_wavelength_arms[arm_name] = {
                "lambda_min": arm_raster[0],
                "lambda_max": arm_raster[-1],
                "lambda_step": mean_pixel_width
            }

        # Print information about the wavelength arms we have identified
        logger.info("Identified 4MOST wavelength arms as follows:")
        for arm_name in sorted(all_wavelength_arms.keys()):
            arm = all_wavelength_arms[arm_name]
            logger.info("{0:10s}: min {1:10.1f} max {2:10.1f} step {3:10.5f}".format(arm_name,
                                                                                     arm["lambda_min"],
                                                                                     arm["lambda_max"],
                                                                                     arm["lambda_step"]))

        self.wavelength_arms = all_wavelength_arms

    def median_filter_flux(self, filename):

        '''
        inputs:
        1. filename of the unnormalised 2dfdr spectrum, be like
        '/priv/miner3/galah/2dfdr_reductions/aaorun/614/150606/1506060029013971.fits' (ccd1)

        outputs:
        1. interpolated, median filtered flux
        2. the common wl grid
        3. s/n
        4. median filtered flux (un-interploated)
        5. original wl grid, without nans
        '''

        if filename.split('.')[-2].endswith('1'):
            grid=wave1
        if filename.split('.')[-2].endswith('2'):
            grid=wave2
        if filename.split('.')[-2].endswith('3'):
            grid=wave3

        #starno=filename.split('_')[-1].split('.')[0]#galahic
        starno=filename.split('.')[:-1]
        hdu=pyfits.open(filename)
        h = pyfits.getheader(filename)
        x = np.arange(h['NAXIS1']) + 1 - h['CRPIX1']
        wave = h['CDELT1']*x + h['CRVAL1']
        d=hdu[0].data
        #getting rid of the nan values at the beginning of the ccd
        nan_end=np.where(np.isnan(np.array(hdu[0].data)[::-1][:2000]))
        nan_beg=np.where(np.isnan(np.array(hdu[0].data)[:2000]))
        try:
            d=hdu[0].data[:-1*(nan_end[0][-1]+1)]
            wave=wave[:-1*(nan_end[0][-1]+1)]
        except:
            pass
        try:
            d=d[nan_beg[0][-1]+1:]
            wave=wave[nan_beg[0][-1]+1:]
        except:
            pass
        nans2=np.where(np.isnan(d))[0] #for nans in the middle. but WHY ARE THERE NANS IN THE MIDDLE!?!!? GARH!
        if len(nans2)!=0:
            d[nans2]=np.mean(d[np.where(np.isnan(d)==False)])
        #median filter smoothing and calculating s/n

        signal = np.median(d)
        noise = 1.4826/np.sqrt(2)*np.median(np.abs(d[1:] - d[:-1]))
        snr = signal / noise
        med_d = sig.medfilt(d,5) #5 pixels window
        w = np.where(abs(med_d - d) > sigma * noise)[0]
        d[w]=med_d[w]
        sn=True

        if filename.split('.')[-2].endswith('4'):
            grid=wave3 #dummy, to be consistent with the function output format.
            #exclude O2 absorption when measuring SNR
            non_tel=np.where(wave>=telluric_region)[0]
            fl=d[non_tel]
            signal = np.median(fl)
            noise = 1.4826/np.sqrt(2)*np.median(np.abs(fl[1:] - fl[:-1]))
            snr = signal / noise

        if snr <3 and filename.split('.')[-2].endswith('4')==False :
            sn=False
            sn_low.append(starno)
        grid_min=min(grid, key=lambda x:abs(x-wave[0]))
        grid_max=min(grid, key=lambda x:abs(x-wave[-1]))
        grid=grid[np.where(grid==grid_min)[0][0]:np.where(grid==grid_max)[0][0]]
        flux=np.interp(grid,wave,d)
        return flux,grid,sn,snr,d,wave

    def find_rv (self, model_ctm,flux,wave,filename):

        '''
        inputs:
        1. model
        2. data, ctm normalised
        3. common wavelength gird
        4. full path to fits

        outputs:
        1. rv
        2. max(correlation coeff) for this particular model
        '''

        def fitfunc2(p,x):
            return p[0]*x**2+p[1]*x+p[2]

        def errfunc2(p,x,y):
            return fitfunc2(p,x)-y

        flux -= np.mean(flux) #get rid of the top hat like thingy
        h=pyfits.getheader(filename)
        if filename.split('.')[-2].endswith('1'):
            cdelt1=wave1[1]-wave1[0]
        if  filename.split('.')[-2].endswith('2'):
            cdelt1=wave2[1]-wave2[0]
        if filename.split('.')[-2].endswith('3'):
            cdelt1=wave3[1]-wave3[0]

        flux /= np.sqrt(np.sum(flux**2)) #normalise the flux
        model_ctm1 = model_ctm- np.mean(model_ctm) #get rid of the top hat like thingy
        model_ctm1 /= np.sqrt(np.sum(model_ctm1**2)) #normalise the flux

        #window to smooth out the ends so it dont break when cross correlating yo
        slope=np.arange(len(wave)*0.10)# creating a ramp using 10% of the spectrum on either side
        window=np.hstack([slope,np.ones(len(wave)-len(slope)*2)*slope[-1],-1*slope+slope[-1]])
        window=window/max(window)
        model_ctm1=window*model_ctm1
        flux=window*flux

        #performing the cross correlation
        coeff=np.correlate(model_ctm1,flux,mode='same')

        max_index=np.where(coeff==max(coeff))[0][0]
        #fitting the cross correlation peak
        x_vals=np.array([max_index-1,max_index,max_index+1])
        d=x_vals[-1]-len(coeff)
        if d>=0:
            x_vals=np.array([max_index-2,max_index-1,max_index])
        try:
            y_vals=coeff[x_vals]
        except IndexError:
            x_vals=np.array([max_index,max_index+1,max_index+2])
        y_vals=coeff[x_vals]
        p0,p1,p2=leastsq(errfunc2,[0,0,0],args=(x_vals,y_vals))[0]
        x=np.arange(min(x_vals),max(x_vals),0.1)
        max_x=x[np.where(p0*x**2+p1*x+p2==max(p0*x**2+p1*x+p2))[0][0]]
        shift=max_x-len(model_ctm1)//2
        dlambda=cdelt1*shift
        velocity=-1*dlambda * 2.99792458e8 / (h['LAMBDAC']*1000)
        #print velocity,max(coeff)
        return velocity,max(coeff)

    def run_rv (self, filename,folder):

        '''
        inputs:
        1. path to file, just one ccd, it finds the rest 2 automatically
        2. folder

        outputs:
        1. 3 rvs for 3 ccds
        '''

        #starno=filename.split('_')[-1].split('.')[0]#galahic
        starno=filename.split('/')[-1][:-6]
        #date=filename.split('/')[-1].split('_')[0][:-1]
        #date=filename[:-1] #not really date
        rvss=[]
        snrs=[]
        #checks sn

        snr1=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'1'))[-3]#folder,date,ccd,galahic
        snr2=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'2'))[-3]#folder,date,ccd,galahic
        snr3=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'3'))[-3]#folder,date,ccd,galahic
        snr4=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'4'))[-3]

        snrs.append(snr1)
        snrs.append(snr2)
        snrs.append(snr3)
        snrs.append(snr4)

        for i in ['1','2','3']: #looping thru the ccds
            filename='%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,i)
            flux,grid,sn,snr,a,b = self.median_filter_flux(filename)
            snrs.append(snr)

            if i =='1':
                wavee=wave1
                mm=m1
                #excluding H regions when fitting continua
                h_s=min(list(grid), key=lambda x:abs(x-4847))
                h_e=min(list(grid), key=lambda x:abs(x-4900))
                h_s_index=np.where(grid==h_s)[0][0]
                h_e_index=np.where(grid==h_e)[0][0]
                s = UnivariateSpline(np.hstack([grid[:h_s_index],grid[h_e_index:]]), np.hstack([flux[:h_s_index],flux[h_e_index:]]), s=1e15) #s=smoothing factor
                flux_ctm=flux/s(grid)

            if i == '2':
                wavee=wave2
                mm=m2
                p0,p1,p2,p3,p4=leastsq(errfunc3,[0,0,0,0,0],args=(grid,flux))[0]
                fit=p0*grid**4+p1*grid**3+p2*grid**2+p3*grid+p4
                flux_ctm=flux/fit#(grid)

            if i == '3':
                wavee=wave3
                mm=m3
                h_s=min(list(grid), key=lambda x:abs(x-6530))
                h_e=min(list(grid), key=lambda x:abs(x-6592))
                h_s_index=np.where(grid==h_s)[0][0]
                h_e_index=np.where(grid==h_e)[0][0]
                s = UnivariateSpline(np.hstack([grid[:h_s_index],grid[h_e_index:]]), np.hstack([flux[:h_s_index],flux[h_e_index:]]), s=1e15)
                flux_ctm=flux/s(grid)

            model_s=np.where(wavee==grid[0])[0][0]
            model_e=np.where(wavee==grid[-1])[0][0]
            models=mm[0].data[:,model_s:model_e+1]

            coeffs=[]
            rvs=[]
            #finding the rvs for all 15 models
            for j in range(len(models))[0:-3]:
                rv,coeff=find_rv(models[j],flux_ctm,grid,filename)
                coeffs.append(coeff)
                rvs.append(rv)
            good_coeff=np.where(np.array(coeffs)>0.3)[0]
            #only take the rvs where the cross correlation is reasonable (ie cross correlation
            #coeff >0.3)
            good_rv=np.array(rvs)[good_coeff]
            #pdb.set_trace()
            weights=np.array(coeffs)[good_coeff]/float(sum(np.array(coeffs)[good_coeff]))
            #final rv for this ccd is the weighted sum of rvs from ~15 models, weighted by ccoeff
            print('weighted rv', sum(weights*good_rv))
            if sum(weights*good_rv) ==0 and i=='2': #incase the cc is crap
                rvss.append(999)
                continue
            if sum(weights*good_rv) ==0 and (i=='1' or i=='3'):
                bad_weights.append(starno)
                return 999,snrs,999,1
            rvss.append(sum(weights*good_rv))

        if rvss != [] and len(rvss)==3 :
            return rvss,snrs,0,0
