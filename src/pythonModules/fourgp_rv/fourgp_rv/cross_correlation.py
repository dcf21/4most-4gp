# -*- coding: utf-8 -*-

"""
This module implements a simple cross-correlation code for determining RVs.

It is a heavily cleaned up version of Jane Lin's GUESS code, as used by GALAH. Original code taken from
<https://github.com/jlin0504/GUESS>.
"""

from math import sqrt

import scipy.signal as sig
import numpy as np
from scipy.optimize import leastsq
import logging
from scipy.interpolate import UnivariateSpline

import fourgp_speclib
from fourgp_degrade.resample import SpectrumResampler

from .templates_resample import logarithmic_raster

logger = logging.getLogger(__name__)


class RvInstanceCrossCorrelation(object):
    """
    A class which is adapted from Jane Lin's GUESS code, as used by GALAH.
    """

    def __init__(self, spectrum_library):
        """
        Instantiate the RV code, and read from disk the library of template spectra used for cross correlation.

        :param spectrum_library:
            A SpectrumLibrary containing the template spectra we use for modelling.

        :type spectrum_library:
            SpectrumLibrary
        """

        assert isinstance(spectrum_library, fourgp_speclib.SpectrumLibrary), \
            "Argument to RvInstanceCrossCorrelation should be a SpectrumLibrary."
        self._spectrum_library = spectrum_library

        # Load template spectra
        spectrum_list = self._spectrum_library.search(continuum_normalised=1)

        # Load library of template spectra
        self._template_spectra = self._spectrum_library.open(ids=[i["specId"] for i in spectrum_list],
                                                             shared_memory=True)

        # Make a list of templates by 4MOST wavelength arm
        self.templates_by_arm = {}
        self.arm_rasters = {}

        # Loop over all the templates we've just loaded, sorting them by arm
        for i in range(len(self._template_spectra)):
            template_metadata = self._template_spectra.get_metadata(i)
            mode = template_metadata['mode']
            arm_name = template_metadata['arm_name']

            if mode not in self.templates_by_arm:
                self.templates_by_arm[mode] = {}
            if arm_name not in self.templates_by_arm[mode]:
                self.templates_by_arm[mode][arm_name] = []

            # Add this template to the list of available templates for this wavelength arm
            self.templates_by_arm[mode][arm_name].append(i)

            # If we haven't already recreated the fixed-log-step wavelength raster for this arm, do it now
            if arm_name not in self.arm_rasters:
                self.arm_rasters[arm_name] = logarithmic_raster(lambda_min=template_metadata['lambda_min'],
                                                           lambda_max=template_metadata['lambda_max'],
                                                           lambda_step=template_metadata['lambda_step'])

    def resample_single_arm(self, input_spectrum, arm_name):
        """
        Resample an input spectrum onto a raster with fixed logarithmic stride, representing a single 4MOST arm. We
        use the same wavelength rasters that the template spectra were sampled onto.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum
        :param arm_name:
            The name of the arm within this 4MOST mode
        :return:
            A Spectrum object containing a single arm of 4MOST data, resampled with a fixed logarithmic stride.
        """

        new_raster = self.arm_rasters[arm_name]

        resampler = SpectrumResampler(input_spectrum=input_spectrum)

        resampled_spectrum = resampler.onto_raster(output_raster=new_raster, resample_errors=True, resample_mask=False)

        return resampled_spectrum

    def estimate_rv_from_single_arm(self, input_spectrum, mode, arm_name):
        """
        Estimate the RV of a spectrum on the basis of data from a single arm. We return a list of RV estimates from
        cross correlation with each of the template spectra, and the chi-squared mismatch of each template spectrum
        which should be used to weight the RV estimates.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum. This should represent one arm of data only, and should
            be sampled onto the same fixed logarithmic wavelength stride as the template spectra.
        :param mode:
            The name of the 4MOST mode this arm is part of -- either LRS or HRS
        :param arm_name:
            The name of the arm within this 4MOST mode
        :return:
            List of [RV value, chi-squared weight]
        """

    def estimate_rv(self, input_spectrum, mode):
        """
        Estimate the RV of a spectrum on the basis of all of the 4MOST arms of either HRS or LRS.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum
        :param mode:
            The name of the 4MOST mode this arm is part of -- either LRS or HRS
        :return:
            [RV, error in RV]
        """

        rv_estimates = []

        arm_names = self.templates_by_arm[mode].keys()

        for arm_name in arm_names:
            new_rv_estimates = self.estimate_rv_from_single_arm(
                input_spectrum=self.resample_single_arm(input_spectrum=input_spectrum, arm_name=arm_name),
                mode=mode,
                arm_name=arm_name
            )
            rv_estimates.extend(new_rv_estimates)

        rv_mean = (sum([rv * weight for rv, weight in rv_estimates]) /
                   sum([weight for rv, weight in rv_estimates]))

        rv_variance = (sum([pow(rv-rv_mean, 2) * weight for rv, weight in rv_estimates]) /
                       sum([weight for rv, weight in rv_estimates]))

        rv_std_dev = sqrt(rv_variance)

        return rv_mean, rv_std_dev

    # --------------------------------------

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
