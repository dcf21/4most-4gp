#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This python package provides a wrapper for functions provided by the 4MOST Facility Simulator, 4FS.
"""

import os
from os import path as os_path
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table


class FourFS:
    def __init__(self,
                 path_to_4fs="/home/dcf21/iwg7_pipeline/OpSys/OpSim"
                 ):
        """
        Instantiate a class for calling 4FS.

        :param path_to_4fs:
            Path to where 4FS binaries can be found.
        """
        self.path_to_4fs = path_to_4fs

        # Create temporary directory
        self.id_string = "4fs_{:d}".format(os.getpid())
        self.tmp_dir = os_path.join("/tmp", self.id_string)
        os.system("mkdir -p {}".format(self.tmp_dir))

    def close(self):
        # Remove temporary directory
        os.system("rm -Rf {}".format(self.tmp_dir))

    @staticmethod
    def make_fits_spectemplate(self,
                               spec_file='test1.spec',
                               output='test1.fits',
                               path='/Users/khawkins/Desktop/BACCHUS/fmost_spec/',
                               resolution=80000,
                               continuum=False
                               ):
        """
        Generate a 4FS readable fits template from the Turbospectrum output.
        :param spec_file:
            Name of the ascii Turbospectrum output
        :param output:
            Name of the output fits file
        :param path:
            Path of the spec file
        |
        :return:
            Fits file written to path + output
        """

        wavelength, continuum_normalised_flux, flux = np.loadtxt(path + spec_file, unpack=True)
        magnitude = 12.0
        lambda_min = min(wavelength)
        delta_lambda = wavelength[1] - wavelength[0]
        if not continuum:
            hdu_1 = fits.PrimaryHDU(flux * 10E-15 / max(flux))
        else:
            hdu_1 = fits.PrimaryHDU((flux * 10E-15 / max(flux)) / continuum_normalised_flux)
        hdu_1.header['CRTYPE1'] = "LINEAR  "
        hdu_1.header['CRPIX1'] = 1.0
        hdu_1.header['CRVAL1'] = lambda_min
        hdu_1.header['CDELT1'] = delta_lambda
        hdu_1.header['CUNIT1'] = "Angstrom"
        hdu_1.header['BUNIT'] = "erg/s/cm2/Angstrom"
        hdu_1.header['ABMAG'] = magnitude
        if resolution is not None:
            fwhm_res = np.mean(wavelength / resolution)
            hdu_1.header['RESOLUTN'] = fwhm_res

        hdu_list = fits.HDUList([hdu_1])
        hdu_list.writeto(path + output)

    def generate_4FS_templatelist(self,
                                  path='/Users/khawkins/Desktop/BACCHUS/fmost_spec/',
                                  resolution=40000,
                                  spectable='/Users/khawkins/Desktop/BACCHUS/fmost_spec/fmost_testsetR40000.tab'
                                  ):
        """
        Generate 4FS template list using the spectra that are listed in the input (spectable) table.
        |
        | INPUT:
        | path = path of the synthetic spectra
        | R = resoulution that will be used for the 4FS
        | spectable = full path of the table with the spectra names, parameters and abundances which is the output of
        | make_synthetic_grid1() function.
        |
        | OUTPUT:
        | text file in the template_list file which contains the input template_list for the 4FS
        --------------------------------------------------------------------------------
        '''
        """

        # Location to write template list to be used for 4FS
        template_list = '/Users/khawkins/Documents/4FS_source/examples/template_list_test.txt'

        allstars = []
        t = Table.read(spectable, format='ascii')
        objname = t['Starname']
        writestr = '#OBJECTNAME           FILENAME                                                 RULESET SIZE REDSHIFT MAG  MAG_RANGE\n'
        for i in range(len(objname)):
            print objname[i]
            # Look for the spectra and its continuum fits files and removes it first just in case an old file is there
            try:
                os.remove(path + '%s.fits' % objname[i])
                os.remove(path + '%s_c.fits' % objname[i])
            except OSError:
                pass

            # Generates the 4FS readable fits files from the output spectra from Turbospectrum
            self.make_fits_spectemplate(spec_file='%s.spec' % objname[i],
                                        output='%s.fits' % objname[i],
                                        path=path,
                                        resolution=resolution)
            self.make_fits_spectemplate(spec_file='%s.spec' % objname[i],
                                        output='%s_c.fits' % objname[i],
                                        path=path,
                                        resolution=resolution,
                                        continuum=True)

            writestr += '%s_SNR05 %s %s %s %s %s %s\n' % (
                objname[i], path + '%s.fits' % objname[i], 'goodSNR05', '0.0', '0.0', '15.0', '15.0')
            writestr += '%s_SNR10 %s %s %s %s %s %s\n' % (
                objname[i], path + '%s.fits' % objname[i], 'goodSNR10', '0.0', '0.0', '15.0', '15.0')
            writestr += '%s_SNR15 %s %s %s %s %s %s\n' % (
                objname[i], path + '%s.fits' % objname[i], 'goodSNR15', '0.0', '0.0', '15.0', '15.0')
            writestr += '%s_SNR20 %s %s %s %s %s %s\n' % (
                objname[i], path + '%s.fits' % objname[i], 'goodSNR20', '0.0', '0.0', '15.0', '15.0')
            writestr += '%s_SNR50 %s %s %s %s %s %s\n' % (
                objname[i], path + '%s.fits' % objname[i], 'goodSNR50', '0.0', '0.0', '15.0', '15.0')
            writestr += '%s_SNR100 %s %s %s %s %s %s\n' % (
                objname[i], path + '%s.fits' % objname[i], 'goodSNR100', '0.0', '0.0', '15.0', '15.0')
            writestr += '%s_SNR250 %s %s %s %s %s %s\n' % (
                objname[i], path + '%s.fits' % objname[i], 'goodSNR250', '0.0', '0.0', '15.0', '15.0')
            writestr += '%s_c %s %s %s %s %s %s\n' % (
                objname[i], path + '%s_c.fits' % objname[i], 'goodSNR250', '0.0', '0.0', '15.0', '15.0')
            allstars.append('%s_SNR05' % objname[i])
            allstars.append('%s_SNR10' % objname[i])
            allstars.append('%s_SNR15' % objname[i])
            allstars.append('%s_SNR20' % objname[i])
            allstars.append('%s_SNR50' % objname[i])
            allstars.append('%s_SNR100' % objname[i])
            allstars.append('%s_SNR250' % objname[i])

        # ---write template list to file---
        f = open(template_list, 'w')
        f.write(writestr)
        f.close()
        print writestr

        return 1

    def combine_spectra(self,
                        stars,
                        path='/Users/khawkins/Documents/4FS_source/examples/outdir_LRS/',
                        setup='LRS',
                        save_output=False
                        ):
        """
        Combine the multiple fits files from the 4FS into a single ascii file with all arms included.

        NOTE: The point at which one are is favoured over another in the stitching
              process is hardcoded should be changed in future verisons!
        |
        |
        | INPUT:
        | stars = which stars to combine from the 4FS output
        | setup = (1) 'LRS' - low res or (2) 'HRS' - high res
        | path = path of 4FS output to be used to combine the spectra
        | plt = True or False to plot the final combined spectra used for diagnostics
        | sav = True (or False) to save the output to file will take the name of star.txt
        |
        :return :
        | text file in the template_list file which contains the input template_list for the 4FS
        """

        # Path to save the final output spectra
        save_path = '/Users/khawkins/Documents/4FS_source/examples/testset/%s/' % setup

        for i in range(len(stars)):
            root_star_name = stars[i].split('_')[0]
            # Load in the three output spectra (blue, green, red arms)---
            d = fits.open(path + 'specout_template_%s_%s_blue.fits' % (stars[i], setup))
            d2 = fits.open(path + 'specout_template_%s_%s_green.fits' % (stars[i], setup))
            d3 = fits.open(path + 'specout_template_%s_%s_red.fits' % (stars[i], setup))

            # Read the data of interest
            lam1 = d[2].data['LAMBDA']
            lam2 = d2[2].data['LAMBDA']
            lam3 = d3[2].data['LAMBDA']
            flux1 = d[2].data['REALISATION'] - d[2].data['SKY']
            flux2 = d2[2].data['REALISATION'] - d2[2].data['SKY']
            flux3 = d3[2].data['REALISATION'] - d3[2].data['SKY']
            snr1 = d[2].data['SNR']
            snr2 = d2[2].data['SNR']
            snr3 = d3[2].data['SNR']

            if setup == 'LRS':
                # ---stitch arms based on where SNR (hardcorded here may need changed in future versions)--
                ind1 = np.where(lam1 <= 5327.7)[0]
                ind2 = np.where((lam2 > 5327.7) & (lam2 <= 7031.7))[0]
                ind3 = np.where(lam3 > 7031.7)[0]
                lam1 = lam1[ind1]
                lam2 = lam2[ind2]
                lamma3 = lam3[ind3]
                flux1 = flux1[ind1]
                flux2 = flux2[ind2]
                flux3 = flux3[ind3]
                snr1 = snr1[ind1]
                snr2 = snr2[ind2]
                snr3 = snr3[ind3]
            lamtot = np.array(list(lam1) + list(lam2) + list(lam3))
            fluxtot = np.array(list(flux1) + list(flux2) + list(flux3))
            SNRtot = np.array(list(snr1) + list(snr2) + list(snr3))

            # Load continuum spectra
            dc = fits.open(path + 'specout_template_%s_c_%s_blue.fits' % (root_star_name, setup))
            d2c = fits.open(path + 'specout_template_%s_c_%s_green.fits' % (root_star_name, setup))
            d3c = fits.open(path + 'specout_template_%s_c_%s_red.fits' % (root_star_name, setup))
            lam1c = d[2].data['LAMBDA']
            lam2c = d2[2].data['LAMBDA']
            lam3c = d3[2].data['LAMBDA']
            flux1c = dc[2].data['FLUENCE']
            flux2c = d2c[2].data['FLUENCE']
            flux3c = d3c[2].data['FLUENCE']
            if setup == 'LRS':
                ind1 = np.where(lam1c <= 5327.7)[0]
                ind2 = np.where((lam2c > 5327.7) & (lam2c <= 7031.7))[0]
                ind3 = np.where(lam3c > 7031.7)[0]
                lam1c = lam1c[ind1]
                lam2c = lam2c[ind2]
                lam3c = lam3c[ind3]
                flux1c = flux1c[ind1]
                flux2c = flux2c[ind2]
                flux3c = flux3c[ind3]

            normflux1 = flux1 / (flux1c * max(d[2].data['FLUENCE']) / max(flux1c))
            normflux2 = flux2 / (flux2c * max(d2[2].data['FLUENCE']) / max(flux2c))
            normflux3 = flux3 / (flux3c * max(d3[2].data['FLUENCE']) / max(flux3c))

            # Combine everything into one set of arrays to be saved
            lamtotc = np.array(list(lam1c) + list(lam2c) + list(lam3c))
            fluxtotc = np.array(list(flux1c * max(d[2].data['FLUENCE']) / max(flux1c)) + list(
                flux2c * max(d2[2].data['FLUENCE']) / max(flux2c)) + list(
                flux3c * max(d3[2].data['FLUENCE']) / max(flux3c)))
            normfluxtot = np.array(list(normflux1) + list(normflux2) + list(normflux3))
            uncertflux = normfluxtot / SNRtot

            # Continuum normalize
            normfluxtot = fluxtot / fluxtotc

            # Remove bad pixels
            normfluxtot[np.where((normfluxtot > 2.0) | (normfluxtot <= 0.00))[
                0]] = 0  # pixel where flux > 2 or flux < 0 get reset done for other downstream codes
            # return lamtotc,normfluxtot

            # Save if wanted
            if save_output:
                np.savetxt(save_path + '%s.txt' % stars[i], np.array([lamtotc, normfluxtot, uncertflux]).T, fmt='%.3f',
                           header='Wavelength\tNormFlux\tNormFluxerr', delimiter='\t')

            return uncertflux, normfluxtot, SNRtot


ETC_input_params_HRS = """
# example parameter file for 4FS_ETC

SIM.CODE_NAME                = 'example_HRS'                                                                          # Human readable codename for this run of the 4FS_TS
SIM.OUTDIR                   = './outdir_HRS'                                                                         # Where should we put output files?

SIM.MODE                     = 'CALC_TEXP'                                                                            # Should we calculate SNR from given TEXP, or TEXP/MAG from given SNR? (CALC_SNR,CALC_TEXP,CALC_MAG)
SIM.OUTPUT                   = 'SUMMARY,SPECTRA_ALL'                                                                  # Which output types to produce? (ADD LIST OF OPTIONS HERE)
SIM.SPECFORMAT               = 'TABLE,NATIVE'                                                                         # Which output spectral formats should be produced? (IMAGE,TABLE;NATIVE,RESAMPLED)
SIM.CLOBBER                  = TRUE                                                                                   # Run in clobber mode? (existing output files will be overwritten)

TEMPLATES.FILENAME           = 'template_list_test.txt'                                                               # Name of file containing the list of spectral templates
RULELIST.FILENAME            = 'rulelist.txt'                                                                         # Name of file containing the list of spectral success rules
RULESETLIST.FILENAME         = 'ruleset.txt'                                                                          # Name of file containing the list of spectral success rulesets

SIM.NUM_FILTERS              = 5                                                                                      # How many filters to read?
SIM.NORM_FILTER.MAGSYS       = 'AB'                                                                                   # Magnitude system for normailsing templates (AB,Vega)
SIM.NORM_FILTER1.NAME        = 'SDSS u'                                                                               # Name of filter bandpass
SIM.NORM_FILTER1.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_u_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER2.NAME        = 'SDSS g'                                                                               # Name of filter bandpass
SIM.NORM_FILTER2.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_g_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER3.NAME        = 'SDSS r'                                                                               # Name of filter bandpass
SIM.NORM_FILTER3.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_r_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER4.NAME        = 'SDSS i'                                                                               # Name of filter bandpass
SIM.NORM_FILTER4.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_i_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER5.NAME        = 'SDSS z'                                                                               # Name of filter bandpass
SIM.NORM_FILTER5.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_z_transmission_curve.fits'               # Name of file containing the normalising filter bandpass

SPECTRO.FIBER_DIAM           = 1.45                                                                            # Fibre diameter (arcsec)
SPECTRO.EFFECTIVE_AREA       = 8.3975                                                                          # Effective collecting area of telescope (m^2)
SPECTRO.SKYSUB_RESIDUAL      = 0.0                                                                             # Fractional uncertaintity on sky subtraction
SPECTRO.NUM_ARMS             = 3                                                                               # Number of spectrograph arms

SPECTRO.ARM1.CODENAME        = 'HRS_blue'                                                                      # Codename for spectrograph arm
SPECTRO.ARM1.RES_FILENAME    = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_resolution_curve_slitlet_10_interp_blue.fits'     # Filename describing spectral resolution
SPECTRO.ARM1.TPUT_FILENAME   = '4FS_ETC_system_model_v0.2/HRS/hrs_blue_material_4fs_efficiency_total.fits'     # Filename describing spectral throughput
SPECTRO.ARM1.APER_SIZE       = 4.0                                                                             # Number of pixels to sum over in the cross-dispersion direction
SPECTRO.ARM1.APER_EEF        = 0.9545                                                                          # Fraction of light in the extraction aperture
SPECTRO.ARM1.PEAK_PIX_FRAC   = 0.3702                                                                          # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)
SPECTRO.ARM1.READ_NOISE      = 2.5                                                                             # CCD read noise (e-/pix)
SPECTRO.ARM1.DARK_CURRENT    = 3.0                                                                             # CCD dark current (e-/hr/pix)
SPECTRO.ARM1.FULL_WELL       = 350000                                                                          # Full well capacity of the CCD (e-/pix)
SPECTRO.ARM1.BINNING.DISP    = 1                                                                               # On-chip binning in dispersion direction
SPECTRO.ARM1.BINNING.CROSS   = 1                                                                               # On-chip binning in cross-dispersion direction
SPECTRO.ARM1.LAMBDA.TYPE     = 'FULLFILE'                                                                      # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE
SPECTRO.ARM1.LAMBDA.FILENAME = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_wavelength_solution_slitlet_10_interp_blue.fits'  # Filename describing wavelength solution

SPECTRO.ARM2.CODENAME        = 'HRS_green'                                                                     # Codename for spectrograph arm
SPECTRO.ARM2.RES_FILENAME    = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_resolution_curve_slitlet_10_interp_green.fits'    # Filename describing spectral resolution
SPECTRO.ARM2.TPUT_FILENAME   = '4FS_ETC_system_model_v0.2/HRS/hrs_green_material_4fs_efficiency_total.fits'    # Filename describing spectral throughput
SPECTRO.ARM2.APER_SIZE       = 4.0                                                                             # Number of pixels to sum over in the cross-dispersion direction
SPECTRO.ARM2.APER_EEF        = 0.9545                                                                          # Fraction of light in the extraction aperture
SPECTRO.ARM2.PEAK_PIX_FRAC   = 0.3702                                                                          # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)
SPECTRO.ARM2.READ_NOISE      = 2.5                                                                             # CCD read noise (e-/pix)
SPECTRO.ARM2.DARK_CURRENT    = 3.0                                                                             # CCD dark current (e-/hr/pix)
SPECTRO.ARM2.FULL_WELL       = 350000                                                                          # Full well capacity of the CCD (e-/pix)
SPECTRO.ARM2.BINNING.DISP    = 1                                                                               # On-chip binning in dispersion direction
SPECTRO.ARM2.BINNING.CROSS   = 1                                                                               # On-chip binning in cross-dispersion direction
SPECTRO.ARM2.LAMBDA.TYPE     = 'FULLFILE'                                                                      # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE
SPECTRO.ARM2.LAMBDA.FILENAME = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_wavelength_solution_slitlet_10_interp_green.fits' # Filename describing wavelength solution

SPECTRO.ARM3.CODENAME        = 'HRS_red'                                                                       # Codename for spectrograph arm
SPECTRO.ARM3.RES_FILENAME    = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_resolution_curve_slitlet_10_interp_red.fits'      # Filename describing spectral resolution
SPECTRO.ARM3.TPUT_FILENAME   = '4FS_ETC_system_model_v0.2/HRS/hrs_red_material_4fs_efficiency_total.fits'      # Filename describing spectral throughput
SPECTRO.ARM3.APER_SIZE       = 4.0                                                                             # Number of pixels to sum over in the cross-dispersion direction
SPECTRO.ARM3.APER_EEF        = 0.9545                                                                          # Fraction of light in the extraction aperture
SPECTRO.ARM3.PEAK_PIX_FRAC   = 0.3702                                                                          # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)
SPECTRO.ARM3.READ_NOISE      = 2.5                                                                             # CCD read noise (e-/pix)
SPECTRO.ARM3.DARK_CURRENT    = 3.0                                                                             # CCD dark current (e-/hr/pix)
SPECTRO.ARM3.FULL_WELL       = 350000                                                                          # Full well capacity of the CCD (e-/pix)
SPECTRO.ARM3.BINNING.DISP    = 1                                                                               # On-chip binning in dispersion direction
SPECTRO.ARM3.BINNING.CROSS   = 1                                                                               # On-chip binning in cross-dispersion direction
SPECTRO.ARM2.LAMBDA.TYPE     = 'FULLFILE'                                                                      # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE
SPECTRO.ARM2.LAMBDA.FILENAME = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_wavelength_solution_slitlet_10_interp_green.fits' # Filename describing wavelength solution

SPECTRO.ARM3.CODENAME        = 'HRS_red'                                                                       # Codename for spectrograph arm
SPECTRO.ARM3.RES_FILENAME    = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_resolution_curve_slitlet_10_interp_red.fits'      # Filename describing spectral resolution
SPECTRO.ARM3.TPUT_FILENAME   = '4FS_ETC_system_model_v0.2/HRS/hrs_red_material_4fs_efficiency_total.fits'      # Filename describing spectral throughput
SPECTRO.ARM3.APER_SIZE       = 4.0                                                                             # Number of pixels to sum over in the cross-dispersion direction
SPECTRO.ARM3.APER_EEF        = 0.9545                                                                          # Fraction of light in the extraction aperture
SPECTRO.ARM3.PEAK_PIX_FRAC   = 0.3702                                                                          # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)
SPECTRO.ARM3.READ_NOISE      = 2.5                                                                             # CCD read noise (e-/pix)
SPECTRO.ARM3.DARK_CURRENT    = 3.0                                                                             # CCD dark current (e-/hr/pix)
SPECTRO.ARM3.FULL_WELL       = 350000                                                                          # Full well capacity of the CCD (e-/pix)
SPECTRO.ARM3.BINNING.DISP    = 1                                                                               # On-chip binning in dispersion direction
SPECTRO.ARM3.BINNING.CROSS   = 1                                                                               # On-chip binning in cross-dispersion direction
SPECTRO.ARM3.LAMBDA.TYPE     = 'FULLFILE'                                                                      # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE
SPECTRO.ARM3.LAMBDA.FILENAME = '4FS_ETC_system_model_v0.2/HRS/4MOST_HRS_wavelength_solution_slitlet_10_interp_red.fits'   # Filename describing wavelength solution

FIBRECOUPLING.TYPE           = 'FILE'                                                                          # Method by which fibre losses are calculated (NONE,FIXED,SEEING,FILE)
FIBRECOUPLING.FILENAME       = '4FS_ETC_system_model_v0.2/fibre_coupling/geometrical_throughput.fits'          # File describing fibre losses
FIBRECOUPLING.FRAC_SKY       = 1.0                                                                             # Fraction of sky light transmitted into fibre

SKY.TRANSMISSION.FILENAME    = '4FS_ETC_system_model_v0.2/sky/paranal_sky_transmission_vectors.fits'           # Name of file containing the sky transmission info
SKY.EMISSION.FILENAME        = '4FS_ETC_system_model_v0.2/sky/paranal_sky_emission_vectors.fits'               # Name of file containing the sky emission info

OBS_PARAMS.INTERP_METHOD     = 'NEAREST'                              # Method to use when interpolating obs params grid: NEAREST,LINEAR,SPLINE
OBS_PARAMS.SKYBRIGHT_TYPE    = 'ZENITH'                               # Is the specified sky brightness to be measured at ZENITH or LOCALly?
OBS_PARAMS.AIRMASS           = "1.3"                                  # List of airmasses to simulate
OBS_PARAMS.IQ                = "1.1"                                  # List of delivered image quality values to simulate (V-band,FWHM,arcsec)
OBS_PARAMS.SKYBRIGHT         = "21.77"                                # List of sky brightnesses to simulate (V-band,ABmag/arcsec2)
OBS_PARAMS.TILT              = "6.0"                                  # List of tilts to simulate (mm)
OBS_PARAMS.MISALIGNMENT      = "0.1"                                  # List of fibre->target misalignments to simulate (arcsec)
OBS_PARAMS.TEXP              = "1200"
OBS_PARAMS.NSUB              = "1"
#OBS_PARAMS.TEXP              = "100,300,500,1000,1200"               # List of total exposure times to simulate (sec)
#OBS_PARAMS.NSUB              = "1,1,1,1,1"                           # List of numbers of sub-exposures to simulate
"""

ETC_input_params_LRS = """
# example parameter file for 4FS_ETC

SIM.CODE_NAME                = 'example_LRS'                                                                          # Human readable codename for this run of the 4FS_TS
SIM.OUTDIR                   = './outdir_LRS'                                                                         # Where should we put output files?

SIM.MODE                     = 'CALC_TEXP'                                                                            # Should we calculate SNR from given TEXP, or TEXP/MAG from given SNR? (CALC_SNR,CALC_TEXP,CALC_MAG)
SIM.OUTPUT                   = 'SUMMARY,SPECTRA_ALL'                                                                  # Which output types to produce? (ADD LIST OF OPTIONS HERE)
SIM.SPECFORMAT               = 'TABLE,NATIVE'                                                                         # Which output spectral formats should be produced? (IMAGE,TABLE;NATIVE,RESAMPLED)
SIM.CLOBBER                  = TRUE                                                                                   # Run in clobber mode? (existing output files will be overwritten)

TEMPLATES.FILENAME           = 'template_list_test.txt'                                                               # Name of file containing the list of spectral templates
RULELIST.FILENAME            = 'rulelist.txt'                                                                         # Name of file containing the list of spectral success rules
RULESETLIST.FILENAME         = 'ruleset.txt'                                                                          # Name of file containing the list of spectral success rulesets

SIM.NUM_FILTERS              = 5                                                                                      # How many filters to read?
SIM.NORM_FILTER.MAGSYS       = 'AB'                                                                                   # Magnitude system for normailsing templates (AB,Vega)
SIM.NORM_FILTER1.NAME        = 'SDSS u'                                                                               # Name of filter bandpass
SIM.NORM_FILTER1.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_u_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER2.NAME        = 'SDSS g'                                                                               # Name of filter bandpass
SIM.NORM_FILTER2.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_g_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER3.NAME        = 'SDSS r'                                                                               # Name of filter bandpass
SIM.NORM_FILTER3.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_r_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER4.NAME        = 'SDSS i'                                                                               # Name of filter bandpass
SIM.NORM_FILTER4.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_i_transmission_curve.fits'               # Name of file containing the normalising filter bandpass
SIM.NORM_FILTER5.NAME        = 'SDSS z'                                                                               # Name of filter bandpass
SIM.NORM_FILTER5.FILENAME    = '4FS_ETC_system_model_v0.2/filter_curves/SDSS_z_transmission_curve.fits'               # Name of file containing the normalising filter bandpass

SPECTRO.FIBER_DIAM           = 1.45                                                                                   # Fibre diameter (arcsec)
SPECTRO.EFFECTIVE_AREA       = 8.3975                                                                                 # Effective collecting area of telescope (m^2)
SPECTRO.SKYSUB_RESIDUAL      = 0.0                                                                                    # Fractional uncertaintity on sky subtraction
SPECTRO.NUM_ARMS             = 3                                                                                      # Number of spectrograph arms

SPECTRO.ARM1.CODENAME        = 'LRS_blue'                                                                             # Codename for spectrograph arm
SPECTRO.ARM1.RES_FILENAME    = '4FS_ETC_system_model_v0.2/LRS/4MOST_LRS_resolution_curve_middle_interp_blue.fits'     # Filename describing spectral resolution
SPECTRO.ARM1.TPUT_FILENAME   = '4FS_ETC_system_model_v0.2/LRS/lrs_blue_material_4fs_efficiency_total.fits'            # Filename describing spectral throughput
SPECTRO.ARM1.APER_SIZE       = 4.0                                                                                    # Number of pixels to sum over in the cross-dispersion direction
SPECTRO.ARM1.APER_EEF        = 0.9545                                                                                 # Fraction of light in the extraction aperture
SPECTRO.ARM1.PEAK_PIX_FRAC   = 0.3702                                                                                 # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)
SPECTRO.ARM1.READ_NOISE      = 2.5                                                                                    # CCD read noise (e-/pix)
SPECTRO.ARM1.DARK_CURRENT    = 3.0                                                                                    # CCD dark current (e-/hr/pix)
SPECTRO.ARM1.FULL_WELL       = 350000                                                                                 # Full well capacity of the CCD (e-/pix)
SPECTRO.ARM1.BINNING.DISP    = 1                                                                                      # On-chip binning in dispersion direction
SPECTRO.ARM1.BINNING.CROSS   = 1                                                                                      # On-chip binning in cross-dispersion direction
SPECTRO.ARM1.LAMBDA.TYPE     = 'FULLFILE'                                                                             # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE
SPECTRO.ARM1.LAMBDA.FILENAME = '4FS_ETC_system_model_v0.2/LRS/4MOST_LRS_wavelength_solution_middle_interp_blue.fits'  # Filename describing wavelength solution

SPECTRO.ARM2.CODENAME        = 'LRS_green'                                                                            # Codename for spectrograph arm
SPECTRO.ARM2.RES_FILENAME    = '4FS_ETC_system_model_v0.2/LRS/4MOST_LRS_resolution_curve_middle_interp_green.fits'    # Filename describing spectral resolution
SPECTRO.ARM2.TPUT_FILENAME   = '4FS_ETC_system_model_v0.2/LRS/lrs_green_material_4fs_efficiency_total.fits'           # Filename describing spectral throughput
SPECTRO.ARM2.APER_SIZE       = 4.0                                                                                    # Number of pixels to sum over in the cross-dispersion direction
SPECTRO.ARM2.APER_EEF        = 0.9545                                                                                 # Fraction of light in the extraction aperture
SPECTRO.ARM2.PEAK_PIX_FRAC   = 0.3702                                                                                 # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)
SPECTRO.ARM2.READ_NOISE      = 2.5                                                                                    # CCD read noise (e-/pix)
SPECTRO.ARM2.DARK_CURRENT    = 3.0                                                                                    # CCD dark current (e-/hr/pix)
SPECTRO.ARM2.FULL_WELL       = 350000                                                                                 # Full well capacity of the CCD (e-/pix)
SPECTRO.ARM2.BINNING.DISP    = 1                                                                                      # On-chip binning in dispersion direction
SPECTRO.ARM2.BINNING.CROSS   = 1                                                                                      # On-chip binning in cross-dispersion direction
SPECTRO.ARM2.LAMBDA.TYPE     = 'FULLFILE'                                                                             # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE
SPECTRO.ARM2.LAMBDA.FILENAME = '4FS_ETC_system_model_v0.2/LRS/4MOST_LRS_wavelength_solution_middle_interp_green.fits' # Filename describing wavelength solution

SPECTRO.ARM3.CODENAME        = 'LRS_red'                                                                              # Codename for spectrograph arm
SPECTRO.ARM3.RES_FILENAME    = '4FS_ETC_system_model_v0.2/LRS/4MOST_LRS_resolution_curve_middle_interp_red.fits'      # Filename describing spectral resolution
SPECTRO.ARM3.TPUT_FILENAME   = '4FS_ETC_system_model_v0.2/LRS/lrs_red_material_4fs_efficiency_total.fits'             # Filename describing spectral throughput
SPECTRO.ARM3.APER_SIZE       = 4.0                                                                                    # Number of pixels to sum over in the cross-dispersion direction
SPECTRO.ARM3.APER_EEF        = 0.9545                                                                                 # Fraction of light in the extraction aperture
SPECTRO.ARM3.PEAK_PIX_FRAC   = 0.3702                                                                                 # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)
SPECTRO.ARM3.READ_NOISE      = 2.5                                                                                    # CCD read noise (e-/pix)
SPECTRO.ARM3.DARK_CURRENT    = 3.0                                                                                    # CCD dark current (e-/hr/pix)
SPECTRO.ARM3.FULL_WELL       = 350000                                                                                 # Full well capacity of the CCD (e-/pix)
SPECTRO.ARM3.BINNING.DISP    = 1                                                                                      # On-chip binning in dispersion direction
SPECTRO.ARM3.BINNING.CROSS   = 1                                                                                      # On-chip binning in cross-dispersion direction
SPECTRO.ARM3.LAMBDA.TYPE     = 'FULLFILE'                                                                             # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE
SPECTRO.ARM3.LAMBDA.FILENAME = '4FS_ETC_system_model_v0.2/LRS/4MOST_LRS_wavelength_solution_middle_interp_red.fits'   # Filename describing wavelength solution

FIBRECOUPLING.TYPE           = 'FILE'                                                                                 # Method by which fibre losses are calculated (NONE,FIXED,SEEING,FILE)
FIBRECOUPLING.FILENAME       = '4FS_ETC_system_model_v0.2/fibre_coupling/geometrical_throughput.fits'                 # File describing fibre losses
FIBRECOUPLING.FRAC_SKY       = 1.0                                                                                    # Fraction of sky light transmitted into fibre

SKY.TRANSMISSION.FILENAME    = '4FS_ETC_system_model_v0.2/sky/paranal_sky_transmission_vectors.fits'                  # Name of file containing the sky transmission info
SKY.EMISSION.FILENAME        = '4FS_ETC_system_model_v0.2/sky/paranal_sky_emission_vectors.fits'                      # Name of file containing the sky emission info

OBS_PARAMS.INTERP_METHOD     = 'NEAREST'                              # Method to use when interpolating obs params grid: NEAREST,LINEAR,SPLINE
OBS_PARAMS.SKYBRIGHT_TYPE    = 'ZENITH'                               # Is the specified sky brightness to be measured at ZENITH or LOCALly?
OBS_PARAMS.AIRMASS           = "1.3"                                  # List of airmasses to simulate
OBS_PARAMS.IQ                = "1.1"                                  # List of delivered image quality values to simulate (V-band,FWHM,arcsec)
OBS_PARAMS.SKYBRIGHT         = "21.77"                                # List of sky brightnesses to simulate (V-band,ABmag/arcsec2)
OBS_PARAMS.TILT              = "6.0"                                  # List of tilts to simulate (mm)
OBS_PARAMS.MISALIGNMENT      = "0.1"                                  # List of fibre->target misalignments to simulate (arcsec)
OBS_PARAMS.TEXP              = "500"                                  # List of total exposure times to simulate (sec)
OBS_PARAMS.NSUB              = "1"                                    # List of numbers of sub-exposures to simulate
"""

rulelist = """
#RULE VARIABLE METRIC OPER VALUE L_MIN L_MAX L_UNIT DELTAL DELTAL_UNIT
SNR60         SNR MEDIAN GE   60.0  522.5 569.0 NM     1.0    PIX
MEDIANSNRLRS  SNR MEDIAN DIV   1.0  560.0 620.0 NM     1.0    PIX
MEDIANSNRHRS  SNR MEDIAN DIV   1.0  520.0 560.0 NM     1.0    PIX
MEDIANSNR     SNR MEDIAN DIV   1.0  618.0 668.0 NM     1.0    PIX
"""
