#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
This python package provides a wrapper for functions provided by the 4MOST Facility Simulator, 4FS.

This file contains the standard input configuration parameters for the Exposure Time Calculator (ETC).
"""

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
OBS_PARAMS.AIRMASS           = "1.4"                                  # List of airmasses to simulate
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
DCFSNR        SNR MEDIAN DIV   1.0  618.0 668.0 NM     1.0    PIX

############################################################################################################
# Title: Interpretation of 4MOST Phase A Spectral Success Criteria for input into "round9" of the 4FS
# Author: T. Dwelly (dwelly@mpe.mpg.de)
# Last updated: 12.09.2016 (rev11)
#
# Scope: This file is the author's interpretation of the Spectral Success Criteria that were used by the seven Design 
#        Reference Surveys (DRSs) that were considered in the 4MOST Phase A. The spectral success criteria have 
#        been translated into the standard format understood by the 4FS_ETC tool. The original (human readable) 
#        definitions are given in the following document (referred to below as the SuccCritDoc):  
#          http://wiki.4most.eu/local--files/4fs-dqct-verification/DRS_SuccessCriteriaFINAL.pdf 
#        In cases where the SuccCritDoc is ambiguous I have referred to the IoA_DQCT python code: 
#          http://wiki.4most.eu/local--files/4fs-dqct-verification/update_exptimes_new.py
#        I have noted where the calculation in update_exptimes_new.py seems to differ 
#        from my interpretation of the SuccCritDoc, and where there are still ambiguities. 
#        Note that the spectral success criteria for GalHaloHR have been updated since PhaseA (see email 
#        exchange, subject "Halo HR survey FoM update" at end of July 2015).
#        See also: rulesetlist_derived_from_phaseA_SuccCritDoc.txt
#
#        History: 
#        rev3              - Copied from rulelist_derived_from_phaseA_SuccCritDoc.txt
#                          - Removed old rest-frame rules for AGN
#        rev4 (20.01.2016) - Updated DRS3 rules and rulesets according to info on IWG2 wiki pages
#        rev5 (11.02.2016) - Renamed "DRSn" to "Sn". 
#                          - Changed S4 rules slightly to require median(SNR) rather than min(SNR).  
#        rev6 (19.02.2016) - Added new rules for 1001MC. 
#        rev7 (02.03.2016) - Added new obs frame rule for Clusters  
#        rev8 (20.05.2016) - Treat GalHaloHR targets differently according to type (Dwarf or Giant)
#        rev9 (25.05.2016) - Change wavelength ranges in Rules for S3 GalDiskLR
#       rev10 (16.08.2016) - Change Rules for S4, move definition of SNR thresholds into Ruleset
#       rev11 (12.09.2016) - Added new Rules for S8
#
############################################################################################################


############################################################################################################
############################################################################################################
## S1 Low resolution Milky Way halo (GalHaloLR)
############################################################################################################
#RULE_NAME         VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN  LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
GalHaloLR_blue          SNR  MEDIAN       GE       10.0      5150.0      5200.0          AA    1.0          AA
GalHaloLR_red           SNR  MEDIAN       GE       10.0      8490.0      8675.0          AA    1.0          AA
#Notes:
#Notes: IoA_DQCT's wavelength range is 8582+-45AA - not 8490.0->8675.0AA
#Notes: SNR is measured per AA but in SuccCritDoc it is not specified.
#Notes: In an email from Mike Irwin, subject "Re: 4MOST - Spectra Quality Criteria for Round 6 (corrections)" 
#Notes: (date 09/01/2013) it is mentioned that the criterion is "median snr > 10 over the wavelength range",
#Notes: rather than the minimum SNR as given in the SuccCritDoc
#Notes:
############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################
## S2 High resolution Milky Way halo (GalHaloHR)
############################################################################################################
#RULE_NAME      VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN  LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
#GalHaloHR_blue      SNR    95PC       GT       50.0      426.1       427.0           NM   0.01          NM 
GalHaloHR_blue       SNR    95PC      DIV        1.0      426.1       427.0           NM   0.01          NM 
#Notes:
#Notes: The above is an updated criterion (May 2016), allowing different rulesets to have different SNR thresholds
#Notes:
############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################
## S3 Low resolution Milky Way disk/bulge (GalDiskLR)
############################################################################################################
#RULE_NAME       VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN  LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
GalDiskLR_CaII        SNR  MEDIAN      DIV        1.0       870.0       875.0          NM      1          AA 
GalDiskLR_Mgb         SNR  MEDIAN      DIV        1.0       570.0       580.0          NM      1          AA 
GalDiskLR_blue        SNR  MEDIAN      DIV        1.0       459.0       464.0          NM      1          AA 
#Notes: 
#Notes: Revised description received 24/05/2016
#Notes: For DRS3 can we ask also for a (hopefully straightforward) test?
#Notes: In our S/N calculations we have given initially quite large regions 
#Notes: where our lines of interest lie. We would like to change these to instead 
#Notes: contain the best nearby continuum regions, so that the lines don't affect 
#Notes: the S/N measurement, as detailed below:
#Notes: Closest continuum regions for the following spectral regions:
#Notes: GalDiskLR_CaII was 835-885 nm ==>  870 - 875
#Notes: GalDiskLR_Mgb was 514-520 nm  == > 570- 580
#Notes: GalDiskLR_blue was 450-470 nm  ==> 459 - 464
#Notes: 
#GalDiskLR_CaII        SNR  MEDIAN      DIV        1.0       835.0       885.0          NM      1          AA 
#GalDiskLR_Mgb         SNR  MEDIAN      DIV        1.0       514.0       520.0          NM      1          AA 
#GalDiskLR_blue        SNR  MEDIAN      DIV        1.0       450.0       470.0          NM      1          AA 
#Notes: Revised rules from IWG2 wiki page, updated 18/01/2015, adapted from description below
#Notes: Successful targets: 
#Notes: ESN:      MEDIAN SN/AA > 30 over range 835-885 nm (Ca II triplet region) and 514-520 nm (Mg triplet region), 
#Notes: DiskDyn:  MEDIAN SN/AA > 30 over range 835-885 nm, 
#Notes: DiskChem: MEDIAN SN/AA > 60 over range 450-470 nm and 514-520 nm, 
#Notes: Bulge:    MEDIAN SN/AA > 60 over range 450-470 nm + 514-520 nm
#Notes: 
############################################################################################################

############################################################################################################
############################################################################################################
## S4 High resolution Milky Way disk/bulge (GalDiskHR)
############################################################################################################
#RULE_NAME      VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN  LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
GalDiskHR_536NM      SNR  MEDIAN      DIV        1.0     535.4         536.1          NM      1          AA
GalDiskHR_545NM      SNR  MEDIAN      DIV        1.0     544.8         546.0          NM      1          AA
GalDiskHR_620NM      SNR  MEDIAN      DIV        1.0     619.0         621.0          NM      1          AA
#Notes:
#Notes: 
#Notes: Update 2016-08-15:
#Notes:   Sub-surveys 1 (bulge), 2 (inner disk), and 3 (outer disk):
#Notes:     median(SNR/AA) >= 100.0 over ranges 535.4 - 536.1 nm and 544.8 - 546.0 nm, and 619.0-621.0 nm
#Notes:   Sub-survey 4 (metallicity-complete local disk):
#Notes:     median(SNR/AA) >= 165.0 over ranges 535.4 - 536.1 nm and 544.8 - 546.0 nm, and 619.0-621.0 nm
#Notes:
#Notes:
#old# #RULE_NAME      VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN  LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
#old# GalDiskHR_blue       SNR  MEDIAN       GE       50.0     419.0         421.0          NM 0.0084          NM
#old# GalDiskHR_red        SNR  MEDIAN       GE       50.0     619.0         621.0          NM 0.0124          NM
#old# #Notes: Wavelength range is not specified, assuming measurements are made over +-10AA from nominal wavelength
#old# #Notes: SNR requirements in SuccCritDoc are: 
#old# #Notes:   "Minimum S/N per pixel requested at resolution 20000 and 2.5 pixel/resolution element"
#old# #Notes: Interpret this requirement by measuring SNR in linear wavelength bins: 
#old# #Notes:  * measure SNR per deltaL=420/(20000*2.5)=0.0084nm @ 420nm and 
#old# #Notes:  * measure SNR per deltaL=620/(20000*2.5)=0.0124nm @ 620nm  
#old# #Notes: IoA_DQCT's criterion for GalDiskHRdwarf and GalDiskHRgiant_MP is max(SNR/pix)>threshold over given 
#old# #Notes: wavelength range, but SuccCritDoc imples that criteria should be min(SNR/pix) >threshold
#old# #Notes:  
#old# #Notes: Slightly updated on 11.02.2016: changed metric from MIN to MEDIAN. 
#old# #Notes:  
############################################################################################################
############################################################################################################


############################################################################################################
############################################################################################################
## S5 X-ray galaxy clusters (Clusters)
############################################################################################################
#RULE_NAME      VARIABLE   METRIC OPERATOR COMP_VALUE  LAMBDA_MIN  LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
ClusEll_SNR_Specgen  SNR   MEDIAN      DIV        1.0      6000.0      6400.0     OBS_AA     1.0      OBS_AA
#Notes: 
#Notes: Clusters rule is now updated to new (March 2016) prescription, based on observer frame SNR
#Notes: Rule based on experiments with SpecGen->4FS_ETC->AutoZ
#Notes: 
#Notes: 
############################################################################################################

############################################################################################################
############################################################################################################
## S6 X-ray AGN (AGN)
############################################################################################################
#RULE_NAME    VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
AGN_SNR_RULEB      SNR  MEDIAN      DIV        1.0      4200.0     5000.0      OBS_AA    1.0      OBS_AA  
AGN_SNR_RULEG      SNR  MEDIAN      DIV        1.0      5500.0     6700.0      OBS_AA    1.0      OBS_AA  
AGN_SNR_RULER      SNR  MEDIAN      DIV        1.0      7200.0     9000.0      OBS_AA    1.0      OBS_AA  
#Notes: 
#Notes: New (Dec 2015) formulation based on BOSS/XMM-XXL measurements
#Notes:
#Notes: - 5% failure rate
#Notes: RULE1: MEDIAN SN/AA>1.1 over the range 4200-5000AA (5% case)
#Notes: RULE2: MEDIAN SN/AA>1.3 over the range 5500-6700AA (5% case)
#Notes: RULE3: MEDIAN SN/AA>1.5 over the range 7200-9000AA (5% case)
#Notes: 
#Notes: - 3% failure rate
#Notes: RULE1: MEDIAN SN/AA>2.1 over the range 4200-5000AA (3% case)
#Notes: RULE2: MEDIAN SN/AA>2.4 over the range 5500-6700AA (3% case)
#Notes: RULE3: MEDIAN SN/AA>2.8 over the range 7200-9000AA (3% case)
#Notes: 
############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################
## S7 Galaxy evolution (WAVES)
############################################################################################################
#RULE_NAME        VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
WAVES_RULE1       SNR       MEDIAN DIV      1.0         6000.0     6400.0     AA          1.0    PIXEL
#Notes:
#Notes: SNR 'per pixel' is deprecated, better to use 'per Angstrom' as it is independent of exact spectral format
#Notes: New formulation, (Dec 2015)
#Notes:
############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################
## S8 Cosmology (Cosmology)
############################################################################################################
#RULE_NAME        VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
COSMOBG_SNR         SNR    MEDIAN    DIV        1.0      5500.0     6700.0      OBS_AA    1.0      OBS_AA 
COSMOLRG_SNR        SNR    MEDIAN    DIV        1.0      6800.0     7300.0      OBS_AA    1.0      OBS_AA
COSMOELG_SNR        SNR    MEDIAN    DIV        1.0      5500.0     6700.0      OBS_AA    1.0      OBS_AA
COSMOQSO_SNR_B      SNR    MEDIAN    DIV        1.0      4200.0     5000.0      OBS_AA    1.0      OBS_AA
COSMOQSO_SNR_G      SNR    MEDIAN    DIV        1.0      5500.0     6700.0      OBS_AA    1.0      OBS_AA
COSMOQSO_SNR_R      SNR    MEDIAN    DIV        1.0      7200.0     9000.0      OBS_AA    1.0      OBS_AA
#Notes:
#Notes: As specified in Johan Richard's submission from 10/09/2016 
#Notes:
#Notes:
#old# Cosmo_ELG_OII          SNR    MEAN       GE        6.0      3726.0     3730.0     REST_AA    1.0      OBS_AA 
#old# Cosmo_ELG_OIII         SNR    MEAN       GE        6.0      5004.0     5010.0     REST_AA    1.0      OBS_AA 
#old# #Notes:
#old# #Notes: IoA_DQCT's wavelength range is 3728AArest +-5 AAobs, not 3726.0->3730.0 AArest as given in SuccCritDoc
#old# #Notes: IoA_DQCT's wavelength range is 5008AArest +-5 AAobs, not 5004.0->5010.0 AArest as given in SuccCritDoc
#old# #Notes: Calculated deltal is assumed to be per observed frame AA
#old# #Notes: IoA_DQCT's calculated measure is MEDIAN but should be MEAN as given in SuccCritDoc
#old# #Notes: Use Clusters rules for Cosmology Ellipticals 
#old# #Notes:
############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################
## ??S9?? Magellanic Clouds (1001MC)
############################################################################################################
#RULE_NAME        VARIABLE  METRIC OPERATOR COMP_VALUE  LAMBDA_MIN LAMBDA_MAX LAMBDA_UNIT DELTAL DELTAL_TYPE
MC_LR_Mg              SNR   MEDIAN      DIV        1.0       514.0      520.0      OBS_NM    1.0      OBS_AA 
MC_LR_Ca              SNR   MEDIAN      DIV        1.0       835.0      885.0      OBS_NM    1.0      OBS_AA 
MC_HR_Mg              SNR   MEDIAN      DIV        1.0       514.0      520.0      OBS_NM    1.0      OBS_AA 
#Notes:
#Notes: LR mode: S/N>30 per AA in the ranges of the Ca (835-885 nm) and Mg (514-520 nm)  triplets.
#Notes: for HR please consider S/N>50 in the range 514-520 nm. That's again a starting point.
#Notes:
############################################################################################################
############################################################################################################

"""


def ruleset(snr_list, snr_definitions):
    output = """
#NAME       REQUIRED_VALUE EXPRESSION
goodSNR250c 250.0          DCFSNR
"""

    for snr_definition in snr_definitions:
        for snr in snr_list:
            output += "goodSNR{0:.1f}_{1:s} {2:13.1f} {1:s}\n".format(snr, snr_definition, snr)

    return output
