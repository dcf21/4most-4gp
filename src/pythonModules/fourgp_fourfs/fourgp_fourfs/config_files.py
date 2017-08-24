#!/usr/bin/env python
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
DCFSNR        SNR MEDIAN DIV   1.0  618.0 668.0 NM     1.0    PIX
"""


def ruleset(snr_list, snr_definitions):
    output = """
#NAME       REQUIRED_VALUE EXPRESSION
goodSNR2500 2500.0         DCFSNR
"""

    for snr_definition in snr_definitions:
        for snr in snr_list:
            output += "goodSNR{0:.1f}_{1:s} {2:13.1f} {1:s}\n".format(snr, snr_definition, snr)

    return output
