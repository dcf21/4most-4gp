#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This python package provides a wrapper for functions provided by the 4MOST Facility Simulator, 4FS.
"""

import os
from os import path as os_path
import numpy as np
import astropy.io.fits as fits
import logging

from fourgp_speclib import Spectrum

import config_files

logger = logging.getLogger(__name__)


class FourFS:
    def __init__(self,
                 path_to_4fs="/home/dcf21/iwg7_pipeline/OpSys/ETC",
                 snr_list=None
                 ):
        """
        Instantiate a class for calling 4FS.

        :param path_to_4fs:
            Path to where 4FS binaries can be found.

        :param snr_list:
            List of the SNRs that we want 4FS to degrade input spectra to.
        """
        self.path_to_4fs = path_to_4fs
        self.template_counter = 0
        self.metadata_store = {}

        if snr_list is None:
            snr_list = ("05", "10", "15", "20", "50", "100", "250")
        self.snr_list = snr_list

        # Create temporary directory
        self.id_string = "4fs_{:d}".format(os.getpid())
        self.tmp_dir = os_path.join("/tmp", self.id_string)
        os.system("mkdir -p {}".format(self.tmp_dir))

        # Put standard 4FS configuration files into temporary working directory
        with open(os_path.join(self.tmp_dir, "ETC_input_params_HRS.txt"), "w") as f:
            f.write(config_files.ETC_input_params_HRS)

        with open(os_path.join(self.tmp_dir, "ETC_input_params_LRS.txt"), "w") as f:
            f.write(config_files.ETC_input_params_LRS)

        with open(os_path.join(self.tmp_dir, "rulelist.txt"), "w") as f:
            f.write(config_files.rulelist)

        with open(os_path.join(self.tmp_dir, "ruleset.txt"), "w") as f:
            f.write(config_files.ruleset(snr_list=self.snr_list))

        # Extract 4MOST telescope description files
        cwd = os.getcwd()
        os.chdir(self.tmp_dir)
        path_to_telescope_description = os_path.split(os_path.abspath(__file__))[0]
        os.system("tar xvfz {}".format(os_path.join(path_to_telescope_description,
                                                    "4FS_ETC_system_model_latest.tar.gz")))
        os.chdir(cwd)

    def close(self):
        # Remove temporary directory
        os.system("rm -Rf {}".format(self.tmp_dir))

    def make_fits_spectral_template(self,
                                    input_spectrum, input_spectrum_continuum_normalised,
                                    output_filename='test1.fits',
                                    resolution=50000,
                                    continuum_only=False
                                    ):
        """
        Generate a 4FS readable fits template from the Turbospectrum output.
        
        :param input_spectrum:
            Spectrum object that we want to pass through 4FS.
            
        :type input_spectrum:
            Spectrum
            
        :param input_spectrum_continuum_normalised:
            Continuum-normalised version of the Spectrum object that we want to pass through 4FS.
            
        :type input_spectrum_continuum_normalised:
            Spectrum
            
        :param output_filename:
            Output filename for FITS file we store in our temporary workspace.
            
        :type output_filename:
            str
            
        :param resolution:
            The spectral resolution of the input spectrum, stored in the FITS headers
            
        :type resolution:
            float
        
        :param continuum_only:
            Select whether to output full spectrum, including spectral lines, or just the continuum outline.
            
        :type continuum_only:
            bool
            
        :return:
            None
        """

        assert input_spectrum.raster_hash == input_spectrum_continuum_normalised.raster_hash, \
            "Continuum-normalised spectrum needs to have the same wavelength raster as the original."

        # Extract data from input spectra
        wavelength_raster = input_spectrum.wavelengths
        flux = input_spectrum.values
        continuum_normalised_flux = input_spectrum_continuum_normalised.values

        magnitude = 12.0
        lambda_min = np.min(wavelength_raster)
        delta_lambda = wavelength_raster[1] - wavelength_raster[0]

        if not continuum_only:
            # Apply some normalisation to input flux levels
            hdu_1 = fits.PrimaryHDU(flux * 10E-15 / max(flux))
        else:
            # If we only want continuum, without lines, then divide by continuum_normalised flux
            hdu_1 = fits.PrimaryHDU((flux * 10E-15 / max(flux)) / continuum_normalised_flux)

        hdu_1.header['CRTYPE1'] = "LINEAR  "
        hdu_1.header['CRPIX1'] = 1.0
        hdu_1.header['CRVAL1'] = lambda_min
        hdu_1.header['CDELT1'] = delta_lambda
        hdu_1.header['CUNIT1'] = "Angstrom"
        hdu_1.header['BUNIT'] = "erg/s/cm2/Angstrom"
        hdu_1.header['ABMAG'] = magnitude
        if resolution is not None:
            fwhm_res = np.mean(wavelength_raster / resolution)
            hdu_1.header['RESOLUTN'] = fwhm_res

        hdu_list = fits.HDUList([hdu_1])
        hdu_list.writeto(os_path.join(self.tmp_dir, output_filename))

    def generate_4fs_template_list(self,
                                   spectra_list,
                                   resolution=50000
                                   ):
        """
        Generate 4FS template list from a list of 4GP Spectrum objects

        :param spectra_list:
            A list of the spectra we should pass through 4FS. Each entry in the list should be a list of tuple with
            two entries: (input_spectrum, input_spectrum_continuum_normalised). These reflect the contents of the
            third and second columns of Turbospectrum's ASCII output respectively.

        :type spectra_list:
            (list, tuple) of (list, tuple) of Spectrum objects

        :param resolution:
            The spectral resolution of the input spectrum, stored in the FITS headers

        :type resolution:
            float

        :return:
            The string contents of the template file we wrote for 4FS.
        """

        # Location to write template list to be used for 4FS
        template_list_filename = 'template_list_test.txt'

        # Start writing template list file
        writestr = """\
#OBJECTNAME           FILENAME                                                 RULESET SIZE REDSHIFT MAG  MAG_RANGE
"""

        # Loop over stars we are to process
        star_list = []
        template_number_first = self.template_counter
        for spectra in spectra_list:
            input_spectrum, input_spectrum_continuum_normalised = spectra
            obj_name = input_spectrum.metadata['Starname']
            logger.info("Working on <{}>".format(obj_name))

            # Keep a copy of the metadata on this object, which we put back into the output from 4FS
            self.metadata_store[self.template_counter] = input_spectrum.metadata.copy()

            # Generates the 4FS readable fits files from the output spectra from Turbospectrum
            self.make_fits_spectral_template(input_spectrum=input_spectrum,
                                             input_spectrum_continuum_normalised=input_spectrum_continuum_normalised,
                                             output_filename='template_{}.fits'.format(self.template_counter),
                                             resolution=(resolution is not None))
            self.make_fits_spectral_template(input_spectrum=input_spectrum,
                                             input_spectrum_continuum_normalised=input_spectrum_continuum_normalised,
                                             output_filename='template_{}_c.fits'.format(self.template_counter),
                                             resolution=(resolution is not None),
                                             continuum_only=True)

            for snr in self.snr_list:
                writestr += 'template_{}_SNR{} {} {} {} {} {} {}\n'.format(
                    self.template_counter,
                    snr,
                    os_path.join(self.tmp_dir, 'template_{}.fits'.format(self.template_counter)),
                    'goodSNR{}'.format(snr),
                    '0.0', '0.0', '15.0', '15.0')
            writestr += 'template_{}_c {} {} {} {} {} {}\n'.format(
                self.template_counter,
                os_path.join(self.tmp_dir, 'template_{}.fits'.format(self.template_counter)),
                'goodSNR250', '0.0', '0.0', '15.0', '15.0')

            for snr in self.snr_list:
                star_list.append('template_{}_SNR{}'.format(self.template_counter, snr))
            star_list.append('template_{}_c'.format(self.template_counter))

            # Increment template counter
            self.template_counter += 1

        # Write template list to file
        with open(os_path.join(self.tmp_dir, template_list_filename), 'w') as f:
            f.write(writestr)

        # Return the string contents of the template list
        template_number_last = self.template_counter - 1

        return {'template_file': writestr,
                'template_numbers': range(template_number_first, template_number_last + 1)
                }

    def combine_spectra(self,
                        template_numbers,
                        path='outdir_LRS/',
                        setup='LRS'
                        ):
        """
        4FS produces three output fits files for both 4MOST's LRS and HRS modes. These contain the three spectra
        arms. However, the Cannon expects all its data in a single unified spectrum. This function stitches the
        three arms into a single data set.

        Combine the multiple fits files from the 4FS into a single ascii file with all arms included.

        NOTE: The point at which one are is favoured over another in the stitching
              process is hardcoded should be changed in future versions.

        :param template_numbers:
            Each spectrum is allocated an ID number by the function `generate_4FS_template_list` above when it is
            prepared for input into 4FS. These ID numbers are also present in the output filenames of the observed
            spectra that 4FS produces. This parameter is a list of the IDs of the spectra we are to post-process.

        :type template_numbers:
            List or tuple of ints

        :param path:
            Set the path to the 4FS output FITS files we are to post-process. This path is relative to `self.tmp_dir`.

        :type path:
            str

        :param setup:
            Set to either (1) 'LRS' - low res or (2) 'HRS' - high res. Specifies the mode in which 4MOST is working.

        :type setup:
            str

        :return :
            A list of 4GP Spectrum objects.
        """

        output = {}
        for i in template_numbers:
            output[i] = {}
            for snr in self.snr_list:
                setup_full = "SNR{}_{}".format(snr, setup)
                # Load in the three output spectra -- blue, green, red arms
                d = fits.open(os_path.join(path, 'specout_template_template_{}_{}_blue.fits'.format(i, setup_full)))
                d2 = fits.open(os_path.join(path, 'specout_template_template_{}_{}_green.fits'.format(i, setup_full)))
                d3 = fits.open(os_path.join(path, 'specout_template_template_{}_{}_red.fits'.format(i, setup_full)))

                # Read the data from the FITS files
                wavelengths_1 = d[2].data['LAMBDA']
                wavelengths_2 = d2[2].data['LAMBDA']
                wavelengths_3 = d3[2].data['LAMBDA']
                fluxes_1 = d[2].data['REALISATION'] - d[2].data['SKY']
                fluxes_2 = d2[2].data['REALISATION'] - d2[2].data['SKY']
                fluxes_3 = d3[2].data['REALISATION'] - d3[2].data['SKY']
                snr_1 = d[2].data['SNR']
                snr_2 = d2[2].data['SNR']
                snr_3 = d3[2].data['SNR']

                if setup == 'LRS':
                    # Stitch arms based on where SNR is best -- hardcoded here may need changed in future versions
                    indices_1 = np.where(wavelengths_1 <= 5327.7)[0]
                    indices_2 = np.where((wavelengths_2 > 5327.7) & (wavelengths_2 <= 7031.7))[0]
                    indices_3 = np.where(wavelengths_3 > 7031.7)[0]
                    wavelengths_1 = wavelengths_1[indices_1]
                    wavelengths_2 = wavelengths_2[indices_2]
                    wavelengths_3 = wavelengths_3[indices_3]
                    fluxes_1 = fluxes_1[indices_1]
                    fluxes_2 = fluxes_2[indices_2]
                    fluxes_3 = fluxes_3[indices_3]
                    snr_1 = snr_1[indices_1]
                    snr_2 = snr_2[indices_2]
                    snr_3 = snr_3[indices_3]

                wavelengths_final = np.array(list(wavelengths_1) + list(wavelengths_2) + list(wavelengths_3))
                fluxes_final = np.array(list(fluxes_1) + list(fluxes_2) + list(fluxes_3))
                snr_final = np.array(list(snr_1) + list(snr_2) + list(snr_3))

                # Load continuum spectra
                dc = fits.open(os_path.join(path, 'specout_template_template_{}_c_{}_blue.fits'.format(i, setup)))
                d2c = fits.open(os_path.join(path, 'specout_template_template_{}_c_{}_green.fits'.format(i, setup)))
                d3c = fits.open(os_path.join(path, 'specout_template_template_{}_c_{}_red.fits'.format(i, setup)))

                # Read the data from the FITS files
                wavelengths_1c = dc[2].data['LAMBDA']
                wavelengths_2c = d2c[2].data['LAMBDA']
                wavelengths_3c = d3c[2].data['LAMBDA']
                fluxes_1c = dc[2].data['FLUENCE']
                fluxes_2c = d2c[2].data['FLUENCE']
                fluxes_3c = d3c[2].data['FLUENCE']

                if setup == 'LRS':
                    indices_1 = np.where(wavelengths_1c <= 5327.7)[0]
                    indices_2 = np.where((wavelengths_2c > 5327.7) & (wavelengths_2c <= 7031.7))[0]
                    indices_3 = np.where(wavelengths_3c > 7031.7)[0]
                    wavelengths_1c = wavelengths_1c[indices_1]
                    wavelengths_2c = wavelengths_2c[indices_2]
                    wavelengths_3c = wavelengths_3c[indices_3]
                    fluxes_1c = fluxes_1c[indices_1]
                    fluxes_2c = fluxes_2c[indices_2]
                    fluxes_3c = fluxes_3c[indices_3]

                normalised_fluxes_1 = fluxes_1 / (fluxes_1c * max(d[2].data['FLUENCE']) / max(fluxes_1c))
                normalised_fluxes_2 = fluxes_2 / (fluxes_2c * max(d2[2].data['FLUENCE']) / max(fluxes_2c))
                normalised_fluxes_3 = fluxes_3 / (fluxes_3c * max(d3[2].data['FLUENCE']) / max(fluxes_3c))

                # Combine everything into one set of arrays to be saved
                wavelengths_final_c = np.array(list(wavelengths_1c) + list(wavelengths_2c) + list(wavelengths_3c))

                fluxes_final_c = np.array(
                    list(fluxes_1c * max(d[2].data['FLUENCE']) / max(fluxes_1c)) +
                    list(fluxes_2c * max(d2[2].data['FLUENCE']) / max(fluxes_2c)) +
                    list(fluxes_3c * max(d3[2].data['FLUENCE']) / max(fluxes_3c)))

                normalised_fluxes_final = np.array(list(normalised_fluxes_1) +
                                                   list(normalised_fluxes_2) +
                                                   list(normalised_fluxes_3))
                uncertainty_fluxes_final = normalised_fluxes_final / snr_final

                # Do continuum normalisation
                normalised_fluxes_final = fluxes_final / fluxes_final_c

                # Remove bad pixels
                # Any pixels where flux > 2 or flux < 0 get reset to zero for other downstream codes
                normalised_fluxes_final[
                    np.where((normalised_fluxes_final > 2.0) | (normalised_fluxes_final <= 0.00))[0]] = 0

                # Turn data into a 4GP Spectrum object
                metadata = self.metadata_store[i]
                metadata['continuum_normalised'] = 0
                metadata['SNR'] = float(snr)
                spectrum = Spectrum(wavelengths=wavelengths_final,
                                    values=fluxes_final,
                                    value_errors=fluxes_final / snr_final,
                                    metadata=metadata)

                metadata['continuum_normalised'] = 1
                spectrum_continuum_normalised = Spectrum(wavelengths=wavelengths_final,
                                                         values=normalised_fluxes_final,
                                                         value_errors=uncertainty_fluxes_final,
                                                         metadata=metadata)

                output[i][snr] = {
                    "spectrum": spectrum,
                    "spectrum_continuum_normalised": spectrum_continuum_normalised
                }

        return output

    def process_spectra(self, spectra_list, resolution=50000):
        # Generate input configuration files for 4FS
        template_information = self.generate_4fs_template_list(spectra_list=spectra_list,
                                                               resolution=resolution)

        # Find the 4FS binary
        fourfs_command = os_path.join(self.path_to_4fs, "4FS_ETC")
        assert os_path.exists(fourfs_command), "Could not find 4FS binary!"

        # Switching into our temporary working directory where 4FS can find all its input files
        cwd = os.getcwd()
        os.chdir(self.tmp_dir)

        # Make sure there aren't any old 4FS outputs lying around
        os.system("rm -Rf outdir_LRS outdir_HRS")

        # Run 4FS
        os.system("{} PARAM_FILENAME=ETC_input_params_LRS.txt > lrs_output.log 2>&1".format(fourfs_command))
        os.system("{} PARAM_FILENAME=ETC_input_params_HRS.txt > hrs_output.log 2>&1".format(fourfs_command))

        # Stitched 4MOST wavebands together
        output = {}
        for mode in ("LRS", "HRS"):
            stitched_spectra = self.combine_spectra(template_numbers=template_information['template_numbers'],
                                                    path="outdir_{}".format(mode),
                                                    setup=mode
                                                    )
            output[mode] = stitched_spectra

        # Switch back into the user's cwd
        os.chdir(cwd)

        # Return output spectra to user
        return output
