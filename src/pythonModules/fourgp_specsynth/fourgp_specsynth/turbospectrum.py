#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A class which wraps Turbospectrum
"""

import numpy as np
import subprocess
import os
from os import path as os_path
import re
import glob
import logging

from solar_abundances import solar_abundances

logger = logging.getLogger(__name__)


class TurboSpectrum:
    def __init__(self,
                 turbospec_path="/home/dcf21/iwg7_pipeline/turbospectrum-15.1/exec-gf-v15.1",
                 interpol_path="/home/dcf21/iwg7_pipeline/interpol_marcs",
                 line_list_paths="/home/dcf21/iwg7_pipeline/fromBengt/line_lists/3700-9500",
                 marcs_grid_path="/home/dcf21/iwg7_pipeline/fromBengt/marcs_grid"):
        """
        Instantiate a class for generating synthetic stellar spectra using Turbospectrum.

        :param turbospec_path:
            Path where the turbospectrum binaries 'babsma' and 'bsyn' can be found.

        :type turbospec_path:
            str

        :param interpol_path:
            Path where the compiled interpol_modeles.f binary can be found.

        :type interpol_path:
            str

        :param line_list_paths:
            Path(s) where line lists for synthetic spectra can be found. Specify as either a string, or a list of
            strings.

        :type line_list_paths:
            list or str

        :param marcs_grid_path:
            Path where a grid of MARCS .mod files can be found. These contain the model atmospheres we use.

        :type marcs_grid_path:
            str
        """

        if not isinstance(line_list_paths, (list, tuple)):
            line_list_paths = [line_list_paths]

        self.turbospec_path = turbospec_path
        self.interpol_path = interpol_path
        self.line_list_paths = line_list_paths
        self.marcs_grid_path = marcs_grid_path

        # Default spectrum parameters
        self.lambda_min = 5100
        self.lambda_max = 5200
        self.lambda_delta = 0.05
        self.metallicity = -1.5
        self.turbulent_velocity = 1.0
        self.free_abundances = None
        self.sphere = None
        self.alpha = None
        self.s_process = 0
        self.r_process = 0
        self.star_name = "anonymous_star"
        self.verbose = True
        self.line_list_files = None
        self.spec_file = "test1.spec"

        # Create temporary directory
        self.id_string = "turbospec_%d" % os.getpid()
        self.tmp_dir = os_path.join("/tmp", self.id_string)
        os.system("mkdir -p %s" % self.tmp_dir)

        # Look up what MARCS models we have
        self.marcs_values = None
        self.marcs_value_keys = []
        self.marcs_models = {}
        self._fetch_marcs_grid()

    def _fetch_marcs_grid(self):
        """
        Get a list of all of the MARCS models we have.

        :return:
            None
        """

        pattern = r"([sp])(\d\d\d\d)_g(....)_m(...)_t(..)_(..)_z(.....)_" \
                  r"a(.....)_c(.....)_n(.....)_o(.....)_r(.....)_s(.....).mod"

        self.marcs_values = {
            "spherical": [], "temperature": [], "log_g": [], "mass": [], "turbulence": [], "model_type": [],
            "metallicity": [], "a": [], "c": [], "n": [], "o": [], "r": [], "s": []
        }

        self.marcs_value_keys = self.marcs_values.keys()
        self.marcs_value_keys.sort()
        self.marcs_models = {}

        marcs_models = glob.glob(os_path.join(self.marcs_grid_path,"*"))
        for item in marcs_models:

            # Extract model parameters from .mod filename
            filename = os_path.split(item)[1]
            re_test = re.match(pattern, filename)
            assert re_test is not None, "Could not parse MARCS model filename <{}>".format(filename)

            try:
                model = {
                    "spherical": re_test.group(1),
                    "temperature": float(re_test.group(2)),
                    "log_g": float(re_test.group(3)),
                    "mass": float(re_test.group(4)),
                    "turbulence": float(re_test.group(5)),
                    "model_type": re_test.group(6),
                    "metallicity": float(re_test.group(7)),
                    "a": float(re_test.group(8)),
                    "c": float(re_test.group(9)),
                    "n": float(re_test.group(10)),
                    "o": float(re_test.group(11)),
                    "r": float(re_test.group(12)),
                    "s": float(re_test.group(13))
                }
            except ValueError:
                logger.error("Could not parse MARCS model filename <{}>".format(filename))
                raise

            # Keep a list of all of the parameter values we've seen
            for parameter, value in model.iteritems():
                if value not in self.marcs_values[parameter]:
                    self.marcs_values[parameter].append(value)

            # Keep a list of all the models we've got in the grid
            dict_iter = self.marcs_models
            for parameter in self.marcs_value_keys:
                value = model[parameter]
                if value not in dict_iter:
                    dict_iter[value] = {}
                dict_iter = dict_iter[value]
            dict_iter["filename"] = item

        # Sort model parameter values into order
        for parameter in self.marcs_value_keys:
            self.marcs_values[parameter].sort()

    def configure(self, lambda_min=None, lambda_max=None, lambda_delta=None, metallicity=None,
                  turbulent_velocity=None, free_abundances=None,
                  sphere=None, spec_file=None, alpha=None, s_process=None, r_process=None, star_name=None,
                  line_list_paths=None, line_list_files=None,
                  verbose=True):
        """
        Set the stellar parameters of the synthetic spectra to generate. This can be called as often as needed
        to generate many synthetic spectra with one class instance. All arguments are optional; any which are not
        supplied will remain unchanged.

        :param lambda_min:
            Short wavelength limit of the synthetic spectra we generate. Unit: A.
        :param lambda_max:
            Long wavelength limit of the synthetic spectra we generate. Unit: A.
        :param lambda_delta:
            Wavelength step of the synthetic spectra we generate. Unit: A.
        :param metallicity:
            Metallicity of the star we're synthesising.
        :param turbulent_velocity:
        :param free_abundances:
            List of elemental abundances to use in stellar model. These are passed to Turbospectrum.
        :param sphere:
            Select whether to use a spherical model (True) or a plane-parallel model (False).
        :param spec_file:
        :param alpha:
            Alpha enhancement to use in stellar model.
        :param s_process:
            S-Process element enhancement to use in stellar model.
        :param r_process:
            R-Process element enhancement to use in stellar model.
        :param star_name:
            The name of the star being modelled.
        :param line_list_paths:
            List of paths where we should search for line lists.
        :param line_list_files:
            List of line list files to use. If not specified, we use all files in `line_list_paths`
        :param verbose:
            Let Turbospectrum print debugging information to terminal?
        :return:
            None
        """

        if lambda_min is not None:
            self.lambda_min = lambda_min
        if lambda_max is not None:
            self.lambda_max = lambda_max
        if lambda_delta is not None:
            self.lambda_delta = lambda_delta
        if metallicity is not None:
            self.metallicity = metallicity
        if star_name is not None:
            self.star_name = star_name
        if turbulent_velocity is not None:
            self.turbulent_velocity = turbulent_velocity
        if free_abundances is not None:
            self.free_abundances = free_abundances
        if sphere is not None:
            self.sphere = sphere
        if spec_file is not None:
            self.spec_file = spec_file
        if alpha is not None:
            self.alpha = alpha
        if s_process is not None:
            self.s_process = s_process
        if r_process is not None:
            self.r_process = r_process
        if line_list_paths is not None:
            if not isinstance(line_list_paths, (list, tuple)):
                line_list_paths = [line_list_paths]
            self.line_list_paths = line_list_paths
        if line_list_files is not None:
            self.line_list_files = line_list_files
        if verbose is not None:
            self.verbose = verbose

        # ------ assume an alpha enhancement ----
        if self.alpha is None:
            if self.metallicity < -1.0:
                self.alpha = 0.4
            elif -1.0 < self.metallicity < 0.0:
                self.alpha = -0.4 * self.metallicity
            else:
                self.alpha = 0

    def _generate_model_atmosphere(self, t_eff, g, f):
        """
        Generates an interpolated model atmosphere from the MARCS grid using the interpol.f routine developed by
        T. Masseron (Masseron, PhD Thesis, 2006). This is a python wrapper for that fortran code.

        :param t_eff:
            Effective temperature of requested model atmosphere.
        :param g:
            Log(g) of requested model atmosphere.
        :param f:
            [Fe/H] metallicity of requested model atmosphere.

        """

        if self.verbose:
            stdout = subprocess.STDOUT
            stderr = subprocess.STDOUT
        else:
            stdout = open('/dev/null', 'w')
            stderr = subprocess.STDOUT

        # ----- defines the point at which Plane-Parallel vs spherical model atmosphere models are used
        if g >= 3.0:
            if self.star_name is not None:
                self.marcs_model_name = '%s_%ig%.2fm0.0z%.2f.int' % (self.star_name, t_eff, g, f)
            else:
                self.marcs_model_name = '%ig%.2fm0.0z%.2f.int' % (t_eff, g, f)

            output = os_path.join(self.mod_path, self.marcs_model_name)

            try:
                p = subprocess.Popen(
                    [self.interp_path + 'interpol_planparallel_in.com', str(t_eff), str(g), str(f), output],
                    stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
                stdout, stderr = p.communicate()
            except subprocess.CalledProcessError:
                raise RuntimeError('Plane-Parallel Model atmosphere interpolation failed ....')
            self.sphere = 'F'

        else:
            if self.star_name is not None:
                self.marcs_model_name = '%s_%ig%.2fm1.0z%.2f.int' % (self.star_name, t_eff, g, f)
            else:
                self.marcs_model_name = '%ig%.2fm1.0z%.2f.int' % (t_eff, g, f)

            output = os_path.join(self.mod_path, self.marcs_model_name)

            try:
                p = subprocess.Popen(
                    [self.interp_path + 'interpol_spherical_in.com', str(t_eff), str(g), str(f), output],
                    stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
                stdout, stderr = p.communicate()
            except subprocess.CalledProcessError:
                raise RuntimeError('Spherical Model atmosphere interpolation failed ....')
            self.sphere = 'T'

        return

    def make_babsmabysn_file(self):
        '''
        --------------------------------------------------------------------------------
        | PURPOSE:
        | Generate the parameter file for both the babsma and bsyn codes for turbospectrum
        |
        --------------------------------------------------------------------------------
        '''

        # Checks that model metallicity and input metallicity are consistent
        modelmetal = float(self.marcs_model_name.split('z')[1].split('.int')[0])
        if (modelmetal > self.metallicity + 0.05) or (modelmetal < self.metallicity - 0.05):
            logger.warn("Atmosphere model ({:.2f}) and input metallicity ({:.2f}) not consistent; proceed with caution".
                        format(modelmetal, self.metallicity))

        # Allow for user input abundances as a dictionary of the form {element: abundance}
        if self.free_abundances is None:
            individual_abundances = "'INDIVIDUAL ABUNDANCES:'   '0'\n"
        else:
            individual_abundances = "'INDIVIDUAL ABUNDANCES:'   '{:d}'\n".format(len(self.free_abundances))

            for element, abundance in self.free_abundances.iteritems():
                assert element in solar_abundances, "Cannot proceed as solar abundance for element <{}> is unknown". \
                    format(element)

                individual_abundances += "%i  {:.2f}\n" % (
                    int(z[elind]), float(solar_abundances[element]) + float(abundance))

        # Make a list of line-list files
        line_lists = "'NFILES   :' '{:d}'\n".format(len(self.line_list_files))
        for item in self.line_list_files:
            line_lists += "%s\n".format(item)

        # Build bysn configuration file
        bsyn_config = """\
'LAMBDA_MIN:'    '{this[lambda_min]:.3f}'
'LAMBDA_MAX:'    '{this[lambda_max]:.3f}'
'LAMBDA_STEP:'    '{this[lambda_delta]:.3f}'
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
'ABFIND        :' '.false.'
'MODELOPAC:' '{this[tmp_dir]}/model_opacity.opac'
'RESULTFILE :' '{this[tmp_dir]}/{this[spec_file]}.spec'
'METALLICITY:'    '{this[metallicity]:.2f}'
'ALPHA/Fe   :'    '{this[alpha]:.2f}'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '{this[r_process]:.2f}'
'S-PROCESS  :'    '{this[s_process]:.2f}'
{individual_abundances}
'ISOTOPES : ' '2'
'3.006  0.075
'3.007  0.925
{line_lists}
'SPHERICAL:'  '{this[sphere]}'
  30
  300.00
  15
  1.30
""".format(this=self.__dict__, individual_abundances=individual_abundances, line_lists=line_lists)

        # Build babsma configuration file
        babsma_config = """\
'LAMBDA_MIN:'    '{this[lambda_min]:.3f}'
'LAMBDA_MAX:'    '{this[lambda_max]:.3f}'
'LAMBDA_STEP:'    '{this[lambda_delta]:.3f}'
'MODELINPUT:' '{this[tmp_dir]}/'
'MARCS-FILE:' '.false.'
'MODELOPAC:' '{this[tmp_dir]}/model_opacity.opac'
'METALLICITY:'    '{this[metallicity]:.2f}'
'ALPHA/Fe   :'    '{this[alpha]:.2f}'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '{this[r_process]:.2f}'
'S-PROCESS  :'    '{this[s_process]:.2f}'
{individual_abundances}
'XIFIX:' 'T'
{this[turbulent_velocity]:.2f}
""".format(this=self.__dict__, individual_abundances=individual_abundances)

        return babsma_config, bsyn_config

    def synthesise(self, stellarpar=None):

        if stellarpar is None:
            babsma_in, bsyn_in = self.make_babsmabysn_file()
        else:
            if len(stellarpar) != 3:
                raise ValueError('The stellarpar varible must be a list/array of length 3 = (Teff, logg, [Fe/H])')
            else:
                T, g, f = stellarpar
                logger.info('Generating model atmosphere with T=%.1f, log g = %.2f, metallicity = %.2f' % (T, g, f))
                self._generate_model_atmosphere(T, g, f)
                babsma_in, bsyn_in = self.make_babsmabysn_file()

        if self.verbose:
            stdout = None
            stderr = subprocess.STDOUT
        else:
            stdout = open('/dev/null', 'w')
            stderr = subprocess.STDOUT

        try:
            pr1 = subprocess.Popen([os_path.join(self.turbospec_path, 'babsma_lu')],
                                   stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
            for line in babsma_in.splitlines():
                pr1.stdin.write(line)
            stdout, stderr = pr1.communicate()
        except subprocess.CalledProcessError:
            raise RuntimeError('babsma failed ....')
        logger.info("%s %s" % (pr1.returncode, stderr))

        if self.verbose:
            stdout = None
            stderr = subprocess.STDOUT
        else:
            stdout = open('/dev/null', 'w')
            stderr = None
        try:
            pr = subprocess.Popen([os_path.join(self.turbospec_path, 'bsyn_lu')],
                                  stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
            for line in bsyn_in:
                pr.stdin.write(line)
            stdout, stderr = pr.communicate()
        except subprocess.CalledProcessError:
            raise RuntimeError('babsma failed ....')
        logger.info("%s %s" % (pr.returncode, stderr))
        return pr.returncode
