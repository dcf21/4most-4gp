#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A class which wraps Turbospectrum
"""

from numpy import *
import subprocess
import os
from os import path as os_path
import glob
import warnings
import logging

logger = logging.getLogger(__name__)


class TurboSpectrum:
    def __init__(self, turbospec_path, interpol_path, line_list_paths, marcs_grid_path):
        """

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


def set_config(self, lambda_min=None, lambda_max=None, lambda_delta=None, metallicity=None,
                 turbulent_velocity=None, free_abund=None,
                 sphere=None, spec_file=None, alpha=None, s_process=None, r_process=None, star_name=None,
                 line_list_paths=None, line_list_files=None,
                 verbose=True):

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
        if free_abund is not None:
            self.free_abundances = free_abund
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

        # ------ assume an alpha enhancment ----
        if self.alpha is None:
            if self.metallicity < -1.0:
                self.alpha = 0.4
            elif -1.0 < self.metallicity < 0.0:
                self.alpha = -0.4 * self.metallicity
            else:
                self.alpha = 0

    def generate_model_atmosphere(self, t_eff, g, f):
        """
        Generates an interpolated model atmosphere from the MARCs grid using the interpol.f routine developed by
        T. Masseron (Masseron, PhD Thesis, 2006). This is a python wrapper for that fortran code.

        INPUT :
        (1) Temperature, (2) Surface gravity, (3) Metallicity of the model atmospshere

        OUTPUT :
        full path of model atmosphere and the type of interpolation done: Plane-parallel (P), or Spherical (S)

        OPTIONAL:
        modpath = the path to write the model atmosphere

        """

        if self.verbose:
            stdout = open('/dev/null', 'w')
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

        # -----checks that model metallicity and input metallicity are consistent
        modelmetal = float(self.marcs_model_name.split('z')[1].split('.int')[0])
        if modelmetal > self.metallicity + 0.05 or modelmetal < self.metallicity - 0.05:
            warnings.warn(
                'Atmopshere model (%.2f) and input metallicity (%.2f) not consistent; proceed with caution' % (
                    modelmetal, self.metallicity), RuntimeWarning)
            logger.warn(
                'WARNING: Atmopshere model (%.2f) and input metallicity (%.2f) not consistent; proceed with caution' % (
                    modelmetal, self.metallicity))
        # -----------

        # -----build bysn----
        s1 = "'LAMBDA_MIN:'    '%.3f'\n" % self.lambda_min
        s1 += "'LAMBDA_MAX:'    '%.3f'\n" % self.lambda_max
        s1 += "'LAMBDA_STEP:'    '%.3f'\n" % self.lambda_delta
        s1 += "'INTENSITY/FLUX:' 'Flux'\n"
        s1 += "'COS(THETA)    :' '1.00'\n"
        s1 += "'ABFIND        :' '.false.'\n"
        s1 += "'MODELOPAC:' '%s%sopac'\n" % (
            self.continuum_opacity_dir, self.marcs_model_name.split('.int')[0] + 't%.2fl%i-%i' % (
                self.turbulent_velocity, self.lambda_min, self.lambda_max))
        s1 += "'RESULTFILE :' '%s'\n" % (self.spec_path + self.spec_file)
        s1 += "'METALLICITY:'    '%.2f'\n" % self.metallicity
        s1 += "'ALPHA/Fe   :'    '%.2f'\n" % self.alpha
        s1 += "'HELIUM     :'    '0.00'\n"
        s1 += "'R-PROCESS  :'    '%.2f'\n" % self.r_process
        s1 += "'S-PROCESS  :'    '%.2f'\n" % self.s_process

        # ---- Allow for user input abundances in the form [Abund, value]
        if self.free_abundances is None:
            individual_abundances = "'INDIVIDUAL ABUNDANCES:'   '0'\n"
        else:
            individual_abundances = "'INDIVIDUAL ABUNDANCES:'   '%i'\n" % (shape(self.free_abundances)[0])
            solar_abundances = loadtxt(self.solar_abundance_file, dtype=str)
            z = solar_abundances[:, 0]
            els = solar_abundances[:, 1]
            log_eps = solar_abundances[:, 2]
            if shape(self.free_abundances)[1] == 2:
                for i in arange(len(self.free_abundances)):
                    elind = where(els == self.free_abundances[i][0])[0]
                    if len(elind) != 1:
                        raise ValueError(
                            'element %s not found in %s' % (self.free_abundances[i][0], self.solar_abundance_file))
                    else:
                        elind = elind[0]
                    individual_abundances += '%i  %.2f\n' % (
                    int(z[elind]), float(log_eps[elind]) + float(self.free_abundances[i][1]))
            else:
                raise ValueError(
                    'shape of free_abundances is not correct... Please try again... free_abund = [[element,log(eps)]]')
        s1 += individual_abundances + '\n'
        s1 += "'ISOTOPES : ' '2'\n"
        s1 += '3.006  0.075\n'
        s1 += '3.007  0.925\n'

        s1 += "'NFILES   :' '%i'\n" % (len(self.line_list_files))
        for i in arange(len(self.line_list_files)):
            s1 += self.line_list_dir + self.line_list_files[i] + '\n'
        s1 += "'SPHERICAL:'  '%s'\n" % self.sphere
        s1 += '  30\n'
        s1 += '  300.00\n'
        s1 += '  15\n'
        s1 += '  1.30\n'

        # ------- babsma----
        s2 = "'LAMBDA_MIN:'    '%.3f'\n" % self.lambda_min
        s2 += "'LAMBDA_MAX:'    '%.3f'\n" % self.lambda_max
        s2 += "'LAMBDA_STEP:'    '%.3f'\n" % self.lambda_delta
        s2 += "'MODELINPUT:' '%s'\n" % (self.mod_path + self.marcs_model_name)
        s2 += "'MARCS-FILE:' '.false.'\n"  # if it is not a marcs model atmo (i.e. it is a interoplated one) needs to be false
        s2 += "'MODELOPAC:' '%s%sopac'\n" % (
            self.continuum_opacity_dir, self.marcs_model_name.split('.int')[0] + 't%.2fl%i-%i' % (
                self.turbulent_velocity, self.lambda_min, self.lambda_max))
        s2 += "'METALLICITY:'    '%.2f'\n" % self.metallicity
        s2 += "'ALPHA/Fe   :'    '%.2f'\n" % self.alpha
        s2 += "'HELIUM     :'    '0.00'\n"
        s2 += "'R-PROCESS  :'    '%.2f'\n" % self.r_process
        s2 += "'S-PROCESS  :'    '%.2f'\n" % self.s_process
        s2 += individual_abundances + '\n'
        s2 += "'XIFIX:' 'T'\n"
        s2 += '%.2f' % self.turbulent_velocity

        return s2, s1

    def Turbosynth(self, stellarpar=None):

        if stellarpar is None:
            babsma_in, bsyn_in = self.make_babsmabysn_file()
        else:
            if len(stellarpar) != 3:
                raise ValueError('The stellarpar varible must be a list/array of length 3 = (Teff, logg, [Fe/H])')
            else:
                T, g, f = stellarpar
                logger.info('Generating model atmosphere with T=%.1f, log g = %.2f, metallicity = %.2f' % (T, g, f))
                self.generate_model_atmosphere(T, g, f)
                babsma_in, bsyn_in = self.make_babsmabysn_file()

        if self.verbose:
            stdout = None
            stderr = subprocess.STDOUT
        else:
            stdout = open('/dev/null', 'w')
            stderr = subprocess.STDOUT

        try:
            pr1 = subprocess.Popen([self.turbo_path + 'babsma_lu'], stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
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
            pr = subprocess.Popen([self.turbo_path + 'bsyn_lu'], stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
            for line in bsyn_in:
                pr.stdin.write(line)
            stdout, stderr = pr.communicate()
        except subprocess.CalledProcessError:
            raise RuntimeError('babsma failed ....')
        logger.info("%s %s" % (pr.returncode, stderr))
        return pr.returncode
