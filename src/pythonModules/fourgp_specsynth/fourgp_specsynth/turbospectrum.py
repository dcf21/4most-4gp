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
from operator import itemgetter

from solar_abundances import solar_abundances, periodic_table

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
        self.stellar_mass = 1
        self.log_g = 0
        self.t_eff = 5000
        self.turbulent_velocity = 1.0
        self.free_abundances = None
        self.sphere = None
        self.alpha = None
        self.s_process = 0
        self.r_process = 0
        self.verbose = True
        self.line_list_files = None

        # Create temporary directory
        self.id_string = "turbospec_{:d}".format(os.getpid())
        self.tmp_dir = os_path.join("/tmp", self.id_string)
        os.system("mkdir -p {}".format(self.tmp_dir))

        # Look up what MARCS models we have
        self.counter_marcs = 0
        self.marcs_model_name = "default"
        self.counter_spectra = 0
        self.marcs_values = None
        self.marcs_value_keys = []
        self.marcs_models = {}
        self._fetch_marcs_grid()

    def close(self):
        # Remove temporary directory
        os.system("rm -Rf {}".format(self.tmp_dir))

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

        marcs_models = glob.glob(os_path.join(self.marcs_grid_path, "*"))
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

    def configure(self, lambda_min=None, lambda_max=None, lambda_delta=None,
                  metallicity=None, log_g=None, t_eff=None, stellar_mass=None,
                  turbulent_velocity=None, free_abundances=None,
                  sphere=None, alpha=None, s_process=None, r_process=None,
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
            Metallicity of the star we're synthesizing.
        :param t_eff:
            Effective temperature of the star we're synthesizing.
        :param log_g:
            Log(gravity) of the star we're synthesizing.
        :param stellar_mass:
            Mass of the star we're synthesizing (solar masses).
        :param turbulent_velocity:
        :param free_abundances:
            List of elemental abundances to use in stellar model. These are passed to Turbospectrum.
        :param sphere:
            Select whether to use a spherical model (True) or a plane-parallel model (False).
        :param alpha:
            Alpha enhancement to use in stellar model.
        :param s_process:
            S-Process element enhancement to use in stellar model.
        :param r_process:
            R-Process element enhancement to use in stellar model.
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
        if t_eff is not None:
            self.t_eff = t_eff
        if log_g is not None:
            self.log_g = log_g
        if stellar_mass is not None:
            self.stellar_mass = stellar_mass
        if turbulent_velocity is not None:
            self.turbulent_velocity = turbulent_velocity
        if free_abundances is not None:
            self.free_abundances = free_abundances
        if sphere is not None:
            self.sphere = sphere
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

    def _generate_model_atmosphere(self):
        """
        Generates an interpolated model atmosphere from the MARCS grid using the interpol.f routine developed by
        T. Masseron (Masseron, PhD Thesis, 2006). This is a python wrapper for that fortran code.
        """
        self.counter_marcs += 1
        self.marcs_model_name = "marcs_{:08d}".format(self.counter_marcs)

        if self.verbose:
            stdout = subprocess.STDOUT
            stderr = subprocess.STDOUT
        else:
            stdout = open('/dev/null', 'w')
            stderr = subprocess.STDOUT

        # Defines default point at which plane-parallel vs spherical model atmosphere models are used
        spherical = self.sphere
        if spherical is None:
            spherical = (self.log_g < 3)

        # Default MARCS model settings
        marcs_parameters = {"spherical": spherical,
                            "mass": 0,
                            "turbulence": 0,
                            "model_type": "st",
                            "a": 0, "c": 0, "n": 0, "o": 0, "r": 0, "s": 0}

        # Pick MARCS settings which bracket requested stellar parameters
        interpolate_parameters = ("metallicity", "log_g", "temperature")

        interpolate_parameters_around = {"temperature": self.t_eff,
                                         "log_g": self.log_g,
                                         "metallicity": self.metallicity
                                         }

        for key in interpolate_parameters:
            value = interpolate_parameters_around[key]
            options = self.marcs_values[key]
            if (value < options[0]) or (value > options[-1]):
                raise ValueError("Value of parameter <{}> needs to be in range {} to {}. You requested {}.".
                                 format(key, options[0], options[-1], value))
            for index in range(len(options)):
                if value <= options[index]:
                    break
            marcs_parameters[key] = [options[index], options[index + 1], index, index+1]

        # Loop over eight vertices of cuboidal cell in parameter space, collecting MARCS models
        marcs_model_list = []
        failures = True
        while failures:
            marcs_model_list = []
            failures = 0
            n_vertices = 2 ** len(interpolate_parameters)
            for vertex in range(n_vertices):  # Loop over 8 vertices
                dict_iter = self.marcs_models  # Navigate through dictionary tree of MARCS models we have
                try:
                    for parameter in self.marcs_value_keys:
                        value = marcs_parameters[parameter]
                        # When we encounter Teff, log_g or metallicity, we get two options, not a single value
                        # Choose which one to use by looking at the binary bits of <vertex> as it counts from 0 to 7
                        if isinstance(value, (list, tuple)):
                            option_number = int(bool(vertex & (2 ** interpolate_parameters.index(parameter))))  # 0 or 1
                            value = value[option_number]
                        if value not in dict_iter:
                            dict_iter[value] = {}
                        dict_iter = dict_iter[value]  # Step to next level of dictionary tree
                    dict_iter = dict_iter['filename']
                except KeyError:
                    dict_iter = None
                    failures += 1
                marcs_model_list.append(dict_iter)

            # If there are MARCS models missing from the corners of the cuboid we tried, see which face had the most
            # corners missing, and move that face out by one grid row
            if failures:
                n_faces = 2 * len(interpolate_parameters)
                failures_per_face = []
                for cuboid_face_no in range(n_faces):  # Loop over 6 faces of cuboid
                    failure_count = 0
                    parameter_no = int(cuboid_face_no / 2)
                    option_no = cuboid_face_no & 1
                    for vertex in range(n_vertices):  # Loop over 8 vertices
                        if marcs_model_list[vertex] is None:
                            failure_option_no = int(bool(vertex & (2 ** parameter_no)))  # This is 0/1
                            if option_no == failure_option_no:
                                failure_count += 1
                    failures_per_face.append([failure_count, parameter_no, option_no])
                failures_per_face.sort(key=itemgetter(0))

                face_to_move = failures_per_face[-1]
                failure_count, parameter_no, option_no = face_to_move
                parameter_to_move = interpolate_parameters[parameter_no]
                options = self.marcs_values[parameter_to_move]
                parameter_descriptor = marcs_parameters[parameter_to_move]

                if option_no == 0:
                    parameter_descriptor[2] -= 1
                    assert parameter_descriptor[2] >= 0, \
                        "Value of parameter <{}> needs to be in range {} to {}. You requested {}, " \
                        "and due to missing models we could not interpolate.". \
                            format(parameter_to_move, options[0], options[-1], value)
                    parameter_descriptor[0] = options[parameter_descriptor[2]]
                else:
                    parameter_descriptor[3] += 1
                    assert parameter_descriptor[3] < len(options), \
                        "Value of parameter <{}> needs to be in range {} to {}. You requested {}, " \
                        "and due to missing models we could not interpolate.". \
                            format(parameter_to_move, options[0], options[-1], value)
                    parameter_descriptor[1] = options[parameter_descriptor[3]]

        # Now we run the FORTRAN model interpolator
        output = os_path.join(self.tmp_dir, self.marcs_model_name)
        model_test = "{}.test".format(output)
        try:
            p = subprocess.Popen([os_path.join(self.interpol_path, 'interpol_modeles')],
                                 stdin=subprocess.PIPE, stdout=stdout, stderr=stderr)
            for line in marcs_model_list:
                p.stdin.write("{}\n".format(line))
            p.stdin.write("{}.interpol\n".format(output))
            p.stdin.write("{}.alt\n".format(output))
            p.stdin.write("{}\n".format(self.t_eff))
            p.stdin.write("{}\n".format(self.log_g))
            p.stdin.write("{}\n".format(self.metallicity))
            p.stdin.write(
                ".true.\n")  # the test option is set to .true. if you want to plot comparison model (model_test)
            p.stdin.write(".false.\n")  # MARCS binary format (.true.) or MARCS ASCII web format (.false.)?
            p.stdin.write("{}\n".format(model_test))

            stdout, stderr = p.communicate()
        except subprocess.CalledProcessError:
            raise RuntimeError('MARCS model atmosphere interpolation failed ....')

        return

    def make_babsma_bysn_file(self):
        """
        Generate the configurations files for both the babsma and bsyn binaries in Turbospectrum.
        """

        # If we've not been given an explicit alpha enhancement value, assume one
        alpha = self.alpha
        if alpha is None:
            if self.metallicity < -1.0:
                alpha = 0.4
            elif -1.0 < self.metallicity < 0.0:
                alpha = -0.4 * self.metallicity
            else:
                alpha = 0

        # Allow for user input abundances as a dictionary of the form {element: abundance}
        if self.free_abundances is None:
            individual_abundances = "'INDIVIDUAL ABUNDANCES:'   '0'\n"
        else:
            individual_abundances = "'INDIVIDUAL ABUNDANCES:'   '{:d}'\n".format(len(self.free_abundances))

            for element, abundance in self.free_abundances.iteritems():
                assert element in solar_abundances, "Cannot proceed as solar abundance for element <{}> is unknown". \
                    format(element)

                atomic_number = periodic_table.index(element)
                individual_abundances += "{:d}  {:.2f}\n".format(int(atomic_number),
                                                                 float(solar_abundances[element]) + float(abundance))

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
'RESULTFILE :' '{this[tmp_dir]}/spectrum_{this[counter_spectra]:08d}.spec'
'METALLICITY:'    '{this[metallicity]:.2f}'
'ALPHA/Fe   :'    '{alpha:.2f}'
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
""".format(this=self.__dict__, alpha=alpha, individual_abundances=individual_abundances, line_lists=line_lists)

        # Build babsma configuration file
        babsma_config = """\
'LAMBDA_MIN:'    '{this[lambda_min]:.3f}'
'LAMBDA_MAX:'    '{this[lambda_max]:.3f}'
'LAMBDA_STEP:'    '{this[lambda_delta]:.3f}'
'MODELINPUT:' '{this[tmp_dir]}/{this[marcs_model_name]}'
'MARCS-FILE:' '.false.'
'MODELOPAC:' '{this[tmp_dir]}/model_opacity.opac'
'METALLICITY:'    '{this[metallicity]:.2f}'
'ALPHA/Fe   :'    '{alpha:.2f}'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '{this[r_process]:.2f}'
'S-PROCESS  :'    '{this[s_process]:.2f}'
{individual_abundances}
'XIFIX:' 'T'
{this[turbulent_velocity]:.2f}
""".format(this=self.__dict__, alpha=alpha, individual_abundances=individual_abundances)

        return babsma_config, bsyn_config

    def synthesise(self):
        """
        Invoke Turbospectrum to synthesise a single spectrum.
        """
        self.counter_spectra += 1

        logger.info("Generating model atmosphere with T={:.1f}, log_g = {:.2f}, metallicity = {:.2f}".
                    format(self.t_eff, self.log_g, self.metallicity))

        self._generate_model_atmosphere()
        babsma_in, bsyn_in = self.make_babsma_bysn_file()

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

        return {
            "return_code": pr.returncode,
            "output_file": os_path.join(self.tmp_dir, "spectrum_{:08d}.spec".format(self.counter_spectra))
        }
