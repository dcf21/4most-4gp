# -*- coding: utf-8 -*-

"""

Here we define implementations of the tasks we perform to analyse the spectra of FGK stars -- e.g.
determining RVs, or continuum normalising spectra.

"""

import gzip
import json
import logging
from os import path as os_path

import numpy as np
from fourgp_cannon.cannon_wrapper_casey_old import CannonInstanceCaseyOld
from fourgp_degrade import SpectrumProperties, SpectrumResampler
from fourgp_rv import RvInstanceCrossCorrelation
from fourgp_speclib import SpectrumLibrarySqlite

from .dummy_processes import ContinuumNormalisationDummy, DecisionMakerDummy
from .pipeline import Pipeline, PipelineTask, PipelineFailure


class TaskCleanUpNans(PipelineTask):
    """
    Task to make sure the spectrum doesn't have any NaNs in it, as these cause numerical interpolation to fail.
    """

    @staticmethod
    def task_name():
        """
        :return:
            Name for this pipeline task.
        """
        return "clean_up_nans"

    def __init__(self):
        """
        No initialisation required.
        """

        super(TaskCleanUpNans, self).__init__()

    def task_implementation(self, input_spectrum, spectrum_analysis):
        """
        Task to make sure the spectrum doesn't have any NaNs in it, as these cause numerical interpolation to fail.

        :param input_spectrum:
            The Spectrum object we are to continuum normalise.
        :type input_spectrum:
            Spectrum
        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            The continuum-normalised spectrum.
        """

        output_spectrum = input_spectrum.copy()

        # Replace errors which are nans with a large value, otherwise they cause numerical failures in the RV code
        output_spectrum.value_errors[np.isnan(output_spectrum.value_errors)] = 1000.

        # Check for NaNs in actual spectrum
        if not np.all(np.isfinite(output_spectrum.values)):
            logging.warning("Warning: NaN values in test spectrum!")
            output_spectrum.value_errors[np.isnan(output_spectrum.values)] = 1000.
            output_spectrum.values[np.isnan(output_spectrum.values)] = 1.

        return {
            'spectrum': output_spectrum,
            'metadata': {}
        }


class TaskContinuumNormalisationFirstPass(PipelineTask):
    """
    First-pass continuum normalisation, done before RV correction.
    """

    @staticmethod
    def task_name():
        """
        :return:
            Name for this pipeline task.
        """
        return "continuum_normalise_pass_1"

    def __init__(self):
        """
        Instantiate the dummy continuum normalisation code.
        """

        super(TaskContinuumNormalisationFirstPass, self).__init__()
        self.normaliser = ContinuumNormalisationDummy()

    def task_implementation(self, input_spectrum, spectrum_analysis):
        """
        Do the first pass continuum normalisation on a spectrum.

        :param input_spectrum:
            The Spectrum object we are to continuum normalise.
        :type input_spectrum:
            Spectrum
        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            The continuum-normalised spectrum.
        """

        return {
            'spectrum': self.normaliser.continuum_normalise(
                spectrum_flux_normalised=input_spectrum,
                spectrum_analysis=spectrum_analysis
            ),
            'metadata': {}
        }


class TaskRVCorrect(PipelineTask):
    """
    Estimate the radial velocity of a star by the cross correlation method.
    """

    @staticmethod
    def task_name():
        """
        :return:
            Name for this pipeline task.
        """
        return "rv_correct"

    def __init__(self,
                 workspace,
                 fourmost_mode,
                 rv_cross_correlation_library="rv_templates_resampled",
                 rv_upsampling=1):
        """
        Initialise the cross-correlation RV code.

        :param workspace:
            Directory where we expect to find spectrum libraries.
        :type workspace:
            str
        :param fourmost_mode:
            The name of the 4MOST mode we are operating, either hrs or lrs
        :type fourmost_mode:
            str
        :param rv_cross_correlation_library:
            The name of the spectrum library we are to get our cross correlation templates from.
        :type rv_cross_correlation_library:
            str
        :param rv_upsampling:
            The upsampling factor to apply to input spectra before cross correlating them with the templates.
        :type rv_upsampling:
            int
        """
        super(TaskRVCorrect, self).__init__()
        self.fourmost_mode = fourmost_mode

        # Open spectrum library containing cross-correlation templates
        template_library = SpectrumLibrarySqlite(
            path=os_path.join(workspace, rv_cross_correlation_library),
            create=False,
        )

        # Instantiate RV code
        self.rv_code = RvInstanceCrossCorrelation(spectrum_library=template_library,
                                                  upsampling=rv_upsampling)

    def task_implementation(self, input_spectrum, spectrum_analysis):
        """
        RV correct a spectrum.

        :param input_spectrum:
            The Spectrum object we are to RV correct.
        :type input_spectrum:
            Spectrum
        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            The RV-corrected spectrum.
        """

        rv_mean, rv_std_dev, stellar_parameters, rv_estimates_by_weight = self.rv_code.estimate_rv(
            input_spectrum=input_spectrum,
            mode=self.fourmost_mode.upper()
        )

        spectrum_object_rest_frame = input_spectrum.correct_radial_velocity(v=rv_mean)
        resampler = SpectrumResampler(input_spectrum=spectrum_object_rest_frame)

        return {
            'spectrum': resampler.match_to_other_spectrum(other=input_spectrum),
            'metadata': {
                'rv': rv_mean,
                'rv_uncertainty': rv_std_dev,
                'stellar_parameter': stellar_parameters,
                'rv_estimates_by_weight': rv_estimates_by_weight
            }
        }


class TaskContinuumNormalisationSecondPass(PipelineTask):
    """
    Second-pass continuum normalisation, done after the RV of the spectrum has been estimated, and this radial
    velocity has been estimated.
    """

    @staticmethod
    def task_name():
        """
        :return:
            Name for this pipeline task.
        """
        return "continuum_normalise_pass_2"

    def __init__(self):
        """
        Instantiate the dummy continuum normalisation code.
        """
        super(TaskContinuumNormalisationSecondPass, self).__init__()
        self.normaliser = ContinuumNormalisationDummy()

    def task_implementation(self, input_spectrum, spectrum_analysis):
        """
        Do the second pass continuum normalisation on a spectrum.

        :param input_spectrum:
            The Spectrum object we are to continuum normalise.
        :type input_spectrum:
            Spectrum
        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            The continuum-normalised spectrum.
        """

        return {
            'spectrum': self.normaliser.continuum_normalise(
                spectrum_flux_normalised=input_spectrum,
                spectrum_analysis=spectrum_analysis
            ),
            'metadata': {}
        }


class TaskCannonAbundanceAnalysis(PipelineTask):
    """
    Run the spectrum through the Cannon to estimate stellar parameters and abundances.
    """

    @staticmethod
    def task_name():
        """
        :return:
            Name for this pipeline task.
        """
        return "cannon_abundance_analysis"

    def __init__(self, workspace, reload_cannon_from_file):
        """
        Instantiate a Cannon, and load training data from disk.

        :param workspace:
            Directory where we expect to find spectrum libraries.
        :type workspace:
            str
        :param reload_cannon_from_file:
            The filename of the output files containing the trained Cannon that we are to reload.
        :type reload_cannon_from_file:
            str
        """

        super(TaskCannonAbundanceAnalysis, self).__init__()

        # Load the JSON data that summarises the Cannon training that we're about to reload
        json_summary_filename = "{}.summary.json.gz".format(reload_cannon_from_file)
        cannon_pickle_filename = "{}.cannon".format(reload_cannon_from_file)

        with gzip.open(json_summary_filename, "rt") as f:
            summary_json = json.loads(f.read())

        raster = np.array(summary_json['wavelength_raster'])
        test_labels = summary_json['labels']

        # If we're doing our own continuum normalisation, we need to treat each wavelength arm separately
        wavelength_arm_breaks = SpectrumProperties(raster).wavelength_arms()['break_points']

        self.model = CannonInstanceCaseyOld(training_set=None,
                                            wavelength_arms=wavelength_arm_breaks,
                                            load_from_file=cannon_pickle_filename,
                                            label_names=test_labels,
                                            censors=None
                                            )

    def task_implementation(self, input_spectrum, spectrum_analysis):
        """
        Do the second pass continuum normalisation on a spectrum.

        :param input_spectrum:
            The Spectrum object we are to continuum normalise.
        :type input_spectrum:
            Spectrum
        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            The continuum-normalised spectrum.
        """

        # Pass spectrum to the Cannon
        labels, cov, meta = self.model.fit_spectrum(spectrum=input_spectrum)

        # Check whether Cannon failed
        if labels is None:
            raise PipelineFailure("The Cannon failed to fit this spectrum.")

        # Return results
        return {
            'spectrum': input_spectrum,
            'metadata': {
                'label_values': labels,
                'covariance_matrix': cov,
                'metadata': meta
            }
        }


class TaskDecisionMaker(PipelineTask):
    """
    Inspect the results of this run of the pipeline, and decide whether the run was a success or a failure.
    """

    @staticmethod
    def task_name():
        """
        :return:
            Name for this pipeline task.
        """
        return "decision_maker"

    def __init__(self):
        super(TaskDecisionMaker, self).__init__()
        self.oracle = DecisionMakerDummy()

    def task_implementation(self, input_spectrum, spectrum_analysis):
        """
        Do the second pass continuum normalisation on a spectrum.

        :param input_spectrum:
            The Spectrum object we are to continuum normalise.
        :type input_spectrum:
            Spectrum
        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            The continuum-normalised spectrum.
        """

        success = self.oracle.consult(spectrum_analysis=spectrum_analysis)

        return {
            'spectrum': input_spectrum,
            'metadata': {
                'success': success
            }
        }


class PipelineFGK(Pipeline):
    """
    Pipeline for FGK stars.
    """

    def __init__(self, workspace, fourmost_mode, reload_cannon_from_file):
        """
        Pipeline for processing spectra of FGK stars, for either 4MOST HRS or LRS.

        :param workspace:
            Directory where we expect to find spectrum libraries.
        :type workspace:
            str
        :param fourmost_mode:
            String containing either "hrs" or "lrs" indicating which 4MOST mode we are analysing.
        :type fourmost_mode:
            str
        :param reload_cannon_from_file:
            The filename of the output files containing the trained Cannon that we are to reload.
        :type reload_cannon_from_file:
            str
        """

        # Check that 4MOST mode is valid
        assert fourmost_mode in ("hrs", "lrs"), \
            "4MOST mode <{}> is not recognised.".format(fourmost_mode)

        # Initialise an empty pipeline
        super(PipelineFGK, self).__init__()

        # Create pipeline, step by step
        self.append_task(task_name="continuum_normalise_pass_1",
                         task_implementation=TaskContinuumNormalisationFirstPass()
                         )
        self.append_task(task_name="clean_up_nans",
                         task_implementation=TaskCleanUpNans()
                         )
        self.append_task(task_name="rv_correct",
                         task_implementation=TaskRVCorrect(
                             workspace=workspace,
                             fourmost_mode=fourmost_mode
                         )
                         )
        self.append_task(task_name="continuum_normalise_pass_2",
                         task_implementation=TaskContinuumNormalisationSecondPass()
                         )
        self.append_task(task_name="cannon_abundance_analysis",
                         task_implementation=TaskCannonAbundanceAnalysis(
                             workspace=workspace,
                             reload_cannon_from_file=reload_cannon_from_file
                         )
                         )
        self.append_task(task_name="decision_maker",
                         task_implementation=TaskDecisionMaker()
                         )
