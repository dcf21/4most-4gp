# -*- coding: utf-8 -*-

"""

Here we define implementations of the tasks we perform to analyse the spectra of FGK stars -- e.g.
determining RVs, or continuum normalising spectra.

"""

import gzip
import json
from os import path as os_path

import numpy as np
from fourgp_cannon.cannon_wrapper_casey_old import CannonInstanceCaseyOld
from fourgp_rv import RvInstanceCrossCorrelation
from fourgp_speclib import SpectrumLibrarySqlite

from .dummy_processes import ContinuumNormalisationDummy, DecisionMakerDummy
from .pipeline import Pipeline, PipelineTask, PipelineFailure


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
            mode=self.fourmost_mode
        )

        return {
            'spectrum': input_spectrum.apply_rv(rv=rv_mean),
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
        # We need the JSON data as well as the pickle that Andy Casey's code saves, as this contains information
        # such as the censoring masks which we need to reproduce. It also makes this code immune to the detail that Andy
        # Casey's various versions of the Cannon save their training data in incompatible formats.
        json_summary_filename = "{}.summary.json.gz".format(reload_cannon_from_file)

        with gzip.open(json_summary_filename, "rt") as f:
            summary_json = json.loads(f.read())

        training_spectrum_library_name = summary_json['train_library']
        label_names = summary_json['labels']
        censoring_masks = summary_json['censoring_mask']

        # Reload training spectra
        spectra = SpectrumLibrarySqlite.open_and_search(
            library_spec=training_spectrum_library_name,
            workspace=workspace,
            extra_constraints={"continuum_normalised": 1}
        )
        training_library, training_library_items = [spectra[i] for i in ("library", "items")]

        # Load training set
        training_library_ids_all = [i["specId"] for i in training_library_items]
        training_spectra_all = training_library.open(ids=training_library_ids_all)

        # JSON serialisation turns censoring masks from numpy arrays into tuples. Turn them back into numpy
        # arrays before passing to the Cannon
        if censoring_masks is None:
            censoring_masks_numpy = None
        else:
            censoring_masks_numpy = dict([(label, np.array(mask))
                                          for label, mask in censoring_masks.items()])

        # Instantiate the wrapper to the Cannon
        self.cannon = CannonInstanceCaseyOld(training_set=training_spectra_all,
                                             label_names=label_names,
                                             censors=censoring_masks_numpy
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
        labels, cov, meta = self.cannon.fit_spectrum(spectrum=input_spectrum)

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
