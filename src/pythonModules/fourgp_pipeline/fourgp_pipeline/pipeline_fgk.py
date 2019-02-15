# -*- coding: utf-8 -*-

"""

Here we define implementations of the tasks we perform to analyse the spectra of FGK stars -- e.g.
determining RVs, or continuum normalising spectra.

"""

from fourgp_rv import RvInstanceCrossCorrelation

from .dummy_processes import ContinuumNormalisationDummy, DecisionMakerDummy
from .pipeline import Pipeline, PipelineTask


class TaskContinuumNormalisationFirstPass(PipelineTask):
    @staticmethod
    def task_name():
        return "continuum_normalise_pass_1"

    def __init__(self):
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

        return self.normaliser.continuum_normalise(
            spectrum_flux_normalised=input_spectrum,
            spectrum_analysis=spectrum_analysis
        )


class TaskRVCorrect(PipelineTask):
    @staticmethod
    def task_name():
        return "rv_correct"

    def __init__(self,
                 rv_cross_correlation_library="rv_templates_resampled",
                 rv_upsampling=1):
        super(TaskRVCorrect, self).__init__()
        self.rv_code = RvInstanceCrossCorrelation(spectrum_library=rv_cross_correlation_library,
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

        return self.rv_code.estimate_rv(input_spectrum=input_spectrum,
                                        mode=spectrum_analysis.fourmost_mode)


class TaskContinuumNormalisationSecondPass(PipelineTask):
    @staticmethod
    def task_name():
        return "continuum_normalise_pass_2"

    def __init__(self):
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

        return self.normaliser.continuum_normalise(
            spectrum_flux_normalised=input_spectrum,
            spectrum_analysis=spectrum_analysis
        )


class TaskCannonAbundanceAnalysis(PipelineTask):
    @staticmethod
    def task_name():
        return "cannon_abundance_analysis"

    def __init__(self):
        super(TaskCannonAbundanceAnalysis, self).__init__()
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

        return self.normaliser.continuum_normalise(
            spectrum_flux_normalised=input_spectrum,
            spectrum_analysis=spectrum_analysis
        )


class TaskDecisionMaker(PipelineTask):
    @staticmethod
    def task_name():
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

        return self.oracle.consult(spectrum_analysis=spectrum_analysis)


class PipelineFGK(Pipeline):
    """
    Pipeline for FGK stars.
    """

    def __init__(self, **kwargs):
        # Initialise an empty pipeline
        super(PipelineFGK, self).__init__()

        # Create pipeline, step by step
        self.append_task(task_name="continuum_normalise_pass_1",
                         task_implementation=TaskContinuumNormalisationFirstPass()
                         )
        self.append_task(task_name="rv_correct",
                         task_implementation=TaskRVCorrect(),
                         )
        self.append_task(task_name="continuum_normalise_pass_2",
                         task_implementation=TaskContinuumNormalisationSecondPass()
                         )
        self.append_task(task_name="cannon_abundance_analysis",
                         task_implementation=TaskCannonAbundanceAnalysis()
                         )
        self.append_task(task_name="decision_maker",
                         task_implementation=TaskDecisionMaker()
                         )
