# -*- coding: utf-8 -*-

"""

The `Pipeline` class represents a pipeline which runs a sequence of tasks for analysing spectra. By defining new
descendents of the PipelineTask class, and appending them to a Pipeline, it is
possible to configure which 4GP classes it uses to perform each task within the
pipeline -- e.g. determining RVs, or continuum normalising spectra.

"""

from .spectrum_analysis import SpectrumAnalysis


class Pipeline:
    """
    Class which represents a pipeline for analysing spectra.
    """

    def __init__(self):
        """
        Base class representing a pipeline which runs a sequence of PipelineTasks in turn.
        """

        self.task_names = {}
        self.task_list = []

    def append_task(self, task_name, task_implementation):
        """
        Appends a new task to the end of the actions performed by this pipeline.

        :param task_name:
            The name of the task we are to carry out.
        :type task_name:
            str
        :param task_implementation:
            The descendent of the class PipelineTask which implements this task.
        :return:
            None
        """

        assert task_name not in self.task_names, "Pipeline has multiple tasks with the name <{}>".format(task_name)

        assert isinstance(task_implementation, PipelineTask), \
            "A Pipeline task must be a descendent of the class PipelineTask."

        self.task_names[task_name] = True
        self.task_list.append({'name': task_name, 'implementation': task_implementation})

    def analyse_spectrum(self, input_spectrum):
        """
        Analyse a spectrum through the 4GP pipeline.

        :param input_spectrum:
            The Spectrum object we are to analyse.
        :type input_spectrum:
            Spectrum
        :return:
            A SpectrumAnalysis object representing the analysis of this spectrum.
        """

        # Create a structure to hold the results from analysing this spectrum
        spectrum_analysis = SpectrumAnalysis(input_spectrum=input_spectrum)

        # Run each task in turn
        for task in self.task_list:
            task['implementation'].run_task(spectrum_analysis=spectrum_analysis)

        # Return result
        return spectrum_analysis


class PipelineTask:
    """
    A class representing a task which needs to be performed by the pipeline, and exposing it via a standard calling API
    """

    def __init__(self, configuration=None):
        if configuration is None:
            configuration = {}

        self.configuration = configuration

    def run_task(self, spectrum_analysis):
        input_spectrum = spectrum_analysis.intermediate_results[-1]
        self.task_implementation(input_spectrum=input_spectrum,
                                 spectrum_analysis=spectrum_analysis
                                 )

    def task_implementation(self, input_spectrum, spectrum_analysis):
        raise NotImplementedError("The task implementation must be specified for each descendent of the "
                                  "PipelineTask class")
