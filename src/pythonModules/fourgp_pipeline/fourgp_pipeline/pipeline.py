# -*- coding: utf-8 -*-

"""

The `Pipeline` class represents a pipeline which runs a sequence of tasks for analysing spectra. By defining new
descendents of the PipelineTask class, and appending them to a Pipeline, it is
possible to configure which 4GP classes it uses to perform each task within the
pipeline -- e.g. determining RVs, or continuum normalising spectra.

"""

import logging

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

    def analyse_spectrum(self, input_spectrum, spectrum_identifier):
        """
        Analyse a spectrum through the 4GP pipeline.

        :param input_spectrum:
            The Spectrum object we are to analyse.
        :type input_spectrum:
            Spectrum
        :param spectrum_identifier:
            Some string name that we can use in logging messages to identify which spectrum we are working on.
        :type spectrum_identifier:
            str
        :return:
            A SpectrumAnalysis object representing the analysis of this spectrum.
        """

        logging.info("Working on spectrum <{}>".format(spectrum_identifier))

        # Create a structure to hold the results from analysing this spectrum
        spectrum_analysis = SpectrumAnalysis(input_spectrum=input_spectrum)

        spectrum_analysis.store_result(task_name="initial_input",
                                       output_spectrum=input_spectrum,
                                       output_metadata=input_spectrum.metadata
                                       )

        # Run each task in turn. Stop running tasks if something breaks.
        for task in self.task_list:
            if not spectrum_analysis.failure:
                task['implementation'].run_task(spectrum_analysis=spectrum_analysis)

        # Return result
        return spectrum_analysis


class PipelineTask:
    """
    A class representing a task which needs to be performed by the pipeline, and exposing it via a standard calling API.
    """

    def __init__(self, configuration=None):
        """
        A class representing a task which needs to be performed by the pipeline, and exposing it via a standard
        calling API.

        :param configuration:
            Optional dictionary, containing configuration parameters to pass to this task.
        """

        if configuration is None:
            configuration = {}

        self.configuration = configuration

    @staticmethod
    def task_name():
        """
        All pipeline tasks must have a defined name.
        :return:
            string name
        """

        raise NotImplementedError("Descendents of the class PipelineTask must define a name for themselves")

    def run_task(self, spectrum_analysis):
        """
        Run this pipeline task, as the next step in the analysis of a spectrum, whose analysis hitherto is summarised
        in the structure spectrum_analysis.

        :param spectrum_analysis:
            The analysis hitherto of the spectrum we are to perform a pipeline task upon.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            None
        """

        task_name = self.task_name()

        # Fetch the most recent intermediate result from the analysis of this spectrum
        input_spectrum = spectrum_analysis.intermediate_results[-1]

        # Run this task on that intermediate result
        try:
            logging.info("Running task <{task_name}>".format(task_name=task_name))
            result = self.task_implementation(
                input_spectrum=input_spectrum,
                spectrum_analysis=spectrum_analysis
            )

            # Store the output of this step
            spectrum_analysis.store_result(
                task_name=task_name,
                output_spectrum=result['spectrum'],
                output_metadata=result['metadata']
            )
        except PipelineFailure:
            spectrum_analysis.report_failure(task_name=self.task_name())

    def task_implementation(self, input_spectrum, spectrum_analysis):
        """
        Desecent classes must insert code here to perform some task on the spectrum <input_spectrum>, or to analyse
        the analysis done so far, as described in spectrum_analysis.

        :param input_spectrum:
            The latest intermediate result produced by this pipeline.
        :type input_spectrum:
            Spectrum
        :param spectrum_analysis:
            The complete analysis of this spectrum done so far by this pipeline.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            None
        """

        raise NotImplementedError("The task implementation must be specified for each descendent of the "
                                  "PipelineTask class")


class PipelineFailure(Exception):
    """
    An exception that PipelineTasks should raise when shit happens.
    """

    pass
