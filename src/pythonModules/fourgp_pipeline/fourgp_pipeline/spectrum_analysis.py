# -*- coding: utf-8 -*-

"""

The `SpectrumAnalysis` class is created by the `Pipeline` class, and represents
the analysis of a particular spectrum.

"""

import logging


class SpectrumAnalysis:
    """
    Class which is created by the `Pipeline` class, and represents
    the analysis of a particular spectrum.
    """

    def __init__(self, input_spectrum):
        """
        Class which represents the analysis of a particular spectrum.

        :param input_spectrum:
            The Spectrum object we are to analyse.
        :type input_spectrum:
            Spectrum
        """

        self.input_spectrum = input_spectrum
        self.failure = False
        self.tasks_which_failed = []
        self.intermediate_results = []
        self.spectrum_by_task_name = {}
        self.metadata_by_task_name = {}

    def store_result(self, task_name, output_spectrum, output_metadata):
        """
        Store the output from some pipeline process.

        :param task_name:
            The name of the pipeline task we have just executed.
        :type task_name:
            str
        :param output_spectrum:
            The spectrum output from this process. If the process does not produce a spectrum as an output, then
            it should pass through the input spectrum that it worked on, so that that is accessible to the next
            pipeline step.
        :type output_spectrum:
            Spectrum
        :param output_metadata:
            Dictionary of metadata produced by this step of the pipeline.
        :type output_metadata:
            dict

        :return:
            None
        """

        self.intermediate_results.append(output_spectrum)
        self.spectrum_by_task_name[task_name] = output_spectrum
        self.metadata_by_task_name[task_name] = output_metadata

    def report_failure(self, task_name):
        """
        Method we call to flag up that some particular analysis step failed.

        :param task_name:
            The name of the pipeline task which failed.
        :return:
            None
        """

        logging.info("Task <{}> failed.".format(task_name))
        self.failure = True
        self.tasks_which_failed.append(task_name)
