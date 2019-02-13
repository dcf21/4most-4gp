# -*- coding: utf-8 -*-

"""

The `Pipeline` class represents a pipeline for analysing spectra. It is
possible to configure which 4GP classes it uses to perform each task within the
pipeline -- e.g. determining RVs, or continuum normalising spectra.

"""

from .dummy_processes import ContinuumNormalisationDummy


class Pipeline:
    """
    Class which represents a pipeline for analysing spectra.
    """

    def __init__(self, **kwargs):
        """
        Create a pipeline which uses to default implementations for each task in the analysis of spectra.

        By specifying keyword arguments to this constructor, it is possible to override any of the defaults in the
        dictionary below

        :param kwargs:
            Override any of the default classes used for each task in the analysis of spectra.
        """

        self.task_implementations = {
            "continuum_normalise_pass_1": ContinuumNormalisationDummy(),
            "continuum_normalise_pass_2": ContinuumNormalisationDummy()
        }

        for item_name, item_implementation in kwargs.items():
            self.set_task_implementation(task_name=item_name, task_implementation=item_implementation)

    def set_task_implementation(self, task_name, task_implementation):
        """
        Set a new implementation for one of the tasks involved in the analysis of spectra.

        :param task_name:
            The name of the task we are to set a new implementation for. This must be one of the keys of the
            dictionary self.task_implementations.
        :type task_name:
            str
        :param task_implementation:
            The class which implements this task.
        :return:
            None
        """

        assert task_name in self.task_implementations, "Unknown task <{}>".format(task_name)

        if not isinstance(task_implementation, (list, tuple)):
            task_implementation = [task_implementation]

        self.task_implementations[task_name] = task_implementation

    def append_task_implementation(self, task_name, task_implementation):
        """
        Add an additional new implementation for one of the tasks involved in the analysis of spectra.

        :param task_name:
            The name of the task we are to set a new implementation for. This must be one of the keys of the
            dictionary self.task_implementations.
        :type task_name:
            str
        :param task_implementation:
            The class which implements this task.
        :return:
            None
        """

        assert task_name in self.task_implementations, "Unknown task <{}>".format(task_name)

        self.task_implementations[task_name].append(task_implementation)

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
