# -*- coding: utf-8 -*-

"""

The `PipelineManager` class fetches work to do, processes spectra through a
pipeline, and posts results back to some destination. You will probably want to
create your own subclass which communicates with the specific servers (e.g.
4OR) that you want to fetch work from.

"""

from .pipeline import Pipeline


class PipelineManager:
    """
    Class which fetches work to do, processes spectra through a
    pipeline, and posts results back to some destination.
    """

    def __init__(self, pipeline):
        """
        Class which fetches work to do, processes spectra through a pipeline, and posts results back to some
        destination.

        :param pipeline:
            The Pipeline we are to run spectra through.
        :type pipeline:
            Pipeline
        """

        # Instantiate the pipeline that we're going to use to analyse spectra
        self.pipeline = pipeline

    def fetch_work(self):
        """
        Check to see if we have any spectra to analyse. If yes, return the next Spectrum object which needs
        analysing. If we have nothing to do, return None.

        :return:
            Spectrum object to be analysed
        """

        raise NotImplementedError("The fetch work method needs to be implemented to talk to whatever database you "
                                  "wish to store your work in.")

    def post_result(self, spectrum_analysis):
        """
        Post the results from analysing a spectrum back to whatever database you are using to store the pipeline's
        output

        :param spectrum_analysis:
            A SpectrumAnalysis object containing the results of the pipeline's analysis of the spectrum.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            True for success; False for failure.
        """

        raise NotImplementedError("The fetch work method needs to be implemented to talk to whatever database you "
                                  "wish to store your work in.")

    def handle_failure(self, spectrum_analysis):
        """
        Define the actions we follow if we fail to analyse a spectrum. In this base implementation, we do nothing.

        :param spectrum_analysis:
            A SpectrumAnalysis object containing the failed analysis of a spectrum.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            None
        """

        return None

    def poll_server_for_work(self):
        """
        See if we have any spectra to analyse. If so, analyse a spectrum, and post the results back. Otherwise,
        return without doing anything.

        :return:
            The number of spectra we analysed.
        """

        # Query whether we have any work to do
        job_description = self.fetch_work()

        # If we have no work to do, exit immediately
        if job_description is None:
            return 0

        input_spectrum = job_description['spectrum']
        spectrum_identifier = job_description['spectrum_identifier']

        # If we got a spectrum, analyse it now
        analysis = self.pipeline.analyse_spectrum(input_spectrum=input_spectrum,
                                                  spectrum_identifier=spectrum_identifier)

        # Post results back
        self.post_result(spectrum_analysis=analysis)

        # Check for failure
        if analysis.failure:
            self.handle_failure(spectrum_analysis=analysis)

        # Finished
        return 1

    def set_pipeline(self, new_pipeline):
        """
        Define an entirely new pipeline which we should use for analysing spectra.

        :param new_pipeline:
            An instance of the Pipeline class, representing the pipeline we should use to analyse spectra.
        :type new_pipeline:
            Pipeline
        :return:
            None
        """

        assert isinstance(new_pipeline, Pipeline), \
            "An analysis pipeline must be an instance of the base class Pipeline."

        self.pipeline = new_pipeline

    def set_task_implementation(self, task_name, task_implementation):
        """
        Set a new implementation for one of the tasks involved in the analysis of spectra.

        :param task_name:
            The name of the task we are to set a new implementation for. This must be one of the keys of the
            dictionary self.pipeline.task_implementations.
        :type task_name:
            str
        :param task_implementation:
            The class which implements this task.
        :return:
            None
        """

        self.pipeline.set_task_implementation(task_name=task_name, task_implementation=task_implementation)
