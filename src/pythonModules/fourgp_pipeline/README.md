# fourgp_pipeline

This python package defines classes for automatically processing spectra through an analysis pipeline.

The `PipelineManager` class fetches work to do, processes spectra through a pipeline, and posts results back to some destination. You will probably want to create your own subclass which communicates with the specific servers (e.g. 4OR) that you want to fetch work from.

The `Pipeline` class represents a pipeline for analysing spectra. It is possible to configure which 4GP classes it uses to perform each task within the pipeline -- e.g. determining RVs, or continuum normalising spectra.

The `SpectrumAnalysis` class is created by the `Pipeline` class, and represents the analysis of a particular spectrum.

# Contact details
This code is maintained by:

Dominic Ford  
Lund Observatory  
Box 43  
SE-221 00 Lund  
Sweden
