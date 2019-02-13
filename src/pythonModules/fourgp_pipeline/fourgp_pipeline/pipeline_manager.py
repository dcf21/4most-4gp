# -*- coding: utf-8 -*-

"""

The `PipelineManager` class fetches work to do, processes spectra through a
pipeline, and posts results back to some destination. You will probably want to
create your own subclass which communicates with the specific servers (e.g.
4OR) that you want to fetch work from.

"""


class PipelineManager:
    """
    Class which fetches work to do, processes spectra through a
    pipeline, and posts results back to some destination.
    """
