#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convenience wrapper functions for the emcee package.
"""

import numpy as np
import acor


class SamplerData(object):
    """
    Emcee sampler wrapper that keeps the data in a similar manner
    but keeps the initialization values

    Attributes
    ----------

    Initialization with a Sampler

    sampler: sampler instance
        emcee sampler

    theta0: ndarray
        initial guess
    """
    def __init__(self, sampler, theta0):
        self.theta0 = theta0
        self.chain = sampler.chain
        self.lnprobability = sampler.lnprobability
        self.naccepted = sampler.naccepted
        self.iterations = sampler.iterations
        self.attrs = {}

    @property
    def dim(self):
        return len(self.theta0)

    @property
    def acceptance_fraction(self):
        """
        The fraction of proposed steps that were accepted.
        """
        return self.naccepted / self.iterations

    @property
    def flatlnprobability(self):
        """
        A shortcut to return the equivalent of ``lnprobability`` but aligned
        to ``flatchain`` rather than ``chain``.
        """
        return self.lnprobability.flatten()

    @property
    def flatchain(self):
        """
        A shortcut for accessing chain flattened along the zeroth (walker)
        axis.
        """
        s = self.chain.shape
        return self.chain.reshape(s[0] * s[1], s[2])

    @property
    def acor(self):
        """
        The autocorrelation time of each parameter in the chain (length:
        ``dim``) as estimated by the ``acor`` module.

        """
        s = self.dim
        t = np.zeros(s)
        for i in range(s):
            t[i] = acor.acor(self.chain[:, :, i])[0]
        return t


def run_mcmc(sampler, pos0, N, rstate0=None, lnprob0=None, **kwargs):
    """
    Iterate :func:`sample` for ``N`` iterations and return the result

    Parameters
    ----------
    pos0:
        The initial position vector.  Can also be None to resume from
        where :func:``run_mcmc`` left off the last time it executed.

    N:
        The number of steps to run.

    lnprob0: (optional)
        The log posterior probability at position ``p0``. If ``lnprob``
        is not provided, the initial value is calculated.

    rstate0: (optional)
        The state of the random number generator. See the
        :func:`random_state` property for details.

    kwargs: (optional)
        Other parameters that are directly passed to :func:`sample`.

    Returns
    -------
    t: tuple
        This returns the results of the final sample in whatever form
        :func:`sample` yields.  Usually, that's: ``pos``, ``lnprob``,
        ``rstate``, ``blobs`` (blobs optional)
    """
    if pos0 is None:
        if sampler._last_run_mcmc_result is None:
            raise ValueError("Cannot have pos0=None if run_mcmc has never been called.")
        pos0 = sampler._last_run_mcmc_result[0]
        if lnprob0 is None:
            rstate0 = sampler._last_run_mcmc_result[1]
        if rstate0 is None:
            rstate0 = sampler._last_run_mcmc_result[2]

    results = None
    for results in sampler.sample(pos0, lnprob0, rstate0, iterations=N, **kwargs):
        sampler._last_run_mcmc_result = results[:3]

    return results
