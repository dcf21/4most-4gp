# -*- coding: utf-8 -*-

"""
Code to generate random radial velocities when performing synthetic tests of the performance of RV codes.
"""

import random


def random_radial_velocity():
    """
    Pick a random radial velocity in km/s, following the approximate probability distribution expected for 4MOST
    galactic stars.

    :return:
        Radial velocity in km/s
    """

    # We first of all have a random switch to choose between two populations of stars
    distribution_selector = random.uniform(a=0, b=100)

    if distribution_selector < 10:
        # Pick 10% of stars from a uniform distribution from -200 to 200
        # This population is added so that we have a small fraction of high velocity stars
        return random.uniform(a=-200, b=200)  # Unit km/s
    else:
        # Pick 90% of stars from a Gaussian distribution with a standard deviation of 25 km/s
        return random.gauss(mu=0, sigma=25)
