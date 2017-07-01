#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This python package is a small collection of data about 4MOST -- e.g. its wavelength coverage and spectral resolution.
"""


class FourMost:
    def __init__(self):
        # Data from <https://www.4most.eu/cms/facility/>
        self.bands = {
            "LRS": {"R": 7500,
                    "lambda_min": 3700,
                    "lambda_max": 9500,
                    "line_lists_edvardsson": "fromBengt/line_lists/3700-9500/"
                    },
            "HRS_1": {"R": 20000,
                      "lambda_min": 3926,
                      "lambda_max": 4355,
                      "line_lists_edvardsson": "fromBengt/line_lists/3926-4355/"
                      },
            "HRS_2": {"R": 20000,
                      "lambda_min": 5160,
                      "lambda_max": 5730,
                      "line_lists_edvardsson": "fromBengt/line_lists/5160-5730/"
                      },
            "HRS_3": {"R": 20000,
                      "lambda_min": 6100,
                      "lambda_max": 6790,
                      "line_lists_edvardsson": "fromBengt/line_lists/6100-6790/"
                      }

        }
