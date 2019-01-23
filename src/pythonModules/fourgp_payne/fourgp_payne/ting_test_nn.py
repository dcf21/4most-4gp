# -*- coding: utf-8 -*-

import os
from multiprocessing import Pool

import numpy as np
from scipy.optimize import curve_fit


def test_nn(payne_status, threads, num_labels, test_spectra, test_spectra_errors):
    # set number of threads per CPU
    os.environ['OMP_NUM_THREADS'] = '{:d}'.format(1)

    # =====================================================================
    # number of processor
    num_CPU = threads

    # ======================================================================
    # testing spectra

    Y_u_all = test_spectra.T
    Y_u_all_err = test_spectra_errors.T

    # ======================================================================
    # load NN results
    w_array_0 = payne_status["w_array_0"]
    w_array_1 = payne_status["w_array_1"]
    w_array_2 = payne_status["w_array_2"]
    b_array_0 = payne_status["b_array_0"]
    b_array_1 = payne_status["b_array_1"]
    b_array_2 = payne_status["b_array_2"]
    x_min = payne_status["x_min"]
    x_max = payne_status["x_max"]

    # -----------------------------------------------------------------------
    # define sigmoid function
    def sigmoid_def(z):
        return 1.0 / (1.0 + np.exp(-z))

    # =======================================================================
    # make spectroscopic mask

    # if desired, bodge some of the errors to be huge in order to create a mask

    # ===========================================================================
    # fit best models
    def fit_func(input_param, *labels):
        predict_flux = w_array_2 * sigmoid_def(np.sum(w_array_1 * (sigmoid_def(np.dot(
            w_array_0, labels[:-1]) + b_array_0)), axis=1) + b_array_1) \
                       + b_array_2

        # perform radial velocity shift
        # f_interp = interpolate.interp1d(wavelength_template, predict_flux,
        #                                 bounds_error=False, kind="linear", fill_value="extrapolate")
        # return f_interp(wavelength_template + labels[-1] * wavelength_template / 10 ** 5)

        return predict_flux

    # ---------------------------------------------------------------------------
    # define function to perform testing step in batch
    def fit_spectrum(spec_no):
        p0_test = np.zeros(num_labels)

        # set bounds
        bounds = np.zeros((num_labels, 2))
        bounds[:, 0] = -0.5
        bounds[:, 1] = 0.5
        bounds[-1, 0] = -5.
        bounds[-1, 1] = 5.

        try:
            popt, pcov = curve_fit(fit_func, [spec_no], Y_u_all[:, spec_no],
                                   p0=p0_test,
                                   sigma=Y_u_all_err[:, spec_no],
                                   absolute_sigma=True, bounds=bounds.T)
        except IOError:  # doctored this because I want to see exceptions for now...
            popt = np.zeros(num_labels) - 9999.

        return popt

    # ============================================================================
    # fit spectra
    pool = Pool(num_CPU)
    recovered_results = np.array(pool.map(fit_spectrum, range(Y_u_all.shape[1]))).T

    # -------------------------------------------------------------------------------
    # initiate chi^2
    chi2 = []

    # loop over all spectra
    for j in range(recovered_results.shape[1]):
        predict_flux = w_array_2 * sigmoid_def(np.sum(w_array_1 * (sigmoid_def(np.dot(
            w_array_0, recovered_results[:-1, j]) + b_array_0)), axis=1) + b_array_1) \
                       + b_array_2

        # radial velocity
        # f_interp = interpolate.interp1d(wavelength_template, predict_flux,
        #                                 bounds_error=False, kind="linear", fill_value="extrapolate")
        # predict_flux = f_interp(wavelength_template \
        #                         + recovered_results[-1, j] * wavelength_template / 10 ** 5)
        chi2.append(np.mean((predict_flux - Y_u_all[:, j]) ** 2 / (Y_u_all_err[:, j] ** 2)))
    chi2 = np.array(chi2)

    # ----------------------------------------------------------------------------
    # rescale back to original values
    for i in range(recovered_results.shape[0] - 1):
        ind_invalid = (recovered_results[i, :] < -100.)
        recovered_results[i, :] = (recovered_results[i, :] + 0.5) * (x_max[i] - x_min[i]) + x_min[i]
        recovered_results[i, ind_invalid] = -999.

    # save array
    return {
        'results': recovered_results.T,
        'chi2': chi2,
        'num_pix': np.sum(Y_u_all_err[:, 0] != 999.)
    }
