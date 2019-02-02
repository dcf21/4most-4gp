# -*- coding: utf-8 -*-

import logging
from multiprocessing import Pool
import numpy as np
import torch
from torch.autograd import Variable


def train_pixel(params):
    pixel_no = params[0]
    dim_in = params[1]
    x = params[2]
    x_valid = params[3]
    y_row = params[4]
    y_valid_row = params[5]

    logger = logging.getLogger(__name__)
    logger.info("Training pixel {:6d}".format(pixel_no))

    # define neural network
    neuron_count = 3

    model = torch.nn.Sequential(
        torch.nn.Linear(dim_in, neuron_count),
        torch.nn.Sigmoid(),
        torch.nn.Linear(neuron_count, 1),
        torch.nn.Sigmoid(),
        torch.nn.Linear(1, 1)
    )

    # define optimizer
    learning_rate = 0.01  # Yuan-Sen set this to 0.001. Doctored by dcf for speed.
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    # ==============================================================================
    # convergence counter
    current_loss = np.inf
    count = 0
    t = 0

    # -----------------------------------------------------------------------------
    # train the neural network
    while count < 5:  # Yuan-Sen set this to 20

        # training
        y_pred = model(x)[:, 0]
        loss = ((y_pred - y_row).pow(2) / (0.01 ** 2)).mean()

        # validation
        y_pred_valid = model(x_valid)[:, 0]
        loss_valid = (((y_pred_valid - y_valid_row).pow(2)
                       / (0.01 ** 2)).mean()).item()

        # -----------------------------------------------------------------------------
        # check convergence

        # Set number of iterations of optimizer to run between checking progress. Yuan-Sen set to 10,000
        if t % 400 == 0:
            if loss_valid >= current_loss:
                count += 1
            else:
                count = 0
                # record the best loss
                current_loss = loss_valid

                # record the best parameters
                model_numpy = []
                for param in model.parameters():
                    model_numpy.append(param.data.numpy())

        # -----------------------------------------------------------------------------
        # optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        t += 1

    # -----------------------------------------------------------------------------
    # return parameters
    return model_numpy

    # =============================================================================


def train_nn(threads, batch_number, batch_count, labelled_set, normalized_flux, normalized_ivar, dispersion):
    """
    Train the neural network

    :param batch_number:
        If training pixels in multiple batches on different machines, then this is the number of the batch of pixels
        we are to train. It should be in the range 0 .. batch_count-1 inclusive.

    :param batch_count:
        If training pixels in multiple batches on different machines, then this is the number of batches.

    :param labelled_set:
        2D numpy array containing the label values for the training set.
        labelled_set[spectrum_num][label_num] = label_value

    :param normalized_flux:
        2D numpy array containing the continuum-normalised fluxes for the training stars.
        normalized_flux[spectrum_num][pixel_num] = flux

    :param normalized_ivar:
        2D numpy array containing the uncertainty in the normalized_flux.
        normalized_ivar[spectrum_num][pixel_num] = uncertainty

    :param dispersion:
        1D numpy array containing the wavelength (in Angstrom) associated with each pixel.

    :return:
        The optimized neural network weights.
        output['w_array_0'][pixel_number] = list of w0 weights for that pixel
    """

    # set number of threads per CPU
    # mkl.set_num_threads(1)
    torch.set_num_threads(1)

    # ------------------------------------------------------------------------------
    # number of CPUs for parallel computing
    num_CPU = threads

    # ------------------------------------------------------------------------------
    # Rearrange storage of training data in memory to be [pixel_n][spectrum_id]
    # This means that training each pixel involves data stored sequentially in memory!
    normalized_flux = normalized_flux.T.copy().T

    # ==============================================================================
    # restore training spectra
    x = labelled_set  # [spectrum_id][label_n]
    y = normalized_flux  # [spectrum_id][pixel_n]

    # and validation spectra (fudge for now)
    x_valid = x
    y_valid = y

    # scale the labels
    x_max = np.max(x, axis=0)
    x_min = np.min(x, axis=0)
    x = (x - x_min) / (x_max - x_min) - 0.5
    x_valid = (x_valid - x_min) / (x_max - x_min) - 0.5

    # -----------------------------------------------------------------------------
    # dimension of the input
    dim_in = x.shape[1]
    num_pix = y.shape[1]

    # make pytorch variables
    x = Variable(torch.from_numpy(x)).type(torch.FloatTensor)
    y = Variable(torch.from_numpy(y), requires_grad=False).type(torch.FloatTensor)
    x_valid = Variable(torch.from_numpy(x_valid)).type(torch.FloatTensor)
    y_valid = Variable(torch.from_numpy(y_valid),
                       requires_grad=False).type(torch.FloatTensor)

    # =============================================================================
    # loop over all pixels

    # Work out which batch of pixels we are to work on
    pixel_start = (num_pix // batch_count + 1) * batch_number
    pixel_end = min(num_pix, (num_pix // batch_count + 1) * (batch_number + 1))


    # train in parallel
    with Pool(num_CPU) as pool:
      net_array = pool.map(train_pixel, [[i, dim_in, x, x_valid, y[:, i], y_valid[:, i]]
                                         for i in range(pixel_start, pixel_end)])

    # train in serial mode
    # net_array = []
    # for i in range(pixel_start, pixel_end):
    #     net_array.append(train_pixel([i, dim_in, x, x_valid, y[:, i], y_valid[:, i]]))

    # extract parameters
    w_array_0 = np.array([net_array[i][0] for i in range(len(net_array))])
    b_array_0 = np.array([net_array[i][1] for i in range(len(net_array))])
    w_array_1 = np.array([net_array[i][2][0] for i in range(len(net_array))])
    b_array_1 = np.array([net_array[i][3][0] for i in range(len(net_array))])
    w_array_2 = np.array([net_array[i][4][0][0] for i in range(len(net_array))])
    b_array_2 = np.array([net_array[i][5][0] for i in range(len(net_array))])

    # save parameters and remember how we scale the labels
    return {
        'w_array_0': w_array_0,
        'w_array_1': w_array_1,
        'w_array_2': w_array_2,
        'b_array_0': b_array_0,
        'b_array_1': b_array_1,
        'b_array_2': b_array_2,
        'x_max': x_max,
        'x_min': x_min
    }
