# -*- coding: utf-8 -*-

from multiprocessing import Pool

# import mkl
import numpy as np
import logging
import torch
from torch.autograd import Variable

# This is a bit of a fudge, but prevents pytorch from failing with...
# OSError: [Errno 24] Too many open files: '/tmp/pymp-kucaes4t' 
# import torch.multiprocessing
# torch.multiprocessing.set_sharing_strategy('file_system')

def train_pixel(params):
    pixel_no = params[0]
    dim_in = params[1]
    x = params[2]
    x_valid = params[3]
    y_row = params[4]
    y_valid_row = params[5]

    # define neural network
    neuron_count = 1

    model = torch.nn.Sequential(
        torch.nn.Linear(dim_in, neuron_count),
        torch.nn.Sigmoid(),
        torch.nn.Linear(neuron_count, 1),
        torch.nn.Sigmoid(),
        torch.nn.Linear(1, 1)
    )

    # define optimizer
    learning_rate = 0.001
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    # ==============================================================================
    # convergence counter
    current_loss = np.inf
    count = 0
    t = 0

    # -----------------------------------------------------------------------------
    # train the neural network
    while count < 10:

        # training
        y_pred = model(x)[:, 0]
        loss = ((y_pred - y_row).pow(2) / (0.01 ** 2)).mean()

        # validation
        y_pred_valid = model(x_valid)[:, 0]
        loss_valid = (((y_pred_valid - y_valid_row).pow(2)
                       / (0.01 ** 2)).mean()).item()

        # -----------------------------------------------------------------------------
        # check convergence
        if t % 1000 == 0:
            if loss_valid >= current_loss:
                count += 1
            else:
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


def train_nn(threads, labelled_set, normalized_flux, normalized_ivar, dispersion):
    logger = logging.getLogger(__name__)

    # set number of threads per CPU
    # mkl.set_num_threads(1)

    # ------------------------------------------------------------------------------
    # number of CPUs for parallel computing
    num_CPU = threads

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

    # train in parallel
    # pool = Pool(num_CPU)
    # net_array = pool.map(train_pixel, [[i, dim_in, x, x_valid, y[:, i], y_valid[:, i]]
    #                                    for i in range(num_pix)])

    for i in range(num_pix):
        logger.info("Training pixel {:6d}/{:6d}".format(i, num_pix))
        train_pixel([i, dim_in, x, x_valid, y[:, i], y_valid[:, i]])

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
