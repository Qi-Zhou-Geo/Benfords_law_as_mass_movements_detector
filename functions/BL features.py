#!/usr/bin/python
# -*- coding: UTF-8 -*-


#__modification time__ = 2024-02-03
#__author__ = Qi Zhou, Helmholtz Centre Potsdam - GFZ German Research Centre for Geosciences
#__find me__ = qi.zhou@gfz-potsdam.de, qi.zhou.geo@gmail.com, https://github.com/Nedasd
# Please do not distribute this code without the author's permission


import numpy as np
import pandas as pd
from scipy.stats import iqr, ks_2samp, mannwhitneyu


# Benford's Law feature is based on the following references
# Sambridge, Malcolm, et al. "Benford's law in the natural sciences." Geophysical research letters 37.22 (2010).
# Zhou, Qi, et al. "Benford's law as mass movement detector in seismic signals." (2023).
def calBL_feature(data, ruler=100):
    '''
    Parameters
    ----------
    data: np array data after dtrend, dmean, or raw data
    ruler: shit the data from 0 to ruler 100. If you use the raw data, please set "ruler" as 0

    Returns: numpy array with shape 1*16
    -------
    '''
    data = np.abs(data)
    # BL theoretical distribution value
    BL_distribution = [0.301, 0.176, 0.125, 0.097, 0.079, 0.067, 0.058, 0.051, 0.046]
    dataSelected = data[data >= 100]

    # <editor-fold desc="iq, max, min">
    iq = float("{:.2f}".format(iqr(dataSelected)))
    max_amp = float("{:.2f}".format(np.max(dataSelected - ruler)))
    min_amp = float("{:.2f}".format(np.min(dataSelected - ruler)))
    # </editor-fold>

    # <editor-fold desc="ovserved first digit frequency">
    # if you have more efficient approach for firstDigit_frequency, please let me know
    amp_data = pd.DataFrame(dataSelected)
    amp_data = amp_data.astype(str)

    d = (amp_data.iloc[:, 0]).str[0: 1]
    d = list(d)

    digit_count = np.empty((0, 9))
    for digit in range(1, 10):
        first_digit = d.count(str(digit))
        digit_count = np.append(digit_count, first_digit)

    firstDigit_frequency = digit_count / np.sum(digit_count)
    firstDigit_frequency = [float('{:.3f}'.format(i)) for i in firstDigit_frequency]
    # </editor-fold>

    # <editor-fold desc="goodness, KS, Mann Whitne U, alpha">
    frequency = np.empty((0, 9))
    for a in range(0, 9):
        first_firstDigit_frequency = pow((firstDigit_frequency[a] - BL_distribution[a]), 2) / BL_distribution[a]
        frequency = np.append(frequency, first_firstDigit_frequency)
    goodness = (1 - pow(sum(frequency), 0.5)) * 100
    goodness = float("{:.3f}".format(goodness))

    statistic, pvalue = ks_2samp(BL_distribution, firstDigit_frequency, alternative='two-sided', method='exact')
    ks = float("{:.3f}".format(pvalue))  # pvalue


    statistic, pvalue = mannwhitneyu(BL_distribution, firstDigit_frequency, alternative='two-sided', method='exact')
    MannWhitneU = float("{:.3f}".format(pvalue))  # pvalue

    if ks >= 0.95 and MannWhitneU >= 0.95:
        follow = 1 # follow BL
    else:
        follow = 0 # do not follow BL

    sum_d = []
    y_min = np.min(dataSelected)
    for s in range(0, len(dataSelected)):
        i = np.log(dataSelected[s] / y_min)
        sum_d.append(i)
    alpha = 1 + len(dataSelected) / np.sum(sum_d)
    alpha = float("{:.4f}".format(alpha))
    # </editor-fold>

    output = np.array([max_amp, min_amp, iq, goodness, alpha, ks, MannWhitneU, follow], dtype=float)
    output = np.append(firstDigit_frequency, output)
    return output
