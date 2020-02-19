import numpy as np
import sys
import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar, exp
from sklearn.preprocessing import minmax_scale
from extensions.utils.data_utils import lineprof, model_lineprof


def fit_yds(env):
    """
    Fits a sinc**2 function to the far-field diffraction pattern from a 
    double slit arrangement
    
    :param env: envelope of far field diffraction line profile
    """

    x = np.linspace(0, 10, len(env))

    k = env.argmax()

    def func(x, a, b, c, d):
        return (a * np.sinc(b * x - c) + d) ** 2

    popt, pcov = curve_fit(func, x, env, p0=[env[k], 1, x[k], np.min(env)])

    fit = func(x, *popt)

    return fit
