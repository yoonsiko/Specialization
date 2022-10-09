from enthalpy import const, enthalpy
from Keq import K_smr, K_wgsr
import scipy.optimize as opt
import autograd.numpy as np
from autograd import grad, jacobian


def postATR_temp(input, output):
    n8CH4, n8H2O, n8H2, n8CO, n8CO2, \
    n9CH4, n9H2O, n9H2, n9CO, n9CO2, T8, Q3 = input

    T9 = output

    stream8 = np.array([n8CH4, n8H2O, n8H2, n8CO, n8CO2])
    stream9 = np.array([n9CH4, n9H2O, n9H2, n9CO, n9CO2])

    eq1 = np.dot(stream8, enthalpy(const, T8)) - \
           np.dot(stream9, enthalpy(const, T9)) + Q3

    return eq1

