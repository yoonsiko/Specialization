from enthalpy import const, enthalpy
from Keq import K_smr, K_wgsr
import scipy.optimize as opt
import autograd.numpy as np



def itsr(input, output):
    n9CH4, n9H2O, n9H2, n9CO, n9CO2, T9 = input
    n10CH4, n10H2O, n10H2, n10CO, n10CO2, Q4 = output
    n10 = n10CH4 + n10H2O + n10H2 + n10CO + n10CO2
    stream9 = np.array([n9CH4, n9H2O, n9H2, n9CO, n9CO2])
    stream10 = np.array([n10CH4, n10H2O, n10H2, n10CO, n10CO2])


    T10 = 523  # Setting the reactor temperature and solving for the heat required instead
    itsr_ksi1 = n9H2O - n10H2O

    eq1 = n10CH4 - n9CH4
    eq2 = n10H2 - n9H2 - itsr_ksi1
    eq3 = n10CO - n9CO + itsr_ksi1
    eq4 = n10CO2 - n9CO2 - itsr_ksi1
    eq5 = K_wgsr(T10) * ((n10CO / n10) * (n10H2O / n10)) - (((n10CO2 / n10) * (n10H2 / n10)))
    eq6 = np.dot(stream9, enthalpy(const, T9)) - \
           np.dot(stream10, enthalpy(const, T10)) + Q4

    return np.array([eq1,eq2,eq3,eq4,eq5,eq6])