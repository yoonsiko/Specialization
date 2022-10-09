import scipy.optimize as opt
import numpy as np
from enthalpy import const, enthalpy, enthalpy_heavy, heavy_const

n2H2O = 8440.46
T2 = 423.0
n3stream = 2676.81 + 8442.14 + 45.52 + 208.86 + 229.41 + 84.92 + 48.28 + 126.69
inletstream = np.array([78.24, 0.0, 0.0, 0.0, 1.34, 6.10, 6.7, 2.48, 1.41, 3.7]) * 3424 / 100
other = np.array([30.0, 311.0,n2H2O,T2])

inletstream = np.append(inletstream, other)


def mixingTemp(input,guess):
    n1CH4, n1H2O, n1H2, n1CO, n1CO2, n1C2H6, n1C3H8, n1i_C4H10, n1n_C4H10, n1_C5, P1, T1, n2H2O, T2 = input
    n3stream = input[:10]+np.array([0, n2H2O, 0, 0, 0, 0, 0, 0, 0, 0])
    T3 = guess[0]

    eq1 = np.dot(input[:10], enthalpy_heavy(heavy_const, T1)) + \
          np.dot(np.array([0, n2H2O, 0, 0, 0, 0, 0, 0, 0, 0]), enthalpy_heavy(heavy_const, T2)) - \
          np.dot(n3stream, enthalpy_heavy(heavy_const, T3))

    return eq1

def T3function(outlet):
    return mixingTemp(inletstream,outlet)

#x0 = np.array([400])
#sol = opt.fsolve(T3function,x0)
#print(sol)