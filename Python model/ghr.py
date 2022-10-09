# Author: Yoonsik Oh
# Project: SUBPRO Blue hydrogen

#import numpy as np
from enthalpy import const, enthalpy
from Keq import K_smr, K_wgsr
import scipy.optimize as opt
import autograd.numpy as np
from autograd import grad, jacobian


# inlet/outlet = [molar flow, T, P]
# molar flow = CH4, H2O, H2, CO, CO2
Mm = np.array([16, 18, 2, 28, 44])

def GHR_Reactor(inlet, outlet):
    n0CH4, n0H2O, n0H2, n0CO, n0CO2, T0 = inlet
    nCH4, nH2O, nH2, nCO, nCO2, Q = outlet
    ksi1 = n0CH4 - nCH4
    ksi2 = nCO2 - n0CO2
    ntot = nCH4 + nH2O + nH2 + nCO + nCO2
    T = 973

    eq1 = K_smr(T)*((nCH4/ntot) * (nH2O/ntot )) - (((nCO/ntot) * (nH2/ntot) ** 3))

    eq2 = K_wgsr(T)*((nCO/ntot) * (nH2O/ntot)) - (((nCO2/ntot) * (nH2/ntot)))

    eq3 = nH2O - n0H2O + ksi1 + ksi2

    eq4 = nH2 - n0H2 - 3 * ksi1 - ksi2

    eq5 = nCO - n0CO - ksi1 + ksi2

    #MB = np.array([eq1, eq2, eq3, eq4, eq5])

    eq6 = np.dot(inlet[:5],enthalpy(const, T0)) - np.dot(outlet[:5],enthalpy(const, T)) + Q

    #EB = np.array([eq6])


    balances = np.array([eq1,eq2,eq3,eq4,eq5,eq6])

    return balances


inletstream = np.array([9.9508860e+02, 2.2681504e+03, 1.4973120e+02, 3.4660000e-01,
                  5.2683200e+01, 753])
outletstream = np.array([4000.0, 6000.0, 100.0, 100.0, 100.0, 28000.0,])


def function(outlet):
    return GHR_Reactor(inletstream,outlet)

def solve(f,guess,grad_function):
    sol = opt.fsolve(f,guess,fprime=grad_function)
    print(sol)

#plswork = outletstream

#grad_function = jacobian(function)
#grad_function(plswork)



#sol = opt.fsolve(function,plswork,fprime=grad_function)
#Mm @ inletstream[:5]

#print(function(sol))
#solve(function,plswork,grad_function)

#n7CH4, n7H2O, n7H2, n7CO, n7CO2, P7, Q3 = sol