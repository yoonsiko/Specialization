import autograd.numpy as np

# Order of coefficients: A B C
# Order of components: C2H6, C3H8, n-C4H10, i-C4H10, C5+

cp_const = np.array([[1.702, 9.081, -2.164],
                     [3.470, 1.450, 0.000],
                     [3.249, 0.422, 0.000],
                     [3.376, 0.557, 0.000],
                     [5.457, 0.557, 0.000],
                     [1.131, 19.225, -5.561],
                     [1.213, 28.785, -8.824],
                     [1.935, 36.915, -11.402],
                     [1.935, 36.915, -11.402],
                     [2.464, 45.351, -14.111]])

H2Ocoeff = np.array([3.470, 1.450, 0.000])

def cp(v, T):  # J/mol*K
    cp = v[:,0] + v[:,1]*T*10**(-3) + v[:,2]*T**2*10**(-6)
    return cp

def cp_print(f,v,T):
    cp = f(v,T)
    print("CH4:",cp[0],"\nH2O:",cp[1],"\nH2:",cp[2],"\nCO:",cp[3],"\nCO2:",cp[4],"\nC2H6:",cp[5],"\nC3H8:",cp[6],"\nn-C4H10:",cp[7],"\ni-C4H10:", cp[8],"\nC5+:",cp[9])


#print(cp_print(cp, cp_const, 298.15))

def single_cp(coeff, T):
    cp = coeff[0] + coeff[1]*T*10**(-3) + coeff[2]*T**2*10**(-6)
    return cp