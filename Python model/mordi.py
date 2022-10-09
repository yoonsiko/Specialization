#import numpy as np
from enthalpy import const, enthalpy
import autograd.numpy as np
test = []
for i in range(117):
    test.append(i)

print(test)

print(test[7:12])
A = np.array([1,2,3,4,5])
B = np.array([0,0,0,0,0])

C = np.append(A,B)
print(C)

inletstream = np.array([78.24,0,0,0,1.34,6.10,6.7,2.48,1.41,3.7])*3424/100
other = np.array([30,311])
inletstream = np.append(inletstream,other)
print(inletstream)

a,c,v = np.array([1,2,3])
