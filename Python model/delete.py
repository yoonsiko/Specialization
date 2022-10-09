import numpy as np

satisfaction = np.array([92,92,89,91,79,75,80,84,84,81,77,91,89,87])
abandon = np.array([14,13,11,13,11,27,15,5,7,17,7,3,6,4])

average = 0
new_satisfaction = np.zeros(len(satisfaction))
for i in range(len(satisfaction)):
    new_satisfaction[i] = (satisfaction[i] + abandon[i]*0.7)

call = np.array([25,45,67,56,47,87,94,63,53,89,94,34,45,35])

for i in range(len(new_satisfaction)):
    if call[i] > 80:
        new_satisfaction[i] = (satisfaction[i] - 0.4*(call[i]-80))
    else:
        new_satisfaction[i] = satisfaction[i]

print(new_satisfaction)

print((sum(new_satisfaction)/len(new_satisfaction))-(sum(satisfaction)/len(satisfaction)))

print((0.4*(7+14+9+14))/4)