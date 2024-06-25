import numpy as np

data = np.loadtxt("./NumVortices_dr400.txt")

# merge data of the second column (summing them) each 4 rows
# like this
# x0 y0
# x1 y1
# x2 y2
# x3 y3
# x4 y4
# x5 y5
# x6 y6
# x7 y7
# x8 y8
# ...

# to this
# x3 y0+y1+y2+y3
# x7 y4+y5+y6+y7
# x11 y8+y9+y10+y11
# ...

aux = 0
new_data = []
divider = 8
for i in range(len(data)):
    if i % divider == divider - 1:
        new_data.append([data[i, 0], aux + data[i, 1]])
        aux = 0
    else:
        aux += data[i, 1]

new_data = np.array(new_data)

# save the result in 'NumVortices.txt'
np.savetxt("NumVortices.txt", new_data)
