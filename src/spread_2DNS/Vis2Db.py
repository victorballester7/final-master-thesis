import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import pylab


def readslice(inputfilename, nx, ny):
    f = open(inputfilename, 'rb')
    f.seek(4)
    field = np.fromfile(f, dtype='d', count=int(nx * ny))
    f.close()
    return field


def get_num_files(folder_path):
    num_files = len([f for f in os.listdir(folder_path)
                     if os.path.isfile(os.path.join(folder_path, f))])
    return num_files


def get_dim(dim_dir):
    f = open(dim_dir, 'r')
    dim = f.readlines()
    f.close()
    return int(dim[0])


dir = './data/output/'
output_dir = './images/'
STR = 'ww.'
num_procs = 40
outnum_nd = get_num_files(dir) // (num_procs * 4)  # we have ww, ps, fw and fp
color = 'RdBu_r'

dim_dir = './data/dim.txt'
reso = get_dim(dim_dir)
nx = reso
ny = int(reso / num_procs)
cut = reso - ny * num_procs
dim_dir = './data/dim.txt'
print('reso=', reso)
print('reso/nprocs=', reso / num_procs)
print('ny1=', ny + 1)
print('ny2=', ny)
print('cut=', cut)
print(reso, '=', cut, '*', ny + 1, '+', num_procs - cut, '*', ny)

data = np.zeros((outnum_nd, reso, reso))


for file in range(outnum_nd):
    data1 = []
    data2 = np.zeros((reso, reso))
    for islice in range(cut):
        filename = dir + 'hd2D' + STR + \
            str("%03d" % islice) + '.' + str("%03d" % file) + '.out'
        f = open(filename, 'rb')
        f.seek(4)
        data1 = np.fromfile(f, dtype='d', count=int(nx * (ny + 1)))
        f.close()
        for r in range(nx * (ny + 1)):
            i = r - int(r / nx) * nx
            j = int(r / nx) + (ny + 1) * islice
            data2[i, j] = data1[r]
    for islice in range(cut, num_procs):
        filename = dir + 'hd2D' + STR + \
            str("%03d" % islice) + '.' + str("%03d" % file) + '.out'
        f = open(filename, 'rb')
        f.seek(4)
        data1 = np.fromfile(f, dtype='d', count=int(nx * (ny)))
        f.close()
        for r in range(nx * ny):
            i = r - int(r / nx) * nx
            j = int(r / nx) + ny * (islice - cut) + (ny + 1) * cut
            data2[i, j] = data1[r]
    data[file, :, :] = data2

zmax = np.max(data)
zmin = np.min(data)

zmax /= 2
zmin /= 2

for file in range(outnum_nd):
    print('****************************')
    print('#', 'hd2D' + STR +
          str("%03d" % file) + '.out', pylab.size(data1))
    fig = plt.figure()
    ax = fig.add_subplot()
    data2 = data[file, :, :]
    myplot = ax.imshow(data2, cmap=color, vmin=zmin, vmax=zmax, origin='lower')
    # add colorbar
    cbar = plt.colorbar(myplot)
    filename = output_dir + 'FlowD_' + STR + \
        str("%03d" % file) + '.png'
    plt.savefig(filename)
print('****************************')
