import scipy
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys
from pylab import *


def readslice(inputfilename, nx, ny):
    f = open(inputfilename, 'rb')
    f.seek(4)
    field = np.fromfile(f, dtype='d', count=int(nx * ny))
    f.close()
    return field


dir = './data/output/'
output_dir = './images/'
STR = 'ww.'
num_procs = 40
outnum_st = 1
outnum_nd = 20
outnum = 6
reso = 512
nx = reso
ny = int(reso / num_procs)
cut = reso - ny * num_procs
print('reso=', reso)
print('reso/nprocs=', reso / num_procs)
print('ny1=', ny + 1)
print('ny2=', ny)
print('cut=', cut)
print(reso, '=', cut, '*', ny + 1, '+', num_procs - cut, '*', ny)


# outnum = 0
# data1 = []
# data2 = np.zeros((reso, reso))
# for islice in range(cut):
#     print('#', islice, 'hd2D' + 'fw.' + str("%03d" %
#                                             islice) + '.' + str("%03d" %
#                                                                 outnum) + '.out', size(data1))
#     f = open(dir + 'hd2D' + 'fw.' + str("%03d" %
#                                         islice) + '.' + str("%03d" %
#                                                             outnum) + '.out', 'rb')
#     f.seek(4)
#     data1 = np.fromfile(f, dtype='d', count=int(nx * (ny + 1)))
#     f.close()
#     for r in range(nx * (ny + 1)):
#         i = r - int(r / nx) * nx
#         j = int(r / nx) + (ny + 1) * islice
#         data2[i, j] = data1[r]
#         # print(r,i,j,data2[i,j])
#     print('****************************')
# for islice in range(cut, num_procs):
#     print('#', islice, 'hd2D' + 'fw.' + str("%03d" %
#                                             islice) + '.' + str("%03d" %
#                                                                 outnum) + '.out', size(data1))
#     f = open(dir + 'hd2D' + 'fw.' + str("%03d" %
#                                         islice) + '.' + str("%03d" %
#                                                             outnum) + '.out', 'rb')
#     f.seek(4)
#     data1 = np.fromfile(f, dtype='d', count=int(nx * (ny)))
#     f.close()
#     for r in range(nx * ny):
#         i = r - int(r / nx) * nx
#         j = int(r / nx) + ny * (islice - cut) + (ny + 1) * cut
#         data2[i, j] = data1[r]
# figure(1)
# im1 = plt.imshow(data2, cmap=cm.cool)
# cbar1 = plt.colorbar(im1)
# # plt.savefig('FlowD_'+STR+str("%03d" % outnum)+'png' ) #% (dtemp,stemp,qname))
# plt.show()

for outnum in range(outnum_st, outnum_nd):
    data1 = []
    data2 = np.zeros((reso, reso))
    for islice in range(cut):
        print('#', islice, 'hd2D' + STR + str("%03d" %
                                              islice) + '.' + str("%03d" %
                                                                  outnum) + '.out', size(data1))
        f = open(dir + 'hd2D' + STR + str("%03d" %
                                          islice) + '.' + str("%03d" %
                                                              outnum) + '.out', 'rb')
        f.seek(4)
        data1 = np.fromfile(f, dtype='d', count=int(nx * (ny + 1)))
        f.close()
        for r in range(nx * (ny + 1)):
            i = r - int(r / nx) * nx
            j = int(r / nx) + (ny + 1) * islice
            data2[i, j] = data1[r]
            # print(r,i,j,data2[i,j])
    print('****************************')
    for islice in range(cut, num_procs):
        print('#', islice, 'hd2D' + STR + str("%03d" %
                                              islice) + '.' + str("%03d" %
                                                                  outnum) + '.out', size(data1))
        f = open(dir + 'hd2D' + STR + str("%03d" %
                                          islice) + '.' + str("%03d" %
                                                              outnum) + '.out', 'rb')
        f.seek(4)
        data1 = np.fromfile(f, dtype='d', count=int(nx * (ny)))
        f.close()
        for r in range(nx * ny):
            i = r - int(r / nx) * nx
            j = int(r / nx) + ny * (islice - cut) + (ny + 1) * cut
            data2[i, j] = data1[r]
            # print(r,i,j,data2[i,j])
    figure(2)
    im1 = plt.imshow(data2, cmap=cm.hot)
    cbar1 = plt.colorbar(im1)
    plt.savefig(output_dir + 'FlowD_' + STR +
                str("%03d" % outnum) + '.png')  # %
    # (dtemp,stemp,qname))
    # plt.show()
