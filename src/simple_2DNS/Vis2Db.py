import numpy as np
import matplotlib.pyplot as plt
import os
import pylab


def readslice(inputfilename, nx, ny):
    f = open(inputfilename, "rb")
    f.seek(4)
    field = np.fromfile(f, dtype="d", count=int(nx * ny))
    f.close()
    return field


def get_num_files(folder_path):
    num_files = len(
        [
            f
            for f in os.listdir(folder_path)
            if os.path.isfile(os.path.join(folder_path, f))
        ]
    )
    return num_files


def get_dim(dim_dir):
    f = open(dim_dir, "r")
    dim = f.readlines()
    f.close()
    return int(dim[0])


input_dir = "data/simple_2DNS/output/"
output_dir = "images/simple_2DNS/"
STR = "fw."
outnum_nd = get_num_files(input_dir) // 4  # we have ww, ps, fw and fp
homogeneous = False  # True if we want the same colors for all images

dim_dir = "data/simple_2DNS/dim.txt"
reso = get_dim(dim_dir)
nx = reso
ny = nx
color = "RdBu_r"
print("reso=", reso)
script_dir = os.path.dirname(__file__)

data = np.zeros((outnum_nd, reso, reso))

for file in range(outnum_nd):
    data1 = []
    data2 = np.zeros((reso, reso))
    filename = os.path.join(
        script_dir, "../../" + input_dir + "hd2D" + STR + str("%03d" % file) + ".out"
    )
    f = open(filename, "rb")
    f.seek(4)
    data1 = np.fromfile(f, dtype="d", count=int(nx * (ny)))
    f.close()
    for r in range(nx * ny):
        i = r - int(r / nx) * nx
        j = int(r / nx)
        data2[i, j] = data1[r]
    data[file, :, :] = data2

if homogeneous:
    zmax = np.max(data)
    zmin = np.min(data)

    zmax /= 2
    zmin /= 2

for file in range(outnum_nd):
    print("****************************")
    print("#", "hd2D" + STR + str("%03d" % file) + ".out", pylab.size(data[file, :, :]))
    fig = plt.figure()
    ax = fig.add_subplot()
    data2 = data[file, :, :]
    if homogeneous:
        myplot = ax.imshow(data2, cmap=color, vmin=zmin, vmax=zmax, origin="lower")
    else:
        myplot = ax.imshow(data2, cmap=color, origin="lower")
    # add colorbar
    cbar = plt.colorbar(myplot)
    filename = os.path.join(
        script_dir, "../../" + output_dir + "FlowD_" + STR + str("%03d" % file) + ".png"
    )
    plt.savefig(filename)
    plt.close()
print("****************************")
