import numpy as np
import matplotlib.pyplot as plt
import os


def readslice(inputfilename, nx, ny):
    f = open(inputfilename, "rb")
    f.seek(4)
    field = np.fromfile(f, dtype="d", count=int(nx * ny))
    f.close()
    return field


def get_num_files(folder_path, STR):
    num_files = [
        f
        for f in os.listdir(folder_path)
        if os.path.isfile(os.path.join(folder_path, f))
    ]
    num_files = len([f for f in num_files if STR in f])
    return num_files


def get_dim(dim_dir):
    f = open(dim_dir, "r")
    dim_str = f.readlines()
    f.close()
    # split the line and return the two values
    dim, num_procs = dim_str[0].split()
    return int(dim), int(num_procs)


input_dir = "./data/output/"
output_dir = "./images/"
STR = "ww."
core = 0
dim_dir = "./data/dim.txt"
homogeneous = False  # True if we want the same colors for all images
first_file = 1


reso, num_procs = get_dim(dim_dir)
outnum_nd = get_num_files(input_dir, STR) // num_procs  # we have ww, ps, fw and fp
nx = reso
ny = nx
color = "RdBu_r"
print("reso =", reso)
print("num_procs =", num_procs)
print("core =", core)
print("outnum_nd=", outnum_nd)
script_dir = os.path.dirname(__file__)

data = np.zeros((outnum_nd, reso, reso))


for file in range(outnum_nd):
    data1 = []
    data2 = np.zeros((reso, reso))
    filename = (
        input_dir
        + "hd2D"
        + STR
        + str("%03d" % core)
        + "."
        + str("%03d" % (file + first_file))
        + ".out"
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

    print("file=", file)

if homogeneous:
    zmax = np.max(data)
    zmin = np.min(data)

    zmax /= 2
    zmin /= 2

for file in range(outnum_nd):
    fig = plt.figure()
    ax = fig.add_subplot()
    data2 = data[file, :, :]
    if homogeneous:
        myplot = ax.imshow(data2, cmap=color, vmin=zmin, vmax=zmax, origin="lower")
    else:
        myplot = ax.imshow(data2, cmap=color, origin="lower")
    # add colorbar
    cbar = plt.colorbar(myplot)
    filename = (
        output_dir
        + "FlowD_"
        + STR
        + str("%03d" % core)
        + "."
        + str("%03d" % (file + first_file))
        + ".png"
    )
    plt.savefig(filename)
    plt.close()
