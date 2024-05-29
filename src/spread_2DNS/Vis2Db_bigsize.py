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


dir = "./data/output/"
output_dir = "./images/"
STR = "ww."
color = "RdBu_r"
homogeneous = False  # True if we want the same colors for all images
first_file = 1

dim_dir = "./data/dim.txt"
reso, num_procs = get_dim(dim_dir)
outnum_nd = get_num_files(dir, STR) // num_procs  # we have ww, ps, fw and fp
nx = reso
ny = int(reso / num_procs)
cut = reso - ny * num_procs
dim_dir = "./data/dim.txt"
print("reso=", reso)
print("num_procs=", num_procs)
print("outnum_nd=", outnum_nd)
print("reso/nprocs=", reso / num_procs)
print("ny1=", ny + 1)
print("ny2=", ny)
print("cut=", cut)
print(reso, "=", cut, "*", ny + 1, "+", num_procs - cut, "*", ny)

zmax, zmin = 0, 0

for file in range(outnum_nd + 1 - first_file):
    data1 = []
    data2 = np.zeros((reso, reso))
    for islice in range(cut):
        filename = (
            dir
            + "hd2D"
            + STR
            + str("%03d" % islice)
            + "."
            + str("%03d" % (file + first_file))
            + ".out"
        )
        f = open(filename, "rb")
        f.seek(4)
        data1 = np.fromfile(f, dtype="d", count=int(nx * (ny + 1)))
        f.close()
        for r in range(nx * (ny + 1)):
            i = r - int(r / nx) * nx
            j = int(r / nx) + (ny + 1) * islice
            data2[i, j] = data1[r]
    for islice in range(cut, num_procs):
        filename = (
            dir
            + "hd2D"
            + STR
            + str("%03d" % islice)
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
            j = int(r / nx) + ny * (islice - cut) + (ny + 1) * cut
            data2[i, j] = data1[r]

    if file == 0:
        zmax = np.max(data2)
        zmin = np.min(data2)
        zmax = max(abs(zmax), abs(zmin))
        zmax = zmax / 4
        zmin = -zmax

    fig = plt.figure()
    ax = fig.add_subplot()
    myplot = ax.imshow(data2, cmap=color, vmin=zmin, vmax=zmax, origin="lower")
    # add colorbar
    cbar = plt.colorbar(myplot)
    filename = output_dir + "FlowD_" + STR + str("%03d" % (file + first_file)) + ".png"
    plt.savefig(filename, dpi=300)
    plt.close()

    print("file=", file)
