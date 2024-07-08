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


# use tex for the labels
plt.rc("text", usetex=True)

dir = "./data/output/"
output_dir = "./images/"
STR = "ww."
color = "RdBu_r"
homogeneous = True  # True if we want the same colors for all images
only_one_file = True
myfile = 35
first_file = 1
FONTSIZE = 16

dim_dir = "./data/dim.txt"
reso, num_procs = get_dim(dim_dir)
outnum_nd = get_num_files(dir, STR) // num_procs  # we have ww, ps, fw and fp
if only_one_file:
    outnum_nd = 1
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

data = np.zeros((outnum_nd, reso, reso))

for file in range(outnum_nd + 1 - first_file):
    data1 = []
    data2 = np.zeros((reso, reso))
    for islice in range(cut):
        if only_one_file:
            filename = (
                dir
                + "hd2D"
                + STR
                + str("%03d" % islice)
                + "."
                + str("%03d" % myfile)
                + ".out"
            )
        else:
            filename = (
                dir
                + "hd2D"
                + STR
                + str("%03d" % islice)
                + "."
                + str("%03d" % (myfile + first_file))
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
        if only_one_file:
            filename = (
                dir
                + "hd2D"
                + STR
                + str("%03d" % islice)
                + "."
                + str("%03d" % myfile)
                + ".out"
            )
        else:
            filename = (
                dir
                + "hd2D"
                + STR
                + str("%03d" % islice)
                + "."
                + str("%03d" % (myfile + first_file))
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
    data[file, :, :] = data2

    print("file=", file)


if homogeneous:
    zmax = np.max(data)
    zmin = np.min(data)

    zmax = max(abs(zmax), abs(zmin))

    zmax /= 2
    zmin = -zmax

for file in range(outnum_nd + 1 - first_file):
    fig = plt.figure()
    ax = fig.add_subplot()
    data2 = data[file, :, :]
    if homogeneous:
        myplot = ax.imshow(
            data2,
            cmap=color,
            vmin=zmin,
            vmax=zmax,
            origin="lower",
            interpolation="nearest",
        )
    else:
        myplot = ax.imshow(data2, cmap=color, origin="lower", interpolation="nearest")

    # add circle to represent the annulus
    circle1 = plt.Circle(
        (reso / 2, reso / 2),
        reso / 4,
        color="black",
        fill=False,
        linestyle="--",
        linewidth=1,
    )
    circle2 = plt.Circle(
        (reso / 2, reso / 2),
        reso / 4 + 70,
        color="black",
        fill=False,
        linestyle="--",
        linewidth=1,
    )
    ax.add_artist(circle1)
    ax.add_artist(circle2)

    # rescale the domain to [-pi, pi] x [-pi, pi]
    ax.set_xlim([0, reso])
    ax.set_ylim([0, reso])
    ax.set_xticks([0, reso / 4, reso / 2, 3 * reso / 4, reso])
    ax.set_yticks([0, reso / 4, reso / 2, 3 * reso / 4, reso])
    ax.set_xticklabels(
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"], fontsize=FONTSIZE
    )
    ax.set_yticklabels(
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"], fontsize=FONTSIZE
    )

    ax.set_xlabel(r"$x$", fontsize=FONTSIZE)
    ax.set_ylabel(r"$y$", fontsize=FONTSIZE)

    # add colorbar
    fig.tight_layout()
    # cbar = plt.colorbar(myplot)
    if only_one_file:
        filename = output_dir + "FlowD_" + STR + str("%03d" % myfile) + "_annulus.pdf"
    else:
        filename = (
            output_dir + "FlowD_" + STR + str("%03d" % (file + first_file)) + ".pdf"
        )
    plt.savefig(filename, bbox_inches="tight")
    plt.close()
