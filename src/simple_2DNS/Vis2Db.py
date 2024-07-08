import numpy as np
import matplotlib.pyplot as plt
import os
import pylab
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


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
    dim = f.readlines()
    f.close()
    return int(dim[0])


# use latex
plt.rc("text", usetex=True)

input_dir = "data/simple_2DNS/output/"
output_dir = "images/simple_2DNS/"
STR = "fw."
outnum_nd = get_num_files(input_dir, STR)  # we have ww, ps, fw and fp
homogeneous = True  # True if we want the same colors for all images
first_file = 1
paint_arrows = True

FONTSIZE = 16

dim_dir = "data/simple_2DNS/dim.txt"
reso = get_dim(dim_dir)
nx = reso
ny = nx
color = "RdBu_r"
print("reso=", reso)
print("outnum_nd=", outnum_nd)

script_dir = os.path.dirname(__file__)

data = np.zeros((outnum_nd, reso, reso))

for file in range(outnum_nd + 1 - first_file):
    data1 = []
    data2 = np.zeros((reso, reso))
    filename = os.path.join(
        script_dir,
        "../../"
        + input_dir
        + "hd2D"
        + STR
        + str("%03d" % (file + first_file))
        + ".out",
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
    zmax = max(abs(zmax), abs(zmin))
    zmax /= 2
    zmin = -zmax

for file in range(outnum_nd + 1 - first_file):
    print("****************************")
    print(
        "#",
        "hd2D" + STR + str("%03d" % (file + first_file)) + ".out",
        pylab.size(data[file - 1, :, :]),
    )
    fig = plt.figure()
    ax = fig.add_subplot()
    data2 = data[file, :, :]
    if homogeneous:
        myplot = ax.imshow(data2, cmap=color, vmin=zmin, vmax=zmax, origin="lower")
    else:
        myplot = ax.imshow(data2, cmap=color, origin="lower")
    # change scale axis 0 -> -pi, 2024 -> pi
    ax.set_xticks([0, 512, 1024, 1536, 2048])
    ax.set_xticklabels(["$-\pi$", "$-\pi/2$", "$0$", "$\pi/2$", "$\pi$"])
    ax.set_yticks([0, 512, 1024, 1536, 2048])
    ax.set_yticklabels(["$-\pi$", "$-\pi/2$", "$0$", "$\pi/2$", "$\pi$"])

    # increase font size
    plt.xticks(fontsize=FONTSIZE)
    plt.yticks(fontsize=FONTSIZE)

    ax.set_xlabel(r"$x$", fontsize=FONTSIZE)
    ax.set_ylabel(r"$y$", fontsize=FONTSIZE)

    if paint_arrows:
        # add arrows with 2 heads
        ax.arrow(
            reso / 2 - reso / 16,
            reso / 2 - reso / 12,
            reso / 8,
            0,
            head_width=30,
            head_length=40,
            fc="k",
            ec="k",
        )

        ax.arrow(
            reso / 2 + reso / 16,
            reso / 2 - reso / 12,
            -reso / 8,
            0,
            head_width=30,
            head_length=40,
            fc="k",
            ec="k",
        )

        # add label under the arrows
        ax.text(
            reso / 2,
            reso / 2 - reso / 6,
            r"$2\pi/k_r$",
            fontsize=FONTSIZE,
            ha="center",
        )

    # axins = zoomed_inset_axes(ax, 10, loc=1)  # zoom = 6

    # # sub region of the original image
    # ll = 25
    # x1, x2, y1, y2 = -ll, ll, -ll, ll
    # x1 += 1024
    # x2 += 1024
    # y1 += 1024
    # y2 += 1024

    # axins.imshow(
    #     data2,
    #     cmap=color,
    #     vmin=zmin,
    #     vmax=zmax,
    #     origin="lower",
    # )
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)

    #

    # # now set lim

    # # hide ticks
    # plt.xticks(visible=False)
    # plt.yticks(visible=False)

    # # draw a bbox of the region of the inset axes in the parent axes and
    # # connecting lines between the bbox and the inset axes area
    # mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    # add colorbar
    # cbar = plt.colorbar(myplot)
    # tight margins
    if paint_arrows:
        filename = os.path.join(
            script_dir,
            "../../"
            + output_dir
            + "FlowD_"
            + STR
            + str("%03d" % (file + first_file))
            + "_arrow.pdf",
        )
    else:
        filename = os.path.join(
            script_dir,
            "../../"
            + output_dir
            + "FlowD_"
            + STR
            + str("%03d" % (file + first_file))
            + ".pdf",
        )
    plt.savefig(filename, bbox_inches="tight")
    plt.close()
print("****************************")
