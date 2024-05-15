import numpy as np
import matplotlib.pyplot as plt
import os
import pylab
from matplotlib.animation import FuncAnimation


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


input_dir = "data/simple_2DNS/output/"
output_dir = "videos/simple_2DNS/"
STR = "ww."
outnum_nd = get_num_files(input_dir, STR)   # we have ww, ps, fw and fp
homogeneous = True  # True if we want the same colors for all images
save = True # True if we want to save the animationm, False if we want to show it

dim_dir = "data/simple_2DNS/dim.txt"
reso = get_dim(dim_dir)
nx = reso
ny = nx
color = "RdBu_r"
print("reso=", reso)
print("outnum_nd=", outnum_nd)

script_dir = os.path.dirname(__file__)

data = np.zeros((outnum_nd, reso, reso))

for file in range(1, outnum_nd+1):
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
    data[file-1, :, :] = data2

    print("file=", file)

fig, ax = plt.subplots()

if homogeneous:
    zmax = np.max(data)
    zmin = np.min(data)

    zmax /= 2
    zmin /= 2
    plot = ax.imshow(data[0], cmap=color, vmin=zmin, vmax=zmax, origin='lower')
else:
    plot = ax.imshow(data[0], cmap=color, origin='lower')


cbar = plt.colorbar(plot)

ax.set_aspect('equal')

def update(frame):

    # Update data

    plot.set_data(data[frame])

    return plot

# Create the animation
ani = FuncAnimation(fig, update, frames=outnum_nd,
                    blit=False)

if save:
    # Save the animation
    filename = output_dir + 'simple_2DNS.mp4'
    ani.save(filename, writer='ffmpeg', dpi=300)
else:
    # Show the animation
    plt.show()
