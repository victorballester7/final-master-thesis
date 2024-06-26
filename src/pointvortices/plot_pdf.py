import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from read_data import get_num_frames, set_data, get_misc
import os

matplotlib.use("Agg")  # for not having memory leak

plt.rc("text", usetex=True)

# frame = 440
# folder = "/R4"
frame = 550
folder = "/R32"
FONTSIZE = 16

# get script path
script_dir = os.path.dirname(__file__)
script_path = os.path.abspath(script_dir)

folder_path = script_path + "/../../data/pointvortices"

cmap = plt.get_cmap("RdBu")

# Get the colors for the minimum and maximum values
blue_max = cmap(30)[:3]  # Lowest value (blue)
red_max = cmap(220)[:3]  # Highest value (red)

print("blue_max: ", blue_max)
print("red_max: ", red_max)

# Read data
if len(sys.argv) < 2:
    print("Usage: python animate.py <type>")
    print("type: disk, dipoles")
    sys.exit(1)

mytype = sys.argv[1]

if mytype == "disk":
    program = "/disk"
elif mytype == "dipoles":
    program = "/dipoles"
else:
    print("Unknown type, exiting...")
    sys.exit(1)

folder_path += program + folder + "/positions"
L, output_pos = get_misc(folder_path + "/../misc.txt")

# L/32 is just a small margin added to saw better the plots
L = L + L / 32
# L = 10
xmin, xmax = -L, L
ymin, ymax = -L, L


# Read data
num_frames = get_num_frames(folder_path)

# set data
data_blocks_x, data_blocks_y, circulation = set_data(folder_path)

num_vortices = [len(data_blocks_x[i]) for i in range(num_frames)]


# Print some information
# print("Length of data: ", num_frames)
# print("Number of frames: ", num_frames)
# print("Interval: ", interval)
# print("Expected duration: ", duration)

# Axes

# Create a figure
fig, ax = plt.subplots()
ax.set_xlabel("$x$", fontsize=FONTSIZE)
ax.set_ylabel("$y$", fontsize=FONTSIZE)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# set ticks -L to L ---> -pi to pi
ax.set_xticks([-L, -L / 2, 0, L / 2, L])
ax.set_xticklabels(["$-\pi$", "$-\pi/2$", "$0$", "$\pi/2$", "$\pi$"], fontsize=FONTSIZE)
ax.set_yticks([-L, -L / 2, 0, L / 2, L])
ax.set_yticklabels(["$-\pi$", "$-\pi/2$", "$0$", "$\pi/2$", "$\pi$"], fontsize=FONTSIZE)

ax.set_aspect("equal")

# Update data
x = np.array(data_blocks_x[frame])
y = np.array(data_blocks_y[frame])
c = np.array(circulation[frame])
# colors = [red_max if c[i] else blue_max for i in range(num_vortices[frame])]

for i in range(num_vortices[frame]):
    if c[i]:
        ax.plot(x[i], y[i], "o", markersize=2, color=red_max)
    else:
        ax.plot(x[i], y[i], "o", markersize=2, color=blue_max)
    # ax.plot(x[i], y[i], "o", markersize=2, color=colors[i])


# Save the animation
filename = (
    script_dir
    + "/../../images/pointvortices"
    + program
    + "/pointvortices."
    + folder[1:]
    + "."
    + str("%05d" % frame)
    + ".pdf"
)

# save the figure
fig.savefig(
    filename,
    bbox_inches="tight",
)
