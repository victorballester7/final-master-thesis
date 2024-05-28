import matplotlib.pyplot as plt
import numpy as np
import sys
import time
from read_data import get_num_frames, set_data, get_misc
import os

# count time
UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()

duration = 10  # in seconds

# get script path
script_dir = os.path.dirname(__file__)
script_path = os.path.abspath(script_dir)

folder_path = script_path + "/../../data/pointvortices"

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

folder_path += program + "/positions"
L, output_pos = get_misc(folder_path + "/../misc.txt")
L = L + 1
L = 10
xmin, xmax = -L, L
ymin, ymax = -L, L


# Read data
num_frames = get_num_frames(folder_path)

# set data
data_blocks_x, data_blocks_y, circulation = set_data(folder_path)

num_vortices = [len(data_blocks_x[i]) for i in range(num_frames)]

FPS = int(num_frames / duration) + 1  # +1 to ensure positivity
interval = 1000 / FPS  # interval between frames in milliseconds

# Print some information
# print("Length of data: ", num_frames)
# print("Number of frames: ", num_frames)
# print("Interval: ", interval)
# print("Expected duration: ", duration)

# Axes

for frame in range(num_frames):
    # Create a figure
    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")

    # Update data
    x = np.array(data_blocks_x[frame])
    y = np.array(data_blocks_y[frame])
    c = np.array(circulation[frame])
    colors = ["red" if c[i] else "blue" for i in range(num_vortices[frame])]

    colors = np.array(colors)

    for i in range(num_vortices[frame]):
        ax.plot(x[i], y[i], "o", markersize=2, color=colors[i])

    # add time text
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    time_text.set_text(
        "Time = %d %.3d" % (frame * output_pos // 1000, frame * output_pos % 1000)
    )

    # Save the animation
    filename = (
        script_dir
        + "/../../images/pointvortices"
        + program
        + "/pointvortices."
        + str("%09d" % frame)
        + ".jpg"
    )

    plt.savefig(filename)
    plt.close()

    if frame % 10 == 0:
        print("Frame", frame, "/", num_frames)

end_time = time.time()
