import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import sys
import time
from read_data import *

duration = 10  # in seconds


def animate(folder_path: str, save: bool = False):
    # read L from misc.txt

    L, output_pos = get_misc(folder_path + "/../misc.txt")
    L = L + 1
    xmin, xmax = -L, L
    ymin, ymax = -L, L

    # Create a figure
    fig, ax = plt.subplots()

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

    # Create data
    vortices = []
    for i in range(num_vortices[0]):
        (aux,) = ax.plot(
            [], [], "o", markersize=2, color="red" if circulation[i] == True else "blue"
        )
        vortices.append(aux)

    # add time text
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    # Axes
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")

    def update(frame):
        # Update data
        x = data_blocks_x[frame]
        y = data_blocks_y[frame]

        for j in range(num_vortices):
            vortices[j].set_data([x[j]], [y[j]])

        time_text.set_text(
            "Time = %d %.3d" % (frame * output_pos // 1000, frame * output_pos % 1000)
        )

        return [vortices, time_text]

    # Create the animation
    ani = FuncAnimation(fig, update, frames=num_frames, interval=interval, blit=False)

    end_time = time.time()
    print(
        "Total time for animating: ",
        int(UNIT_TIME * (end_time - start_time)),
        LABEL_TIME,
    )

    if save:
        # Save the animation
        filename = (
            "videos/pointvortices/pointvortices." + "N=" + str(num_vortices) + ".mp4"
        )
        ani.save(filename, writer="ffmpeg", dpi=300, fps=FPS)
    else:
        # Show the animation
        plt.show()


# count time
UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
folder_results = "data/pointvortices/positions"

# Read data

# if there is an argument, then save the animation
if len(sys.argv) > 1:
    save = True
    animate(folder_results, save)
else:
    animate(folder_results)
