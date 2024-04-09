import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import sys
import time
from read_data import *

duration = 10  # in seconds
xmin, xmax = -5, 5
ymin, ymax = -5, 5


def animate(folder_path: str, save: bool = False):
    # Create a figure
    fig, ax = plt.subplots()

    # Read data
    num_frames = get_num_frames(folder_path)

    # set data
    data_blocks_x, data_blocks_y, circulation = set_data(folder_path)

    num_vortices = len(data_blocks_x[0])

    FPS = int(num_frames / duration) + 1  # +1 to ensure positivity
    interval = 1000 / FPS  # interval between frames in milliseconds

    # Print some information
    print("Length of data: ", num_frames)
    print("Number of frames: ", num_frames)
    print("Interval: ", interval)
    print("Expected duration: ", duration)

    # Create data
    # vortices = [ax.plot([], [], 'o', lw=2, color='red' if circulation[i] == True else 'blue')
    # for i in range(num_vortices)]
    print(circulation)
    vortices = []
    for i in range(num_vortices):
        aux, = ax.plot(
            [], [], 'o', markersize=2, color='red' if circulation[i] == True else 'blue')
        vortices.append(aux)

    # Axes
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')

    def update(frame):

        # Update data
        x = data_blocks_x[frame]
        y = data_blocks_y[frame]

        for j in range(num_vortices):
            vortices[j].set_data([x[j]], [y[j]])

        return vortices

    # Create the animation
    ani = FuncAnimation(fig, update, frames=num_frames,
                        interval=interval, blit=False)

    end_time = time.time()
    print("Total time for animating: ", int(
        UNIT_TIME * (end_time - start_time)), LABEL_TIME)

    # if save:
    #     # Save the animation
    #     filename = 'data/videos/' + plot_var + '_Pr=' + str(Pr) + '_Le=' + str(
    #         Le) + '_Ra_T=' + str(Ra_T) + '_R_rho=' + str(R_rho) + '.mp4'
    #     ani.save(filename, writer='ffmpeg', dpi=300, fps=FPS)
    # else:
    # Show the animation
    plt.show()


# count time
UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
folder_results = 'data/pointvortices'

# Read data

# if there is an argument, then save the animation
if len(sys.argv) > 1:
    save = True
    animate(folder_results, save)
else:
    animate(folder_results)
