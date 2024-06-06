import numpy as np
import os


def read_data_pos(file_path):
    """Read data from file. The file is of the form:
    All are numbers expect the penultimate column which is a string
    x_00 x_01 ... x_0m
    x_10 x_11 ... x_1m
    ...
    x_n0 x_n1 ... x_nm
    """
    # first column stores in R_rho
    # second column stores in Le
    # third column stores in state
    # fourth column stores in times
    with open(file_path, "r") as file:
        lines = file.readlines()

    x = []
    y = []
    circulation = []

    # forget the first line
    # lines = lines[1:]
    for line in lines:
        values = line.split()
        x.append(float(values[0]))
        y.append(float(values[1]))
        circulation.append(float(values[2]))

    x = np.array(x)
    x = np.array(x)
    circulation = np.array(circulation)
    circulation = np.where(circulation > 0, True, False)

    return x, y, circulation


def get_num_frames(folder_path):
    num_frames = len(
        [
            f
            for f in os.listdir(folder_path)
            if os.path.isfile(os.path.join(folder_path, f))
        ]
    )
    return num_frames


def set_data(folder_path):
    num_frames = get_num_frames(folder_path)

    times = []
    data_blocks_x = []
    data_blocks_y = []
    data_blocks_circulation = []
    for i in range(num_frames):
        file_path = f"{folder_path}/positions.{i:05d}.txt"
        x, y, circulation = read_data_pos(file_path)
        data_blocks_x.append(x)
        data_blocks_y.append(y)
        data_blocks_circulation.append(circulation)

    # data_blocks_x = np.array(data_blocks_x)
    # data_blocks_y = np.array(data_blocks_y)

    return data_blocks_x, data_blocks_y, data_blocks_circulation


def get_misc(dim_dir):
    f = open(dim_dir, "r")
    dim_str = f.readlines()
    f.close()
    # split the line and return the two values
    radius, output_pos = dim_str[0].split()
    return float(radius), int(output_pos)
