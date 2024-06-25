# FOR NONSTATIONARY VORTEX SIMULATIONS
import numpy as np
import os

first_file = 0
num_zeros_str = 5

program = "/disk"


def get_num_files(folder_path):
    num_files = len(
        [
            f
            for f in os.listdir(folder_path)
            if os.path.isfile(os.path.join(folder_path, f))
        ]
    )
    return num_files


def read_data(input_dir, file_name, ext):
    # let node be 000, 001, 002, ..., 999
    data = []
    N = get_num_files(input_dir)
    file_name = input_dir + file_name
    for n in range(first_file, N):
        # print('core', n)
        node = str(n).zfill(num_zeros_str)
        with open(file_name + node + ext, "r") as f:
            # split by columns and convert to float
            data_node = f.readlines()
            data_node = [row.split() for row in data_node]
            # print(data_node)
            data_node = np.array(data_node, dtype=float)
            data.append(data_node)
    # each data_node may have different number of rows
    # so we take the minimum number of rows from all data_node and truncate
    min_rows = min([data_node.shape[0] for data_node in data])
    data = [data_node[:min_rows] for data_node in data]

    data = np.array(data)

    return data


def read_energy_times(file_name):
    # read values of the form:
    # str0  time0
    # str1  time1
    # ...
    # strN  timeN
    with open(file_name, "r") as f:
        data = f.readlines()
        data = [row.split() for row in data]
        data = np.array(data)
        # check if there are rows whose first column is the same
        # if so, remove the first first occurrence of the whole row
        # this is to avoid duplicates
        _, index = np.unique(data[:, 0], return_index=True)
        data = data[index]
        data_times_str = data[:, 0]
        data_times = data[:, 1].astype(float)
    return data_times_str, data_times


def getAveragedRadius(input_dir, output_file, file_name):
    ext = ".txt"
    _, times = read_energy_times(input_dir + "../energy_times.txt")
    data = read_data(input_dir, file_name, ext)

    # Do it only for the Energy (second column)
    data = data[:, :, 1]

    # linspace from 0 to 40 including end and expluding 0 with dr = 0.08
    radius = np.linspace(0.08, 40, 500)

    meanRadius = np.zeros((data.shape[0], 2))
    meanRadius[:, 0] = times[: len(meanRadius)]

    # multiply each row by the radius (do it vectorized)
    f = radius**2
    meanRadius[:, 1] = np.sqrt(np.sum(data * f, axis=1) / np.sum(data, axis=1))

    with open(output_file, "w") as f:
        for row in meanRadius:
            for elem in row:
                f.write(str(elem) + " ")
            f.write("\n")


def average_data(root):
    extra = "/NoRemovalVortices"

    # NumVortices

    input_dir = "data/pointvortices" + program + extra + "/NumVortices/"
    file_name = "NumVortices."
    output_file = "data/pointvortices" + program + extra + "/NumVorticesMeanRadius.txt"

    input_dir = root + input_dir
    output_file = root + output_file

    print("Merging data from NumVortices...")
    getAveragedRadius(input_dir, output_file, file_name)


if __name__ == "__main__":
    # check if there is a folder called average in data

    # get script directory
    script_dir = os.path.dirname(__file__)
    script_path = os.path.abspath(script_dir)

    root = script_path + "/../../"

    average_data(root)
