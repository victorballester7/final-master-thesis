import numpy as np
import os

num_zeros_str = 4


def get_num_files(folder_path):
    num_files = len(
        [
            f
            for f in os.listdir(folder_path)
            if os.path.isfile(os.path.join(folder_path, f))
        ]
    )
    return num_files


def read_data(input_dir, file_name, times_str, ext):
    # let node be 000, 001, 002, ..., 999
    data = []
    N = get_num_files(input_dir)
    file_name = input_dir + file_name
    for n in range(N):
        # print('core', n)
        with open(file_name + times_str[n] + ext, "r") as f:
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


def read_spectra_times(file_name):
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


def getAveragedRadius(input_dir, output_file, file_name, average=False):
    ext = ".txt"
    if average:
        times_str, times = read_spectra_times(input_dir + "../../spectra_times.txt")
    else:
        times_str, times = read_spectra_times(input_dir + "../spectra_times.txt")
    data = read_data(input_dir, file_name, times_str, ext)

    # Do it only for the Energy (second column)
    data = data[:, :, 1]

    radius = np.arange(0, data.shape[1]) + 1

    meanRadius = np.zeros((data.shape[0], 2))
    meanRadius[:, 0] = times[: len(meanRadius)]

    # multiply each row by the radius (do it vectorized)
    f = radius**2
    meanRadius[:, 1] = np.sum(data * f, axis=1) / np.sum(data, axis=1)

    with open(output_file, "w") as f:
        for row in meanRadius:
            for elem in row:
                f.write(str(elem) + " ")
            f.write("\n")


def average_data(root):
    # check if there is a folder called average in data
    extra = ""
    average = False
    if os.path.exists(root + "/data/average"):
        extra = "average/"
        average = True
    # EnergyProf

    input_dir = "/data/" + extra + "EnergyProf/"
    file_name = "Energy."
    output_file = "/data/" + extra + "EnergyMeanRadius.txt"

    input_dir = root + input_dir
    output_file = root + output_file

    print("Merging data from EnergyProf...")
    getAveragedRadius(input_dir, output_file, file_name, average)

    # EnstrophyProf

    input_dir = "/data/" + extra + "EnstrophyProf/"
    file_name = "Enstrophy."
    output_file = "/data/" + extra + "EnstrophyMeanRadius.txt"

    input_dir = root + input_dir
    output_file = root + output_file

    print("Merging data from EnstrophyProf...")
    getAveragedRadius(input_dir, output_file, file_name, average)


if __name__ == "__main__":
    # check if there is a folder called average in data

    # get script directory
    script_dir = os.path.dirname(__file__)
    script_path = os.path.abspath(script_dir)

    average_data(script_path)
