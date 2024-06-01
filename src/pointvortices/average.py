import numpy as np
import os

first_file = 100
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


def average_data_onestage(input_dir, output_file, file_name, rate=False):
    ext = ".txt"
    data = read_data(input_dir, file_name, ext)

    if rate:
        new_data = np.zeros((data.shape[0], data.shape[1], 2))
        # for each data file, divide the 2nd column by the 3rd column, while keeping the 1st column

        for i in range(data.shape[0]):
            dr = data[i, 1, 0] - data[i, 0, 0]
            data[i, :, 1] = data[i, :, 1] / data[i, :, 2] / (data[i, :, 0] * dr)
            # remove the 3rd column
            new_data[i] = data[i, :, [0, 1]].T
        data = new_data

    data_avg = np.mean(data, axis=0)

    with open(output_file, "w") as f:
        for row in data_avg:
            for elem in row:
                f.write(str(elem) + " ")
            f.write("\n")


def average_data(root):
    # EnergyFlux

    input_dir = "data/pointvortices" + program + "/EnergyFlux/"
    file_name = "Energy."
    output_file = "data/pointvortices" + program + "/average/EnergyFlux.txt"

    input_dir = root + input_dir
    output_file = root + output_file

    print("Merging data from EnergyFlux...")
    average_data_onestage(input_dir, output_file, file_name)

    # EnergyProf

    input_dir = "data/pointvortices" + program + "/EnergyProf/"
    file_name = "Energy."
    output_file = "data/pointvortices" + program + "/average/EnergyProf.txt"

    input_dir = root + input_dir
    output_file = root + output_file

    print("Merging data from EnergyProf...")
    average_data_onestage(input_dir, output_file, file_name)

    # NumVortices

    input_dir = "data/pointvortices" + program + "/NumVortices/"
    file_name = "NumVortices."
    output_file = "data/pointvortices" + program + "/average/NumVortices.txt"

    input_dir = root + input_dir
    output_file = root + output_file

    print("Merging data from NumVortices...")
    average_data_onestage(input_dir, output_file, file_name, True)


if __name__ == "__main__":
    # check if there is a folder called average in data

    # get script directory
    script_dir = os.path.dirname(__file__)
    script_path = os.path.abspath(script_dir)

    root = script_path + "/../../"

    if not os.path.exists(root + "data/pointvortices/disk/average/"):
        os.makedirs(root + "data/pointvortices/disk/average/")
    average_data(root)
