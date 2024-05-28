import numpy as np
import os
import fnmatch

N_CORES = 48


# def get_num_files(folder_path):
#     num_files = len(
#         [
#             f
#             for f in os.listdir(folder_path)
#             if os.path.isfile(os.path.join(folder_path, f))
#         ]
#     )
#     return num_files


def get_num_files(folder_path, node):
    pattern = f"*.{node}.*.txt"
    num_files = len(
        [
            f
            for f in os.listdir(folder_path)
            if os.path.isfile(os.path.join(folder_path, f))
            and fnmatch.fnmatch(f, pattern)
        ]
    )
    return num_files


def read_data(file_name, ext):
    # let node be 000, 001, 002, ..., 999
    data = []
    # I want data to be an array first splited by node, then by rows and then
    # by columns (if necessary)
    for n in range(N_CORES):
        # print('core', n)
        node = str(n).zfill(3)
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


def average_data_onestage(input_file, output_file):
    ext = ".txt"
    data = read_data(input_file, ext)
    data_avg = np.mean(data, axis=0)

    with open(output_file, "w") as f:
        for row in data_avg:
            for elem in row:
                f.write(str(elem) + " ")
            f.write("\n")


def get_num_files_multistages(folder_path):
    # get the minimum number of slices which is the least among all the cores.
    # So I want to count the files in a directory of the form *.node.*.txt, for node = 000, 001, 002, ..., N_CORES and take the minimum
    num_files = []
    for i in range(N_CORES):
        node = str(i).zfill(3)
        num_files.append(get_num_files(folder_path, node))
    return min(num_files)


def average_data_multistages(input_dir, input_file, output_file, multicols):
    N = get_num_files_multistages(input_dir)
    print("Number of files: ", N)
    for i in range(N):
        if i % 10 == 0:
            print("file ", i, " / ", N)
        ext = "." + str(i).zfill(3) + ".txt"
        data = read_data(input_file, ext)

        # do the Lp norm for each row, where p is the index of the column
        # ((1/N_CORES) * sum_cores data^p)^(1/p)
        if multicols:
            data_avg = np.zeros(data[0].shape)
            for col in range(data_avg.shape[1]):
                p = col + 1
                data_avg[:, col] = (np.sum(data[:, :, col] ** p, axis=0) / N_CORES) ** (
                    1.00 / p
                )
        else:
            data_avg = np.mean(data, axis=0)
        with open(output_file + ext, "w") as f:
            for row in data_avg:
                for elem in row:
                    f.write(str(elem) + " ")
                f.write("\n")


def average_data():
    # energy_bal.node.txt

    input_file = "./data/energy_bal."
    output_file = "./data/average/energy_bal.txt"

    print("Merging data from energy_bal...")
    average_data_onestage(input_file, output_file)

    # ----------------------------------------
    # enstrophy_bal.node.txt

    input_file = "./data/enstrophy_bal."
    output_file = "./data/average/enstrophy_bal.txt"

    print("Merging data from enstrophy_bal...")
    average_data_onestage(input_file, output_file)
    # ----------------------------------------
    # kspectrum.node.ext.txt

    input_dir = "./data/kspectrum"
    input_file = input_dir + "/kspectrum."
    output_file = "./data/average/kspectrum/kspectrum"

    print("Merging data from kspectrum...")
    average_data_multistages(input_dir, input_file, output_file, False)

    # ----------------------------------------
    # EnergyProf.node.ext.txt

    input_dir = "./data/EnergyProf"
    input_file = "./data/EnergyProf/Energy."
    output_file = "./data/average/EnergyProf/Energy"

    print("Merging data from EnergyProf...")
    average_data_multistages(input_dir, input_file, output_file, True)

    # ----------------------------------------
    # EnstrophyProf.node.ext.txt

    input_dir = "./data/EnstrophyProf"
    input_file = "./data/EnstrophyProf/Enstrophy."
    output_file = "./data/average/EnstrophyProf/Enstrophy"

    print("Merging data from EnstrophyProf...")
    average_data_multistages(input_dir, input_file, output_file, True)


if __name__ == "__main__":
    # check if there is a folder called average in data
    if not os.path.exists("./data/average"):
        os.makedirs("./data/average")
    if not os.path.exists("./data/average/kspectrum"):
        os.makedirs("./data/average/kspectrum")
    if not os.path.exists("./data/average/EnergyProf"):
        os.makedirs("./data/average/EnergyProf")
    if not os.path.exists("./data/average/EnstrophyProf"):
        os.makedirs("./data/average/EnstrophyProf")

    average_data()
