import numpy as np
import os

N_CORES = 40


def get_num_files(folder_path):
    num_files = len([f for f in os.listdir(folder_path)
                     if os.path.isfile(os.path.join(folder_path, f))])
    return num_files


def read_data(file_name, ext):
    # let node be 000, 001, 002, ..., 999
    data = []
    # I want data to be an array first splited by node, then by rows and then
    # by columns (if necessary)
    for n in range(N_CORES):
        # print('core', n)
        node = str(n).zfill(3)
        with open(file_name + node + ext, 'r') as f:
            # split by columns and convert to float
            data_node = f.readlines()
            data_node = [row.split() for row in data_node]
            # print(data_node)
            data_node = np.array(data_node, dtype=float)
            data.append(data_node)
    data = np.array(data)

    return data


def average_data_onestage(input_file, output_file):
    ext = '.txt'
    data = read_data(input_file, ext)
    data_avg = np.mean(data, axis=0)

    with open(output_file, 'w') as f:
        for row in data_avg:
            for elem in row:
                f.write(str(elem) + ' ')
            f.write('\n')


def average_data_multistages(input_dir, input_file, output_file):
    # num of ext stages
    N = get_num_files(input_dir) // N_CORES
    for i in range(N):
        ext = '.' + str(i).zfill(3) + '.txt'
        data = read_data(input_file, ext)
        data_avg = np.mean(data, axis=0)

        with open(output_file + ext, 'w') as f:
            for row in data_avg:
                for elem in row:
                    f.write(str(elem) + ' ')
                f.write('\n')


def average_data():
    # energy_bal.node.txt

    input_file = './data/energy_bal.'
    output_file = './data/average/energy_bal.txt'

    average_data_onestage(input_file, output_file)

    # ----------------------------------------
    # enstrophy_bal.node.txt

    input_file = './data/enstrophy_bal.'
    output_file = './data/average/enstrophy_bal.txt'

    average_data_onestage(input_file, output_file)
    # ----------------------------------------
    # kspectrum.node.ext.txt

    input_dir = './data/kspectrum'
    input_file = input_dir + '/kspectrum.'
    output_file = './data/average/kspectrum/kspectrum'

    average_data_multistages(input_dir, input_file, output_file)

    # ----------------------------------------
    # EnergyProf.node.ext.txt

    input_dir = './data/EnergyProf'
    input_file = './data/EnergyProf/Energy.'
    output_file = './data/average/EnergyProf/Energy'

    average_data_multistages(input_dir, input_file, output_file)

    # ----------------------------------------
    # EnstrophyProf.node.ext.txt

    input_dir = './data/EnstrophyProf'
    input_file = './data/EnstrophyProf/Enstrophy.'
    output_file = './data/average/EnstrophyProf/Enstrophy'

    average_data_multistages(input_dir, input_file, output_file)


if __name__ == '__main__':
    average_data()
