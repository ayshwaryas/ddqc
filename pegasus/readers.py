import os

import pegasus as pg

from config import DATA_DIR


def read_tm(task_id, tasks_per_tiss):
    tissue = ("Bladder", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Marrow",
              "Spleen", "Thymus", "Tongue", "Trachea")[task_id // tasks_per_tiss]  # select tissue based on task id
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "tabula_muris/droplet/"  # path to the mtx files of this dataset
    filename = "read_info_{}.csv".format(task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if directory.startswith(tissue):  # if directory matches the tissue
            p = data_path + directory + "/"  # path to the mtx
            read_info.write("{},{},{},\n".format(directory, p + "matrix.mtx", "GRCm38"))  # add the file info to csv
    read_info.close()
    adata = pg.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv
    return tissue, is_human, adata


def read_other_10X(task_id, tasks_per_tiss):
    tissue = ("adipose", "ASD_snRNAseq", "liver", "skin")[task_id // tasks_per_tiss]  # select tissue based on task id
    is_human = True  # this is mouse data
    data_path = DATA_DIR + "other/" + tissue + "/"  # path to the mtx files of this dataset
    filename = "read_info_{}.csv".format(task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        p = data_path + directory + "/"  # path to the mtx
        if os.path.isdir(p):
            read_info.write("{},{},{},\n".format(directory, p + "matrix.mtx", "GRCh37"))  # add the file info to csv
    read_info.close()
    adata = pg.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv
    return tissue, is_human, adata


def auto_reader(dataset, task_id, tasks_per_tiss):  # find the reading function based on project name
    if dataset == "mc_tm" or dataset == "tm":
        return read_tm(task_id, tasks_per_tiss)
    if dataset == "mc_other_10X" or dataset == "other_10X":
        return read_other_10X(task_id, tasks_per_tiss)
