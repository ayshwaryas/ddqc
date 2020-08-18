import os

import pandas as pd
import pegasusio as io

from config import DATA_DIR


def read_tm(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_tm", task_id, tasks_per_tiss)
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
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    annotations = pd.read_csv(data_path + "annotations_droplet.csv")
    annotations_cell_type = annotations["cell_ontology_class"]
    annotations_cell_type.index = [tissue + "-" + "-".join(t.rsplit("_", 1)) for t in annotations["cell"]]
    annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = annotations_cell_type

    return tissue, is_human, adata


def read_other_10x(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_other_10X", task_id, tasks_per_tiss)
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
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv
    return tissue, is_human, adata


def read_hca(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_hca", task_id, tasks_per_tiss)
    is_human = True  # this is mouse data
    data_path = DATA_DIR + "HCA/"  # path to the mtx files of this dataset
    filename = "read_info_{}.csv".format(task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if tissue in directory:  # if directory matches the tissue
            p = data_path + directory  # path to the txt
            read_info.write("{},{},{},\n".format(directory.split("_")[1], p, "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    def reformat_ann_names(name):
        name = "-".join(name.split(".", 1)).replace("_", "")
        if name[-1].isdigit():
            name = name[:-1]
        return name

    def reformat_cell_names(name):
        name = name.rsplit("-", 1)
        if name[0][-2] == "-":
            name[0] = name[0][:-2]
        return name[0].replace("-", "") + "-" + name[1]

    annotations = pd.read_csv(data_path + "annotations.csv")
    annotations["cellnames"] = [reformat_ann_names(t) for t in annotations["cellnames"]]
    annotations = annotations.drop_duplicates("cellnames", keep='first')
    annotations.index = annotations["cellnames"]
    annotations_cell_type = annotations.reindex([reformat_cell_names(t) for t in adata.obs.index])["celltype"]
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = list(annotations_cell_type)

    return tissue, is_human, adata


def get_tissue_by_task_id(dataset, task_id, tasks_per_tiss):
    if dataset == "mc_tm" or dataset == "tm":
        tissues = ("Bladder", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Marrow",
              "Spleen", "Thymus", "Tongue", "Trachea")
    elif dataset == "mc_other_10X" or dataset == "other_10X":
        tissues = ("adipose", "ASD_snRNAseq", "liver", "skin")
    elif dataset == "mc_hca" or dataset == "hca":
        tissues = ( "Adipose", "Adrenal-Gland", "Artery", "Ascending-Colon", "Bladder", "Bone-Marrow", "Cerebellum",
                    "Cervix", "Duodenum", "Epityphlon", "Esophagus", "Fallopian-Tube", "Gall-Bladder", "Heart", "Ileum",
                    "Jejunum", "Kidney", "Liver", "Lung", "Muscle", "Omentum", "Pancreas", "Peripheral-Blood", "Pleura",
                    "Prostate", "Rectum", "Sigmoid-Colon", "Spleen", "Stomach", "Temporal-Lobe", "Thyroid", "Trachea",
                    "TransverseColon", "Ureter", "Uterus")
    else:
        return None
    return tissues[task_id // tasks_per_tiss]


def auto_reader(dataset, task_id, tasks_per_tiss):  # find the reading function based on project name
    if dataset == "mc_tm" or dataset == "tm":
        return read_tm(task_id, tasks_per_tiss)
    if dataset == "mc_other_10X" or dataset == "other_10X":
        return read_other_10x(task_id, tasks_per_tiss)
    if dataset == "mc_hca" or dataset == "hca":
        return read_hca(task_id, tasks_per_tiss)
