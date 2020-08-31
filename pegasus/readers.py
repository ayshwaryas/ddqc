import os

import pandas as pd
import pegasusio as io

from config import DATA_DIR


def read_tm(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_tm", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "tabula_muris/droplet/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("tm", task_id)  # filename of csv used by aggregate_matrices
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


def read_ebi_tm(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_ebi_tm", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "ebi_tm/" + tissue + "/" # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("ebi_tm", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    read_info.write("{},{},{},\n".format(tissue, data_path + "matrix.mtx", "GRCm38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    annotations = pd.read_csv(data_path + "annotations.csv")
    annotations_cell_type = annotations["celltype"]
    annotations_cell_type.index = [tissue + "-" + t for t in annotations["Unnamed: 0"]]
    annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = annotations_cell_type

    return tissue, is_human, adata


def read_mca(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_mca", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "mca/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("mca", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if directory.startswith(tissue):  # if directory matches the tissue
            p = data_path + directory  # path to the mtx
            read_info.write("{},{},{},\n".format(directory, p, "GRCm38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    return tissue, is_human, adata


def read_other_10x(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_other_10X", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "other/" + tissue + "/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("other_10X", task_id)  # filename of csv used by aggregate_matrices
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


def read_other(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_other", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "other/" + tissue + "/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("other", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    read_info.write("{},{},{},\n".format(tissue, data_path + "matrix.mtx", "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    annotations = pd.read_csv(data_path + "annotations.csv")
    annotations_cell_type = annotations["celltype"]
    annotations_cell_type.index = [tissue + "-" + t for t in annotations["Unnamed: 0"]]
    annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = annotations_cell_type

    return tissue, is_human, adata


def read_hca(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_hca", task_id, tasks_per_tiss)
    is_human = True  # this is mouse data
    data_path = DATA_DIR + "HCA/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("hca", task_id)  # filename of csv used by aggregate_matrices
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


def read_manton(task_id, tasks_per_tiss):
    tissue = "Manton"  # get_tissue_by_task_id("mc_manton", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = "/broad/hptmp/subraman/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("manton", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if directory.startswith(tissue):  # if directory matches the tissue
            p = data_path + directory + "/"  # path to the mtx
            read_info.write("{},{},{},\n".format(directory, p + "raw_gene_bc_matrices_h5.h5", "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    adata.obs["annotations"] = "Unknown"
    return tissue, is_human, adata


def get_tissue_by_task_id(dataset, task_id, tasks_per_tiss):
    if dataset == "mc_tm" or dataset == "tm":
        tissues = ("Bladder", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Marrow",
              "Spleen", "Thymus", "Tongue", "Trachea")
    elif dataset == "mc_ebi_tm" or dataset == "ebi_tm":
        tissues = ("Adipose", "Bladder", "Bone_Marrow", "Cerebellum", "Cerebral_Cortex", "Colon", "Diaphragm", "Fat",
                   "Heart_and_Aorta", "Hippocampus", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland",
                   "Pancreas", "Skin", "Spleen", "Striatum", "Thymus", "Tongue", "Trachea")
    elif dataset == "mc_mca" or dataset == "mca":
        tissues = ("Bladder", "BoneMarrow", "Brain", "Kidney", "Liver", "Lung", "MammaryGland.Involution",
                   "MammaryGland.Lactation", "MammaryGland.Pregnancy", "MammaryGland.Virgin", "MesenchymalStemCells",
                   "Muscle", "Ovary", "Pancreas", "PeripheralBlood", "Placenta", "Prostate", "SmallIntestine", "Spleen",
                   "Stomach", "Testis", "Thymus", "TrophoblastStemCells", "Uterus")
    elif dataset == "mc_other_10X" or dataset == "other_10X":
        tissues = ("adipose", "ASD_snRNAseq", "liver", "skin")
    elif dataset == "mc_other" or dataset == "other":
        tissues = ("colon-epi_human", "colon-fib_human", "colon-imm_human", "kidney_human", "retina_human")
    elif dataset == "mc_hca" or dataset == "hca":
        tissues = ("Adipose", "Adrenal-Gland", "Artery", "Ascending-Colon", "Bladder", "Bone-Marrow", "Cerebellum",
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
    if dataset == "mc_ebi_tm" or dataset == "ebi_tm":
        return read_tm(task_id, tasks_per_tiss)
    if dataset == "mc_other" or dataset == "other":
        return read_other(task_id, tasks_per_tiss)
    if dataset == "mc_other_10X" or dataset == "other_10X":
        return read_other_10x(task_id, tasks_per_tiss)
    if dataset == "mc_hca" or dataset == "hca":
        return read_hca(task_id, tasks_per_tiss)
    if dataset == "mc_manton" or dataset == "manton":
        return read_hca(task_id, tasks_per_tiss)
