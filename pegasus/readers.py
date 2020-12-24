import os

import pandas as pd
import pegasusio as io
from projects_info import *

from config import DATA_DIR


def read_tm(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_tm", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "mouse/tabula_muris/" + tissue + "/"
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

    annotations = pd.read_csv(data_path + "../annotations_droplet.csv")
    annotations_cell_type = annotations["cell_ontology_class"]
    annotations_cell_type.index = [tissue + "-" + "-".join(t.rsplit("_", 1)) for t in annotations["cell"]]
    annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = annotations_cell_type

    return tissue, is_human, adata


def read_ts24(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_ts24", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "mouse/tabula_senis/24_month/" + tissue + "/"
    filename = "read_info_{}_{}.csv".format("ts24", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if tissue.upper() in directory:  # if directory matches the tissue
            p = data_path + directory + "/"  # path to the mtx
            read_info.write("{},{},{},\n".format(directory[:-3], p, "GRCm38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    annotations = pd.read_csv(data_path + "../annotations.csv")
    annotations_cell_type = pd.Series([t if t != "nan" else "Unknown" for t in annotations["cell_ontology_class_reannotated"]])
    annotations_cell_type.index = ["-".join(t.rsplit("_", 1)) for t in annotations["cell"]]
    annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = annotations_cell_type

    return tissue, is_human, adata


def read_ts30(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_ts30", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "mouse/tabula_senis/30_month/" + tissue + "/"
    filename = "read_info_{}_{}.csv".format("ts30", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if directory.startswith(tissue):  # if directory matches the tissue
            p = data_path + directory + "/"  # path to the mtx
            read_info.write("{},{},{},\n".format(directory[:-3], p, "GRCm38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    annotations = pd.read_csv(data_path + "../annotations.csv")
    annotations_cell_type = pd.Series([t if t != "nan" else "Unknown" for t in annotations["cell_ontology_class_reannotated"]])
    annotations_cell_type.index = [tissue + "-" + "-".join(t.rsplit("_", 1)) for t in annotations["cell"]]
    annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = annotations_cell_type

    return tissue, is_human, adata


def read_ebi_tm(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_ebi_tm", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "mouse/ebi_tm/" + tissue + "/" + tissue  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("ebi_tm", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    read_info.write("{},{},{},\n".format(tissue, data_path + ".h5ad", "GRCm38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    return tissue, is_human, adata


def read_ebi(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_ebi", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "human/ebi/" + tissue + "/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("ebi", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    read_info.write("{},{},{},\n".format(tissue, data_path + "matrix.mtx", "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    adata.obs["annotations"] = "Unknown"
    return tissue, is_human, adata


def read_mca(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_mca", task_id, tasks_per_tiss)
    is_human = False  # this is mouse data
    data_path = DATA_DIR + "mouse/mca/" + tissue + "/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("mca", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if directory.startswith(tissue):  # if directory matches the tissue
            p = data_path + directory  # path to the mtx
            read_info.write("{},{},{},\n".format(directory, p, "GRCm38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename, append_sample_name=False)  # read data
    os.remove(filename)  # remove the info csv

    annotations = pd.read_csv(data_path + "../annotations.csv")
    annotations_cell_type = annotations["Annotation"]
    annotations_cell_type.index = [tissue + "-" + t for t in annotations["Cell.name"]]
    annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = annotations_cell_type

    return tissue, is_human, adata


def read_other_10x(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_other_10X", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "human/other/" + tissue + "/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("other_10X", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        p = data_path + directory + "/"  # path to the mtx
        if os.path.isdir(p):
            read_info.write("{},{},{},\n".format(directory, p + "matrix.mtx", "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    adata.obs["annotations"] = "Unknown"
    return tissue, is_human, adata


def read_other(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_other", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "human/other/" + tissue + "/"  # path to the mtx files of this dataset
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
    data_path = DATA_DIR + "human/HCA/" + tissue + "/" # path to the mtx files of this dataset
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

    annotations = pd.read_csv(data_path + "../annotations.csv")
    annotations["cellnames"] = [reformat_ann_names(t) for t in annotations["cellnames"]]
    annotations = annotations.drop_duplicates("cellnames", keep='first')
    annotations.index = annotations["cellnames"]
    annotations_cell_type = annotations.reindex([reformat_cell_names(t) for t in adata.obs.index])["celltype"]
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = list(annotations_cell_type)

    return tissue, is_human, adata


def read_panglaodb(task_id, tasks_per_tiss):
    tissue = get_tissue_by_task_id("mc_PanglaoDB", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "human/PanglaoDB/" + tissue + "/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("PanglaoDB", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    read_info.write("{},{},{},\n".format(tissue, data_path + "matrix.mtx", "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    adata.obs["annotations"] = "Unknown"
    return tissue, is_human, adata


def read_human_blood(task_id, tasks_per_tiss):
    tissue = "Blood"  # get_tissue_by_task_id("mc_manton", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = "/broad/hptmp/subraman/"  # path to the mtx files of this dataset
    filename = "read_info_{}_{}.csv".format("blood", task_id)  # filename of csv used by aggregate_matrices
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


def read_human_brain(task_id, tasks_per_tiss):
    tissue = "Brain"  # get_tissue_by_task_id("mc_manton", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "human/other/brain/brain.h5"
    filename = "read_info_{}_{}.csv".format("brain", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    read_info.write("{},{},{},\n".format(tissue, data_path, "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    adata.obs["annotations"] = "Unknown"
    return tissue, is_human, adata


def read_human_brain_olfactory(task_id, tasks_per_tiss):
    tissue = "Brain_olfactory"  # get_tissue_by_task_id("mc_manton", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "human/other/brain-olfactory/"
    filename = "read_info_{}_{}.csv".format("brain_olfactory", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    for directory in os.listdir(data_path):  # each directory is one mtx file + genes and barcodes
        if directory.startswith("patient"):  # if directory matches the tissue
            p = data_path + directory + "/"  # path to the mtx
            read_info.write("{},{},{},\n".format(directory, p + "matrix.mtx", "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    adata.obs["annotations"] = "Unknown"
    return tissue, is_human, adata


def read_heart_circulation(task_id, tasks_per_tiss):
    tissue = "Heart_Circulation"  # get_tissue_by_task_id("mc_manton", task_id, tasks_per_tiss)
    is_human = True  # this is human data
    data_path = DATA_DIR + "human/other/heart/heart_circulation.h5ad"
    filename = "read_info_{}_{}.csv".format("heart_circulation", task_id)  # filename of csv used by aggregate_matrices
    read_info = open(filename, "w")  # csv for aggregate_matrices
    read_info.write("Sample,Location,Reference,\n")
    read_info.write("{},{},{},\n".format(tissue, data_path, "GRCh38"))  # add the file info to csv
    read_info.close()
    adata = io.aggregate_matrices(filename)  # read data
    os.remove(filename)  # remove the info csv

    annotations = pd.read_csv(DATA_DIR + "other/heart_circulation_annotations.csv")
    annotations["annotation"] = [t[4:] for t in annotations["Cluster"]]
    annotations.index = annotations["index"]
    new_cell_names = [t[18:].replace("::", "_").replace("_Reseq", "").replace("_Rerun", "").replace("LV_01723_1", "LV_1723_1").replace("LV_P01681_1", "LV_1681_1") for t in adata.obs.index]
    annotations_cell_type = annotations.reindex(new_cell_names)["annotation"]
    annotations_cell_type = annotations_cell_type.fillna("Unknown")
    adata.obs["annotations"] = list(annotations_cell_type)

    return tissue, is_human, adata


def auto_reader(dataset, task_id, tasks_per_tiss):  # find the reading function based on project name
    if dataset == "mc_tm" or dataset == "tm":
        return read_tm(task_id, tasks_per_tiss)
    if dataset == "mc_ts24" or dataset == "ts24":
        return read_ts24(task_id, tasks_per_tiss)
    if dataset == "mc_ts30" or dataset == "ts30":
        return read_ts30(task_id, tasks_per_tiss)
    if dataset == "mc_ebi_tm" or dataset == "ebi_tm":
        return read_ebi_tm(task_id, tasks_per_tiss)
    if dataset == "mc_ebi" or dataset == "ebi":
        return read_ebi(task_id, tasks_per_tiss)
    if dataset == "mc_mca" or dataset == "mca":
        return read_mca(task_id, tasks_per_tiss)
    if dataset == "mc_other" or dataset == "other":
        return read_other(task_id, tasks_per_tiss)
    if dataset == "mc_other_10X" or dataset == "other_10X":
        return read_other_10x(task_id, tasks_per_tiss)
    if dataset == "mc_hca" or dataset == "hca":
        return read_hca(task_id, tasks_per_tiss)
    if dataset == "mc_PanglaoDB" or dataset == "PanglaoDB":
        return read_panglaodb(task_id, tasks_per_tiss)
    if dataset == "mc_blood" or dataset == "blood":
        return read_human_blood(task_id, tasks_per_tiss)
    if dataset == "mc_brain" or dataset == "brain":
        return read_human_brain(task_id, tasks_per_tiss)
    if dataset == "mc_brain_olfactory" or dataset == "brain_olfactory":
        return read_human_brain_olfactory(task_id, tasks_per_tiss)
    if dataset == "mc_heart_circulation" or dataset == "heart_circulation":
        return read_heart_circulation(task_id, tasks_per_tiss)

