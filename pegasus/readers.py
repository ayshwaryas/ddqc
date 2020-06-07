import os

import pegasus as pg
from config import DATA_DIR


def read_tm(task_id, tasks_per_tiss):
    tissue = ("Bladder", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Marrow",
              "Spleen", "Thymus", "Tongue", "Trachea")[task_id // tasks_per_tiss]
    is_human = False
    data_path = DATA_DIR + "tabula_muris/droplet/"
    filename = "read_info_{}.csv".format(task_id)
    read_info = open(filename, "w")
    read_info.write("Sample,Location,\n")
    for directory in os.listdir(data_path):
        if directory.startswith(tissue):
            read_info.write("{},{},\n".format(directory, data_path + directory + "/"))  # + "raw_gene_bc_matrices_h5.h5"))
    read_info.close()
    adata = pg.aggregate_matrices(filename)
    os.remove(filename)
    return tissue, is_human, adata


def auto_reader(dataset, task_id, tasks_per_tiss):
    if dataset == "mc_tm" or dataset == "tm":
        return read_tm(task_id, tasks_per_tiss)
