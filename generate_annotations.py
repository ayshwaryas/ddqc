import shutil

import pandas as pd

from config.config import DATA_DIR, RIBO_PREFIXES, MITO_PREFIXES
from filtering import initial_qc
from projects_info import *
from readers_old import auto_reader
from utils import safe_mkdir

path = "config/read_info/"
renaming_dict = {
    "mc_tm": "tabula_muris",
    "mc_ts24": "tabula_senis_24m",
    "mc_ts30": "tabula_senis_30m",
    "mc_ebi_tm": "tabula_muris_smartseq2",
    "mc_ebi": "ebi_sc_experiment",
    "mc_mca": "mouse_cell_atlas",
    "mc_other_10X": "human_other",
    "mc_other": "human_other",
    "mc_hca": "human_tissue_atlas",
    "mc_PanglaoDB": "PanglaoDB",
    "mc_brain": "human_other",
    "mc_blood": "human_other",
    "mc_olfactory_epithelium": "human_other",
    "mc_heart_circulation": "human_other",
    "mc_krasnow_lung": "human_other",
}

with open(path + "projects.csv", "w") as fout:
    fout.write("project,tissue,is_human,annotations\n")

for project in N_TISSUES_PROJ.keys():
    new_project = renaming_dict[project]
    safe_mkdir(path + new_project)
    for task_id in range(N_TISSUES_PROJ[project]):
        tissue, is_human, adata, path_to_annotations = auto_reader(project, task_id, 1)
        if path_to_annotations != "N/A":
            path_to_annotations = path_to_annotations.split(DATA_DIR)[1] + "annotations_pg.csv"
        else:
            path_to_annotations = "Unknown"

        with open(path + "projects.csv", "a") as fout:
            fout.write(",".join([new_project, tissue, str(is_human), path_to_annotations]) + "\n")

        dataset_list = pd.read_csv(tissue + ".csv")
        dataset_list = dataset_list[["Sample", "Location", "Reference"]]
        dataset_list['Location'] = [t.split(DATA_DIR)[1] for t in
                                    dataset_list['Location']]  # update location with relevant directory prefix
        with open(tissue + ".csv", "w") as fout:  # write modified csv
            fout.write(dataset_list.to_csv(index=False))
        shutil.move(tissue + ".csv", path + new_project + "/" + tissue + ".csv")

        if adata is not None:
            mito_prefix = MITO_PREFIXES["human"] if is_human else MITO_PREFIXES["mouse"]
            ribo_prefix = RIBO_PREFIXES["human"] if is_human else RIBO_PREFIXES["mouse"]
            adata = initial_qc(adata, 10, 100, mito_prefix, ribo_prefix)
            adata = adata.obs["annotations"]
            with open(DATA_DIR + path_to_annotations, "w") as fout:
                fout.write(adata.to_csv())
