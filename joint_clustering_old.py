import subprocess
import sys

import pandas as pd
import pegasus as pg

from config.config import *
from filtering import initial_qc
from reading import get_project_info, read_tissue
from utils import safe_mkdir, cluster_data, add_cd_scores, marker_dict_to_df, assign_cell_types, save_to_csv


# function that creates all the relevant directories
def create_joint_dir(project, tissue, res):
    # directory with unfiltered clustering
    task_directory = "{}-joint_clustering_old".format(res)  # name of the directory for this task
    task_name = tissue + "-" + task_directory  # task name to put on plots
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + task_directory + "/"  # directory for saving output

    # create all the directories, if they dont exist
    print("Creating Joint Directories")
    safe_mkdir(results_dir)

    return task_directory, task_name, results_dir


# assign colors for joint clustering plot based on which methods retained which cells
def assign_colors(adata, project, tissue, res):
    mad_cells = pd.read_csv(OUTPUT_DIR + project + "/" + tissue + "/" + str(res) + "-mad-2/!cells.csv")
    cutoff_cells = pd.read_csv(OUTPUT_DIR + project + "/" + tissue + "/" + str(res) + "-cutoff-10/!cells.csv")
    adata.obs["color"] = "Neither"
    adata.obs.loc[mad_cells["barcodekey"], "color"] = "MAD2 only"
    adata.obs.loc[cutoff_cells["barcodekey"], "color"] = "Cutoff only"
    adata.obs.loc[list(set(cutoff_cells["barcodekey"]).intersection(set(mad_cells["barcodekey"]))), "color"] = "All"

    adata.obs["passed_qc"] = (adata.obs.color != "Neither")
    pg.filter_data(adata)
    return adata


def joint_main(project, task_id, tissue=None):
    if tissue is None:
        tissue, is_human, annotations = get_project_info(project, task_id=task_id)
    else:
        is_human, annotations = get_project_info(project, tissue=tissue)

    print("joint clustering old task.id:{} - tissue:{}, res:{}, project:{}".format(task_id, tissue, resolution, project))

    adata = read_tissue(project, tissue, annotations)

    task_directory, task_name, results_dir = create_joint_dir(project, tissue, resolution)
    assign_colors(adata, project, tissue, resolution)

    # filtering
    mito_prefix = MITO_PREFIXES["human"] if is_human else MITO_PREFIXES["mouse"]
    ribo_prefix = RIBO_PREFIXES["human"] if is_human else RIBO_PREFIXES["mouse"]
    adata = initial_qc(adata, basic_genes_filter, basic_mito_filter, mito_prefix,
                       ribo_prefix)  # perform initial qc

    adata, marker_dict = cluster_data(adata, compute_markers=True, compute_reductions=True, resolution=resolution)
    adata = add_cd_scores(adata, is_human)  # add cell death scores
    markers = marker_dict_to_df(marker_dict)
    clusters = assign_cell_types(adata, markers, tissue)

    print("Writing results")
    save_to_csv(adata, results_dir)
    with open(results_dir + "!clusters.csv", "w") as fout:
        fout.write(clusters.to_csv())
    with open(results_dir + "!markers.csv", "w") as fout:
        fout.write(markers.to_csv())

    print(
        subprocess.check_output("Rscript plotting.R {} {} {}".format(task_name, results_dir, "joint"),
                                shell=True).decode('UTF-8'))


if __name__ == '__main__':
    if local:  # for debug outside of cluster
        proj = input("Project: ").strip()
        t_id = int(input("Task ID (-1 to input specific tissue): ").strip())
        if t_id == -1:
            tiss = input("Tissue: ").strip()
            joint_main(proj, 0, tissue=tiss)
        else:
            joint_main(proj, t_id)
    else:  # project and task id are provided as commandline args
        proj = sys.argv[1]
        t_id = int(sys.argv[2]) - 1
        joint_main(proj, t_id)
