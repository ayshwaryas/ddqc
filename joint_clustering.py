import os
import shutil
import subprocess
import sys

import pandas as pd

from config.config import OUTPUT_DIR, local, resolution
from reading import get_project_info
from utils import safe_mkdir


# function that creates all the relevant directories
def prepare_joint_dir(project, tissue, res):
    # directory with unfiltered clustering
    source_dir = OUTPUT_DIR + project + "/" + tissue + "/" + str(res) + "-none-0/"
    task_directory = "{}-joint_clustering".format(res)  # name of the directory for this task
    task_name = tissue + "-" + task_directory  # task name to put on plots
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + task_directory + "/"  # directory for saving output

    # create all the directories, if they dont exist
    print("Creating Joint Directories")
    safe_mkdir(results_dir)

    # copy files
    for file in ("!cells.csv", "!clusters.csv", "!markers.csv"):
        assert os.path.isfile(source_dir + "!cells.csv")  # check if file exists
        shutil.copyfile(source_dir + file, results_dir + file)

    return task_directory, task_name, results_dir


# assign colors for joint clustering plot based on which methods retained which cells
def assign_colors(adata, project, tissue, res):
    mad_cells = pd.read_csv(OUTPUT_DIR + project + "/" + tissue + "/" + str(res) + "-mad-2/!cells.csv")
    cutoff_cells = pd.read_csv(OUTPUT_DIR + project + "/" + tissue + "/" + str(res) + "-cutoff-10/!cells.csv")
    adata["color"] = "Neither"
    adata["color"][mad_cells["barcodekey"]] = "MAD2 only"
    adata["color"][cutoff_cells["barcodekey"]] = "Cutoff only"
    adata["color"][list(set(cutoff_cells["barcodekey"]).intersection(set(mad_cells["barcodekey"])))] = "All"
    return adata


def joint_main(project, task_id, tissue=None):
    if tissue is None:
        tissue = get_project_info(project, task_id=task_id)[0]

    print("joint clustering task.id:{} - tissue:{}, res:{}, project:{}".format(task_id, tissue, resolution, project))

    task_directory, task_name, results_dir = prepare_joint_dir(project, tissue, resolution)
    adata = pd.read_csv(results_dir + "!cells.csv")
    adata.set_index(adata["barcodekey"], inplace=True)
    assign_colors(adata, project, tissue, resolution)

    print("Writing results")
    with open(results_dir + "!cells.csv", "w") as fout:
        fout.write(adata.to_csv())

    print(
        subprocess.check_output("Rscript plots/JC_plots.R {} {}".format(task_name, results_dir), shell=True).decode(
            'UTF-8'))


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
