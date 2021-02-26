import subprocess
import shutil
import sys

import pandas as pd
import pegasus as pg

import paths
from config import OUTPUT_DIR, SOURCE_DIR_PREFIX
from filters import initial_qc
from local_config import local
from mc import write_markers, save_to_csv
from readers import auto_reader
from utils import cluster_data, safe_mkdir, add_cd_scores
from projects_info import *

TASKS_PER_TISS = 1  # how many different methods per one tissue. Used to determine method and param from task id


# function that creates all the relevant directories
def prepare_fc_dirs(project, tissue, res, method):
    source_dir = OUTPUT_DIR + project + "/" + tissue + "/" + str(res) + "-none-0/"
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + "filtered_cells_plots/" + method + "/"  # directory for saving output

    # set results dir globally, so other functions can access it
    paths.set_results_dir(results_dir)

    # create all the directories, if they dont exist
    print("Creating Output Directories")
    safe_mkdir(OUTPUT_DIR + project + "/" + tissue + "/" + "filtered_cells_plots/")
    safe_mkdir(results_dir)

    shutil.copyfile(source_dir + "!cells.csv", results_dir + "!cells.csv")
    shutil.copyfile(source_dir + "!clusters.csv", results_dir + "!clusters.csv",)
    shutil.copyfile(source_dir + "!markers.csv", results_dir + "!markers.csv")

    return results_dir


def assign_colors(adata, project, tissue, res):
    mad_cells = pd.read_csv(SOURCE_DIR_PREFIX + project + "/" + tissue + "/" + str(res) + "-mad-2/!cells.csv")
    cutoff_cells = pd.read_csv(SOURCE_DIR_PREFIX + project + "/" + tissue + "/" + str(res) + "-cutoff-10/!cells.csv")

    adata["color"] = "Neither"
    adata["color"][mad_cells["barcodekey"]] = "MAD2 only"
    adata["color"][cutoff_cells["barcodekey"]] = "Cutoff only"
    adata["color"][list(set(cutoff_cells["barcodekey"]).intersection(set(mad_cells["barcodekey"])))] = "All"

    return adata


def main(project, task_id):
    tissue = get_tissue_by_task_id(project, task_id, TASKS_PER_TISS)
    res = 1.4
    # determine the method and param based on task id
    method = "joint"

    print("task.id:{} - tissue:{}, res:{}, project:{}, method:{}".format(task_id, tissue, res, project, method))

    prepare_fc_dirs(project, tissue, res, method)
    adata = pd.read_csv(paths.results_dir + "!cells.csv")
    adata.set_index(adata["barcodekey"], inplace=True)
    assign_colors(adata, project, tissue, res)

    with open(paths.results_dir + "!cells.csv", "w") as fout:
        fout.write(adata.to_csv())

    # launch seurat plot script
    #print(subprocess.check_output("Rscript r_fc_plots.R {} {} {} {} {}".format(project, tissue, res, method, "y"),
                                  #shell=True).decode('UTF-8'))


if __name__ == '__main__':
    if local:  # for debug outside of cluster
        proj = "mc_tm"
        for t_id in [5]:
            main(proj, t_id)
    else:  # project and task id are provided as commandline args
        proj = sys.argv[1]
        t_id = int(sys.argv[2]) - 1
        main(proj, t_id)
