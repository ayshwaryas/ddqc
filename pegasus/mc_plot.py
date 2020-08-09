import subprocess
import sys

import pandas as pd

import paths
from config import OUTPUT_DIR, SOURCE_DIR_PREFIX
from filters import initial_qc
from local_config import local
from mc import write_markers, save_to_csv
from readers import auto_reader
from utils import cluster_data, safe_mkdir, add_cd_scores

TASKS_PER_TISS = 2  # how many different methods per one tissue. Used to determine method and param from task id


# function that creates all the relevant directories
def create_fc_dirs(tissue, method):
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + "filtered_cells_plots/" + method + "/"  # directory for saving output

    # set results dir globally, so other functions can access it
    paths.set_results_dir(results_dir)

    # create all the directories, if they dont exist
    print("Creating Output Directories")
    safe_mkdir(OUTPUT_DIR + project + "/" + tissue + "/" + "filtered_cells_plots/")
    safe_mkdir(results_dir)

    return results_dir


def filter_cells_by_csv(adata, tissue, res, method):
    mad_cells = pd.read_csv(SOURCE_DIR_PREFIX + project + "/" + tissue + "/" + str(res) + "-mad-2/!cells.csv")
    cutoff_cells = pd.read_csv(SOURCE_DIR_PREFIX + project + "/" + tissue + "/" + str(res) + "-cutoff-10/!cells.csv")
    outlier_cells = pd.read_csv(SOURCE_DIR_PREFIX + project + "/" + tissue + "/" + str(res) + "-outlier-0/!cells.csv")

    adata.obs["color"] = "Did Not Pass"
    adata.obs["color"][mad_cells["Unnamed: 0"]] = "MAD2 only"
    adata.obs["color"][cutoff_cells["Unnamed: 0"]] = "C10 only"

    if method == "no_outlier":
        adata.obs["color"][list(set(cutoff_cells["Unnamed: 0"]).intersection(set(mad_cells["Unnamed: 0"])))] = "All"
    elif method == "all":
        adata.obs["color"][outlier_cells["Unnamed: 0"]] = "Outlier only"
        adata.obs["color"][
            list(set(cutoff_cells["Unnamed: 0"]).intersection(set(mad_cells["Unnamed: 0"])))] = "MAD2 and C10"
        adata.obs["color"][
            list(set(cutoff_cells["Unnamed: 0"]).intersection(set(outlier_cells["Unnamed: 0"])))] = "Outlier and C10"
        adata.obs["color"][
            list(set(outlier_cells["Unnamed: 0"]).intersection(set(mad_cells["Unnamed: 0"])))] = "MAD2 and Outlier"

        adata.obs["color"][list(set(cutoff_cells["Unnamed: 0"]).intersection(set(mad_cells["Unnamed: 0"])).intersection(
            set(outlier_cells["Unnamed: 0"])))] = "All"

    adata = adata[adata.obs.color != "Did Not Pass"]

    return adata


def main():
    tissue, is_human, adata = auto_reader(project, task_id, TASKS_PER_TISS)  # read the data for current task id
    res = 1.4  # this resolution gives results closest to seurat
    # determine the method and param based on task id
    method = ("all", "no_outlier")[task_id % TASKS_PER_TISS]

    print("task.id:{} - tissue:{}, res:{}, project:{}, method:{}".format(task_id, tissue, res, project, method))

    create_fc_dirs(tissue, method)
    adata = filter_cells_by_csv(adata, tissue, res, method)
    adata = initial_qc(adata, 100, 3, is_human)
    adata, marker_dict = cluster_data(adata, compute_markers=True, compute_reductions=True, resolution=res)
    adata = add_cd_scores(adata, is_human)

    # write the results
    write_markers(marker_dict)
    save_to_csv(adata)

    # launch seurat plot script
    print(subprocess.check_output("1 {} {} {} {} {}".format(project, task_id, tissue, res, method),
                                  shell=True).decode('UTF-8'))


if __name__ == '__main__':
    if local:  # for debug outside of cluster
        project = "mc_tm"
        for task_id in [4]:
            main()
    else:  # project and task id are provided as commandline args
        project = sys.argv[1]
        task_id = int(sys.argv[2]) - 1
        main()

