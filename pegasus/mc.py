import subprocess
import sys

import pandas as pd
import pegasus as pg

import paths
from config import do_counts, do_genes, do_mito, do_ribo, OUTPUT_DIR
from filters import filter_cells
from local_config import local
from readers import auto_reader
from utils import cluster_data, safe_mkdir

TASKS_PER_TISS = 4  # how many different methods per one tissue. Used to determine method and param from task id


# function that creates all the relevant directories
def create_dirs(tissue, res, method, param):
    res = str(res)
    param = str(param)
    task_directory = res + "-" + method + "-" + param  # name of the directory for this task
    task_name = tissue + "-" + task_directory  # task name to put on plots
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + task_directory + "/"  # directory for saving output

    # set results dir globally, so other functions can access it
    paths.set_results_dir(results_dir)

    # create all the directories, if they dont exist
    print("Creating Output Directories")
    safe_mkdir(OUTPUT_DIR)
    safe_mkdir(OUTPUT_DIR + project + "/")
    safe_mkdir(OUTPUT_DIR + project + "/" + tissue)
    safe_mkdir(results_dir)

    return task_directory, task_name, results_dir


# save all relevant cell info to csv for plots in seurat
def save_to_csv(adata):
    df = adata.obs

    # add dimensional reductions to data frame
    df["pca1"] = [t[0] for t in list(adata.obsm["X_pca"])]
    df["pca2"] = [t[1] for t in list(adata.obsm["X_pca"])]

    # tsne is turned off for now
    df["tsne1"] = [0] * len(df["pca1"])  # [t[0] for t in list(adata.obsm["X_fitsne"])]
    df["tsne2"] = [0] * len(df["pca2"])  # [t[1] for t in list(adata.obsm["X_fitsne"])]
    df["umap1"] = [t[0] for t in list(adata.obsm["X_umap"])]
    df["umap2"] = [t[1] for t in list(adata.obsm["X_umap"])]

    with open(paths.results_dir + "!cells.csv", "w") as fout:
        fout.write(df.to_csv())  # write df to csv


# write markers to csv
def write_markers(marker_dict, min_log_fc=0.25, min_pct=25):
    frames = []
    # iterate through all keys in the markers dict
    for cl in marker_dict.keys():
        for d in marker_dict[cl].keys():
            df = marker_dict[cl][d]
            df['cluster'] = cl
            df['up/down'] = d

            # filter markers based on log_fc and pct
            df = df[(df["mean_logExpr"] >= min_log_fc) & (
                    (df["percentage"] >= min_pct) | (df["percentage_other"] >= min_pct))]

            frames.append(df)
    result = pd.concat(frames) # merge all marker data frames together
    with open(paths.results_dir + "!markers.csv", "w") as fout:
        fout.write(result.to_csv())


def main():
    tissue, is_human, adata = auto_reader(project, task_id, TASKS_PER_TISS)  # read the data for current task id
    res = 1.4  # this resolution gives results closest to seurat
    # determine the method and param based on task id
    # none - no additional filtering; cutoff - min 200 genes, max 10% mito; outlier and mad - data driven methods
    method, param = (("none", 0), ("cutoff", 10), ("outlier", 0), ("mad", 2))[task_id % TASKS_PER_TISS]
    print(
        "task.id:{} - tissue:{}, res:{}, method:{}, param:{}, project:{}, do.counts:{}, do.genes:{}, do.mito:{}, do.ribo:{}".format(
            task_id, tissue, res, method, param, project, do_counts, do_genes, do_mito, do_ribo))

    task_directory, task_name, results_dir = create_dirs(tissue, res, method, param)
    adata = filter_cells(adata, res, method, param, do_counts, do_genes, do_mito, do_ribo)  # perform filtering
    adata, marker_dict = cluster_data(adata, compute_markers=True, compute_reductions=True, resolution=res)

    # write the results
    write_markers(marker_dict)
    save_to_csv(adata)
    pg.write_output(adata, results_dir + task_name)

    # launch seurat plot script
    print(subprocess.check_output("Rscript r_plots.R {} {} {} {}".format(project, task_id, tissue, res),
                                      shell=True).decode('UTF-8'))


if __name__ == '__main__':
    if local:  # for debug outside of cluster
        project = "mc_tm"
        task_id = 8
    else:  # project and task id are provided as commandline args
        project = sys.argv[1]
        task_id = int(sys.argv[2]) - 1
    main()
