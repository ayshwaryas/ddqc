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

TASKS_PER_TISS = 4


def create_dirs(tissue, res, method, param):
    res = str(res)
    param = str(param)
    task_directory = res + "-" + method + "-" + param
    task_name = tissue + "-" + task_directory
    # robjs_dir = OUTPUT_DIR + "robjs/" + project + "/"  # directory for saving R objects
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + task_directory + "/"  # directory for saving other output

    paths.set_results_dir(results_dir)

    print("Creating Output Directories")
    safe_mkdir(OUTPUT_DIR)
    safe_mkdir(OUTPUT_DIR + project + "/")
    safe_mkdir(OUTPUT_DIR + project + "/" + tissue)
    safe_mkdir(results_dir)

    return task_directory, task_name, results_dir


def save_to_csv(adata):
    df = adata.obs
    df["pca1"] = [t[0] for t in list(adata.obsm["X_pca"])]
    df["pca2"] = [t[1] for t in list(adata.obsm["X_pca"])]
    df["tsne1"] = [t[0] for t in list(adata.obsm["X_fitsne"])]
    df["tsne2"] = [t[1] for t in list(adata.obsm["X_fitsne"])]
    df["umap1"] = [t[0] for t in list(adata.obsm["X_umap"])]
    df["umap2"] = [t[1] for t in list(adata.obsm["X_umap"])]

    with open(paths.results_dir + "!cells.csv", "w") as fout:
        fout.write(df.to_csv())


def write_markers(marker_dict, min_log_fc=0.25, min_pct=25):
    frames = []
    for cl in marker_dict.keys():
        for d in marker_dict[cl].keys():
            df = marker_dict[cl][d]
            df['cluster'] = cl
            df['up/down'] = d

            df = df[(df["mean_logExpr"] >= min_log_fc) & (
                    (df["percentage"] >= min_pct) | (df["percentage_other"] >= min_pct))]

            frames.append(df)
    result = pd.concat(frames)
    with open(paths.results_dir + "!markers.csv", "w") as fout:
        fout.write(result.to_csv())


def main():
    tissue, is_human, adata = auto_reader(project, task_id, TASKS_PER_TISS)
    res = 1.4
    method, param = (("none", 0), ("cutoff", 10), ("outlier", 0), ("mad", 2))[task_id % TASKS_PER_TISS]
    print(
        "task.id:{} - tissue:{}, res:{}, method:{}, param:{}, project:{}, do.counts:{}, do.genes:{}, do.mito:{}, do.ribo:{}".format(
            task_id, tissue, res, method, param, project, do_counts, do_genes, do_mito, do_ribo))

    task_directory, task_name, results_dir = create_dirs(tissue, res, method, param)

    adata = filter_cells(adata, res, method, param, do_counts, do_genes, do_mito, do_ribo)

    # pg.qc_metrics(adata, mito_prefix="mt-")
    # pg.filter_data(adata)

    adata, marker_dict = cluster_data(adata, compute_markers=True, compute_reductions=True, resolution=res)

    write_markers(marker_dict)
    save_to_csv(adata)
    pg.write_output(adata, results_dir + task_name)

    print(subprocess.check_output("Rscript r_plots.R {} {} {} {}".format(project, task_id, tissue, res),
                                      shell=True).decode('UTF-8'))


if __name__ == '__main__':
    if local:
        project = "mc_tm"
        task_id = 8
    else:
        project = sys.argv[1]
        task_id = int(sys.argv[2]) - 1
    main()
