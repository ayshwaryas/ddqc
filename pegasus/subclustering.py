import subprocess

import pandas as pd
import pegasus as pg

import paths
from config import OUTPUT_DIR
from filters import initial_qc
from mc import write_markers, save_to_csv
from readers import auto_reader
from utils import safe_mkdir, cluster_data, add_cd_scores

project = "mc_tm"
method = "mad"
task_id = 5 * 3 + 2
cluster = 30
res = 0.2

results_dir_prefix = OUTPUT_DIR + project + "/!subclustering/"

if method == "mad":
    TASKS_PER_TISS = 3
    method, param = (("none", 0), ("cutoff", 10), ("mad", 2))[task_id % TASKS_PER_TISS]
elif method == "joint":
    TASKS_PER_TISS = 1
else:
    exit()

tissue, is_human, adata = auto_reader(project, task_id, TASKS_PER_TISS)  # read the data for current task id

if method == "mad":
    results_dir = results_dir_prefix + tissue + "_" + str(res) + "-" + method + "-" + str(param) + "-c" + str(cluster) + "/"
    source_dir = OUTPUT_DIR + project + "/" + tissue + "/" + str(1.4) + "-" + method + "-" + str(param) + "/"
elif method == "joint":
    results_dir = results_dir_prefix + tissue + "_" + str(res) + "-" + method + "-c" + str(cluster) + "/"
    source_dir = OUTPUT_DIR + project + "/" + tissue + "/filtered_cells_plots/no_outlier/"
safe_mkdir(results_dir_prefix)
safe_mkdir(results_dir)
paths.set_results_dir(results_dir)

adata = initial_qc(adata, 100, 3, is_human)
cells = pd.read_csv(source_dir + "!cells.csv")
cluster_cells = cells[cells.louvain_labels == (cluster + 1)]
adata.obs["passed_qc"] = False
adata.obs["passed_qc"][cluster_cells["barcodekey"]] = True
pg.filter_data(adata)

adata, marker_dict = cluster_data(adata, compute_markers=True, compute_reductions=True, resolution=res)
adata = add_cd_scores(adata, is_human)

# write the results
write_markers(marker_dict, min_log_fc=0, min_pct=0)
save_to_csv(adata)

print(subprocess.check_output("Rscript r_subclustering_plots.R {} {} {} {} {} {}".format(project, tissue, res, method, param, results_dir),
                                  shell=True).decode('UTF-8'))
