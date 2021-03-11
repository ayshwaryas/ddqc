import os

from config.config import MITO_PREFIXES, RIBO_PREFIXES, MC_METHODS, resolution, \
    basic_genes_filter, basic_mito_filter, OUTPUT_DIR
from filtering import filter_cells
from method_comparison import create_dirs
from reading import get_project_info, read_tissue


def check_finished(prj, tiss):
    results_dir = OUTPUT_DIR + prj + "/" + tiss + "/"
    return os.path.isdir(results_dir)


def run_htapp():
    for i, row in get_project_info(project="htapp").iterrows():
        project = row["project"]
        tissue = row["tissue"]

        if check_finished(project, tissue):
            print("Skipping {} {}".format(project, tissue))
            continue

        is_human, annotations = get_project_info(project, tissue=tissue)
        adata = read_tissue(project, tissue, annotations)
        method, param = MC_METHODS[0]
        task_directory, task_name, results_dir = create_dirs(project, tissue, resolution, method, param)
        # filtering
        mito_prefix = MITO_PREFIXES["human"] if is_human else MITO_PREFIXES["mouse"]
        ribo_prefix = RIBO_PREFIXES["human"] if is_human else RIBO_PREFIXES["mouse"]
        adata = filter_cells(adata, resolution, method, param, basic_n_genes=basic_genes_filter,
                             basic_percent_mito=basic_mito_filter, mito_prefix=mito_prefix, ribo_prefix=ribo_prefix)
        with open(results_dir + "!cells.csv", "w") as fout:
            fout.write(adata.obs.to_csv())  # write df to csv


run_htapp()
