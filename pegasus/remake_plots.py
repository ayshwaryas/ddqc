import subprocess

from readers import get_tissue_by_task_id

N_TISSUES_PROJ = {"mc_tm": 12, "mc_other_10X": 4, "mc_hca": 35, "mc_other": 5, "mc_manton": 1, "mc_PanglaoDB": 5,
                  "mc_blood": 1, "mc_brain_olfactory": 1, "mc_brain": 1, "mc_heart_circulation": 1}

task = "mc"  # mc or mc_plot
project = "mc_other"

if task == "mc":
    from mc import TASKS_PER_TISS

    script_name = "r_plots.R"
elif task == "mc_plot":
    from mc_plot import TASKS_PER_TISS

    script_name = "r_fc_plots.R"
else:
    exit()

max_tid = TASKS_PER_TISS * N_TISSUES_PROJ[project]

for tid in range(max_tid):
    res = 1.4
    tissue = get_tissue_by_task_id(project, tid, TASKS_PER_TISS)
    if task == "mc":
        method, param = (("none", 0), ("cutoff", 10), ("mad", 2))[tid % TASKS_PER_TISS]
        print(subprocess.check_output(
            "Rscript r_plots.R {} {} {} {} {}".format(project, tissue, res, method, param),
            shell=True).decode('UTF-8'))
    elif task == "mc_plot":
        method = "no_outlier"  # ("all", "no_outlier")[tid % TASKS_PER_TISS]
        print(subprocess.check_output("Rscript r_fc_plots.R {} {} {} {} ".format(project, tissue, res, method),
                                      shell=True).decode('UTF-8'))
