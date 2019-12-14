import subprocess
import os
import sys
from time import gmtime, strftime


def submit_command():
    command = (COMMAND_PREFIX + SCRIPT_PATH[script] + " {}").format(max_task_id, script + "_" + tissue, " ".join(args), project)
    print(command)
    result = subprocess.check_output(command, shell=True).decode('UTF-8')
    print(result)
    result = result.split()[2:4]
    result = result[0].split(".") + [result[1][2:-2]]
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()), "|", "|".join(result), file=open("submissions.txt", "a"))
    return result[0]


TISSUE_COUNT = {
    "ebi": 13,
    "mca": 24,
    "tm": 12,
    "ts24": 13,
    "ts30": 10,
    "other": 7,
}
SCRIPT_PATH = {
    "mc": "job_scripts/mc_init.sh",
    "mc_plot": "job_scripts/mc_plot_init.sh"
}

TASKS_PER_TISS_MC = 16
TASKS_PER_TISS_MC_PLOT = 4
COMMAND_PREFIX = "qsub -t 1-{} -N {} {} "

script = sys.argv[1]
tissue = sys.argv[2]
if len(sys.argv) > 3:
    args = sys.argv[3:]
else:
    args = []

if tissue not in TISSUE_COUNT:
    print("invalid tissue: ", sys.argv[2])
    exit()

if script == "mc":
    max_task_id = TISSUE_COUNT[tissue] * TASKS_PER_TISS_MC
    project = ("mc_" + tissue)
    submit_command()
elif script == "mc_plot":
    max_task_id = TISSUE_COUNT[tissue] * TASKS_PER_TISS_MC_PLOT
    project = ("mc_" + tissue)
    submit_command()
elif script == "mc+mc_plot":
    max_task_id = TISSUE_COUNT[tissue] * TASKS_PER_TISS_MC
    script = "mc"
    project = ("mc_" + tissue)
    job_id = submit_command()
    
    args.append("-hold_jid {}".format(job_id))
    max_task_id = TISSUE_COUNT[tissue] * TASKS_PER_TISS_MC_PLOT
    script = "mc_plot"
    project = ("mc_" + tissue)
    submit_command()
else:
    print("invalid task type: ", sys.argv[1])
    exit()