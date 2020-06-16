import subprocess
import os
import sys
from time import gmtime, strftime


def submit_command():
    param = SCRIPT_PARAMETERS[script]
    command = (COMMAND_PREFIX).format(max_task_id, param[0], 
        param[1], param[1], script + "_" + tissue, " ".join(args), 
        param[2], project)
    print(command)
    result = subprocess.check_output(command, shell=True).decode('UTF-8')
    print(result)
    result = result.split()[2:4]
    result = result[0].split(".") + [result[1][2:-2]]
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()), "|", "|".join(result), file=open("submissions.txt", "a"))
    return result[0]


TISSUE_COUNT = {
    "ebi": 13,
    "ebi_tm": 22,
    "mca": 24,
    "tm": 12,
    "ts24": 13,
    "ts30": 10,
    "other": 7,
    "other_10X": 4,
    "PanglaoDB": 5,
    
}

SCRIPT_PARAMETERS = {
    "mc": [150, 1, "scripts/mc.R", "job_scripts/init.sh"],
    "mc_plot": [150, 1, "scripts/mc_plot.R", "job_scripts/init.sh"],
    "pg_mc": [150, 1, "scripts/pegasus/mc.py", "job_scripts/pg_init.sh"],
    "pg_mc_plot": [150, 1, "scripts/pegasus/mc_plot.py", "job_scripts/pg_init.sh"],
    "custom": [8, 1, "scripts/custom.R", "job_scripts/init.sh"]
}

TASKS_PER_TISS_MC = 4
TASKS_PER_TISS_MC_PLOT = 1
dir_path = os.path.dirname(os.path.realpath(__file__))[:os.path.dirname(os.path.realpath(__file__)).rfind("/")]
COMMAND_PREFIX = "qsub -t 1-{} -l h_vmem={}G -pe smp {} -binding linear:{} -l h_rt=72:00:00 -j y -o logs/ -N {} {} job_scripts/init.sh " + dir_path + " {} {}"


script = sys.argv[1]
tissue = sys.argv[2]
if len(sys.argv) > 3:
    args = sys.argv[3:]
else:
    args = []

if tissue not in TISSUE_COUNT:
    print("invalid tissue: ", sys.argv[2])
    exit()

if script == "mc" or script == "custom":
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
