from mc import main as mc_main
from mc_plot import main as mc_plot_main
import traceback

project = "mc_other_10X"
for task_id in range(9, 12):
    try:
        mc_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
        exit()
for task_id in range(3, 4):
    try:
        mc_plot_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
        exit()
