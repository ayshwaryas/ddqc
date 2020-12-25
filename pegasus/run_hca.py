from mc import main as mc_main
from mc_plot import main as mc_plot_main
import traceback

project = "mc_ts30"
for task_id in [15, 16, 17, 24, 25, 26]:
    try:
        mc_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
for task_id in [5, 8]:
    try:
        mc_plot_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())

project = "mc_ts24"
for task_id in range(6, 9):
    try:
        mc_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
for task_id in range(3, 4):
    try:
        mc_plot_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())