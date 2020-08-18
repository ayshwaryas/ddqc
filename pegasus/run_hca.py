from mc_plot import main
import traceback

project = "mc_tm"
for task_id in range(12 * 4 ):
    try:
        main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
