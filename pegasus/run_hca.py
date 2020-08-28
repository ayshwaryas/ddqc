from mc import main
import traceback

project = "mc_other_10X"
for task_id in range(5 * 3):
    try:
        main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
