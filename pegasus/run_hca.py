from mc import main
import traceback

project = "mc_hca"
for task_id in range(140):
    try:
        main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
