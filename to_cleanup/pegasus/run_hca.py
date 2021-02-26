from mc import main as mc_main
from mc_plot import main as mc_plot_main
from mc_new_joint import main as mc_new_joint_main
import traceback
from projects_info import N_TISSUES_PROJ


for project in ["mc_heart_circulation"]:#N_TISSUES_PROJ.keys():
    for task_id in range(N_TISSUES_PROJ[project]):
        print(project, task_id)
        try:
            mc_new_joint_main(project, task_id)
        except Exception as e:
            print("========= FAIL TASK_ID %s =========" % task_id)
            print(traceback.format_exc())
