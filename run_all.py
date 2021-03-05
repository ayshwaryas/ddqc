import traceback

from method_comparison import mc_main
from joint_clustering import joint_main
from config.config import MC_TASKS_PER_TISSUE
from reading import get_project_info

log = open("log.txt", "w+")
for project in set(get_project_info()["project"]):
    for tissue in get_project_info(project)["tissue"]:
        for method in range(MC_TASKS_PER_TISSUE):
            try:
                mc_main(project, task_id=method, tissue=tissue)
            except Exception as e:
                print("FAILED MC {} {} method {}:\n".format(project, tissue, method, traceback.format_exc()), file=log)
        try:
            joint_main(project, 0, tissue=tissue)
        except Exception as e:
            print("FAILED JC {} {}:\n".format(project, tissue, traceback.format_exc(), file=log))

log.close()
