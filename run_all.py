import os
import sys
import traceback

from method_comparison import mc_main
from joint_clustering import joint_main
from config.config import MC_TASKS_PER_TISSUE, OUTPUT_DIR
from reading import get_project_info


def check_finished(prj, tiss):
    results_dir = OUTPUT_DIR + prj + "/" + tiss + "/"
    return os.path.isdir(results_dir)


def run_all(organism=None, project=None, tissue=None):
    log = open("log.txt", "w+")
    for i, row in get_project_info(project=project, tissue=tissue).iterrows():
        project = row["project"]
        tissue = row["tissue"]
        is_human = row["is_human"]

        if (organism == "mouse" and is_human) or (organism == "human" and not is_human):
            continue

        if check_finished(project, tissue):
            print("Skipping {} {}".format(project, tissue))
            continue

        for method in range(MC_TASKS_PER_TISSUE):
            try:
                mc_main(project, task_id=method, tissue=tissue)
            except Exception as e:
                print("FAILED MC {} {} method {}:\n{}".format(project, tissue, method, traceback.format_exc()), file=log)
        try:
            joint_main(project, 0, tissue=tissue)
        except Exception as e:
            print("FAILED JC {} {}:\n{}".format(project, tissue, traceback.format_exc()), file=log)

        print("Finished {} {}".format(project, tissue))

    log.close()


run_all(sys.argv[1], sys.argv[2], sys.argv[3])
