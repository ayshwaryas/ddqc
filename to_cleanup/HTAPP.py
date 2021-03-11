import os

from config.config import DATA_DIR
import pandas as pd

path = DATA_DIR + "human/HTAPP/"
read_info_dir = "../config/read_info/htapp/"
project_csv = open("../config/read_info/projects.csv", "a")
project = "htapp"
is_human = True
genome = "GRCh38"

for tissue in sorted(os.listdir(path)):
    if os.path.isdir(path + tissue):
        ann_path = path.split(DATA_DIR)[1] + tissue + "/annotations_pg.csv"
        project_csv.write("{},{},{},{}\n".format(project, tissue, str(is_human), ann_path))
        with open(read_info_dir + tissue + ".csv", "w") as fout:
            fout.write("Sample,Location,Reference\n")
            for file in sorted(os.listdir(path + tissue + "/")):
                if file.endswith(".h5"):
                    fout.write("{},{},{}\n".format(file.split("_raw_feature")[0], (path + tissue + "/" + file).split(DATA_DIR)[1], genome))

    # if file.endswith(".h5"):
    #     tissue = (file.split("_channel")[0]).split("_", 1)[1]
    #     ann_path = path.split(DATA_DIR)[1] + tissue + "_annotations_pg.csv"
    #     project_csv.write("{},{},{},{}\n".format(project, tissue, str(is_human), ann_path))
    #     with open(read_info_dir + tissue + ".csv", "w") as fout:
    #         fout.write("Sample,Location,Reference\n")
    #         fout.write("{},{},{}\n".format(tissue, (path + file).split(DATA_DIR)[1], genome))
    # if file.endswith(".csv"):
    #     metadata = pd.read_csv(path + file)
    #     tissue = (file.split("_channel")[0]).split("metadata_")[1]
    #     barcodes = [t.rsplit("-", 1)[1] for t in list(metadata["Unnamed: 0"])]
    #     annotations = pd.DataFrame({"barcodekey": barcodes, "annotations": metadata["annotate"]})
    #     with open(path + tissue + "_annotations_pg.csv", "w") as fout:
    #         fout.write(annotations.to_csv())

project_csv.close()
