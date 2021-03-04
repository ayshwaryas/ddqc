import os

import pandas as pd
import pegasusio as io

from config.config import DATA_DIR


def get_project_info(project, task_id=None, tissue=None):  # function that parses projects.csv and returns relevant info
    projects = pd.read_csv("config/read_info/projects.csv")
    assert project in set(projects['project'])  # check if project exists in project list
    project_tissues = projects[projects['project'] == project]  # subset tissues for the project
    if task_id:  # get info based on task_id
        assert task_id < len(project_tissues['project'])
        tissue = project_tissues['tissue'][task_id]
        is_human = project_tissues['is_human'][task_id]
        annotations = project_tissues['annotations'][task_id]
        return tissue, is_human, annotations
    elif tissue:  # get info based on tissue
        assert tissue in set(project_tissues['tissue'])
        tissue_info = project_tissues[project_tissues['tissue'] == tissue]
        is_human = tissue_info['is_human'][1]
        annotations = tissue_info['annotations'][1]
        return is_human, annotations
    else:  # if no tissue or task_id provided return pandas project with project info
        return project_tissues


def read_tissue(project, tissue, annotations="Unknown"):  # function that reads and aggregates one tissue
    dataset_list_path = "config/read_info/{}/{}.csv".format(project, tissue)  # path to tissue read info
    read_info_filename = "read_info_{}_{}.csv".format(project, tissue)  # filename of a current read info copy
    assert os.path.isfile(dataset_list_path)  # check if tissue read info exists
    dataset_list = pd.read_csv(dataset_list_path)
    dataset_list['Location'] = [DATA_DIR + t for t in dataset_list['Location']]  # update location with relevant directory prefix
    with open(read_info_filename, "w") as fout:  # write modified csv
        fout.write(dataset_list.to_csv())
    adata = io.aggregate_matrices(read_info_filename)
    os.remove(read_info_filename)  # remove current read info copy

    # add annotations to adata
    if annotations != 'Unknown':  # TODO: Add Annotations Support
        ann_df = pd.read_csv(DATA_DIR + annotations)
        ann_df = ann_df.reindex(adata.obs.index)["annotations"]
        ann_df = ann_df.fillna("Unknown")
        adata.obs["annotations"] = ann_df
    else:
        adata.obs['annotations'] = 'Unknown'
    return adata
