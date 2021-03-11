import sys

from joint_clustering_old import joint_main

project = sys.argv[1]
tissue = sys.argv[2]
joint_main(project, 0, tissue=tissue)
