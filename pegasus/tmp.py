import pandas as pd


markers = pd.read_csv("~/Downloads/markers3.csv")
markers = markers[(markers.p_val < 0.05) & (markers.avg_logFC > 1)]
cell_types = dict()

for i in range(len(markers.gene)):
    ct = list(markers.cluster)[i]
    gene = list(markers.gene)[i]
    avg_log_fc = list(markers.avg_logFC)[i]
    p_val = list(markers.p_val)[i]
    if p_val < 0.05 and avg_log_fc > 0.7:
        if ct not in cell_types:
            cell_types[ct] = []
        cell_types[ct].append(gene)

markers_list = []
for ct in sorted(list(cell_types.keys())):
    markers_list.append([ct + "/HHA"] + cell_types[ct][:50])

j = 0
f = True
while f:
    f = False
    for i in range(len(markers_list)):
        if j >= len(markers_list[i]):
            print(end=",")
        else:
            print(markers_list[i][j], end=",")
            f = True
    print()
    j += 1
