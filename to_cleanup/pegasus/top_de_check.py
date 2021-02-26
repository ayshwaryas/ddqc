import re

import pandas as pd


dir_prefix = "/Users/michaelalperovich/Documents/primes_storage/"
results_dir = dir_prefix + "output_annotated/mc_ebi_tm/Cerebellum/filtered_cells_plots/no_outlier/"


markers = pd.read_csv(results_dir + "!markers.csv")
clusters = pd.read_csv(results_dir + "!clusters.csv").to_dict()

fout = open(results_dir + "!low_percent_DE.csv", "w")
fout2 = open(results_dir + "!high_percent_DE.csv", "w")
fout3 = open(results_dir + "!high_percent_DE_by_cluster.csv", "w")
fout.write("cluster,rank,gene,avgLogFC,percentage,p_val\n")
fout2.write("cluster,rank,gene,avgLogFC,percentage,p_val\n")
fout3.write("cluster,genes\n")

for i in range(len(clusters["Unnamed: 0"])):
    genes = clusters["markers"][i].split("; ")
    j = 1
    high_genes = []
    for gene in genes[:min(25, len(genes))]:
        if re.match("^Rp[sl]\d", gene, flags=re.IGNORECASE) or re.match("^MT-", gene, flags=re.IGNORECASE):
            continue
        marker = markers[(markers.feature == gene) & (markers.cluster == i + 1) & (markers["up/down"] == "up")]
        log_fc = list(marker.log2FC)[0]
        p_val = list(marker.t_pval)[0]
        percentage = list(marker.percentage)[0]
        if percentage < 50 and log_fc >= 1 and p_val <= 0.05:
            fout.write("{},{},{},{},{},{}\n".format(str(clusters["cluster"][i]), j, gene, log_fc, percentage, p_val))
        if percentage > 65 and log_fc >= 1 and p_val <= 0.05:
            fout2.write("{},{},{},{},{},{}\n".format(str(clusters["cluster"][i]), j, gene, log_fc, percentage, p_val))
            high_genes.append(gene)
        j += 1
    fout3.write(str(clusters["cluster"][i]) + "," + "; ".join(high_genes) + "\n")

fout.close()
fout2.close()