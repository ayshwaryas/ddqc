import csv
import os
import sys
from PyPDF2 import PdfFileReader, PdfFileWriter
#output_path = "/broad/hptmp/malperov/output_combined"
#data_path = "/broad/hptmp/malperov/output/"
output_path = "/Users/michaelalperovich/Documents/primes/results/"
data_path = "/Users/michaelalperovich/Documents/primes/results/"

def list_only_dirs(thedir):
    return [name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name))]


def pdf_cat(input_files, output_stream, res):
    #os.system("ulimit -n 4096")
    input_streams = []
    try:
        for input_file in input_files:
            if res in input_file:
                input_streams.append(open(input_file, 'rb'))
        if not input_streams:
            return
        writer = PdfFileWriter()
        for reader in map(PdfFileReader, input_streams):
            for n in range(reader.getNumPages()):
                writer.addPage(reader.getPage(n))
        writer.write(open(output_stream, "wb"))
    except OSError:
        print(input_files[0])
    finally:
        for f in input_streams:
            f.close()


def merge_tissue(tissue):
    fout = open(output_path + "summaries/" + tissue + ".csv", "w")
    writer = csv.writer(fout)
    header = ["Tissue", "Res/Mehtod", "Cluster", "Annotation", "Annotation2", "%Annotation1", "%Annotation2",
              "Cell Type", "Cell Type2", "Sum Avg_LogFC", "Delta Score", "Cells", "Mean Genes", "Median Genes", 
              "Mean Percent.mt", "Median Percent.mt"] + ["Marker " + str(t) for t in range(1, 11)] + ["Score Gene " + str(t) for t in range(1, 11)]
    writer.writerow(header)
    for d in sorted(list_only_dirs(data_path + "/" + tissue + "/")):
            fin = open(data_path + "/" + tissue + "/" + d + "/" + "!summary.csv")
            reader = csv.reader(fin)
            rows = [[tissue, d] + row[1:15] + row[15].split("; ")[:10] + row[16].split("; ")[:10] for row in reader][1:]
            writer.writerows(rows)
    fout.close()


data_name = sys.argv[1]
data_path += data_name + "/"
os.makedirs(output_path + "combined/")
output_path += "combined/" + data_name + "_combined/"
if not os.path.exists(output_path):
    os.makedirs(output_path)
    os.makedirs(output_path + "plots/")
    os.makedirs(output_path + "summaries/")
fout1 = open(output_path + "score_summary.csv", "w")
writer = csv.writer(fout1)
writer.writerows([["Tissue/Mehtod", "Avg_LogFC", "Delta Score", "Cells", "Mean Genes", "Mean Percent.mt"]])
for t in list_only_dirs(data_path):
    files = []
    try:
        merge_tissue(t)
    except FileNotFoundError as e:
        print(e)
    for d in sorted(list_only_dirs(data_path + t + "/")):
        if ".csv" not in d:
            for f in sorted(os.listdir(data_path + t + "/" + d)):
                if ".pdf" in f and f != "marker_plots":
                    files.append(data_path + t + "/" + d + "/" + f)
    for res in ["0.5-", "1-", "1.5-", "2-"]:
        pdf_cat(files, output_path + "plots/" + t + "-" + res[:-1] + ".pdf", res)
    try:
        fin = open(data_path + t + "/" "score_summary.csv")
        reader = csv.reader(fin)
        writer.writerows([[t, r[0][r[0].find("-") + 1:]] + r[1:] for r in reader])
    except FileNotFoundError:
        print(t)

        
    