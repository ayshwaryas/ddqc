import os

directory_in = "/Volumes/easystore/primes_storage/data/mca/"
directory_out = "/Volumes/easystore/primes_storage/data/mca_pg/"

done = """
.DS_Store
"""


for file in os.listdir(directory_in):
    if file in done:
        continue
    print(file, "started")
    fin = open(directory_in + file, "r")
    fout = open(directory_out + file.replace(".txt", ".csv"), "w")
    for line in fin.readlines():
        line = line.replace("\"", "")
        line = line.replace(" ", ",")
        fout.write(line)
    fin.close()
    fout.close()
    print(file, "complete")
