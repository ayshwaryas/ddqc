csv_files = ["/Users/michaelalperovich/Documents/primes_storage/data/ct_markers.csv",]
             #"/Users/michaelalperovich/Documents/primes_storage/data/marker_genes/lung.csv"]

search_dict = dict()
for filename in csv_files:
    with open(filename, "r") as fin:
        header = []
        for line in fin.readlines():
            if not header:
                for ct in line.strip().split(","):
                    header.append(ct)
            else:
                l = line.strip().split(",")
                for i in range(len(l)):
                    if l[i]:
                        if l[i] not in search_dict:
                            search_dict[l[i]] = []
                        search_dict[l[i]].append(header[i])

while True:
    genes = input().split("; ")
    if not genes:
        exit()
    cell_type_dict = dict()
    for gene in genes:
        for ct in search_dict.get(gene.upper(), "-"):
            if ct == "-":
                continue
            if ct not in cell_type_dict:
                cell_type_dict[ct] = set()
            cell_type_dict[ct].add(gene)
        print(gene + ": " + "; ".join(search_dict.get(gene.upper(), "-")))
    for ct in cell_type_dict.keys():
        print(ct + ": " + "; ".join(sorted(list(cell_type_dict[ct]))))
    print("\n\n============\n")
