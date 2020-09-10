import pegasus as pg


adata = pg.read_input("~/Downloads/ebi_tm.h5ad").to_anndata()
for tissue in set(adata.obs["tissue"]):
    tiss = adata[adata.obs.tissue == tissue]
    tiss.write_h5ad("/Users/michaelalperovich/Downloads/" + tissue + ".h5ad")
    print(set(tiss.obs.tissue))
