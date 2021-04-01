import pegasus as pg
from filtering import filter_cells

adata = pg.read_input("~/Downloads/MantonBM_nonmix_subset.zarr.zip")
adata = filter_cells(adata)