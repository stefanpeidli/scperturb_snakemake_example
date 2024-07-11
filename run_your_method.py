import numpy as np
import scanpy as sc
import pandas as pd
from itertools import product
from tqdm import tqdm
from argparse import ArgumentParser

# Setup CLI parser
parser = ArgumentParser(
                    prog='your_method',
                    description='Example program to run your method',
                    epilog='Text at the bottom of help :^)')
parser.add_argument('-i', '--input')
parser.add_argument('-o', '--output')
# access the arguments
args = parser.parse_args()

### HERE COULD BE THE CODE TO RUN YOUR METHOD ###
adata = sc.read_h5ad(args.input)

def pseudo_bulk(adata, keys, layer=None, min_cells_per_group=10):
    X = []
    Y = []
    for gs in tqdm(product(*[pd.unique(adata.obs[key]) for key in keys])):
        mask = np.logical_and.reduce([adata.obs[key]==g for g, key in zip(gs, keys)])
        ncells = sum(mask)
        if ncells < min_cells_per_group: continue
        Y.append(list(gs)+[ncells])
        X_ = adata[mask].layers[layer] if layer!=None else adata[mask].X
        X.append(np.array(np.sum(X_, axis=0), dtype=int)[0])
    obs = pd.DataFrame(Y, columns=list(keys)+['ncells'])
    return sc.AnnData(np.array(X), obs=obs, var=adata.var)

pseudobulk_data = pseudo_bulk(adata, ['perturbation'], min_cells_per_group=10)

pseudobulk_data.write(args.output)

print("Done!")