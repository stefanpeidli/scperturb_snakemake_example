import decoupler as dc
import scanpy as sc

# Setup CLI parser
from argparse import ArgumentParser
parser = ArgumentParser(
                    prog='your_method',
                    description='Example program to run your method',
                    epilog='Text at the bottom of help :^)')
parser.add_argument('-i', '--input')
parser.add_argument('-o', '--output')
# access the arguments
args = parser.parse_args()

### HERE WOULD BE THE CODE TO RUN YOUR METHOD ###
adata = sc.read_h5ad(args.input)

pseudobulk_data = dc.get_pseudobulk(
    adata,
    sample_col='individual',
    groups_col='perturbation',
    layer='counts',
    mode='sum',
    min_cells=20,
    min_counts=10
)

pseudobulk_data.write(args.output)

print("Done!")