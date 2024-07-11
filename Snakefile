"""
Author: Stefan Peidli
Aim: Example snakemake workflow for scperturb analysis.
Date: 11.07.2024
Run it as: snakemake
DAG: snakemake --forceall --dag | dot -Tpdf > snake_dag.pdf
Rulegraph: snakemake --forceall --rulegraph | dot -Tpdf > snake_rulegraph.pdf
"""

import pandas as pd
import numpy as np
import sys
import yaml
from pathlib import Path

# Load configs for this snakefile
configfile: "config.yaml"
DATA_DIR = Path(config["DATA_DIR"])
OUT_DIR = Path(config["OUT_DIR"])
DATASETS = config["DATASETS"]
DATASET_TABLE = config["DATASET_TABLE"]

# Load the dataset table, make selection of datasets (see config.yaml!) and check if they are available
dataset_table = pd.read_excel(DATASET_TABLE, index_col=0)
available_datasets = dataset_table.index[dataset_table["File Extension"]=='h5ad'].tolist()
if DATASETS == 'all':
    datasets = available_datasets
else:
    datasets = [d for d in DATASETS if d in available_datasets]
wildcard_constraints:
    dataset='|'.join(datasets)

# This here defines the files that you want this pipeline to generate.
# The expand function is basically a list comprehension that returns a list of 
# all possible output files.
rule all:
    input:
        expand(OUT_DIR / '{dataset}_your_method_output.h5ad', dataset=datasets)

# No need to change this rule
rule download_dataset:
	"Downloads an RNA datasets of scperturb from zenodo via wget."
	output: DATA_DIR / '{dataset}.h5ad'
	resources:
		time='4:00:00',
		mem_mb=4000,
		disk_mb=4000
	params:
		sleeptime=np.random.uniform(1,20),
		zenodo_url=config["ZENODO_RNA_URL"]
	shell:
		'''
		sleep {params.sleeptime}
		wget {params.zenodo_url}/files/{wildcards.dataset}.h5ad -O {output}
		'''

# Change this rule to use your method and specify the input and output files as needed.
rule run_your_method:
    "Runs your method on the downloaded dataset."
    input: 
        input_1 = DATA_DIR / '{dataset}.h5ad'
    output: 
        output_1 = OUT_DIR / '{dataset}_your_method_output.h5ad'
    resources:
        time='4:00:00',
        mem_mb=16000,
        disk_mb=16000
    conda:  # Optional, if you want to use a conda environment for your method
        "conda_env_for_your_method.yaml"
    shell:
        "python run_your_method.py -i {input.input_1} -o {output.output_1}"