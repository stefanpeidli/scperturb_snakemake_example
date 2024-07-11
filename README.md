# scperturb_snakemake_example
This repository shows how you can use snakemake to download and use the datasets
in [scperturb](https://www.sanderlab.org/scPerturb/datavzrd/scPerturb_vzrd_v2/dataset_info/index_1.html).

Let's say we are developing a new method for single-cell perturbation analysis. We
want to evaluate our method on some of the datasets in scperturb. For the sake of this example,
we call our method "your_method". Let us realistically assume that your_method is
horribly slow. We want to run your_method on the datasets in scperturb, ideally in parallel.

## Installation
- If you don't have conda installed, you can install it by following the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). I recommend installing using Mamba instead as it is faster.
- Install snakemake with `conda install -c bioconda snakemake`
- Clone this repository with `git clone`
- Configure this example by editing the `config.yaml` file:
    - You can change the `datasets` variable to include the datasets you want to download. The datasets are specified by their names in scperturb. You can find the names of the datasets in the scperturb repository.
    - You should change the directories in the `config.yaml` file to match your system.

## Usage
- You can run the pipeline on your local machine e.g. via `snakemake --cores 4`
- If you want to run the pipeline on SLURM, please configure snakemake. See [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).