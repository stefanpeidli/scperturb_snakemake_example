# scperturb_snakemake_example
This repository shows how you can use snakemake to download and use the datasets
in [scperturb](https://www.sanderlab.org/scPerturb/datavzrd/scPerturb_vzrd_v2/dataset_info/index_1.html).

Let's say we are developing a new method for single-cell perturbation analysis. We
want to evaluate our method on some of the datasets in scperturb. For the sake of this example,
we call our method "your_method". Let us realistically assume that your_method is
horribly slow. We want to run your_method on the datasets in scperturb, ideally in parallel.
This is where snakemake comes in. Snakemake is a workflow management system that allows you to
define a pipeline in a human-readable way.

In this example, your_method is a simple script that reads in a dataset, pseudo-bulks each perturbation
and saves the results. You can just replace the script with your actual method.

![DAG of the example snakemake workflow](snake_dag.pdf "DAG of the example snakemake workflow")


## Installation
- If you don't have conda installed, you can install it by following the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). I recommend using Mamba instead of conda as it is much faster.
- Install snakemake with `conda install bioconda::snakemake pandas numpy pyyaml`
- Clone this repository with `git clone URL_OF_THIS_REPOSITORY`
- Configure the example to your needs by editing the `config.yaml` file:
    - You can change the `datasets` variable to include the datasets you want to download. The datasets are specified by their names in scperturb. You can find the names of the datasets in the scperturb repository.
    Be aware that if you put "all" in there, it will download over 40 datasets with a total size of over 27GB.
    However, these files are compressed, so actually reading them in will take up much more space. I recommend starting with a few datasets first.
    - You should change the directories in the `config.yaml` file to match your system.

## Usage
- You can run the pipeline on your local machine e.g. via `snakemake --cores 4`
- If you want to run the pipeline on SLURM, please configure snakemake. See [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
