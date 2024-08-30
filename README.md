# structural diversity of microbial genomes

This repository contains a snakemake pipeline to analyze the structural diversity of microbial genomes using [PanGraph](https://github.com/neherlab/pangraph).
The pipeline is a simplified version of [the one used for the analysis of _E. coli_ ST131 genomes](https://github.com/mmolari/ecoliST131-structural-evo) in [this preprint](https://www.biorxiv.org/content/10.1101/2024.07.08.602537v1).

The input consists of a set of genbank files, containing each a single record (the bacterial chromosome) forming a collection of elements that can be meaningfully structurally compared.

## setup

### conda environment creation

The pipeline requires a working installation of [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/#). To run, it requires an environment with snakemake (v.7+). This can be created with:

```sh
conda create -n snakemake -c conda-forge -c bioconda snakemake=7
```

### input data

A new dataset can be added with the following three steps:

- all input data files can be placed in the `data/gbk/{acc}.gbk` directory, where `{acc}` is an id of the record (e.g. accession number).
- in addition, a file `datasets/{dataset_name}/acc_nums.txt` should contain a list of ids for all entries of the dataset. Here `{dataset_name}` is the name of the dataset.
- finally, update the `config.yaml` file by adding the new dataset name under the `datasets` key.
  ```yaml
  datasets:
    {dataset_name}:
        guide-strain: "{acc}"
  ```
  In addition you must specify the id of a `guide-strain` from the dataset. This will be used as a reference for the structural comparison.


## running the pipeline

After activating the environment with:
```sh
conda activate snakemake
```

You can run the pipeline on a cluster with SLURM workload manager with:
```sh
snakemake all --profile cluster
```

Or alternatively on your local machine with:

```sh
snakemake -c <n. cores> --use-conda all
```