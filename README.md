
# Apollo-reference: making snakemake apollo-*** pipelines supporting multiple reference genomes

!!IMPORTANT!! This repository hasn't been configured for execution using Snakemake (yet).

## It has 3 purposes:

### [1] connects multi-reference indices, species selection and elegiable reference sequences (.fa) to other apollo-**** pipelines

- based on these files (currently hard-coded and within this repo):
    - supported-reference-species.tsv
    - reference_assembly_data.tsv
    - accession2SRR.tsv
    - see "workflow-update-reference-assemblies.sh"
- provides re-useable argparse additions for other "run_pipeline.py" scripts
    - TODO: currently these are still in apollo-mapping/apollo_mapping.py!
- by default, configured to RIVM-internal usage on HPC environment
    - see config/referencedata.yaml
- which can be easily tweaked into your own local calculation environment

### [2] automated workflow to download stated references

- implemented in/by workflow-update-reference-assemblies.sh
- currently hard-coded to "Apollo_clinical_Candida_species.xlsx", within this repo
- creates supported-reference-species.tsv (by xlsx2csv copy)
- creates reference_assembly_data.tsv (for internal usage)
- workflow creates as well a reference-selection-index database

### [3] partially automated workflow to identify regions with low ANI in multi-clade references

- current use-case Candida auris Clades I-II-III-IV-V
- All clades are mapped (both .fa as corresponding PE .fastq) to a "central main" reference
- mainly works on identification of predominantly softclipped areas and obvious deleted areas
- eventually creates a/some BED file(s) with regions to ignore in variant calling

## For a visual on what both workflows do (downloading & "softclipping" determination"):

- files/apollo-reference.drawio
- files/apollo-reference.drawio.png
- files/apollo-multiclade-reference-masking.drawio
- files/apollo-multiclade-reference-masking.drawio.png

## The event of a redefinition/expansion of multi-species scope will be rarely undertaken.

- (future) fungal reference species selection was fast-forwarded in 2026 Q1/Q2.
- momentarily multi-species scope has been defined to this species dataset alone
- it's rare enough to postpone full automation into a Snakemake pipeine
- it's frequent enough to prepare the full automation of it (in bash)


## Installation & usage

Depending on your purpose with apollo-reference (read on [1], [2], and [3] above), this is how to use it:

### [1a] install it (stand-alone python code, e.g. to use match-ref.py)

There are several options here

1. install it yourself the traditional way
- the very few dependancies are stated in requirements.txt
- this assumes you have standard workhorses minimap2 and samtools already installed
- and .. ready!

2. use the co-packaged conda environment (with conda, mamba or micromamba)
- ./workflow/envs/download_reference_dataset.yaml

3. install the repo locally
```
git clone --depth 1 https://github.com/RIVM-bioinformatics/apollo-reference.git .
pip install --user apollo-reference
```

### [1b] install it (as installable requirement in another repo) 
```
see [apollo-mapping repository](https://github.com/RIVM-bioinformatics/apollo-mapping/pyproject.toml) for an example how to this (assuming hatch)
```

### [2] download stated references

Below instuctions will use micromamba on the conda environment and (incrementally!) downloads your reference dataset.

```
## install micromamba if not installed yet
#curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | sudo tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba
micromamba create -n download_reference_dataset -f workflow/envs/download_reference_dataset
eval "$(micromamba shell hook --shell bash)"
micromamba activate download_reference_dataset
./workflow/workflow-update-reference-assemblies.sh --out /your/desired/output/directory
```

Some important considerations:
- Typically, you would want to use [apollo-mapping](https://github.com/RIVM-bioinformatics/apollo-mapping) and thus need apollo-reference's data
- apollo_mapping.py (without params) expects the reference data at /mnt/db/apollo/reference
- Adjust this path in apollo_reference/config/referencedata.yaml
- Alternatively mv /your/desired/output/directory to /mnt/db/apollo/reference
- Momentarily, apollo-reference is "hard-coded" to one particular excel sheet (data/Apollo_clininal_Candida_species.xlsx)
- In case you want your own set of species, carefully (!) edit this excel sheet

### [3] multi-clade reference low ANI region identification 
```
needs instructions ...
```



