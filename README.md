
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
# here, probably you want a pyhton venv, which I assume you'll take for yourself
git clone --depth 1 https://github.com/RIVM-bioinformatics/apollo-reference.git .
pip install --user apollo-reference
```

### [1b] install it (as installable requirement in another repo) 
```
see [apollo-mapping repository](https://github.com/RIVM-bioinformatics/apollo-mapping/pyproject.toml) for an example how to do this (assuming hatch)
```

### [2] download stated references

**For RIVM-idsbioinfo colleagues wanting to update the reference dataset, please read [2.4] first before you kick-off!**

Below instuctions will use micromamba on the conda environment and (incrementally!) downloads your reference dataset.
It's important to realize that the workflow has been designed to **incrementally** download references. So:
- adding one extra entry will only result in new data being downloaded
- adding/removing/changing any accession will result in a new "dataset" definition
- and result in a new minimap2 index for species identification
This allows (earlier) versions of the multi-reference dataset to persist on the filesystem,
which is a requirement for (upcoming) ISO-validation of the [apollo-mapping](https://github.com/RIVM-bioinformatics/apollo-mapping) pipeline

```
## [2.0] install micromamba if not installed yet
#curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | sudo tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

## [2.1] creare environment with micromamba (or with conda if you prefer)
micromamba create -n download_reference_dataset -f workflow/envs/download_reference_dataset.yaml

## [2.2] activate environment
# works at RIVM-HPC headnode (and probably on any compute environments)
env_path=$(conda info --envs | grep download_reference_dataset | awk '{ print $1 }' | tee /dev/stderr)
conda activate $env_path

## [2.3] execute the workflow that will download the multi-reference dataset ylocally
# Please have some coffee meanwhile, since:
# - dozens of assemblies and taxrecords
# - hundreds of mitochondria
# Will get downloaded and pre-processed.
# At this moment in time [2026/05/22], when downloading de novo, it will claim ~4GB of disk space.
# Mind, that since it's incrementally generating (unique) species identification databases,
# each accession change will trigger an additional ~3GB of disk space.
# Realize disk space used will inflate once the supported reference dataset will grow considerably.
./workflow/workflow-update-reference-assemblies.sh --out /your/desired/output/directory

## [2.4] !!important!! for RIVM-idsbioinfo colleagues that need to update the supported reference dataset
# --------------------------------------------------------------------------------------------------------------
# In case you've downloaded multi-references for external useage, you can safely skip this step
# The apollo-reference repository contains a .xlsx (to stay biologist-fiendly) and generated *.tsv files connected to it.
# Together, these represent the *state* of the multi-reference dataset.
# Any addition/removal/changing of an accession (in .xlsx) changes the files and thereby the state.
# So, after an supported reference dataset update, information should get commited back into the repo.
# - assuming there's indeed any update (check git status!)
# - !important! this should become a new version bumb of the repo (allowing pinned version ISO-verification)
# - please commit the updates
# --------------------------------------------------------------------------------------------------------------

# Assuming you're in the cloned repo's directory.
# Make sure you:
# a) have the excel file in your HPC environment
# b) alternatively, branch to a $bumped_version already now and commit the *.xlxs file from your laptop
cp ~/Downloads/Apollo_clininal_Candida_species.xlsx data/
# now steps [2.1], [2.2], [2.3] 
./workflow/workflow-update-reference-assemblies.sh --out </path/on/hpc/compute, see config/referencedata.yaml>
git status
bumped_version="<one-notched-up-version-bumb-here>"
git branch -m $bumped_version
git add data/*.tsv
git add data/*.xlsx # might have been added by you already
git commit -m "xxxx: multi-reference dataset upgraded to $bumped_version"
git push

## [2.2a] activate environment using micromamba
# works at a desktop/laptop (Ubuntu 24.04), but **NOT** on RIVM-HPC headnode
# severe issue with shell bash-fork bomb
#eval "$(micromamba shell hook --shell bash)"
#micromamba activate download_reference_dataset

## [2.3b] Alternative solution to execute the workflow
# In case activating the environment doesn't work in your hands (e.g. microconda on HPC and bound to micromamba),
# work with the scripts directly IN the actual environment's path

## [2.3b.1] Define the absolute path to your specific environment folder
ENV_PATH="/home/username/micromamba/envs/jouw_omgevings_naam"

## [2.3b.2] Forcefully prepend the environment's binaries to your system PATH
export PATH="$ENV_PATH/bin:$PATH"

## [2.3b.3] Explicitly export environment roots
export CONDA_PREFIX="$ENV_PATH"

## [2.3b.4] Launch the workflow on its exact location
bash /home/username/micromamba/envs/jouw_omgevings_naam/workflow/workflow-update-reference-assemblies.sh

# {"error":"API rate limit exceeded","api-key":"131.224.251.101","count":"4","limit":"3"}
# sleep 2->5 in workflow/scripts/download-mito-accessions.sh
# TODO: make function for "download it", check for {"error:", and in case it happes, sleep (long) and increase sleeptime
# error: no such file: /mnt/db/apollo/reference/reference_assembly_data.tsv
# error: no such file: /mnt/db/apollo/reference/reference_assembly_data.tsv
# TODO: remove symlink
# TODO: discard dated versions
# TODO:     - remind that the user should commit reference_assembly_data.tsv (in apollo-mapping!)
# TODO:     - build in check if 2th tsv is in sync with 1th tsv
# TODO: is /mnt/db/apollo/reference backuped !??!


```

Some important considerations:
- Typically, you would want to use [apollo-mapping](https://github.com/RIVM-bioinformatics/apollo-mapping) and thus need apollo-reference's data
- apollo_mapping.py (without params) expects the reference data at /mnt/db/apollo/reference
- Adjust this path in apollo_reference/config/referencedata.yaml
- Alternatively mv /your/desired/downloaded/multirefdata/directory to /mnt/db/apollo/reference
- Momentarily, apollo-reference is "hard-coded" connected to a particular excel sheet (data/Apollo_clininal_Candida_species.xlsx)
- In case you want your own set of species, carefully (!) edit this excel sheet

### [3] multi-clade reference low ANI region identification 
```
This section is only partially documented here. Mostly because it's at this moment still a work in progress.
A multi-clade reference in apollo-reference (and apollo-mapping) serves the purpose of:
- when comparing towards the "central" whole-genome assembly (WGA) reference 
- knowing based on WGA to WGA mapping which regions corresponds to deletions (in a/each clade)
- knowing based on PE Illumina mapping (of clade.fastq to central reference.fa) where excess soft-clipping of reads occurs
- these regions taken together are a black-list of called variants in apollo-mapping (called to the central reference)
- and the multi-clade concept triggers apollo-mapping to duplicate mapping workflow (to the central and to the clade reference)

In order for this to work, scenario/workflow *[2] download stated references* needs to habe been executed as:
- ./workflow/workflow-update-reference-assemblies.sh --out  </path/on/hpc/compute> --download-fastq
- this will have downloaded the required (mind: if available!) PE fastq data

needs further instructions ...

```



