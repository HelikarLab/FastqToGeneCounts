---
title: Setting up Conda Environment
sidebar: sidebar
permalink: fastq_setup_conda.html
summary: This is an overview of how to set up the Conda Environment
last_updated: October 11, 2022
---

## Installation
Conda is reuired to install and use FastqToGeneCounts. To install Conda, [follow the instructions here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

{% include note.html content="If you are using HCC, conda is already installed. You can skip this step." %}

## Creating a Conda Environment
In most cluster environments (i.e., HCC), you must activate the `conda` module. This can be done as follows:
```bash
module load mamba
```

| Component |                          Description                          | 
|:---------:|:-------------------------------------------------------------:|
| `module`  | The module command is used to load, unload, and list modules. |
|  `load`   |          The load command is used to load a module.           |
|  `mamba`  |  The conda module is used to activate the conda environment.  |

### Create the Environment
Once this is done, we can create a new conda environment with the name "snakemake". This can be done as follows:
```bash
mamba create --name=snakemake
```

|     Component      |                                     Description                                     |
|:------------------:|:-----------------------------------------------------------------------------------:|
|      `mamba`       |     The conda command is used to create, update, and remove conda environments.     |
|      `create`      |            The create command is used to create a new conda environment.            |
| `--name=snakemake` | The --name flag is used to specify the name of the conda environment (`snakemake`). |

### Activate the Environment
Once you have created the Snakemake environment, you should activate it with the following command:
```bash
mamba activate snakemake
```

### Include Additional Channels
From the [Conda website](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html) on channels:

> Channels are the locations where packages are stored. They serve as the base for hosting and managing packages. Conda packages are downloaded from remote channels, which are URLs to directories containing conda packages. The conda [or mamba, in our case] command seraches a default set of channels and packages are automatically downloaded and updated from [these channels]. 

Because we have additional software required by Snakemake, we must add several channels. These commands only need to be entered once, and you don't have to be in the conda environment to do so. The following commands will add the channels to your conda configuration file:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels r
```

If you get an error during this part, try running the following command:
```bash
module load anaconda
```

Then rerun the `conda config...` commands. Once this is done, re-load mamba (as it was unloaded during `module load anaconda`)
```bash
module load mamba
```

## Installing software
### Install Snakemake
Snakemake is required to run the pipeline.
```bash
mamba install --name snakemake --channel bioconda snakemake==6.15.5
```

|      Component       |                      Description                      |
|:--------------------:|:-----------------------------------------------------:|
|       `mamba`        | Use mamba to install additional software more quickly |
|      `install`       |         The mamba command to install software         |
|  `--name snakemake`  |       The environment to install software into        |
| `--channel bioconda` |         The channel to install software from          |
| `snakemake==6.15.5`  |          The software and version to install          |

### Install CookieCutter
CookieCutter is used to install Snakemake profiles, which makes it much easier for us to submit our jobs to a cluster.
```bash
mamba install --name snakemake --channel conda-forge cookiecutter
```

|        Component        |                      Description                      |
|:-----------------------:|:-----------------------------------------------------:|
|         `mamba`         | Use mamba to install additional software more quickly |
|        `install`        |         The mamba command to install software         |
|   `--name snakemake`    |       The environment to install software into        |
| `--channel conda-forge` |         The channel to install software from          |
|     `cookiecutter`      |                The software to install                |

## Test Installations
```bash
snakemake --version
# Returns "6.15.5"

cookiecutter --version
# Returns Cookiecutter `VERSION` from `INSTALLATION_PATH` (Python `VERSION`)
```
