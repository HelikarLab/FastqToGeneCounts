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
module load conda
```
| Component |                          Description                          | 
|:---------:|:-------------------------------------------------------------:|
| `module`  | The module command is used to load, unload, and list modules. |
|  `load`   |          The load command is used to load a module.           |
|  `conda`  |  The conda module is used to activate the conda environment.  |

Once this is done, we can create a new conda environment with the name "snakemake". This can be done as follows:
```bash
conda create --name=snakemake
```

|     Component      |                                     Description                                     |
|:------------------:|:-----------------------------------------------------------------------------------:|
|      `conda`       |     The conda command is used to create, update, and remove conda environments.     |
|      `create`      |            The create command is used to create a new conda environment.            |
| `--name=snakemake` | The --name flag is used to specify the name of the conda environment (`snakemake`). |

## Installing software
### Install Mamba
{% include note.html content="Mamba is [recommended by SnakeMake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), as it is much faster than Conda.<br>Additionally, Conda can have issues installing SnakeMake dependencies." %}
```bash
conda install --name snakemake --channel conda-forge mamba
```

|        Component        |                                    Description                                     | 
|:-----------------------:|:----------------------------------------------------------------------------------:|
|         `conda`         |                                 The conda command                                  |
|        `install`        |                  The install command is used to install software                   |
|   `--name snakemake`    | The --name flag is used to specify the name of the conda environment (`snakemake`) |
| `--channel conda-forge` |   The --channel flag is used to specify the channel to install the software from   |
|         `mamba`         |                        The name of the software to install                         |

### Install SnakeMake
SnakeMake is required to run the pipeline.
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
CookieCutter is used to install SnakeMake profiles, which makes it much easier for us to submit our jobs to a cluster.
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
