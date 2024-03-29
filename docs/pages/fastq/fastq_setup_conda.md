---
title: Setting up Conda Environment
sidebar: sidebar
permalink: fastq_setup_conda.html
summary: This is an overview of how to set up the Conda Environment
last_updated: October 11, 2022
---

## Installation
Conda is reuired to install and use FastqToGeneCounts. To install Conda, [follow the instructions here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

{% include note.html content="If you are using HCC, mamba is already installed. You can skip this step." %}

## Creating a Conda Environment
In most cluster environments (i.e., HCC), you must activate the `mamba` module. This can be done in a terminal as follows:
```bash
module load mamba
```

| Component |                          Description                          | 
|:---------:|:-------------------------------------------------------------:|
| `module`  | The module command is used to load, unload, and list modules. |
|  `load`   |          The load command is used to load a module.           |
|  `mamba`  |  The mamba module is used to activate the mamba environment.  |

### Create the Environment
Once this is done, we can create a new conda environment with the name "snakemake". This can be done as follows. In a terminal:
```bash
mamba create --name=snakemake
```

|     Component      |                                     Description                                     |
|:------------------:|:-----------------------------------------------------------------------------------:|
|      `mamba`       |     The mamba command is used to create, update, and remove mamba environments.     |
|      `create`      |            The create command is used to create a new conda environment.            |
| `--name=snakemake` | The --name flag is used to specify the name of the conda environment (`snakemake`). |

### Activate the Environment
Once you have created the Snakemake environment, you should activate it with the following command:
```bash
mamba activate snakemake
```

#### Initializing Mamba

{% include important.html content="If you get a warning saying you must run `mamba init`, follow the following steps" %}

```bash
mamba init
source ~/.bashrc  # If you are using bash, otherwise use the appropriate shell
```

|     Component      |                                 Description                                 |
|:------------------:|:---------------------------------------------------------------------------:|
|      `mamba`       | The mamba command is used to create, update, and remove mamba environments. |
|       `init`       |                The init command is used to initialize mamba.                |
| `source ~/.bashrc` |   The `source` command is used to run the commands in the `.bashrc` file.   |

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

----

If you get an error during this part, try running the following command:
```bash
module load anaconda
```

Then rerun the `conda config...` commands. Once this is done, re-load mamba (as it was unloaded during `module load anaconda`)
```bash
module load mamba
```

## Installing software
### Install Snakemake and Benchmarking Requirements
Snakemake is required to run the pipeline.
```bash
mamba install --name snakemake --channel bioconda --channel plotly snakemake==6.15.5 tabulate==0.8.10 plotly==5.11.0
```

We must install tabulate version `0.8.10` as anything under the `0.9.*` release causes issues for our current version of Snakemake

|      Component       |                      Description                      |
|:--------------------:|:-----------------------------------------------------:|
|       `mamba`        | Use mamba to install additional software more quickly |
|      `install`       |         The mamba command to install software         |
|  `--name snakemake`  |       The environment to install software into        |
| `--channel bioconda` |         The channel to install software from          |
|  `--channel plotly`  |            Install from the plotly channel            |
| `snakemake==6.15.5`  |          The software and version to install          |
|  `tabulate==0.8.10`  |          The software and version to install          |
|   `plotly==5.11.0`   |           Install plotly for graph creation           |

## Test Installations
```bash
snakemake --version
# Returns "7.25.0"
```
