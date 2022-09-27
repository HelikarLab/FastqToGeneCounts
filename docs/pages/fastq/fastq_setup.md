---
title: Setting up FastqToGeneCounts
sidebar: sidebar
permalink: fastq_setup.html
summary: This is an overview of how to set up the pipeline
last_updated: Sept 22, 2022
---

## Overview
Following the overview is a step-by-step guide on how to set up the pipeline. A brief description is as follows:

1. Create a conda environment containing:
   1. Mamba
   2. Snakemake ([version 6.15.5](https://snakemake.readthedocs.io/en/stable/project_info/history.html#id112))
2. Create a CookieCutter template/profile based on the [Slurm profile](https://github.com/Snakemake-Profiles/slurm)
3. Modification of configuration variables
4. (Optional) Set up [`screen`](https://linux.die.net/man/1/screen) to keep our pipeline running if we close the terminal

## Conda Environment
It should be noted that some of these steps can take a considerable length of time, especially depending on your hardware.

Choose a name for your conda environment. This tutorial will be using `snakemake` for its environment name




### 1. Create the Environment
```bash
conda create -n snakemake
```

### 2. Activate the Environment
```bash
conda activate snakemake
```

### 3. Install Mamba
{% include note.html content="Mamba is [recommended by SnakeMake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), as it is much faster than Conda. Additionally, Conda can have issues installing SnakeMake dependencies." %}
```bash
conda install -n snakemake -c conda-forge mamba
```
   
|    Component     |                          Description                          |
|:----------------:|:-------------------------------------------------------------:|
|     `conda`      |                   The command to run Conda                    |
|    `install`     |               The command to install a package                |
|  `-n snakemake`  | The name of the conda environment to install into (snakemake) |
| `-c conda-forge` |           The channel to install from (Conda Forge)           |
|     `mamba`      |                    The package to install                     |


### 4. Install CookieCutter
```bash
conda install -n snakemake -c conda-forge cookiecutter
```
   
|    Component     |                          Description                          |
|:----------------:|:-------------------------------------------------------------:|
|     `conda`      |                   The command to run Conda                    |
|    `install`     |               The command to install a package                |
|  `-n snakemake`  | The name of the conda environment to install into (snakemake) |
| `-c conda-forge` |           The channel to install from (Conda Forge)           |
|  `cookiecutter`  |                    The package to install                     |

### 5. Test the SnakeMake and CookieCutter installation
```bash
snakemake --version
# Returns `6.15.5`

cookiecutter --version
# Returns Cookiecutter `VERSION` from `INSTALLATION_PATH` (Python `VERSION`)
```

No errors should occur during this installation; if you encounter errors, please [open an issue](https://github.com/HelikarLab/FastqToGeneCounts/issues)

## Creating a CookieCutter Template
This section will assume you are setting up profiles for [Slurm](https://slurm.schedmd.com/documentation.html). If you are not, please reference the [SnakeMake Profile GitHub Page](https://github.com/Snakemake-Profiles/doc) and select your cluster's scheduler.

### CookieCutter Profile Benefits
SnakeMake Profiles provide multiple benefits:
1. No need to write Slurm scripts
2. Jobs are automatically submitted to the scheduler, and killed if they exceed our resources
3. SnakeMake will automatically retry jobs that fail
4. SnakeMake will schedule jobs independently, so that they can run in parallel

Perhaps the most important part of Profiles is Point 4. Instead of requesting, for example, 40 cores, 50GB of RAM, and multiple hours for the entire SnakeMake workflow (as this is the maximum resources we require), Profiles will request small amounts of resources for rules that do not require them. This is especially important for rules that are not CPU intensive, as they take far less time to run. This means we can run more jobs at once, and our jobs will finish faster.

### 1. Creating directories
Create the snakemake directory to hold the profile(s)
```bash
mkdir -p ~/.config/snakemake
```

|       Component       |                Description                 |
|:---------------------:|:------------------------------------------:|
|        `mkdir`        |     The command to create a directory      |
|         `-p`          | Create any parent directories, if required |
| `~/.config/snakemake` |          The directory to create           |

### 2. Creating the Profile
We must first change to the directory we just created, the  we can create the profile
```bash
cd ~/.config/snakemake
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
```

During this section, no details are required. Simply press `Enter` until the profile is created.

### 3. Modifying the Profile
Several details of the profile must be modified to work with our pipeline. These are as follows:
1. `config.yaml` - The configuration file for the profile
2. `cluster_config.yaml` - The configuration file for the cluster

#### 3.1. `config.yaml`
The `config.yaml` is located at `~/.config/snakemake/slurm/config.yaml`. This file contains the configuration for the profile.
Open the file on the cluster, delete its contents, and paste the following:
```yaml
# Default Values
restart-times: 3
jobscript: "slurm-jobscript.sh"
cluster: "slurm-submit.py"
cluster-status: "slurm-status.py"
max-status-checks-per-second: 10
local-cores: 1
latency-wait: 60

# User-modified settings
jobs: 100
printshellcmds: True
max-jobs-per-second: 10  # Default is 1 job per second

# DO NOT MODIFY THESE VALUES
# Snakemake will break if these are changed from their current setting
use-conda: True
conda-frontend: mamba
```

|              Component              |                                                                           Description                                                                            |
|:-----------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|         `restart-times: 3`          |                                                            The number of times to retry a failed job                                                             |
|  `jobscript: "slurm-jobscript.sh"`  |                                                                       The jobscript to use                                                                       |
|    `cluster: "slurm-submit.py"`     |                                                               The cluster submission script to use                                                               |
| `cluster-status: "slurm-status.py"` |                                                                 The cluster status script to use                                                                 |
| `max-status-checks-per-second: 10`  |                                                          The maximum number of status checks per second                                                          |
|          `local-cores: 1`           |                                                            The number of cores to use for local jobs                                                             |
|         `latency-wait: 60`          |                                                The number of seconds to wait before checking the status of a job                                                 |
|             `jobs: 100`             |                                                            The maximum number of jobs to run at once                                                             |
|       `printshellcmds: True`        |                                                            Whether or not to print the shell commands                                                            |
|      `max-jobs-per-second: 10`      |                                                         The maximum number of jobs to submit per second                                                          |
|          `use-conda: True`          |        <span style="color: blue">DO NOT MODIFY</span><br>Whether or not to use Conda to manage dependencies.<br>Modifying this value will break SnakeMake        |
|       `conda-frontend: mamba`       | <span style="color: blue">DO NOT MODIFY</span><br>Use mamba when SnakeMake installs environments.<br>Modifying this value will significantly slow down SnakeMake |


#### 3.2. `cluster_config.yaml`
The `cluster_config.yaml` is located at `~/.config/snakemake/slurm/cluster_config.yaml`. This file contains default configuration values for Slurm to use when submitting jobs to the cluster.
Open the file on the cluster, delete its contents, and paste the following:
```yaml
__default__ :
   job-name  : "{rule}.{wildcards}"
   ntasks    : "1"
   cpus-per-task : "{threads}"
   nodes     : "1"
   output : "logs/{rule}/{rule}.{wildcards.tissue_name}.{wildcards.tag}.output"
   error  : "logs/{rule}/{rule}.{wildcards.tissue_name}.{wildcards.tag}.output"
```

|    Component    |                                     Description                                      |
|:---------------:|:------------------------------------------------------------------------------------:|
|   `job-name`    |                             The default name of the job                              |
|    `ntasks`     |                          The default number of tasks to run                          |
| `cpus-per-task` |       The default number of cores to use.<br>It is set by each SnakeMake rule        |
|     `nodes`     | The number of cluster nodes to request.<br>This should most likely never be changed. |
|    `output`     |                         How (and where) to save output logs                          |
|     `error`     |                          How (and where) to save error logs                          |


