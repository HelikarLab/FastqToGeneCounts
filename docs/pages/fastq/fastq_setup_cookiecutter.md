---
title: Setting up CookieCutter Profiles
sidebar: sidebar
permalink: fastq_setup_cookiecutter.html
summary: This is how to set up the cookiecutter profile
last_updated: Oct 11, 2022
---

## Overview
This section will assume you are setting up profiles for [Slurm](https://slurm.schedmd.com/documentation.html). If you are not, please reference the [SnakeMake Profile GitHub Page](https://github.com/Snakemake-Profiles/doc) and select your cluster's scheduler.

## CookieCutter Profile Benefits
SnakeMake Profiles have multiple benefits:
1. No need to write Slurm scripts 
2. Jobs are automatically submitted to the scheduler, and killed if they exceed our resources 
3. SnakeMake will automatically retry jobs that fail 
4. SnakeMake will schedule jobs independently, so that they can run in parallel

Perhaps the most important part of Profiles is Point 4. Instead of requesting, for example, 40 cores, 50GB of RAM, and multiple hours for the entire SnakeMake workflow (as this is the maximum resources we require), Profiles will request small amounts of resources for rules that do not require them. This is especially important for rules that are not CPU intensive, as they take far less time to run. This means we can run more jobs at once, and our jobs will finish faster.

## Setup
### Creating Directories
Create the snakemake directory to hold the profile(s)
```bash
profile_dir="${HOME}/.config/snakemake"
mkdir -p "$profile_dir"
```

|         Component         |                                               Description                                               |
|:-------------------------:|:-------------------------------------------------------------------------------------------------------:|
|       `profile_dir`       | The directory to hold the profile(s). `"${HOME}` is an environment variable and does not need to be set |
| `mkdir -p "$profile_dir"` |                                   Create the `$profile_dir` directory                                   |

### Creating the Profile
Using CookieCutter arguments, we can define what profile we would like to install, and where it should be installed. `"$profile_dir` was defined in the previous code block, which is set to our configuration directory (where SnakeMake will look for profiles).

When installing the profile, simply press `Enter` until setup is complete. We will be modifying values manually, later.

```bash
template="gh:Snakemake-Profiles/slurm"
cookiecutter --output-dir "$profile_dir" "$template"
```

|      Component       |                                     Description                                      |
|:--------------------:|:------------------------------------------------------------------------------------:|
|   `template=. . .`   | The template to use when creating the profile. This is the profile we will be using. |
| `cookiecutter . . .` |                          The command to create the profile.                          |

### Modify the Profile
Several details of the profile must be modified to work with our pipeline. THese are as follows:
1. `config.yaml`: THe configuration file for the profile
2. `cluster_config.yaml`: The configuration file for the cluster.

#### `config.yaml`
The `config.yaml` is located at `~/.config/snakemake/slurm/config.yaml` (assuming you used `slurm`, and not a different scheduler). This file contains configuration details for the profile. Open the file, delete its contents, and paste the following

{% include note.html content="If you are on the command line, use `nano`, for its ease of use." %}

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


#### `cluster_config.yaml`
This file contains default configuration values for Slurm to use when submitting jobs to the cluster. We have set two rules to have specific values, and all other rules to use `__default__` values. `get_screen_genomes` and `generate_genome` use specific values because they do not contain wildcards in their output files, and as a result, using `__default__` values will cause an error. To set up these values, create a new file under the `~/.config/snakemake/slurm` directory called `cluster_config.yaml`.

Paste the following contents:  
```yaml
get_screen_genomes:
  job-name: "get_screen_genomes"
  output: "logs/{rule}/{rule}.output"
  error: "logs/{rule}/{rule}.output"

generate_genome:
  job-name: "generate_genome"
  output: "logs/{rule}/{rule}.output"
  error: "logs/{rule}/{rule}.output"

__default__ :
   job-name: "{rule}.{wildcards}"
   ntasks: "1"
   cpus-per-task: "{threads}"
   nodes: "1"
   output: "logs/{rule}/{rule}.{wildcards.tissue_name}.{wildcards.tag}.output"
   error: "logs/{rule}/{rule}.{wildcards.tissue_name}.{wildcards.tag}.output"
```


{% include note.html content="Again, if on the command line, use `nano` for its ease of use" %}

|    Component    |                                     Description                                      |
|:---------------:|:------------------------------------------------------------------------------------:|
|   `job-name`    |                             The default name of the job                              |
|    `ntasks`     |                          The default number of tasks to run                          |
| `cpus-per-task` |       The default number of cores to use.<br>It is set by each SnakeMake rule        |
|     `nodes`     | The number of cluster nodes to request.<br>This should most likely never be changed. |
|    `output`     |                         How (and where) to save output logs                          |
|     `error`     |                          How (and where) to save error logs                          |
