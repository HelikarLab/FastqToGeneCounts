---
title: Running FastqToGeneCounts
sidebar: sidebar
permalink: fastq_running.html
summary: This is an overview of how to run the pipeline
last_updated: Sept 22, 2022
---

## Overview
This section goes over how to execute the workflow<br>
The following topics will be covered:
1. (Optional) Setting up Screen
2. Using SnakeMake's dry run
3. Executing the workflow

## (Optional) Using Screen
Unfortunately, Snakemake does not offer a method of closing the terminal while keeping the jobs running. This makes sense, as the main `snakemake --profile cluster` command is tied directly to the main terminal process. To overcome this, we will simply start a Screen session. This allows us to close the main terminal window, while keeping our SSH connection/instance alive.

[Read more about Screen here](https://stackoverflow.com/questions/40527629/)

Alternatively, you can run snakemake in a bash script submitted to SLURM, as explained in [this Google Doc](https://docs.google.com/presentation/d/1gxlxbIObhxitgrPLp7lByYFwrFhEdvEm4mILmygAATY/edit#slide=id.g11366b6085b_0_0)

First, set a large scrollback for Screen, so we can view more lines after we have detached from the terminal. Execute the following:
```bash
echo "defscrollback 10000" >> ~/.screenrc
```

Once this is done, we can start a screen session with the following command:
```bash
screen -S snakemake
```

To leave the screen session while keeping it running, do the following:
1. Press and hold the `control` key
2. Press `a`. **Continue holding control``
3. Press `d`
4. The session will exit. Verify the session is alive, but detached by executing `screen -ls`
   1. It should say `(Detacted)` next to the session name

To re-enter a screen session, execute the following:
```bash
# View all screen sessions
screen -ls

# This will show the following output (if a screen session is running)
> There is a screen on:
>	184700.snakemake	(Detached)
> 1 Socket in /run/screen/S-joshl.

# Pick the session you would like to enter (we are going to re-enter the `snakemake` session)
screen -r snakemake
```


## SnakeMake Dry Run
It is **highly** recommended to run the workflow in dry run mode first to ensure that the workflow will run as expected. A dry run does several things:
1. It checks the syntax of the Snakefile
2. Allows you to see what steps in the workflow will be executed
3. Ensures preliminary configuration is set up properly

A dry run does not truly execute any components of the pipeline. No results will be generated

Execute the following to perform a dry run
```bash
# Activate our conda environment
module load mamba
mamba activate snakemake

# Change to the FastqToGeneCounts directory
cd /work/helikarlab/joshl/FastqToGeneCounts

# Perfom a dry run
snakemake --profile cluster --dry-run
```

{% include note.html content="If you did renamed the `cluster` directory to something else, replace the `--profile cluster` with the name of your directory" %}
{% include note.html content="If you receive an error when running `snakemake --profile cluster --dry-run`, replcae `slurm` with `./cluster`" %}

After several seconds, many lines should move through the terminal.<br>
It should end with `This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.`

If this is not the case, an error has occured, and it will need to be investigated before continuing. If you are having troubles, please [Open an Issue](https://github.com/HelikarLab/FastqToGeneCounts/issues)

## Execution
Once you have confirmed that a dry-run will execute successfully, it is time to start a real run of the workflow.<br>

{{site.data.alerts.note}}
If you have started a screen session, now is the time to re-enter the session<br><br>
To see what sessions are available:
<pre>screen -ls</pre>

To re-enter a session:
<pre>screen -r SESSION_NAME</pre>

{{site.data.alerts.end}}

The following steps will start the workflow:
```bash
# Activate the snakemake environment
module load mamba
mamba activate snakemake

# Make sure you are in the FastqToGeneCounts directory!
cd /work/helikarlab/joshl/FastqToGeneCounts

# Start the workflow
snakemake --profile cluster
```

{% include note.html content="If you started a session with screen, exit the session with `CTRL+a`, `d`" %}

<br>
Any log files will be found in the `logs` directory of the project directory.<br>
Each rule has its own output folder, with output files containing the information they are running on (tissue name, run number, etc.)
