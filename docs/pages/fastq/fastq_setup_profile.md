---
title: Setting up Profiles
sidebar: sidebar
permalink: fastq_setup_profile.html
summary: This is how to set up a slurm profile
last_updated: Oct 11, 2022
---

## Overview
This section will assume you are setting up profiles for [Slurm](https://slurm.schedmd.com/documentation.html). If you are not, please reference the [Snakemake Profile GitHub Page](https://github.com/Snakemake-Profiles/doc) and select your cluster's scheduler.

## Profile Benefits
Snakemake Profiles have multiple benefits:
1. No need to write Slurm scripts 
2. Jobs are automatically submitted to the scheduler, and killed if they exceed our resources 
3. Snakemake will automatically retry jobs that fail 
4. Snakemake will schedule jobs independently, so that they can run in parallel

Perhaps the most important part of Profiles is Point 4. Instead of requesting, for example, 40 cores, 50GB of RAM, and multiple hours for the entire Snakemake workflow (as this is the maximum resources we require), Profiles will request small amounts of resources for rules that do not require them. This is especially important for rules that are not CPU intensive, as they take far less time to run. This means we can run more jobs at once, and our jobs will finish faster.

## Setup
A default profile was downloaded with this repository under the `cluster` directory. If you are not a part of our wonderful Helikar Lab, you must edit the `--account=` line to match your slurm account. This information can be found with the command listed below. In this example, the value listed under `Def Acct` should be included, like this: `--account=helikarlab`
```bash
> sacctmgr show user $USER accounts

      User   Def Acct  Def WCKey     Admin
---------- ---------- ---------- ---------
     joshl helikarlab                 None
```
