---
title: Downloading FastqToGeneCounts
sidebar: sidebar
permalink: fastq_download.html
summary: This is an overview of how to download the pipeline
last_updated: Sept 22, 2022
---

# Overview
Downloading the pipeline is as simple as using Git, or using the GitHub website.<br>
Several things will be downloaded during this section:
1. The pipeline itself
2. [Conda](https://conda.io/) (and [Mamba](https://github.com/mamba-org/mamba))

## Pipeline Download
### First-time Download
#### Using Git
1. In a terminal, log in to the cluster (if you are using one)
2. Navigate to the directory where you would like to download the pipeline (i.e., `cd {{ site.data.terms.work_dir }}`)
3. Execute the following command to download the pipeline:
```bash
git clone https://github.com/HelikarLab/FastqToGeneCounts.git
```

{{ site.data.alerts.note }}
<p>Optionally, define the destination folder with:</p> 
<pre>
git clone https://github.com/HelikarLab/FastqToGeneCounts.git DESTINATION_FOLDER
</pre>
{{ site.data.alerts.end }}

You can now navigate into the downloaded directory and run the pipeline.

#### Using GitHub
1. Navigate to the [GitHub repository](https://github.com/HelikarLab/FastqToGeneCounts)
2. Click the green "Code" button, which looks like this:
    ![GitHub Code Button](/images/code_icon.png)
3. Click `Download ZIP`, and follow the prompts (if any) to download the pipeline.
4. Unzip the downloaded file
5. Once this is done, you must upload this directory onto the cluster (if you are using one)
   1. It is recommended to use a cluster, as the [STAR aligner](https://github.com/alexdobin/STAR) has high memory requirements.
   2. You should upload into the appropriate "work" directory (i.e., `{{ site.data.terms.work_dir }}`)
6. You can now navigate into the directory (on the cluster) and run the pipeline

### Updating the Pipeline
1. In a terminal, navigate to the directory where you downloaded the pipeline.
2. Execute the following command to update the pipeline:
```bash
git pull
```

Unfortunately, there is no easy method of updating the pipeline strictly using the GitHub website. It is recommended to use Git to update.

If any errors occur during the update process, please [open an issue](https://github.com/HelikarLab/FastqToGeneCounts/issues) on our GitHub page for assistance.

## Conda Download
Nearly every cluster has conda installed. If you are running this in an environment that does not have conda installed (such as a local computer), please follow the instructions here to download and install conda for your system: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

This installation will use MiniConda, which is a slim version of Conda, capable of installing the necessary software for the pipeline.
