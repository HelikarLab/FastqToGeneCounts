---
title: Setting up FastqToGeneCounts
sidebar: sidebar
permalink: fastq_setup.html
summary: This is an overview of how to set up the pipeline
last_updated: Sept 22, 2022
published: false
---



## Modifying Configuration Variables
The file `snakemake_config.yaml` contains configuration variables for the SnakeMake pipeline. They define various settings such as input CSV files, where to save results, and what rules to perform. The content of the file is as follows:

{% for line in site.data.snakemake_config %}
{% if line[0] == "generation" %}
{% include yaml_key.html content="generation" %}
{% endif %}
{% endfor %}


```yaml
{% for line in site.data.snakemake_config %}
{% if line[0] == "generation" %} 
{{- line[0] }}:
   {% for element in site.data.snakemake_config.generation %}
   {{- element[0] -}}: {{ element[1] }}
   {% endfor %}
{% else %}
{{- line[0] -}}: {{ line[1] }}
{% endif %}
{% endfor %}
```
