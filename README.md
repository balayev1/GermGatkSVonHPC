# GermGatkSVonHPC
This is HPC version of GATK-SV pipeline

## Setup
- Clone gatk-sv github directory using:
``` bash
git clone https://github.com/broadinstitute/gatk-sv.git
```
This directory contains all the necessary WDL files to execute the pipeline.

- Zip all wdl scripts from gatk-sv directory:
```bash
cd gatk-sv
zip -r deps.zip wdl/ scripts/
cd ~
```

_ Clone this github repository:
``` bash
git clone https://github.com/balayev1/GermGatkSVonHPC.git
cd GermGatkSVonHPC
```

- Since we are using HPC/Slurm to execute the pipeline, we have to configure Cromwell workflow engine where WDL files are executed to Slurm. For this, we generated `cromwell.config` file to make Cromwell configured to slurm. 

- You also need to download cromwell executable jar file (see https://github.com/broadinstitute/cromwell/releases)

## Implementation
- First module of the pipeline is GatherSampleEvidence to gather per-sample structural variant evidence. Prior to running the script,

I) prepare sample information file (see example in data/sample.tsv). 

II) Adjust the paths in `GatherSampleEvidence_inputs.json` file for all the parameters (see https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json to download files). 
* Don't specify sample_id, bam_or_cram_file and bam_or_cram_index slots because they will be automatically read from your sample information file.  

III) Adjust all the variables commented by the arrow in `execute_GatherSampleEvidence_sbatch.sh` file.

IV) Run:
``` bash 
bash execute_GatherSampleEvidence_sbatch.sh
```