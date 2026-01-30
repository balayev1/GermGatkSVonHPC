# GermGatkSVonHPC
An optimized Nextflow implementation of the **GATK-SV pipeline** for HPC clusters utilizing **Slurm**

---

## 📋 Prerequisites

You will need following modules in your environment:
* **Nextflow** (>= 25.10.2)
* **Java** (OpenJDK 17.0.2+)
* **Singularity / Apptainer**
* **Cromwell JAR** (>= v91) — [Download via GitHub](https://github.com/broadinstitute/cromwell/releases)

---

## ⚙️ Setup & Installation

### 1. Prepare GATK-SV Dependencies
Clone the official Broad GATK-SV repository and package the WDL scripts into a zip file. Cromwell requires this `deps.zip` to resolve sub-workflow imports.

```bash
# Clone the official Broad repository
git clone [https://github.com/broadinstitute/gatk-sv.git](https://github.com/broadinstitute/gatk-sv.git)
cd gatk-sv

# Package WDL and script dependencies
zip -r deps.zip wdl/ scripts/
DEPS_ZIP_PATH="$(realpath deps.zip)"
cd ..

# Clone GermGatkSVonHPC github repository
git clone [https://github.com/balayev1/GermGatkSVonHPC.git](https://github.com/balayev1/GermGatkSVonHPC.git)
cd GermGatkSVonHPC
cp "$DEPS_ZIP_PATH" .
```

Since we are using HPC/Slurm to execute the nextflow GATK-SV pipeline, we have to configure Cromwell workflow engine where WDL files are executed using Slurm. For this, we generated `cromwell.config` file to make Cromwell configured to Slurm. 

## 🛠️ Configuration Workflow

The pipeline is organized into three sequential modules. Each requires specific resource allocations in `nextflow.config` and path configurations in the corresponding JSON templates located in the `data/` directory.



### 1. Module: CollectInsertSizeMetrics (Picard) 
* **Manifest Preparation:** Create a tab-separated manifest (see `data/sample.tsv`) containing:
    * `Sample ID`, `BAM/CRAM path`, `BAI/CRAI path`, and `Gender`.
* **Nextflow Config:** Update the parameters for `COLLECT_INSERT_SIZE` block.

### 2. Module: GatherSampleEvidence
* **JSON Template:** Edit `data/GatherSampleEvidence_inputs.json`. 
    * Update paths for required resources (reference genome, contig list and et.c.). To find some of resources, see [link](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json)
    * **Note:** Leave `sample_id`, `bam_or_cram_file`, and `bam_or_cram_index` empty; the pipeline injects these per sample.
* **Nextflow Config:** Adjust `GATHER_SAMPLE_EVIDENCE` block parameters.

### 3. Module: EvidenceQC
* **JSON Template:** Edit `data/EvidenceQC_inputs.json`. This module estimates ploidy and assigns sex across the batch.
    * Update paths for required resources (reference genome, masks, inputs from GatherSampleEvidence and et.c.).
* **Nextflow Config:** Adjust `EVIDENCE_QC` block parameters (typically requires higher memory for joint analysis).

---

## 🚀 Execution

The pipeline is designed to be launched via Slurm script `slurm/nextflow.slurm`.

### 1. Environment Preparation
Ensure to specify required modules and singularity cache directory in your Slurm script (`slurm/nextflow.slurm`):
```bash
module load nextflow
module load java/openjdk-17.0.2
module load singularity

export NXF_SINGULARITY_CACHEDIR="path/to/cachedir"
```

### 2. Adjust the job parameters

### 3. Execute Pipeline
``` bash 
sbatch slurm/nextflow.slurm
```

## 📂 Output Directory Structure
Upon completion, results are moved from the Nextflow `work/` directory to your specified `${params.outdir}`. The structure is organized as follows:
```text
results/
├── CollectInsertSizeMetrics/
│   └── [Sample_ID]/
│
├── GatherSampleEvidence_out/
│   └── [Sample_ID]/
│
└── EvidenceQC_out/
    ├── [results]